# compare_sampling_GFS_ERA5.py
#
# Compare actual sampling + interpolation behavior between GFS and ERA5 classes,
# while allowing different native time resolutions.
#
# Usage:
#   python -m tests.compare_sampling_GFS_ERA5
#
# Notes:
# - Avoids cftime calendar subtraction errors by converting timestamps to plain
#   Python datetimes for comparison math.
# - Preserves native timestamps for each class when calling getNewCoord().
# - Compares actual elapsed real hours instead of requiring the same native
#   time index value between datasets with different cadence.
# - Tracks warning counts and prints a summary at the end.

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, asdict
from datetime import datetime
import traceback
from typing import Any, Dict, Tuple

import numpy as np

from GFS import GFS
from ERA5 import ERA5
import config_earth


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DT_FOR_GETNEWCOORD = 60.0  # seconds for one-step trajectory compare

SAMPLE_LATITUDES = [26.0, 30.0, 35.0, 40.0, 44.0]
SAMPLE_LONGITUDES = [-136.0, -120.0, -100.0, -80.0, -77.0]
SAMPLE_ALTITUDES = [100.0, 1000.0, 5000.0, 10000.0, 20000.0]

NUM_TIME_SAMPLES = 5

UV_WARN_ATOL = 0.15
COORD_WARN_ATOL_DEG = 1e-4 # ~11 meters when converted from latitude
ELAPSED_HOURS_WARN_ATOL = 1e-9  # shared timestamps should align exactly in real time


# -----------------------------------------------------------------------------
# Global warning counters
# -----------------------------------------------------------------------------
WARNING_COUNTS = defaultdict(int)
TOTAL_CHECKS = defaultdict(int)


# -----------------------------------------------------------------------------
# Dataclasses
# -----------------------------------------------------------------------------
@dataclass
class SamplePoint:
    timestamp_py: Any
    timestamp_gfs: Any
    timestamp_era: Any
    lat: float
    lon: float
    alt: float


@dataclass
class SamplerDebug:
    source: str
    timestamp: Any

    # real elapsed time from that source's model start
    elapsed_hours: float

    # native index space for that source
    native_index: float
    t0_idx: int
    t1_idx: int

    lat_idx: int
    lon_idx: int
    chosen_lat: float
    chosen_lon: float
    z_nearest_idx: int
    z_nearest_value: float

    u0_pt: float
    v0_pt: float
    u1_pt: float
    v1_pt: float

    u_old: float
    v_old: float
    u_new: float
    v_new: float


@dataclass
class OneStepDebug:
    source: str
    lat_new: float
    lon_new: float
    x_wind_vel: float
    y_wind_vel: float
    x_wind_vel_old: float
    y_wind_vel_old: float
    bearing: float
    closest_lat: float
    closest_lon: float
    closest_alt_value: float


# -----------------------------------------------------------------------------
# Time helpers
# -----------------------------------------------------------------------------
def to_py_datetime(t):
    """
    Convert cftime / netCDF4 / python datetime-like objects to a plain python datetime.
    """
    if isinstance(t, datetime):
        return t.replace(tzinfo=None)

    return datetime(
        int(t.year),
        int(t.month),
        int(t.day),
        int(t.hour),
        int(t.minute),
        int(t.second),
        int(getattr(t, "microsecond", 0)),
    )


def seconds_since(start, end):
    """
    Return elapsed seconds between two datetime-like objects after converting
    both to plain python datetimes.
    """
    start_dt = to_py_datetime(start)
    end_dt = to_py_datetime(end)
    return (end_dt - start_dt).total_seconds()


def build_common_times(gfs: GFS, era: ERA5, n: int):
    """
    Return common times as triples:
        (python_datetime, gfs_native_time, era_native_time)
    """
    gfs_pairs = [(to_py_datetime(t), t) for t in gfs.time_convert]
    era_pairs = [(to_py_datetime(t), t) for t in era.time_convert]

    era_lookup = {py_t: native_t for py_t, native_t in era_pairs}

    common = []
    for py_t, gfs_native in gfs_pairs:
        if py_t in era_lookup:
            common.append((py_t, gfs_native, era_lookup[py_t]))

    common.sort(key=lambda x: x[0])

    if not common:
        raise RuntimeError("No shared timestamps found between GFS and ERA5.")

    return common[:min(n, len(common))]


# -----------------------------------------------------------------------------
# Array/profile helpers
# -----------------------------------------------------------------------------
def fill_missing_1d(data) -> np.ndarray:
    """
    Make a plain float ndarray from a masked or regular 1D array and fill NaNs linearly.
    """
    arr = np.asanyarray(data)

    if np.ma.isMaskedArray(arr):
        arr = arr.filled(np.nan)

    arr = np.asarray(arr, dtype=float)

    if np.isnan(arr).any():
        nans = np.isnan(arr)
        good = ~nans
        if good.sum() == 0:
            raise ValueError("All values are NaN in profile.")
        if good.sum() == 1:
            arr[nans] = arr[good][0]
        else:
            idx = np.arange(arr.size)
            arr[nans] = np.interp(idx[nans], idx[good], arr[good])

    return arr


def sorted_unique_profile(h, u, v) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Sort profile by height ascending and drop duplicate heights.
    """
    h = fill_missing_1d(h)
    u = fill_missing_1d(u)
    v = fill_missing_1d(v)

    order = np.argsort(h)
    h = h[order]
    u = u[order]
    v = v[order]

    h_unique, unique_idx = np.unique(h, return_index=True)
    u_unique = u[unique_idx]
    v_unique = v[unique_idx]

    return h_unique, u_unique, v_unique


def safe_interp_uv(h, u, v, alt: float) -> Tuple[float, float, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Robust altitude interpolation on a vertical profile.
    """
    h, u, v = sorted_unique_profile(h, u, v)

    alt_clamped = float(np.clip(alt, h[0], h[-1]))
    u_pt = float(np.interp(alt_clamped, h, u))
    v_pt = float(np.interp(alt_clamped, h, v))
    return u_pt, v_pt, h, u, v, alt_clamped


def normalize_lon_for_compare(lon: float) -> float:
    """
    Normalize longitude to [-180, 180).
    """
    return ((lon + 180.0) % 360.0) - 180.0


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def fmt_float(x: Any, digits: int = 6) -> str:
    try:
        return f"{float(x):.{digits}f}"
    except Exception:
        return str(x)


# -----------------------------------------------------------------------------
# GFS sampler debug
# -----------------------------------------------------------------------------
def debug_sample_gfs(gfs: GFS, sample: SamplePoint) -> SamplerDebug:
    """
    Reproduce the actual GFS sampling logic as closely as possible.
    """
    elapsed_sec = seconds_since(gfs.gfs_time, sample.timestamp_py)
    elapsed_hours = elapsed_sec / 3600.0

    # Native index space for current GFS code path
    native_index = elapsed_hours / 3.0
    t0 = int(native_index)
    t1 = t0 + 1

    lat_idx = gfs.getNearestLat(sample.lat, gfs.LAT_LOW, gfs.LAT_HIGH)
    lon_idx = gfs.getNearestLon(sample.lon, gfs.LON_LOW, gfs.LON_HIGH)

    z_idx = gfs.getNearestAlt(native_index, sample.lat, sample.lon, sample.alt)

    v_0 = gfs.vgdrps0[t0, :, lat_idx, lon_idx]
    u_0 = gfs.ugdrps0[t0, :, lat_idx, lon_idx]
    h0 = gfs.hgtprs[t0, :, lat_idx, lon_idx]
    u0_pt, v0_pt, h0s, u0s, v0s, alt0 = safe_interp_uv(h0, u_0, v_0, sample.alt)

    v_1 = gfs.vgdrps0[t1, :, lat_idx, lon_idx]
    u_1 = gfs.ugdrps0[t1, :, lat_idx, lon_idx]
    h1 = gfs.hgtprs[t1, :, lat_idx, lon_idx]
    u1_pt, v1_pt, h1s, u1s, v1s, alt1 = safe_interp_uv(h1, u_1, v_1, sample.alt)

    fp = [t0, t1]
    u_old = float(np.interp(native_index, fp, [u0_pt, u1_pt]))
    v_old = float(np.interp(native_index, fp, [v0_pt, v1_pt]))

    bearing_t0, speed_t0 = gfs.interpolateBearing(h0s, u0s, v0s, alt0)
    bearing_t1, speed_t1 = gfs.interpolateBearing(h1s, u1s, v1s, alt1)
    bearing_interpolated, speed_interpolated = gfs.interpolateBearingTime(
        bearing_t0, speed_t0, bearing_t1, speed_t1, native_index
    )
    u_new = float(speed_interpolated * np.cos(np.radians(bearing_interpolated)))
    v_new = float(speed_interpolated * np.sin(np.radians(bearing_interpolated)))

    z_profile = fill_missing_1d(gfs.hgtprs[t0, :, lat_idx, lon_idx])
    z_nearest_value = float(z_profile[z_idx])

    return SamplerDebug(
        source="GFS",
        timestamp=sample.timestamp_py,
        elapsed_hours=float(elapsed_hours),
        native_index=float(native_index),
        t0_idx=t0,
        t1_idx=t1,
        lat_idx=int(lat_idx),
        lon_idx=int(lon_idx),
        chosen_lat=float(gfs.lat[lat_idx]),
        chosen_lon=float(normalize_lon_for_compare(gfs.lon[lon_idx])),
        z_nearest_idx=int(z_idx),
        z_nearest_value=z_nearest_value,
        u0_pt=u0_pt,
        v0_pt=v0_pt,
        u1_pt=u1_pt,
        v1_pt=v1_pt,
        u_old=u_old,
        v_old=v_old,
        u_new=u_new,
        v_new=v_new,
    )


# -----------------------------------------------------------------------------
# ERA5 sampler debug
# -----------------------------------------------------------------------------
def debug_sample_era5(era: ERA5, sample: SamplePoint) -> SamplerDebug:
    """
    Reproduce the actual ERA5 sampling logic as closely as possible.
    """
    elapsed_sec = seconds_since(era.model_start_datetime, sample.timestamp_py)
    elapsed_hours = elapsed_sec / 3600.0

    native_index = elapsed_hours / era.resolution_hr
    t0 = int(native_index)
    t1 = t0 + 1

    lat_idx = era.getNearestLatIdx(sample.lat, era.lat_max_idx, era.lat_min_idx)
    lon_idx = era.getNearestLonIdx(sample.lon, era.lon_min_idx, era.lon_max_idx)
    z_idx = era.getNearestAltbyIndex(t0, lat_idx, lon_idx, sample.alt)

    v_0 = era.vgdrps0[t0, :, lat_idx, lon_idx]
    u_0 = era.ugdrps0[t0, :, lat_idx, lon_idx]
    h0 = era.hgtprs[t0, :, lat_idx, lon_idx]
    u0_pt, v0_pt, h0s, u0s, v0s, alt0 = safe_interp_uv(h0, u_0, v_0, sample.alt)

    v_1 = era.vgdrps0[t1, :, lat_idx, lon_idx]
    u_1 = era.ugdrps0[t1, :, lat_idx, lon_idx]
    h1 = era.hgtprs[t1, :, lat_idx, lon_idx]
    u1_pt, v1_pt, h1s, u1s, v1s, alt1 = safe_interp_uv(h1, u_1, v_1, sample.alt)

    fp = [t0, t1]
    u_old = float(np.interp(native_index, fp, [u0_pt, u1_pt]))
    v_old = float(np.interp(native_index, fp, [v0_pt, v1_pt]))

    bearing_t0, speed_t0 = era.interpolateBearing(h0s, u0s, v0s, alt0)
    bearing_t1, speed_t1 = era.interpolateBearing(h1s, u1s, v1s, alt1)
    bearing_interpolated, speed_interpolated = era.interpolateBearingTime(
        bearing_t0, speed_t0, bearing_t1, speed_t1, native_index
    )
    u_new = float(speed_interpolated * np.cos(np.radians(bearing_interpolated)))
    v_new = float(speed_interpolated * np.sin(np.radians(bearing_interpolated)))

    z_profile = fill_missing_1d(era.hgtprs[t0, :, lat_idx, lon_idx])
    z_nearest_value = float(z_profile[z_idx])

    return SamplerDebug(
        source="ERA5",
        timestamp=sample.timestamp_py,
        elapsed_hours=float(elapsed_hours),
        native_index=float(native_index),
        t0_idx=t0,
        t1_idx=t1,
        lat_idx=int(lat_idx),
        lon_idx=int(lon_idx),
        chosen_lat=float(era.lat[lat_idx]),
        chosen_lon=float(normalize_lon_for_compare(era.lon[lon_idx])),
        z_nearest_idx=int(z_idx),
        z_nearest_value=z_nearest_value,
        u0_pt=u0_pt,
        v0_pt=v0_pt,
        u1_pt=u1_pt,
        v1_pt=v1_pt,
        u_old=u_old,
        v_old=v_old,
        u_new=u_new,
        v_new=v_new,
    )


# -----------------------------------------------------------------------------
# One-step getNewCoord debug
# -----------------------------------------------------------------------------
def debug_one_step_gfs(gfs: GFS, sample: SamplePoint, dt: float) -> OneStepDebug:
    out = gfs.getNewCoord(
        {"timestamp": sample.timestamp_gfs, "lat": sample.lat, "lon": sample.lon, "alt": sample.alt},
        dt,
    )
    return OneStepDebug(
        source="GFS",
        lat_new=float(out[0]),
        lon_new=float(normalize_lon_for_compare(out[1])),
        x_wind_vel=float(out[2]),
        y_wind_vel=float(out[3]),
        x_wind_vel_old=float(out[4]),
        y_wind_vel_old=float(out[5]),
        bearing=float(out[6]),
        closest_lat=float(out[7]),
        closest_lon=float(normalize_lon_for_compare(out[8])),
        closest_alt_value=float(out[9]),
    )


def debug_one_step_era5(era: ERA5, sample: SamplePoint, dt: float) -> OneStepDebug:
    out = era.getNewCoord(
        {"timestamp": sample.timestamp_era, "lat": sample.lat, "lon": sample.lon, "alt": sample.alt},
        dt,
    )
    return OneStepDebug(
        source="ERA5",
        lat_new=float(out[0]),
        lon_new=float(normalize_lon_for_compare(out[1])),
        x_wind_vel=float(out[2]),
        y_wind_vel=float(out[3]),
        x_wind_vel_old=float(out[4]),
        y_wind_vel_old=float(out[5]),
        bearing=float(out[6]),
        closest_lat=float(out[7]),
        closest_lon=float(normalize_lon_for_compare(out[8])),
        closest_alt_value=float(out[9]),
    )


# -----------------------------------------------------------------------------
# Comparisons
# -----------------------------------------------------------------------------
def compare_sampler_debug(a: SamplerDebug, b: SamplerDebug) -> Dict[str, float]:
    return {
        # apples-to-apples time comparison
        "elapsed_hours_abs_diff": abs(a.elapsed_hours - b.elapsed_hours),

        # native index spaces may legitimately differ
        "native_index_abs_diff": abs(a.native_index - b.native_index),

        "chosen_lat_abs_diff": abs(a.chosen_lat - b.chosen_lat),
        "chosen_lon_abs_diff": abs(normalize_lon_for_compare(a.chosen_lon - b.chosen_lon)),
        "z_nearest_value_abs_diff": abs(a.z_nearest_value - b.z_nearest_value),

        "u0_pt_abs_diff": abs(a.u0_pt - b.u0_pt),
        "v0_pt_abs_diff": abs(a.v0_pt - b.v0_pt),
        "u1_pt_abs_diff": abs(a.u1_pt - b.u1_pt),
        "v1_pt_abs_diff": abs(a.v1_pt - b.v1_pt),

        "u_old_abs_diff": abs(a.u_old - b.u_old),
        "v_old_abs_diff": abs(a.v_old - b.v_old),
        "u_new_abs_diff": abs(a.u_new - b.u_new),
        "v_new_abs_diff": abs(a.v_new - b.v_new),
    }


def compare_one_step(a: OneStepDebug, b: OneStepDebug) -> Dict[str, float]:
    return {
        "lat_new_abs_diff": abs(a.lat_new - b.lat_new),
        "lon_new_abs_diff": abs(normalize_lon_for_compare(a.lon_new - b.lon_new)),
        "x_wind_vel_abs_diff": abs(a.x_wind_vel - b.x_wind_vel),
        "y_wind_vel_abs_diff": abs(a.y_wind_vel - b.y_wind_vel),
        "x_wind_vel_old_abs_diff": abs(a.x_wind_vel_old - b.x_wind_vel_old),
        "y_wind_vel_old_abs_diff": abs(a.y_wind_vel_old - b.y_wind_vel_old),
        "bearing_abs_diff": abs(a.bearing - b.bearing),
        "closest_lat_abs_diff": abs(a.closest_lat - b.closest_lat),
        "closest_lon_abs_diff": abs(normalize_lon_for_compare(a.closest_lon - b.closest_lon)),
        "closest_alt_value_abs_diff": abs(a.closest_alt_value - b.closest_alt_value),
    }


# -----------------------------------------------------------------------------
# Reporting
# -----------------------------------------------------------------------------
def print_sampler_report(sample: SamplePoint, gfs_dbg: SamplerDebug, era_dbg: SamplerDebug) -> None:
    diffs = compare_sampler_debug(gfs_dbg, era_dbg)

    print_section(
        f"SAMPLE: time={sample.timestamp_py}, lat={sample.lat}, lon={sample.lon}, alt={sample.alt}"
    )

    print("GFS sampler debug:")
    for k, v in asdict(gfs_dbg).items():
        print(f"  {k}: {v}")

    print("\nERA5 sampler debug:")
    for k, v in asdict(era_dbg).items():
        print(f"  {k}: {v}")

    print("\nSampler diffs:")
    for k, v in diffs.items():
        TOTAL_CHECKS[k] += 1
        flag = ""

        if k == "elapsed_hours_abs_diff":
            if v > ELAPSED_HOURS_WARN_ATOL:
                flag = "  <-- WARN"
                WARNING_COUNTS[k] += 1

        elif k == "native_index_abs_diff":
            # informational only
            pass

        elif ("u_" in k or "v_" in k):
            if v > UV_WARN_ATOL:
                flag = "  <-- WARN"
                WARNING_COUNTS[k] += 1

        elif ("lat" in k or "lon" in k):
            if v > COORD_WARN_ATOL_DEG:
                flag = "  <-- WARN"
                WARNING_COUNTS[k] += 1

        print(f"  {k}: {fmt_float(v, 8)}{flag}")


def print_one_step_report(sample: SamplePoint, gfs_step: OneStepDebug, era_step: OneStepDebug) -> None:
    diffs = compare_one_step(gfs_step, era_step)

    print("\nOne-step getNewCoord debug:")
    print("GFS:")
    for k, v in asdict(gfs_step).items():
        print(f"  {k}: {v}")

    print("\nERA5:")
    for k, v in asdict(era_step).items():
        print(f"  {k}: {v}")

    print("\nOne-step diffs:")
    for k, v in diffs.items():
        key = f"step_{k}"
        TOTAL_CHECKS[key] += 1

        flag = ""
        if ("wind" in k and v > UV_WARN_ATOL) or (("lat" in k or "lon" in k) and v > COORD_WARN_ATOL_DEG):
            flag = "  <-- WARN"
            WARNING_COUNTS[key] += 1

        print(f"  {k}: {fmt_float(v, 8)}{flag}")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    print_section("INSTANTIATING GFS / ERA5")
    gfs = GFS(config_earth.simulation["start_coord"])
    era = ERA5(config_earth.simulation["start_coord"])

    common_times = build_common_times(gfs, era, NUM_TIME_SAMPLES)

    print("\nCommon times being tested:")
    for timestamp_py, _, _ in common_times:
        print(" ", timestamp_py)

    summary_rows = []
    total_samples = 0
    failed_samples = 0

    for timestamp_py, timestamp_gfs, timestamp_era in common_times:
        for lat in SAMPLE_LATITUDES:
            for lon in SAMPLE_LONGITUDES:
                for alt in SAMPLE_ALTITUDES:
                    total_samples += 1
                    sample = SamplePoint(
                        timestamp_py=timestamp_py,
                        timestamp_gfs=timestamp_gfs,
                        timestamp_era=timestamp_era,
                        lat=lat,
                        lon=lon,
                        alt=alt,
                    )

                    try:
                        gfs_dbg = debug_sample_gfs(gfs, sample)
                        era_dbg = debug_sample_era5(era, sample)

                        gfs_step = debug_one_step_gfs(gfs, sample, DT_FOR_GETNEWCOORD)
                        era_step = debug_one_step_era5(era, sample, DT_FOR_GETNEWCOORD)

                        print_sampler_report(sample, gfs_dbg, era_dbg)
                        print_one_step_report(sample, gfs_step, era_step)

                        sampler_diffs = compare_sampler_debug(gfs_dbg, era_dbg)
                        step_diffs = compare_one_step(gfs_step, era_step)

                        summary_rows.append({
                            "timestamp": timestamp_py,
                            "lat": lat,
                            "lon": lon,
                            "alt": alt,
                            "elapsed_hours_abs_diff": sampler_diffs["elapsed_hours_abs_diff"],
                            "native_index_abs_diff": sampler_diffs["native_index_abs_diff"],
                            "u_new_abs_diff": sampler_diffs["u_new_abs_diff"],
                            "v_new_abs_diff": sampler_diffs["v_new_abs_diff"],
                            "u_old_abs_diff": sampler_diffs["u_old_abs_diff"],
                            "v_old_abs_diff": sampler_diffs["v_old_abs_diff"],
                            "lat_new_abs_diff": step_diffs["lat_new_abs_diff"],
                            "lon_new_abs_diff": step_diffs["lon_new_abs_diff"],
                        })

                    except Exception as exc:
                        failed_samples += 1
                        print_section(
                            f"FAILED SAMPLE: time={timestamp_py}, lat={lat}, lon={lon}, alt={alt}"
                        )
                        print("Exception:", exc)
                        traceback.print_exc()

    print_section("SUMMARY")
    print(f"Total samples attempted: {total_samples}")
    print(f"Failed samples: {failed_samples}")
    print(f"Successful samples: {len(summary_rows)}")

    if summary_rows:
        elapsed_hours_max = max(r["elapsed_hours_abs_diff"] for r in summary_rows)
        native_index_max = max(r["native_index_abs_diff"] for r in summary_rows)
        u_new_max = max(r["u_new_abs_diff"] for r in summary_rows)
        v_new_max = max(r["v_new_abs_diff"] for r in summary_rows)
        u_old_max = max(r["u_old_abs_diff"] for r in summary_rows)
        v_old_max = max(r["v_old_abs_diff"] for r in summary_rows)
        lat_new_max = max(r["lat_new_abs_diff"] for r in summary_rows)
        lon_new_max = max(r["lon_new_abs_diff"] for r in summary_rows)

        print("\nWorst-case diffs across all successful samples:")
        print(f"  max elapsed_hours abs diff: {elapsed_hours_max:.12f}")
        print(f"  max native_index abs diff: {native_index_max:.8f}  (informational)")
        print(f"  max u_new abs diff: {u_new_max:.8f}")
        print(f"  max v_new abs diff: {v_new_max:.8f}")
        print(f"  max u_old abs diff: {u_old_max:.8f}")
        print(f"  max v_old abs diff: {v_old_max:.8f}")
        print(f"  max lat_new abs diff: {lat_new_max:.10f}")
        print(f"  max lon_new abs diff: {lon_new_max:.10f}")

        scored = sorted(
            summary_rows,
            key=lambda r: max(
                r["u_new_abs_diff"],
                r["v_new_abs_diff"],
                r["lat_new_abs_diff"] * 1e5,
                r["lon_new_abs_diff"] * 1e5,
            ),
            reverse=True,
        )

        print("\nTop 10 worst samples:")
        for row in scored[:10]:
            print(row)

    print_section("WARNING SUMMARY")

    if TOTAL_CHECKS:
        keys = sorted(TOTAL_CHECKS.keys())
        for k in keys:
            total = TOTAL_CHECKS[k]
            warns = WARNING_COUNTS.get(k, 0)
            pct = (warns / total * 100.0) if total > 0 else 0.0
            print(f"{k:30s} : {warns:6d} / {total:6d}  ({pct:6.2f}%)")
    else:
        print("No checks recorded.")


if __name__ == "__main__":
    main()