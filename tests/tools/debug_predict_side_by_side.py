# debug_predict_side_by_side.py
#
# Run the prediction loop side-by-side with both GFS and ERA5 and compare:
# - getNewCoord() outputs
# - one-step propagated lat/lon
# - running trajectory separation
#
# This mirrors predict.py, but instantiates both forecast models and keeps
# the vertical state shared between them.
#
# Usage:
#   python -m tests.tools.debug_predict_side_by_side
#
# Optional knobs near the top:
#   CONTINUE_WITH = "GFS"   # or "ERA5"
#   MAX_STEPS = None        # or an integer for shorter debug runs
#
# Notes:
# - This is a debug script, not a replacement for predict.py
# - It uses the same vertical solver and same update cadence as predict.py
# - It feeds BOTH models the same coord/time/alt input before either one
#   can diverge from the other

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
import math
import traceback
from typing import Any, Dict, List, Optional

import pandas as pd
from geographiclib.geodesic import Geodesic
from termcolor import colored
import fluids

import config_earth
import solve_states
import GFS
import ERA5


# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
CONTINUE_WITH = "GFS"   # "GFS" or "ERA5"
MAX_STEPS = None        # e.g. 500 for short debug runs; None for full sim

# Warning thresholds
UV_WARN_ATOL = 0.15               # m/s
COORD_WARN_ATOL_DEG = 1e-5        # ~1.1 m latitude
STEP_SEPARATION_WARN_M = 2.0      # one-step endpoint separation
RUNNING_SEPARATION_WARN_M = 10.0  # running traj separation

PRINT_EVERY_HORIZONTAL_UPDATE = True
PRINT_ONLY_WARNINGS = False


# -----------------------------------------------------------------------------
# Counters
# -----------------------------------------------------------------------------
WARNING_COUNTS = defaultdict(int)
TOTAL_CHECKS = defaultdict(int)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def to_py_datetime(t):
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


def normalize_lon_for_compare(lon: float) -> float:
    return ((lon + 180.0) % 360.0) - 180.0


def geodesic_distance_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    g = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
    return float(g["s12"])


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def fmt(x: Any, ndigits: int = 6) -> str:
    try:
        return f"{float(x):.{ndigits}f}"
    except Exception:
        return str(x)


@dataclass
class StepCompare:
    step_idx: int
    timestamp: Any
    input_lat: float
    input_lon: float
    input_alt: float

    gfs_lat_new: float
    gfs_lon_new: float
    gfs_x: float
    gfs_y: float
    gfs_x_old: float
    gfs_y_old: float
    gfs_bearing: float
    gfs_nearest_lat: float
    gfs_nearest_lon: float
    gfs_nearest_alt: float

    era_lat_new: float
    era_lon_new: float
    era_x: float
    era_y: float
    era_x_old: float
    era_y_old: float
    era_bearing: float
    era_nearest_lat: float
    era_nearest_lon: float
    era_nearest_alt: float

    lat_new_abs_diff: float
    lon_new_abs_diff: float
    x_abs_diff: float
    y_abs_diff: float
    x_old_abs_diff: float
    y_old_abs_diff: float
    bearing_abs_diff: float
    nearest_lat_abs_diff: float
    nearest_lon_abs_diff: float
    nearest_alt_abs_diff: float
    step_endpoint_sep_m: float
    running_sep_m: float


def count_warning(key: str, value: float, threshold: float) -> bool:
    TOTAL_CHECKS[key] += 1
    warn = value > threshold
    if warn:
        WARNING_COUNTS[key] += 1
    return warn


def compare_getnewcoord_outputs(step_idx: int, timestamp: Any, coord_in: Dict[str, Any],
                                out_gfs: List[Any], out_era: List[Any],
                                running_sep_m: float) -> StepCompare:
    gfs_lat_new, gfs_lon_new, gfs_x, gfs_y, gfs_x_old, gfs_y_old, gfs_bearing, gfs_nearest_lat, gfs_nearest_lon, gfs_nearest_alt = out_gfs
    era_lat_new, era_lon_new, era_x, era_y, era_x_old, era_y_old, era_bearing, era_nearest_lat, era_nearest_lon, era_nearest_alt = out_era

    gfs_lon_new = normalize_lon_for_compare(float(gfs_lon_new))
    era_lon_new = normalize_lon_for_compare(float(era_lon_new))
    gfs_nearest_lon = normalize_lon_for_compare(float(gfs_nearest_lon))
    era_nearest_lon = normalize_lon_for_compare(float(era_nearest_lon))

    lat_new_abs_diff = abs(float(gfs_lat_new) - float(era_lat_new))
    lon_new_abs_diff = abs(normalize_lon_for_compare(float(gfs_lon_new) - float(era_lon_new)))
    x_abs_diff = abs(float(gfs_x) - float(era_x))
    y_abs_diff = abs(float(gfs_y) - float(era_y))
    x_old_abs_diff = abs(float(gfs_x_old) - float(era_x_old))
    y_old_abs_diff = abs(float(gfs_y_old) - float(era_y_old))
    bearing_abs_diff = abs(float(gfs_bearing) - float(era_bearing))
    nearest_lat_abs_diff = abs(float(gfs_nearest_lat) - float(era_nearest_lat))
    nearest_lon_abs_diff = abs(normalize_lon_for_compare(float(gfs_nearest_lon) - float(era_nearest_lon)))
    nearest_alt_abs_diff = abs(float(gfs_nearest_alt) - float(era_nearest_alt))

    step_endpoint_sep_m = geodesic_distance_m(
        float(gfs_lat_new), gfs_lon_new,
        float(era_lat_new), era_lon_new
    )

    return StepCompare(
        step_idx=step_idx,
        timestamp=to_py_datetime(timestamp),
        input_lat=float(coord_in["lat"]),
        input_lon=normalize_lon_for_compare(float(coord_in["lon"])),
        input_alt=float(coord_in["alt"]),

        gfs_lat_new=float(gfs_lat_new),
        gfs_lon_new=gfs_lon_new,
        gfs_x=float(gfs_x),
        gfs_y=float(gfs_y),
        gfs_x_old=float(gfs_x_old),
        gfs_y_old=float(gfs_y_old),
        gfs_bearing=float(gfs_bearing),
        gfs_nearest_lat=float(gfs_nearest_lat),
        gfs_nearest_lon=gfs_nearest_lon,
        gfs_nearest_alt=float(gfs_nearest_alt),

        era_lat_new=float(era_lat_new),
        era_lon_new=era_lon_new,
        era_x=float(era_x),
        era_y=float(era_y),
        era_x_old=float(era_x_old),
        era_y_old=float(era_y_old),
        era_bearing=float(era_bearing),
        era_nearest_lat=float(era_nearest_lat),
        era_nearest_lon=era_nearest_lon,
        era_nearest_alt=float(era_nearest_alt),

        lat_new_abs_diff=lat_new_abs_diff,
        lon_new_abs_diff=lon_new_abs_diff,
        x_abs_diff=x_abs_diff,
        y_abs_diff=y_abs_diff,
        x_old_abs_diff=x_old_abs_diff,
        y_old_abs_diff=y_old_abs_diff,
        bearing_abs_diff=bearing_abs_diff,
        nearest_lat_abs_diff=nearest_lat_abs_diff,
        nearest_lon_abs_diff=nearest_lon_abs_diff,
        nearest_alt_abs_diff=nearest_alt_abs_diff,
        step_endpoint_sep_m=step_endpoint_sep_m,
        running_sep_m=running_sep_m,
    )


def print_step_compare(sc: StepCompare) -> None:
    warn_flags = {
        "lat_new_abs_diff": count_warning("lat_new_abs_diff", sc.lat_new_abs_diff, COORD_WARN_ATOL_DEG),
        "lon_new_abs_diff": count_warning("lon_new_abs_diff", sc.lon_new_abs_diff, COORD_WARN_ATOL_DEG),
        "x_abs_diff": count_warning("x_abs_diff", sc.x_abs_diff, UV_WARN_ATOL),
        "y_abs_diff": count_warning("y_abs_diff", sc.y_abs_diff, UV_WARN_ATOL),
        "x_old_abs_diff": count_warning("x_old_abs_diff", sc.x_old_abs_diff, UV_WARN_ATOL),
        "y_old_abs_diff": count_warning("y_old_abs_diff", sc.y_old_abs_diff, UV_WARN_ATOL),
        "nearest_lat_abs_diff": count_warning("nearest_lat_abs_diff", sc.nearest_lat_abs_diff, COORD_WARN_ATOL_DEG),
        "nearest_lon_abs_diff": count_warning("nearest_lon_abs_diff", sc.nearest_lon_abs_diff, COORD_WARN_ATOL_DEG),
        "nearest_alt_abs_diff": count_warning("nearest_alt_abs_diff", sc.nearest_alt_abs_diff, 1e-6),
        "step_endpoint_sep_m": count_warning("step_endpoint_sep_m", sc.step_endpoint_sep_m, STEP_SEPARATION_WARN_M),
        "running_sep_m": count_warning("running_sep_m", sc.running_sep_m, RUNNING_SEPARATION_WARN_M),
    }

    any_warn = any(warn_flags.values())
    if PRINT_ONLY_WARNINGS and not any_warn:
        return

    print_section(f"STEP {sc.step_idx}  time={sc.timestamp}")

    print(f"Input state: lat={fmt(sc.input_lat, 8)} lon={fmt(sc.input_lon, 8)} alt={fmt(sc.input_alt, 3)}")

    print("\nGFS:")
    print(f"  lat/lon new      : {fmt(sc.gfs_lat_new, 8)}, {fmt(sc.gfs_lon_new, 8)}")
    print(f"  wind new (u,v)   : {fmt(sc.gfs_x)}, {fmt(sc.gfs_y)}")
    print(f"  wind old (u,v)   : {fmt(sc.gfs_x_old)}, {fmt(sc.gfs_y_old)}")
    print(f"  bearing          : {fmt(sc.gfs_bearing)}")
    print(f"  nearest lat/lon  : {fmt(sc.gfs_nearest_lat, 8)}, {fmt(sc.gfs_nearest_lon, 8)}")
    print(f"  nearest alt      : {fmt(sc.gfs_nearest_alt, 3)}")

    print("\nERA5:")
    print(f"  lat/lon new      : {fmt(sc.era_lat_new, 8)}, {fmt(sc.era_lon_new, 8)}")
    print(f"  wind new (u,v)   : {fmt(sc.era_x)}, {fmt(sc.era_y)}")
    print(f"  wind old (u,v)   : {fmt(sc.era_x_old)}, {fmt(sc.era_y_old)}")
    print(f"  bearing          : {fmt(sc.era_bearing)}")
    print(f"  nearest lat/lon  : {fmt(sc.era_nearest_lat, 8)}, {fmt(sc.era_nearest_lon, 8)}")
    print(f"  nearest alt      : {fmt(sc.era_nearest_alt, 3)}")

    def flag(name: str) -> str:
        return "  <-- WARN" if warn_flags.get(name, False) else ""

    print("\nDiffs:")
    print(f"  lat_new_abs_diff      : {fmt(sc.lat_new_abs_diff, 8)}{flag('lat_new_abs_diff')}")
    print(f"  lon_new_abs_diff      : {fmt(sc.lon_new_abs_diff, 8)}{flag('lon_new_abs_diff')}")
    print(f"  x_abs_diff            : {fmt(sc.x_abs_diff, 8)}{flag('x_abs_diff')}")
    print(f"  y_abs_diff            : {fmt(sc.y_abs_diff, 8)}{flag('y_abs_diff')}")
    print(f"  x_old_abs_diff        : {fmt(sc.x_old_abs_diff, 8)}{flag('x_old_abs_diff')}")
    print(f"  y_old_abs_diff        : {fmt(sc.y_old_abs_diff, 8)}{flag('y_old_abs_diff')}")
    print(f"  bearing_abs_diff      : {fmt(sc.bearing_abs_diff, 8)}")
    print(f"  nearest_lat_abs_diff  : {fmt(sc.nearest_lat_abs_diff, 8)}{flag('nearest_lat_abs_diff')}")
    print(f"  nearest_lon_abs_diff  : {fmt(sc.nearest_lon_abs_diff, 8)}{flag('nearest_lon_abs_diff')}")
    print(f"  nearest_alt_abs_diff  : {fmt(sc.nearest_alt_abs_diff, 8)}{flag('nearest_alt_abs_diff')}")
    print(f"  step_endpoint_sep_m   : {fmt(sc.step_endpoint_sep_m, 6)}{flag('step_endpoint_sep_m')}")
    print(f"  running_sep_m         : {fmt(sc.running_sep_m, 6)}{flag('running_sep_m')}")



# -----------------------------------------------------------------------------
# Main debug loop
# -----------------------------------------------------------------------------
def main() -> None:
    print_section("INITIALIZING SIDE-BY-SIDE DEBUG RUN")

    coord0 = dict(config_earth.simulation["start_coord"])

    gfs_model = GFS.GFS(coord0)
    era_model = ERA5.ERA5(coord0)

    # Shared vertical simulation setup, mirroring predict.py
    GMT = 7
    t = config_earth.simulation["start_time"]
    start = t
    min_alt = config_earth.simulation["min_alt"]
    alt_sp = config_earth.simulation["alt_sp"]
    v_sp = config_earth.simulation["v_sp"]
    dt = config_earth.simulation["dt"]
    atm = fluids.atmosphere.ATMOSPHERE_1976(min_alt)
    GFSrate = config_earth.forecast["GFSrate"]

    T_s = [atm.T]
    T_i = [atm.T]
    T_atm = [atm.T]
    el = [min_alt]
    v = [0.0]

    # Keep three coordinate tracks:
    # 1) shared input track (used as the common input to both models)
    # 2) GFS output track
    # 3) ERA output track
    coords_shared = [dict(coord0)]
    coords_gfs = [dict(coord0)]
    coords_era = [dict(coord0)]

    simulation_time = config_earth.simulation["sim_time"] * int(3600 * (1 / dt))
    if MAX_STEPS is not None:
        simulation_time = min(simulation_time, MAX_STEPS)

    e = solve_states.SolveStates()
    descent = False

    results: List[StepCompare] = []

    for i in range(simulation_time):
        # Shared vertical state update, same as predict.py
        T_s_new, T_i_new, T_atm_new, el_new, v_new, q_rad, q_surf, q_int = e.solveVerticalTrajectory(
            t, T_s[i], T_i[i], el[i], v[i], coords_shared[i], alt_sp, v_sp
        )

        if v_new < -3.0 and el_new > 15000:
            descent = True

        if descent:
            v_new = -3
            el_new = el[i] + v_new * dt
            if el_new < min_alt:
                el_new = min_alt
                v_new = 0

        T_s.append(T_s_new)
        T_i.append(T_i_new)
        T_atm.append(T_atm_new)
        el.append(el_new)
        v.append(v_new)

        # Match predict.py time update
        t = t + pd.Timedelta(hours=(1 / 3600 * dt))

        # Default: if no horizontal update this step, carry forward previous lon/lat
        lat_gfs_new = coords_gfs[i]["lat"]
        lon_gfs_new = coords_gfs[i]["lon"]
        lat_era_new = coords_era[i]["lat"]
        lon_era_new = coords_era[i]["lon"]

        # Keep old diagnostic fields in case this isn't a horizontal-update step
        out_gfs = None
        out_era = None

        if i % GFSrate == 0:
            shared_coord_in = {
                "lat": coords_shared[i]["lat"],
                "lon": coords_shared[i]["lon"],
                "alt": coords_shared[i]["alt"],
                "timestamp": coords_shared[i]["timestamp"],
            }

            out_gfs = gfs_model.getNewCoord(shared_coord_in, dt * GFSrate)
            out_era = era_model.getNewCoord(shared_coord_in, dt * GFSrate)

            lat_gfs_new = out_gfs[0]
            lon_gfs_new = out_gfs[1]
            lat_era_new = out_era[0]
            lon_era_new = out_era[1]

            running_sep_m = geodesic_distance_m(
                float(lat_gfs_new), normalize_lon_for_compare(float(lon_gfs_new)),
                float(lat_era_new), normalize_lon_for_compare(float(lon_era_new)),
            )

            sc = compare_getnewcoord_outputs(
                step_idx=i,
                timestamp=shared_coord_in["timestamp"],
                coord_in=shared_coord_in,
                out_gfs=out_gfs,
                out_era=out_era,
                running_sep_m=running_sep_m,
            )
            results.append(sc)

            if PRINT_EVERY_HORIZONTAL_UPDATE:
                print_step_compare(sc)

        # Build next coords
        coord_gfs_new = {
            "lat": lat_gfs_new,
            "lon": lon_gfs_new,
            "alt": el_new,
            "timestamp": t,
        }
        coord_era_new = {
            "lat": lat_era_new,
            "lon": lon_era_new,
            "alt": el_new,
            "timestamp": t,
        }

        # Shared coord follows whichever source you choose, so the loop can continue
        if CONTINUE_WITH.upper() == "ERA5":
            coord_shared_new = dict(coord_era_new)
        else:
            coord_shared_new = dict(coord_gfs_new)

        coords_gfs.append(coord_gfs_new)
        coords_era.append(coord_era_new)
        coords_shared.append(coord_shared_new)

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print_section("SUMMARY")
    print(f"Horizontal comparison steps: {len(results)}")
    print(f"Shared integration path followed: {CONTINUE_WITH.upper()}")

    if results:
        max_step_sep = max(r.step_endpoint_sep_m for r in results)
        max_running_sep = max(r.running_sep_m for r in results)
        max_x_diff = max(r.x_abs_diff for r in results)
        max_y_diff = max(r.y_abs_diff for r in results)
        max_lat_diff = max(r.lat_new_abs_diff for r in results)
        max_lon_diff = max(r.lon_new_abs_diff for r in results)

        print(f"max step_endpoint_sep_m : {fmt(max_step_sep, 6)}")
        print(f"max running_sep_m       : {fmt(max_running_sep, 6)}")
        print(f"max x_abs_diff          : {fmt(max_x_diff, 8)}")
        print(f"max y_abs_diff          : {fmt(max_y_diff, 8)}")
        print(f"max lat_new_abs_diff    : {fmt(max_lat_diff, 8)}")
        print(f"max lon_new_abs_diff    : {fmt(max_lon_diff, 8)}")

        worst = sorted(results, key=lambda r: r.running_sep_m, reverse=True)[:10]
        print("\nTop 10 worst running-separation steps:")
        for r in worst:
            print({
                "step_idx": r.step_idx,
                "timestamp": r.timestamp,
                "input_lat": r.input_lat,
                "input_lon": r.input_lon,
                "input_alt": r.input_alt,
                "step_endpoint_sep_m": r.step_endpoint_sep_m,
                "running_sep_m": r.running_sep_m,
                "x_abs_diff": r.x_abs_diff,
                "y_abs_diff": r.y_abs_diff,
            })

    print_section("WARNING SUMMARY")
    if TOTAL_CHECKS:
        for k in sorted(TOTAL_CHECKS.keys()):
            total = TOTAL_CHECKS[k]
            warns = WARNING_COUNTS.get(k, 0)
            pct = 100.0 * warns / total if total else 0.0
            print(f"{k:24s} : {warns:6d} / {total:6d}  ({pct:6.2f}%)")
    else:
        print("No checks recorded.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print_section("UNHANDLED ERROR")
        print(exc)
        traceback.print_exc()
        raise