"""Synthetic netCDF forecast generator for EarthSHAB tests.

Builds tiny in-memory forecasts with known analytical fields and writes them
in either GFS-style or ERA5-style schema. The same physical scenario can be
written in both schemas, so tests can verify that both readers produce the
same answers (which is the consolidation safety net for the planned unified
forecast format).

Used by:
    - pytest session fixture in tests/conftest.py (writes into tmp_path)
    - manual debugging via CLI:
        python -m tests.tools.make_dummy_forecasts --scenario all_static --schema gfs --out /tmp/foo.nc
        python -m tests.tools.make_dummy_forecasts --all --out tests/fixtures/

Two intentional implementation details worth knowing:

(1) Time encoding mismatch between readers.
    GFS.py:70 hardcodes ``netCDF4.num2date(time_arr[:], units="days since 0001-01-01")``
    regardless of what units the file actually stores. So the GFS writer here
    encodes time as numeric days-since-year-0001. The ERA5 writer uses
    CF-standard "hours since 1970-01-01" with calendar=standard, which ERA5.py
    reads from the file attributes.

(2) Masked-array convention.
    netCDF4.Dataset returns numeric variables as MaskedArray with auto_mask=True
    by default, even when no _FillValue is set. The readers call .filled(np.nan)
    on these slices, which works on empty-mask arrays. So nine of the ten
    scenarios write fully-valid float32 data with no _FillValue. Only
    ``nan_at_top`` sets _FillValue and writes masked entries at the upper
    pressure levels, since that scenario exists specifically to exercise
    fill_missing_data against a real mask. As a side effect, the nine
    no-mask scenarios route both readers through the "results == True" branch
    of determineRanges (the under-tested fallback path in ERA5.py:175). The
    branch produces a working instance, so this is a feature: synthetic tests
    cover code that real subset forecasts don't.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4
import numpy as np


# ---------------------------------------------------------------------------
# Grid constants — chosen in Q10(b) of the design grilling.
# ---------------------------------------------------------------------------

GRID_SHAPE = (4, 10, 9, 9)  # (time, level, lat, lon)
N_TIME, N_LEVEL, N_LAT, N_LON = GRID_SHAPE

# Lat/Lon: 9 points each spanning [0.0, 2.0] in 0.25-deg increments.
# Synthetic frame independent of any real flight config.
LAT_VALUES_ASC = np.linspace(0.0, 2.0, N_LAT).astype(np.float64)
LON_VALUES = np.linspace(0.0, 2.0, N_LON).astype(np.float64)

# Canonical altitude column (meters), ASCENDING with index. This is the
# in-memory order used by every scenario builder. The GFS writer stores
# data in this order directly (matching real GFS, where lev is descending
# hPa). The ERA5 writer reverses both the level axis and the data along
# that axis (matching real ERA5, where level is ascending hPa).
ALT_M_GFS_ORDER = np.array(
    [100.0, 800.0, 1500.0, 3200.0, 5800.0, 7400.0, 9700.0, 11000.0, 13800.0, 16500.0],
    dtype=np.float64,
)

# Pressure levels (hPa) aligned with ALT_M_GFS_ORDER ascending altitude →
# GFS lev variable is descending hPa. Values are illustrative; the readers
# don't actually use lev/level for interpolation (they use hgtprs/z).
LEVELS_HPA_GFS = np.array(
    [1000.0, 925.0, 850.0, 700.0, 500.0, 400.0, 300.0, 250.0, 200.0, 150.0],
    dtype=np.float64,
)

# Time axis: 4 steps, 3 hours apart, starting at a synthetic epoch.
TIME_START = datetime(2025, 1, 1, 0, 0, 0)
TIME_STEP_HOURS = 3
TIMES_DT = [TIME_START + timedelta(hours=TIME_STEP_HOURS * i) for i in range(N_TIME)]

# Earth's standard gravity, used to convert m → geopotential (z = h * g).
G_STANDARD = 9.80665

# Synthetic start coord placed in the interior of the grid (not on an edge,
# not on a grid intersection — picks up off-grid behavior). Tests can choose
# to use this or override with their own coords.
SYNTHETIC_START_COORD = {
    "lat": 1.0,
    "lon": 1.0,
    "alt": 5000.0,
    "timestamp": TIMES_DT[1],
}


# ---------------------------------------------------------------------------
# Scenario arrays
# ---------------------------------------------------------------------------

@dataclass
class ScenarioArrays:
    """Bundle of arrays for one synthetic scenario.

    u, v, h are all shape (N_TIME, N_LEVEL, N_LAT, N_LON), GFS index order
    (lev index 0 = highest altitude, lat ascending, lon 0..360 style).
    h is altitude in meters. The ERA5 writer will reverse lev/lat and convert
    h → z = h * G_STANDARD.

    ``mask`` is None for the 9 all-valid scenarios. For nan_at_top it is a
    boolean array of the same shape with True where the cell should be masked.
    """

    u: np.ndarray
    v: np.ndarray
    h: np.ndarray
    mask: np.ndarray | None = None


def _ones4d() -> np.ndarray:
    return np.ones(GRID_SHAPE, dtype=np.float64)


def _broadcast_alt_to_4d() -> np.ndarray:
    """Broadcast ALT_M_GFS_ORDER (shape (N_LEVEL,)) to (N_TIME, N_LEVEL, N_LAT, N_LON)."""
    return np.broadcast_to(
        ALT_M_GFS_ORDER[None, :, None, None], GRID_SHAPE
    ).astype(np.float64).copy()


def _scenario_all_static() -> ScenarioArrays:
    """u = 5, v = 0 everywhere; altitude is the standard column at every cell."""
    u = 5.0 * _ones4d()
    v = 0.0 * _ones4d()
    h = _broadcast_alt_to_4d()
    return ScenarioArrays(u=u, v=v, h=h)


def _scenario_static_by_level() -> ScenarioArrays:
    """u = level_idx (GFS order: 0=highest altitude), v = 0.

    All lat/lon/time cells at a given level share the same u. Tests altitude
    interpolation in isolation.
    """
    u = _ones4d()
    for k in range(N_LEVEL):
        u[:, k, :, :] = float(k)
    v = 0.0 * _ones4d()
    h = _broadcast_alt_to_4d()
    return ScenarioArrays(u=u, v=v, h=h)


def _scenario_linear_ramp_all_dims() -> ScenarioArrays:
    """u = t_idx + lev_idx + lat_idx + lon_idx, v = 0.

    Composes every lookup + interp. Analytical truth is exactly the sum of
    fractional indices for any (t, alt, lat, lon) query.
    """
    u = _ones4d()
    for ti in range(N_TIME):
        for k in range(N_LEVEL):
            for la in range(N_LAT):
                for lo in range(N_LON):
                    u[ti, k, la, lo] = float(ti + k + la + lo)
    v = 0.0 * _ones4d()
    h = _broadcast_alt_to_4d()
    return ScenarioArrays(u=u, v=v, h=h)


def _scenario_time_ramp_only() -> ScenarioArrays:
    """u = t_idx; constant in space. Isolates time interpolation."""
    u = _ones4d()
    for ti in range(N_TIME):
        u[ti, :, :, :] = float(ti)
    v = 0.0 * _ones4d()
    h = _broadcast_alt_to_4d()
    return ScenarioArrays(u=u, v=v, h=h)


def _scenario_lat_ramp_only() -> ScenarioArrays:
    """u = lat_idx (in GFS ascending-lat indexing); constant in others.

    Isolates lat lookup. ERA5 lat is descending — the writer reverses, so
    the same physical lat value yields the same u in both schemas.
    """
    u = _ones4d()
    for la in range(N_LAT):
        u[:, :, la, :] = float(la)
    v = 0.0 * _ones4d()
    h = _broadcast_alt_to_4d()
    return ScenarioArrays(u=u, v=v, h=h)


def _scenario_lon_ramp_only() -> ScenarioArrays:
    """u = lon_idx; constant in others. Isolates lon lookup."""
    u = _ones4d()
    for lo in range(N_LON):
        u[:, :, :, lo] = float(lo)
    v = 0.0 * _ones4d()
    h = _broadcast_alt_to_4d()
    return ScenarioArrays(u=u, v=v, h=h)


def _scenario_altitude_ramp_only() -> ScenarioArrays:
    """u = level_idx, v = 0. Same as static_by_level; kept distinct here for
    clarity in the test catalog and to allow future divergence (e.g., make
    altitude vary in time)."""
    return _scenario_static_by_level()


def _scenario_bearing_wrap() -> ScenarioArrays:
    """u, v chosen so bearing crosses 0°/360° between adjacent levels.

    All cells at a level share the same (u, v); bearing rotates from level
    to level so that two specific adjacent levels straddle 0°. Tests the
    angle-wrap branch in interpolateBearing.

    Concretely: bearing at level k is (-90 + 20*k) degrees, mod 360.
    For some k, bearing crosses 360° → 0°.
    Speed is held constant at 10 m/s so the wrap test isn't confounded by
    speed differences.
    """
    speed = 10.0
    u = _ones4d()
    v = _ones4d()
    for k in range(N_LEVEL):
        bearing_deg = (-90.0 + 20.0 * k) % 360.0
        u[:, k, :, :] = speed * math.cos(math.radians(bearing_deg))
        v[:, k, :, :] = speed * math.sin(math.radians(bearing_deg))
    h = _broadcast_alt_to_4d()
    return ScenarioArrays(u=u, v=v, h=h)


def _scenario_nan_at_top() -> ScenarioArrays:
    """Static-by-level field with one corner column masked at sub-top altitudes.

    Mask design constraint: the readers' determineRanges sampling slices are
    GFS [time, level=0, lat, lon] (the LOWEST altitude in the new ascending
    convention) and ERA5 [time=0, level=0, lat, lon] (the HIGHEST altitude
    after the ERA5 writer's level reversal). If EITHER slice is fully masked,
    determineRanges crashes (GFS .min() on empty array; ERA5 unpacks 3
    values from a 2D nonzero — a real bug in ERA5.py:167 only exercised by
    partial-mask files).

    To avoid both: mask only the column at (lat=0, lon=0) at level indices
    that don't coincide with EITHER reader's level-0 sampling slice. In the
    canonical (ascending) layout the GFS sample is at index 0 (low alt) and
    the ERA5 sample is at index N_LEVEL-1 after the writer's reversal (also
    low alt, since the writer reverses both level axis and data). So both
    samples land on canonical index 0 → safe to mask any indices > 0.
    Mask the top 2 indices in the canonical ascending column (= true
    high-altitude levels). fill_missing_data has work at that column.
    """
    u = _ones4d()
    for k in range(N_LEVEL):
        u[:, k, :, :] = float(k)
    v = 0.0 * _ones4d()
    h = _broadcast_alt_to_4d()
    mask = np.zeros(GRID_SHAPE, dtype=bool)
    # Mask near the top altitudes at the (lat=0, lon=0) corner. Leave the
    # very highest canonical index (N_LEVEL-1) UNMASKED because that maps
    # to ERA5 disk-level 0, the ERA5 reader's determineRanges sample slice.
    # A partial mask there triggers a real bug in ERA5.py:167 (unpacking
    # 3 values from a 2D nonzero). Mask only canonical indices
    # [N_LEVEL-3, N_LEVEL-2].
    mask[:, N_LEVEL - 3:N_LEVEL - 1, 0, 0] = True
    return ScenarioArrays(u=u, v=v, h=h, mask=mask)


def _scenario_geopotential_vs_height() -> ScenarioArrays:
    """Constant non-trivial u, v; altitude profile is the standard column.

    The GFS writer stores h directly (meters). The ERA5 writer stores
    z = h * g. Both readers should resolve the same physical altitude
    after ERA5's /g conversion. This scenario specifically catches a
    silently-wrong g constant or a missing conversion.
    """
    u = 7.5 * _ones4d()
    v = -3.25 * _ones4d()
    h = _broadcast_alt_to_4d()
    return ScenarioArrays(u=u, v=v, h=h)


SCENARIOS = {
    "all_static": _scenario_all_static,
    "static_by_level": _scenario_static_by_level,
    "linear_ramp_all_dims": _scenario_linear_ramp_all_dims,
    "time_ramp_only": _scenario_time_ramp_only,
    "lat_ramp_only": _scenario_lat_ramp_only,
    "lon_ramp_only": _scenario_lon_ramp_only,
    "altitude_ramp_only": _scenario_altitude_ramp_only,
    "bearing_wrap": _scenario_bearing_wrap,
    "nan_at_top": _scenario_nan_at_top,
    "geopotential_vs_height": _scenario_geopotential_vs_height,
}

SCHEMAS = ("gfs", "era5")

# Speed and v-component used by all_static and (partially) by geopotential_vs_height.
ALL_STATIC_U = 5.0
ALL_STATIC_V = 0.0


# ---------------------------------------------------------------------------
# Schema writers
# ---------------------------------------------------------------------------

_FILL = np.float32(9.9692099683868690e36)  # CF-standard fill for float32.

GFS_TIME_UNITS = "days since 0001-01-01"
ERA5_TIME_UNITS = "hours since 1970-01-01 00:00:00"
ERA5_TIME_CALENDAR = "standard"


def _write_gfs_style(arrays: ScenarioArrays, out_path: Path) -> None:
    """Write a GFS-style netCDF: vars ugrdprs, vgrdprs, hgtprs; dims time, lev, lat, lon.

    Time encoded as days-since-year-0001 because GFS.py:70 hardcodes that
    decoder regardless of file metadata.
    Lat ascending. Lon in [0, 360) style (synthetic frame is [0, 2] so it's
    trivially within the GFS convention with no wrap needed).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    use_fill = arrays.mask is not None

    with netCDF4.Dataset(out_path, "w", format="NETCDF4") as ds:
        ds.createDimension("time", N_TIME)
        ds.createDimension("lev", N_LEVEL)
        ds.createDimension("lat", N_LAT)
        ds.createDimension("lon", N_LON)

        v_time = ds.createVariable("time", "f8", ("time",))
        v_lev = ds.createVariable("lev", "f8", ("lev",))
        v_lat = ds.createVariable("lat", "f8", ("lat",))
        v_lon = ds.createVariable("lon", "f8", ("lon",))

        # Time: days since year 0001, has_year_zero=True to match GFS.py:70.
        #
        # date2num + num2date with these units don't round-trip cleanly —
        # netCDF4 shifts by 2 days for our 2025 dates due to the proleptic
        # vs historical Gregorian calendar offset. We measure the offset
        # empirically and pre-shift the stored numerics so the reader's
        # decode lands exactly on TIMES_DT. After this shift, time_convert
        # entries are cftime.DatetimeGregorian objects that compare ==
        # to stdlib datetime instances at the same wall-clock time.
        time_nums = netCDF4.date2num(
            TIMES_DT, units=GFS_TIME_UNITS, has_year_zero=True
        )
        decoded_unshifted = netCDF4.num2date(
            np.asarray(time_nums), units=GFS_TIME_UNITS, has_year_zero=True
        )
        offset_days = (TIMES_DT[0] - decoded_unshifted[0]).days
        v_time[:] = np.asarray(time_nums, dtype=np.float64) + offset_days
        # We intentionally do NOT set units on time. GFS.py ignores them.

        v_lev[:] = LEVELS_HPA_GFS  # ascending hPa = descending altitude
        v_lat[:] = LAT_VALUES_ASC  # ascending lat (GFS convention)
        v_lon[:] = LON_VALUES      # 0..360 style; synthetic frame [0, 2]

        def _mkvar(name: str) -> netCDF4.Variable:
            if use_fill:
                return ds.createVariable(
                    name, "f4", ("time", "lev", "lat", "lon"),
                    fill_value=_FILL,
                )
            return ds.createVariable(
                name, "f4", ("time", "lev", "lat", "lon")
            )

        v_u = _mkvar("ugrdprs")
        v_v = _mkvar("vgrdprs")
        v_h = _mkvar("hgtprs")

        u_arr = arrays.u.astype(np.float32)
        v_arr = arrays.v.astype(np.float32)
        h_arr = arrays.h.astype(np.float32)

        if use_fill:
            u_arr = np.ma.array(u_arr, mask=arrays.mask)
            v_arr = np.ma.array(v_arr, mask=arrays.mask)
            h_arr = np.ma.array(h_arr, mask=arrays.mask)

        v_u[:] = u_arr
        v_v[:] = v_arr
        v_h[:] = h_arr


def _write_era5_style(arrays: ScenarioArrays, out_path: Path) -> None:
    """Write an ERA5-style netCDF: vars u, v, z; dims time, level, latitude, longitude.

    Conventions vs the canonical in-memory layout (ascending altitude,
    ascending lat, GFS-style 0-360 lon-positive):
      - level stored ASCENDING in hPa = DESCENDING in altitude
        (matches real ERA5; opposite of GFS). Writer reverses both the
        level axis values AND the data along the level axis. ERA5.py
        internally re-reverses with [::-1] before interpolation.
      - lat stored DESCENDING (latitude[0] = max). Writer reverses
        lat axis values AND the data along the lat axis.
      - z stored as geopotential = h * g (m^2/s^2). ERA5.__init__ divides
        by g to recover altitude in meters.
      - Lon kept as [0, 2] (synthetic frame). Real ERA5 stores [-180, 180);
        our values are positive so they're valid in both conventions.
      - time uses CF-standard "hours since 1970-01-01" with calendar attr,
        which ERA5.py:108 reads from the file.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    use_fill = arrays.mask is not None

    with netCDF4.Dataset(out_path, "w", format="NETCDF4") as ds:
        ds.createDimension("time", N_TIME)
        ds.createDimension("level", N_LEVEL)
        ds.createDimension("latitude", N_LAT)
        ds.createDimension("longitude", N_LON)

        v_time = ds.createVariable("time", "f8", ("time",))
        v_level = ds.createVariable("level", "f8", ("level",))
        v_lat = ds.createVariable("latitude", "f8", ("latitude",))
        v_lon = ds.createVariable("longitude", "f8", ("longitude",))

        time_nums = netCDF4.date2num(
            TIMES_DT, units=ERA5_TIME_UNITS, calendar=ERA5_TIME_CALENDAR
        )
        v_time[:] = np.asarray(time_nums, dtype=np.float64)
        v_time.units = ERA5_TIME_UNITS
        v_time.calendar = ERA5_TIME_CALENDAR

        v_level[:] = LEVELS_HPA_GFS[::-1]  # ERA5 stores level ascending hPa
        v_lat[:] = LAT_VALUES_ASC[::-1]    # ERA5 stores lat descending
        v_lon[:] = LON_VALUES

        v_level.units = "hPa"
        v_lat.units = "degrees_north"
        v_lon.units = "degrees_east"

        def _mkvar(name: str) -> netCDF4.Variable:
            if use_fill:
                return ds.createVariable(
                    name, "f4", ("time", "level", "latitude", "longitude"),
                    fill_value=_FILL,
                )
            return ds.createVariable(
                name, "f4", ("time", "level", "latitude", "longitude")
            )

        v_u = _mkvar("u")
        v_v = _mkvar("v")
        v_z = _mkvar("z")

        # Reverse both the level axis (so on-disk altitude descends with idx)
        # and the lat axis (so on-disk lat descends with idx).
        u_arr = arrays.u[:, ::-1, ::-1, :].astype(np.float32)
        v_arr = arrays.v[:, ::-1, ::-1, :].astype(np.float32)
        # ERA5 z is geopotential = h * g.
        z_arr = (arrays.h[:, ::-1, ::-1, :] * G_STANDARD).astype(np.float32)

        if use_fill:
            m = arrays.mask[:, ::-1, ::-1, :]
            u_arr = np.ma.array(u_arr, mask=m)
            v_arr = np.ma.array(v_arr, mask=m)
            z_arr = np.ma.array(z_arr, mask=m)

        v_u[:] = u_arr
        v_v[:] = v_arr
        v_z[:] = z_arr


def build_dummy(scenario: str, schema: str, out_path: str | Path) -> Path:
    """Generate one synthetic .nc file. Returns the output path."""
    if scenario not in SCENARIOS:
        raise ValueError(
            f"Unknown scenario {scenario!r}. Known: {sorted(SCENARIOS)}"
        )
    if schema not in SCHEMAS:
        raise ValueError(
            f"Unknown schema {schema!r}. Expected one of {SCHEMAS}."
        )
    arrays = SCENARIOS[scenario]()
    out = Path(out_path)
    if schema == "gfs":
        _write_gfs_style(arrays, out)
    else:
        _write_era5_style(arrays, out)
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cli() -> None:
    p = argparse.ArgumentParser(
        description="Generate synthetic GFS-style or ERA5-style netCDF "
                    "forecasts for EarthSHAB tests."
    )
    p.add_argument("--scenario", choices=sorted(SCENARIOS),
                   help="Single scenario to write. Required unless --all.")
    p.add_argument("--schema", choices=SCHEMAS,
                   help="Schema to write. Required unless --all.")
    p.add_argument("--out", required=True,
                   help="Output file (single scenario) or directory (--all).")
    p.add_argument("--all", action="store_true",
                   help="Write every (scenario, schema) into --out as a directory.")
    args = p.parse_args()

    out = Path(args.out)
    if args.all:
        out.mkdir(parents=True, exist_ok=True)
        for scenario in SCENARIOS:
            for schema in SCHEMAS:
                fname = f"{scenario}_{schema}.nc"
                path = build_dummy(scenario, schema, out / fname)
                print(f"wrote {path}")
        return

    if not args.scenario or not args.schema:
        p.error("--scenario and --schema are required unless --all is given.")
    path = build_dummy(args.scenario, args.schema, out)
    print(f"wrote {path}")


if __name__ == "__main__":
    _cli()
