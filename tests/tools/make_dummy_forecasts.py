"""Synthetic netCDF forecast generator for EarthSHAB tests.

Builds tiny in-memory forecasts with known analytical fields in the v2
canonical schema (see docs/forecast-schema-v2.md). Each scenario maps to
exactly one .nc file — Phase 5 collapsed the GFS and ERA5 readers into a
single Forecast class, so the GFS-vs-ERA5 cross-reader equality tests
(and their dual-schema writers) are no longer needed.

Used by:
    - pytest session fixture in tests/conftest.py (writes into tmp_path)
    - manual debugging via CLI:
        python -m tests.tools.make_dummy_forecasts --scenario all_static --out /tmp/foo.nc
        python -m tests.tools.make_dummy_forecasts --all --out tests/fixtures/

Canonical schema (matches saveNETCDF.py / migrate_v1.py output):
  dims:       valid_time, pressure_level, latitude, longitude
  data vars:  u, v, z (geopotential m**2/s**-2)
  lat:        DESCENDING
  pressure:   DESCENDING hPa  (→ altitude ASCENDING with index)
  lon:        ASCENDING (synthetic frame [0, 2], inside [-180, 180))
  time:       int64 'seconds since 1970-01-01' / proleptic_gregorian
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
# Grid constants
# ---------------------------------------------------------------------------

GRID_SHAPE = (4, 10, 9, 9)  # (time, level, lat, lon)
N_TIME, N_LEVEL, N_LAT, N_LON = GRID_SHAPE

# Lat/Lon: 9 points each spanning [0.0, 2.0] in 0.25-deg increments.
LAT_VALUES_ASC = np.linspace(0.0, 2.0, N_LAT).astype(np.float64)
LON_VALUES = np.linspace(0.0, 2.0, N_LON).astype(np.float64)

# Canonical altitude column (meters), ASCENDING with index. Data arrays are
# built in this index order; the writer stores them as-is because the v2
# canonical layout has pressure_level descending hPa = altitude ascending
# with index.
ALT_M_GFS_ORDER = np.array(
    [100.0, 800.0, 1500.0, 3200.0, 5800.0, 7400.0, 9700.0, 11000.0, 13800.0, 16500.0],
    dtype=np.float64,
)

# Pressure levels (hPa), DESCENDING — matches the canonical layout and
# aligns with ALT_M_GFS_ORDER ascending altitude.
LEVELS_HPA_DESC = np.array(
    [1000.0, 925.0, 850.0, 700.0, 500.0, 400.0, 300.0, 250.0, 200.0, 150.0],
    dtype=np.float64,
)
# Back-compat alias for callers that still import the old name.
LEVELS_HPA_GFS = LEVELS_HPA_DESC

# Time axis: 4 steps, 3 hours apart, starting at a synthetic epoch.
TIME_START = datetime(2025, 1, 1, 0, 0, 0)
TIME_STEP_HOURS = 3
TIMES_DT = [TIME_START + timedelta(hours=TIME_STEP_HOURS * i) for i in range(N_TIME)]

# Earth's standard gravity, used to convert h (m) → z (m²/s²).
G_STANDARD = 9.80665

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

    u, v, h are all shape (N_TIME, N_LEVEL, N_LAT, N_LON), with:
      - level index 0 = LOWEST altitude (= highest hPa)
      - lat index 0 = SOUTHERNMOST lat (LAT_VALUES_ASC[0])
      - lon index 0 = WESTERNMOST lon (LON_VALUES[0])
    h is altitude in meters. The writer converts to z = h * G_STANDARD and
    reverses the lat axis to produce canonical descending-lat layout.

    ``mask`` is None for the all-valid scenarios. For nan_at_top it is a
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
    """u = level_idx (ascending altitude), v = 0. Tests altitude interpolation."""
    u = _ones4d()
    for k in range(N_LEVEL):
        u[:, k, :, :] = float(k)
    v = 0.0 * _ones4d()
    h = _broadcast_alt_to_4d()
    return ScenarioArrays(u=u, v=v, h=h)


def _scenario_linear_ramp_all_dims() -> ScenarioArrays:
    """u = t_idx + lev_idx + lat_idx + lon_idx, v = 0. Composes every axis."""
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
    """u = lat_idx (ascending-lat indexing); constant in others."""
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
    """Same as static_by_level; kept distinct for clarity in the test catalog."""
    return _scenario_static_by_level()


def _scenario_bearing_wrap() -> ScenarioArrays:
    """u, v chosen so bearing crosses 0°/360° between adjacent levels.

    All cells at a level share the same (u, v); bearing rotates from level
    to level so that two specific adjacent levels straddle 0°. Tests the
    angle-wrap branch in interpolateBearing. Speed held constant at 10 m/s.
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
    """Static-by-level field with one corner column masked at the top altitudes.

    The mask is restricted to the (lat=0, lon=0) corner at the two highest
    canonical altitudes (indices N_LEVEL-3 and N_LEVEL-2). The lowest level
    (index 0) and the very topmost (index N_LEVEL-1) stay unmasked so the
    reader's interpolation has anchors at both endpoints. fill_missing_data
    must linearly fill the masked entries before wind_alt_Interpolate2.
    """
    u = _ones4d()
    for k in range(N_LEVEL):
        u[:, k, :, :] = float(k)
    v = 0.0 * _ones4d()
    h = _broadcast_alt_to_4d()
    mask = np.zeros(GRID_SHAPE, dtype=bool)
    mask[:, N_LEVEL - 3:N_LEVEL - 1, 0, 0] = True
    return ScenarioArrays(u=u, v=v, h=h, mask=mask)


def _scenario_geopotential_vs_height() -> ScenarioArrays:
    """Constant non-trivial u, v; altitude profile is the standard column.

    The writer stores z = h * g (m²/s²). The reader's __init__ divides z by
    g to recover altitude in meters. This scenario catches a silently-wrong
    g constant or a missing conversion.
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

ALL_STATIC_U = 5.0
ALL_STATIC_V = 0.0


# ---------------------------------------------------------------------------
# Canonical v2 writer
# ---------------------------------------------------------------------------

_FILL = np.float32(9.9692099683868690e36)  # CF-standard fill for float32.

TIME_UNITS = "seconds since 1970-01-01"
TIME_CALENDAR = "proleptic_gregorian"


def _write_canonical(arrays: ScenarioArrays, out_path: Path) -> None:
    """Write a v2 canonical netCDF.

    Layout per docs/forecast-schema-v2.md:
      - dims:       valid_time, pressure_level, latitude, longitude
      - data vars:  u, v, z(=h*g)
      - lat:        DESCENDING — writer reverses LAT_VALUES_ASC and the data
      - level:      DESCENDING hPa (LEVELS_HPA_DESC as-is)
      - lon:        ASCENDING (LON_VALUES as-is)
      - time:       int64 seconds-since-epoch, proleptic_gregorian
      - global:     Conventions = 'CF-1.7'
                    institution = '' (the Forecast reader's source-detection
                    falls back to the filename, then to 'unknown' — fine for
                    synthetic fixtures)
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    use_fill = arrays.mask is not None

    with netCDF4.Dataset(out_path, "w", format="NETCDF4") as ds:
        ds.createDimension("valid_time", N_TIME)
        ds.createDimension("pressure_level", N_LEVEL)
        ds.createDimension("latitude", N_LAT)
        ds.createDimension("longitude", N_LON)

        v_time = ds.createVariable("valid_time", "i8", ("valid_time",))
        v_level = ds.createVariable("pressure_level", "f8", ("pressure_level",))
        v_lat = ds.createVariable("latitude", "f8", ("latitude",))
        v_lon = ds.createVariable("longitude", "f8", ("longitude",))

        epoch_secs = netCDF4.date2num(
            TIMES_DT, units=TIME_UNITS, calendar=TIME_CALENDAR,
        )
        v_time[:] = np.asarray(epoch_secs, dtype=np.int64)
        v_time.units = TIME_UNITS
        v_time.calendar = TIME_CALENDAR
        v_time.standard_name = "time"
        v_time.long_name = "time"

        v_level[:] = LEVELS_HPA_DESC
        v_level.units = "hPa"
        v_level.long_name = "pressure"
        v_level.standard_name = "air_pressure"
        v_level.positive = "down"
        v_level.stored_direction = "decreasing"

        v_lat[:] = LAT_VALUES_ASC[::-1]
        v_lat.units = "degrees_north"
        v_lat.long_name = "latitude"
        v_lat.standard_name = "latitude"
        v_lat.stored_direction = "decreasing"

        v_lon[:] = LON_VALUES
        v_lon.units = "degrees_east"
        v_lon.long_name = "longitude"
        v_lon.standard_name = "longitude"

        def _mkvar(name: str) -> netCDF4.Variable:
            if use_fill:
                return ds.createVariable(
                    name, "f4",
                    ("valid_time", "pressure_level", "latitude", "longitude"),
                    fill_value=_FILL,
                )
            return ds.createVariable(
                name, "f4",
                ("valid_time", "pressure_level", "latitude", "longitude"),
            )

        v_u = _mkvar("u")
        v_v = _mkvar("v")
        v_z = _mkvar("z")

        # Reverse the lat axis so the on-disk layout is descending lat.
        # Level/lon order in the source arrays already matches the canonical
        # layout (level: 0 = lowest alt = highest hPa → matches descending hPa
        # with index 0 = highest hPa).
        u_arr = arrays.u[:, :, ::-1, :].astype(np.float32)
        v_arr = arrays.v[:, :, ::-1, :].astype(np.float32)
        z_arr = (arrays.h[:, :, ::-1, :] * G_STANDARD).astype(np.float32)

        if use_fill:
            m = arrays.mask[:, :, ::-1, :]
            u_arr = np.ma.array(u_arr, mask=m)
            v_arr = np.ma.array(v_arr, mask=m)
            z_arr = np.ma.array(z_arr, mask=m)

        v_u[:] = u_arr
        v_v[:] = v_arr
        v_z[:] = z_arr

        v_u.units = "m s**-1"
        v_u.long_name = "U component of wind"
        v_u.standard_name = "eastward_wind"
        v_v.units = "m s**-1"
        v_v.long_name = "V component of wind"
        v_v.standard_name = "northward_wind"
        v_z.units = "m**2 s**-2"
        v_z.long_name = "Geopotential"
        v_z.standard_name = "geopotential"

        ds.Conventions = "CF-1.7"
        ds.institution = ""  # synthetic; reader falls back to filename then 'unknown'


def build_dummy(scenario: str, out_path: str | Path) -> Path:
    """Generate one synthetic .nc file. Returns the output path."""
    if scenario not in SCENARIOS:
        raise ValueError(
            f"Unknown scenario {scenario!r}. Known: {sorted(SCENARIOS)}"
        )
    arrays = SCENARIOS[scenario]()
    out = Path(out_path)
    _write_canonical(arrays, out)
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cli() -> None:
    p = argparse.ArgumentParser(
        description="Generate synthetic v2-canonical netCDF forecasts for "
                    "EarthSHAB tests."
    )
    p.add_argument("--scenario", choices=sorted(SCENARIOS),
                   help="Single scenario to write. Required unless --all.")
    p.add_argument("--out", required=True,
                   help="Output file (single scenario) or directory (--all).")
    p.add_argument("--all", action="store_true",
                   help="Write every scenario into --out as a directory.")
    args = p.parse_args()

    out = Path(args.out)
    if args.all:
        out.mkdir(parents=True, exist_ok=True)
        for scenario in SCENARIOS:
            path = build_dummy(scenario, out / f"{scenario}.nc")
            print(f"wrote {path}")
        return

    if not args.scenario:
        p.error("--scenario is required unless --all is given.")
    path = build_dummy(args.scenario, out)
    print(f"wrote {path}")


if __name__ == "__main__":
    _cli()
