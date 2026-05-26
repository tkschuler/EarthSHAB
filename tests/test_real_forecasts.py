"""Smoke + invariant tests against the real forecasts pinned in config_earth.

Per Q6(b), the real-forecast tier runs against ONLY the files that
``config_earth`` currently points at — not every .nc in the forecasts
directory. The two files are:

    GFS:  src/EarthSHAB/forecasts/gfs_0p25_<YYYYMMDD>_<HH>.nc
          (derived from config_earth.forecast['forecast_start_time'])
    ERA5: src/EarthSHAB/forecasts/<config_earth.netcdf_era5['filename']>

These tests are intentionally LOOSE — they do not assert analytical truth
(impossible on real data) and do not assert cross-reader equality (GFS
and ERA5 are different meteorological models). They assert:

  1. The reader loads without raising.
  2. The schema has the expected variables and dimensions.
  3. The configured start_coord lat/lon is within the file's lat/lon bounds.
  4. The configured start_time is within the file's time range.
  5. A sample query at the configured start_coord returns finite values
     within physical bounds.
  6. All three wind_interpolation methods produce finite values within
     wide tolerance of each other at the same query.

These tests will need their physical bounds updated if the config rotates
to a profoundly different region/altitude — but the test structure itself
is config-rotation-stable.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import EarthSHAB.config_earth as config_earth
from EarthSHAB.ERA5 import ERA5
from EarthSHAB.GFS import GFS

from tests.tools.reader_adapter import (
    get_new_coord,
    is_gfs,
    nearest_lat,
    nearest_lon,
    wind_interp,
)


# Loose physical bounds — chosen to catch order-of-magnitude bugs, not
# subtle errors. Real wind speeds rarely exceed 100 m/s in the
# stratosphere; balloons fly up to ~30000 m.
_MAX_WIND_SPEED_MPS = 200.0
_MAX_ALT_M = 60000.0


# ---------------------------------------------------------------------------
# 1. Loads without error
# ---------------------------------------------------------------------------

class TestLoadsWithoutError:
    """Smoke: instantiating each real reader doesn't raise."""

    def test_gfs_reader_loads(self):
        # Sanity-check the file exists before opening — gives a clearer
        # failure than netCDF4's I/O error.
        nc_path = Path(config_earth.netcdf_gfs["nc_file"])
        assert nc_path.exists(), f"missing real GFS file: {nc_path}"
        reader = GFS(config_earth.simulation["start_coord"])
        assert reader.file is not None

    def test_era5_reader_loads(self):
        nc_path = (
            Path("src/EarthSHAB/forecasts") / config_earth.netcdf_era5["filename"]
        )
        assert nc_path.exists(), f"missing real ERA5 file: {nc_path}"
        reader = ERA5(config_earth.simulation["start_coord"])
        assert reader.file is not None


# ---------------------------------------------------------------------------
# 2. Schema
# ---------------------------------------------------------------------------

class TestSchema:
    """The two files conform to the expected per-schema variable + dim layout."""

    def test_gfs_schema(self, real_reader):
        if not is_gfs(real_reader):
            pytest.skip("ERA5 reader has its own schema test")
        ds = real_reader.file
        for var in ("time", "lev", "lat", "lon", "ugrdprs", "vgrdprs", "hgtprs"):
            assert var in ds.variables, f"GFS missing variable {var!r}"
        for dim in ("time", "lev", "lat", "lon"):
            assert dim in ds.dimensions, f"GFS missing dimension {dim!r}"

    def test_era5_schema(self, real_reader):
        if is_gfs(real_reader):
            pytest.skip("GFS reader has its own schema test")
        ds = real_reader.file
        for var in ("time", "level", "latitude", "longitude", "u", "v", "z"):
            assert var in ds.variables, f"ERA5 missing variable {var!r}"
        for dim in ("time", "level", "latitude", "longitude"):
            assert dim in ds.dimensions, f"ERA5 missing dimension {dim!r}"


# ---------------------------------------------------------------------------
# 3. Lat/lon bounds contain the configured start coord
# ---------------------------------------------------------------------------

class TestLatLonBounds:
    """The start_coord configured in config_earth.simulation must fall inside
    the file's lat/lon range — otherwise the simulator would be querying
    outside the downloaded subset on the first step."""

    def test_start_coord_lat_in_range(self, real_reader):
        lat = float(config_earth.simulation["start_coord"]["lat"])
        lat_arr = np.asarray(real_reader.lat, dtype=float)
        assert lat_arr.min() <= lat <= lat_arr.max(), (
            f"start_coord lat {lat} outside file range "
            f"[{lat_arr.min()}, {lat_arr.max()}]"
        )

    def test_start_coord_lon_in_range(self, real_reader):
        lon = float(config_earth.simulation["start_coord"]["lon"])
        lon_arr = np.asarray(real_reader.lon, dtype=float)
        if is_gfs(real_reader):
            # GFS stores lon in [0, 360). Normalize for comparison.
            lon = lon % 360.0
        assert lon_arr.min() <= lon <= lon_arr.max(), (
            f"start_coord lon (normalized={lon}) outside file range "
            f"[{lon_arr.min()}, {lon_arr.max()}]"
        )


# ---------------------------------------------------------------------------
# 4. Start time is within the file's time range
# ---------------------------------------------------------------------------

class TestStartTimeInRange:
    def test_start_time_in_range(self, real_reader):
        start = config_earth.simulation["start_time"]
        times = real_reader.time_convert
        assert times[0] <= start <= times[-1], (
            f"start_time {start} outside file time range "
            f"[{times[0]}, {times[-1]}]"
        )


# ---------------------------------------------------------------------------
# 5. Sample query at start coord returns finite, in-bounds values
# ---------------------------------------------------------------------------

class TestSampleQueryFinite:
    def test_get_new_coord_at_start_finite(self, real_reader):
        coord = dict(config_earth.simulation["start_coord"])
        out = get_new_coord(real_reader, coord, dt=1.0)
        assert len(out) == 10
        lat_new, lon_new, u, v, _, _, bearing, c_lat, c_lon, c_alt = out

        for name, val in [
            ("lat_new", lat_new), ("lon_new", lon_new),
            ("u", u), ("v", v),
            ("c_lat", c_lat), ("c_lon", c_lon), ("c_alt", c_alt),
        ]:
            assert math.isfinite(float(val)), f"{name} not finite: {val}"

        assert abs(float(u)) <= _MAX_WIND_SPEED_MPS, f"u out of bounds: {u}"
        assert abs(float(v)) <= _MAX_WIND_SPEED_MPS, f"v out of bounds: {v}"
        assert 0.0 <= float(c_alt) <= _MAX_ALT_M, f"c_alt out of bounds: {c_alt}"

    def test_wind_interp_at_start_finite(self, real_reader):
        coord = dict(config_earth.simulation["start_coord"])
        u, v, u_diag, v_diag = wind_interp(real_reader, coord)
        for name, val in [("u", u), ("v", v), ("u_diag", u_diag), ("v_diag", v_diag)]:
            assert math.isfinite(float(val)), f"{name} not finite: {val}"


# ---------------------------------------------------------------------------
# 6. All three interpolation methods agree within loose tolerance
# ---------------------------------------------------------------------------

# Loose tolerance: real altitude profiles are non-linear, so different
# interpolation methods CAN diverge — but only by O(1 m/s), not by orders
# of magnitude. Tight enough to catch a "wrong method returns garbage" bug,
# loose enough to permit genuine method differences on real data.
_REAL_METHOD_AGREEMENT_MPS = 5.0


class TestThreeMethodsAgreeLoose:
    def test_all_three_methods_finite_and_close(self, real_reader_type):
        coord = dict(config_earth.simulation["start_coord"])
        results = {}
        for method in ("linear_neighbors", "linear_full", "spline_full"):
            config_earth.forecast["wind_interpolation"] = method
            reader = (
                GFS(config_earth.simulation["start_coord"])
                if real_reader_type == "gfs"
                else ERA5(config_earth.simulation["start_coord"])
            )
            u, v, _, _ = wind_interp(reader, coord)
            assert math.isfinite(u) and math.isfinite(v), (
                f"{method} returned non-finite: u={u}, v={v}"
            )
            results[method] = (u, v)

        speeds = [math.hypot(u, v) for u, v in results.values()]
        spread = max(speeds) - min(speeds)
        assert spread <= _REAL_METHOD_AGREEMENT_MPS, (
            f"interpolation methods disagree by more than "
            f"{_REAL_METHOD_AGREEMENT_MPS} m/s on real data: {results}"
        )
