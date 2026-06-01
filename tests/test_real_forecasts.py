"""Smoke + invariant tests against the real forecast pinned in config_earth.

Per Q6(b), the real-forecast tier runs against ONLY the file that
``config_earth`` currently points at — not every .nc in the forecasts
directory. The file is:

    src/EarthSHAB/forecasts/<config_earth.forecast['file']>
    (typically a GFS file derived from forecast_start_time, or an ERA5 file
     pointed at manually)

Phase 5 collapsed the two readers, so these tests now exercise a single
Forecast instance regardless of whether the file is GFS- or ERA5-sourced.
The Forecast.source attribute records the upstream provenance and is
asserted to be one of the expected values.

These tests are intentionally LOOSE — they do not assert analytical truth
(impossible on real data) and do not pin physical values. They assert:

  1. The reader loads without raising.
  2. The schema has the expected v2 canonical variables and dimensions.
  3. The reader's source is one of the expected provenance labels.
  4. The configured start_coord lat/lon is within the file's lat/lon bounds.
  5. The configured start_time is within the file's time range.
  6. A sample query at the configured start_coord returns finite values
     within physical bounds.
  7. All three wind_interpolation methods produce finite values within
     wide tolerance of each other at the same query.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import EarthSHAB.config_earth as config_earth
from EarthSHAB.Forecast import Forecast


# Skip the whole module on CI / dev boxes that don't have the configured
# forecast file on disk. Current-day GFS predictions aren't committed to
# the repo; the dummy-forecast tests already cover everything that can be
# tested without a real file.
if not Path(config_earth.forecast["file"]).exists():
    pytest.skip(
        f"configured forecast file not present: "
        f"{config_earth.forecast['file']} — set forecast['file'] to a "
        f"file that exists (e.g. src/EarthSHAB/forecasts/"
        f"SHAB14V_ERA5_20220822_20220823.nc) to run this tier",
        allow_module_level=True,
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
    """Smoke: instantiating the reader against the configured file doesn't raise."""

    def test_forecast_reader_loads(self):
        # Sanity-check the file exists before opening — gives a clearer
        # failure than netCDF4's I/O error.
        nc_path = Path(config_earth.forecast["file"])
        assert nc_path.exists(), f"missing real forecast file: {nc_path}"
        reader = Forecast(config_earth.simulation["start_coord"])
        assert reader.file is not None


# ---------------------------------------------------------------------------
# 2. Schema
# ---------------------------------------------------------------------------

class TestSchema:
    """The configured file conforms to the v2 canonical variable + dim layout."""

    def test_v2_canonical_schema(self, real_reader):
        ds = real_reader.file
        for var in ("valid_time", "pressure_level", "latitude", "longitude", "u", "v", "z"):
            assert var in ds.variables, f"missing variable {var!r}"
        for dim in ("valid_time", "pressure_level", "latitude", "longitude"):
            assert dim in ds.dimensions, f"missing dimension {dim!r}"


# ---------------------------------------------------------------------------
# 3. Source provenance
# ---------------------------------------------------------------------------

class TestSource:
    """Forecast.source is populated from the file's `institution` attr (with
    filename fallback). It should be one of the recognised provenance labels.
    """

    def test_source_recognised(self, real_reader):
        assert real_reader.source in ("GFS", "ERA5", "unknown"), (
            f"unexpected Forecast.source: {real_reader.source!r}"
        )


# ---------------------------------------------------------------------------
# 4. Lat/lon bounds contain the configured start coord
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
        # v2 canonical lon is in [-180, 180); no source-specific wrap needed.
        lon = float(config_earth.simulation["start_coord"]["lon"])
        lon_arr = np.asarray(real_reader.lon, dtype=float)
        assert lon_arr.min() <= lon <= lon_arr.max(), (
            f"start_coord lon {lon} outside file range "
            f"[{lon_arr.min()}, {lon_arr.max()}]"
        )


# ---------------------------------------------------------------------------
# 5. Start time is within the file's time range
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
# 6. Sample query at start coord returns finite, in-bounds values
# ---------------------------------------------------------------------------

class TestSampleQueryFinite:
    def test_get_new_coord_at_start_finite(self, real_reader):
        coord = dict(config_earth.simulation["start_coord"])
        out = real_reader.getNewCoord(coord, dt=1.0)
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
        diff_time = coord["timestamp"] - real_reader.model_start_datetime
        lat_idx = real_reader.getNearestLatIdx(coord["lat"])
        lon_idx = real_reader.getNearestLonIdx(coord["lon"])
        u, v, u_diag, v_diag = real_reader.wind_alt_Interpolate2(
            coord["alt"], diff_time, lat_idx, lon_idx
        )
        for name, val in [("u", u), ("v", v), ("u_diag", u_diag), ("v_diag", v_diag)]:
            assert math.isfinite(float(val)), f"{name} not finite: {val}"


# ---------------------------------------------------------------------------
# 7. All three interpolation methods agree within loose tolerance
# ---------------------------------------------------------------------------

# Loose tolerance: real altitude profiles are non-linear, so different
# interpolation methods CAN diverge — but only by O(1 m/s), not by orders
# of magnitude. Tight enough to catch a "wrong method returns garbage" bug,
# loose enough to permit genuine method differences on real data.
_REAL_METHOD_AGREEMENT_MPS = 5.0


class TestThreeMethodsAgreeLoose:
    def test_all_three_methods_finite_and_close(self):
        coord = dict(config_earth.simulation["start_coord"])
        results = {}
        for method in ("linear_neighbors", "linear_full", "spline_full"):
            config_earth.forecast["wind_interpolation"] = method
            reader = Forecast(config_earth.simulation["start_coord"])
            diff_time = coord["timestamp"] - reader.model_start_datetime
            lat_idx = reader.getNearestLatIdx(coord["lat"])
            lon_idx = reader.getNearestLonIdx(coord["lon"])
            u, v, _, _ = reader.wind_alt_Interpolate2(
                coord["alt"], diff_time, lat_idx, lon_idx
            )
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
