"""Composite + smoke tests against synthetic forecasts (Q9(c) "composite tier").

The composite under test is ``wind_alt_Interpolate2``. For every scenario,
we assert:

  (1) on-grid query returns the analytical truth within ON_GRID_EPS,
  (2) off-grid query returns the analytical truth within OFF_GRID_EPS.

v2 collapsed the GFS and ERA5 readers into a single Forecast class, so
the cross-reader equality assertion that used to live here is now redundant
and has been removed.

The smoke tier asserts ``getNewCoord`` returns finite values in physical
bounds for every scenario, without exact-truth assertions (getNewCoord adds
a geodesic step, which is a separate concern from the interpolation we're
testing here).

A small number of scenario-specific behavioral tests live at the end:
  - bearing_wrap: assert ``linear_neighbors`` gives the wrap-aware answer
    and not the naive midpoint-of-degrees.
  - geopotential_vs_height: assert the reader's /g conversion produces the
    expected altitude column.
  - nan_at_top: assert fill_missing_data was actually applied (no NaN
    in the returned u, v even though the dummy has masked top levels).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tests.conftest import (
    OFF_GRID_EPS,
    ON_GRID_EPS,
    _patch_config_for_forecast,
)
from tests.tools.make_dummy_forecasts import (
    ALL_STATIC_U,
    ALL_STATIC_V,
    ALT_M_GFS_ORDER,
    LAT_VALUES_ASC,
    LON_VALUES,
    TIMES_DT,
    SYNTHETIC_START_COORD,
)
from tests.tools.scenario_truth import truth


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _wind_interp(reader, coord):
    """Run the reader's composite altitude+time wind interpolation.

    Translates a GFS-style coord dict (lat/lon/alt/timestamp) into the
    Forecast.wind_alt_Interpolate2 signature (alt_m, diff_time, lat_idx,
    lon_idx). Used to be in tests/tools/reader_adapter.py; inlined here
    after v2 collapsed the two readers.
    """
    diff_time = coord["timestamp"] - reader.model_start_datetime
    lat_idx = reader.getNearestLatIdx(coord["lat"])
    lon_idx = reader.getNearestLonIdx(coord["lon"])
    return reader.wind_alt_Interpolate2(
        coord["alt"], diff_time, lat_idx, lon_idx
    )


def _on_grid_query():
    """A coord that lands exactly on a grid intersection in the synthetic frame.

    Picks middle grid indices so the query is comfortably inside the
    interior (no edge clamping). Time = TIMES_DT[1] — same as the reader's
    fixture start_time, so hour_index resolves to exactly 0.
    """
    return {
        "lat": float(LAT_VALUES_ASC[4]),      # middle of 9 → idx 4
        "lon": float(LON_VALUES[4]),
        "alt": float(ALT_M_GFS_ORDER[4]),     # interior altitude
        "timestamp": TIMES_DT[1],
    }


def _off_grid_query():
    """A coord with on-grid lat/lon but off-grid altitude and time.

    ``wind_alt_Interpolate2`` snaps lat/lon to the nearest stored cell (no
    spatial interpolation) and interpolates only in altitude and time. So
    a meaningfully "off-grid" query for this composite varies only those
    two axes. Lat/lon stay on grid intersections to avoid ambiguous
    equidistant-snap behavior.

    Time = TIMES_DT[1] + 1.5 hours = halfway between TIMES_DT[1] and TIMES_DT[2].
    Alt: halfway (in meters) between two stored altitudes.
    """
    from datetime import timedelta
    return {
        "lat": float(LAT_VALUES_ASC[4]),     # on-grid interior
        "lon": float(LON_VALUES[4]),         # on-grid interior
        "alt": 0.5 * (ALT_M_GFS_ORDER[4] + ALT_M_GFS_ORDER[5]),
        "timestamp": TIMES_DT[1] + timedelta(hours=1.5),
    }


def _query_assert_close(reader, scenario, coord, abs_tol):
    """Run wind_interp; assert the returned (u, v) matches analytical truth."""
    u, v, _u_diag, _v_diag = _wind_interp(reader, coord)
    expected_u, expected_v = truth(
        scenario, coord["timestamp"], coord["alt"], coord["lat"], coord["lon"]
    )
    assert u == pytest.approx(expected_u, abs=abs_tol), (
        f"u mismatch: got {u}, expected {expected_u} "
        f"(scenario={scenario}, coord={coord})"
    )
    assert v == pytest.approx(expected_v, abs=abs_tol), (
        f"v mismatch: got {v}, expected {expected_v} "
        f"(scenario={scenario}, coord={coord})"
    )


# ---------------------------------------------------------------------------
# Composite: on-grid analytical truth
# ---------------------------------------------------------------------------

class TestWindInterpOnGrid:
    """For every scenario, the on-grid query should return the exact
    analytical-truth (u, v) within float32 round-trip precision."""

    def test_on_grid_query_matches_truth(self, reader, scenario):
        if scenario == "nan_at_top":
            # nan_at_top reuses static_by_level truth in the valid column;
            # we pick a coord inside the valid (non-masked) altitude range
            # to avoid depending on fill_missing_data's extrapolation choice.
            coord = _on_grid_query()
            # ALT_M_GFS_ORDER[4] = 5800 m, well inside the valid column.
            _query_assert_close(reader, scenario, coord, ON_GRID_EPS)
            return
        coord = _on_grid_query()
        _query_assert_close(reader, scenario, coord, ON_GRID_EPS)


# ---------------------------------------------------------------------------
# Composite: off-grid analytical truth
# ---------------------------------------------------------------------------

class TestWindInterpOffGrid:
    """Off-grid queries should match analytical truth within a slightly
    looser tolerance (one extra interp arithmetic step over float32 inputs)."""

    def test_off_grid_query_matches_truth(self, reader, scenario):
        if scenario == "nan_at_top":
            coord = _off_grid_query()
            _query_assert_close(reader, scenario, coord, OFF_GRID_EPS)
            return
        coord = _off_grid_query()
        _query_assert_close(reader, scenario, coord, OFF_GRID_EPS)


# ---------------------------------------------------------------------------
# Smoke: getNewCoord doesn't crash, returns finite + in-bounds values
# ---------------------------------------------------------------------------

# Physical sanity bounds for a balloon trajectory step.
_MAX_WIND_SPEED_MPS = 200.0
_MAX_ALT_M = 80000.0


class TestGetNewCoordSmoke:
    """Every scenario, call getNewCoord once. Assert:
      - it returns 10 values (the documented return tuple)
      - lat, lon, u, v, altitude are finite
      - wind speeds are within physical bounds
      - resolved altitude is within physical bounds
    """

    def test_get_new_coord_returns_finite_in_bounds(self, reader, scenario):
        coord = _on_grid_query()
        out = reader.getNewCoord(coord, dt=1.0)
        assert len(out) == 10, f"unexpected return arity {len(out)}: {out}"
        lat_new, lon_new, u, v, u_old, v_old, bearing, c_lat, c_lon, c_alt = out

        # Check finiteness on all scalar returns (skip the bearing if it's
        # NaN; bearing math can NaN for zero wind at altitude exactly matching
        # a grid point in all_static-style scenarios).
        for name, val in [
            ("lat_new", lat_new), ("lon_new", lon_new),
            ("u", u), ("v", v),
            ("c_lat", c_lat), ("c_lon", c_lon), ("c_alt", c_alt),
        ]:
            assert math.isfinite(float(val)), f"{name} not finite: {val}"

        assert abs(float(u)) <= _MAX_WIND_SPEED_MPS
        assert abs(float(v)) <= _MAX_WIND_SPEED_MPS
        assert 0.0 <= float(c_alt) <= _MAX_ALT_M


# ---------------------------------------------------------------------------
# Scenario-specific assertions
# ---------------------------------------------------------------------------

class TestBearingWrap:
    """linear_neighbors should produce the wrap-aware answer, NOT the
    naive midpoint of degrees. Specifically: between bearing 350° and 10°,
    midpoint should be 0° (or 360°), not 180°.
    """

    @pytest.fixture(params=["bearing_wrap"])
    def scenario(self, request):
        return request.param

    def test_wrap_corrected_midpoint(self, reader):
        # Construct a coord at the altitude midpoint of the two adjacent
        # levels where bearing crosses 0°/360°. In _scenario_bearing_wrap,
        # bearing(k) = (-90 + 20k) % 360. Find an adjacent (k, k+1) pair
        # whose bearings are on opposite sides of 0°:
        # k=4: -90 + 80 = -10 → 350°. k=5: -90 + 100 = 10°. Crosses 0°.
        # ALT_M_GFS_ORDER[4] = 5800, [5] = 7400. Midpoint = 6600 m.
        coord = {
            "lat": float(LAT_VALUES_ASC[4]),
            "lon": float(LON_VALUES[4]),
            "alt": 0.5 * (ALT_M_GFS_ORDER[4] + ALT_M_GFS_ORDER[5]),
            "timestamp": TIMES_DT[1],
        }
        u, v, _, _ = _wind_interp(reader, coord)
        bearing_deg = (math.degrees(math.atan2(v, u))) % 360.0
        # Wrap-aware midpoint of 350° and 10° is 0°. Allow tolerance up to a
        # small angular error — the speed is constant 10 m/s so the angle
        # is meaningful.
        diff = min(abs(bearing_deg), abs(360.0 - bearing_deg))
        assert diff == pytest.approx(0.0, abs=1.0)  # 1 degree tolerance

    def test_speed_unchanged_across_wrap(self, reader):
        coord = {
            "lat": float(LAT_VALUES_ASC[4]),
            "lon": float(LON_VALUES[4]),
            "alt": 0.5 * (ALT_M_GFS_ORDER[4] + ALT_M_GFS_ORDER[5]),
            "timestamp": TIMES_DT[1],
        }
        u, v, _, _ = _wind_interp(reader, coord)
        speed = math.hypot(u, v)
        # The scenario uses constant speed = 10 m/s; the interpolated speed
        # should be approximately 10 m/s.
        assert speed == pytest.approx(10.0, abs=OFF_GRID_EPS * 10.0)


class TestGeopotentialVsHeight:
    """For geopotential_vs_height, the reader's internal hgtprs (after the
    /g conversion in __init__) should match the synthetic altitude column.
    """

    @pytest.fixture(params=["geopotential_vs_height"])
    def scenario(self, request):
        return request.param

    def test_altitude_column_matches_after_g_conversion(self, reader):
        """Reader's internal hgtprs at any (t, lat, lon) cell should match
        the synthetic altitude column. Canonical v2 stores pressure_level
        descending in hPa → altitude ASCENDING with level index → reader's
        hgtprs column at index i corresponds to ALT_M_GFS_ORDER[i] directly.
        """
        for i, expected_alt in enumerate(ALT_M_GFS_ORDER):
            stored = float(reader.hgtprs[0, i, 0, 0])
            assert stored == pytest.approx(
                float(expected_alt),
                abs=max(abs(expected_alt) * 1e-5, ON_GRID_EPS),
            ), f"altitude mismatch at level {i}"


class TestNanAtTopFilled:
    """nan_at_top scenario writes masked entries at the topmost altitude
    levels. After fill_missing_data, the reader should return finite
    u, v at any query — even one that maps to the masked altitude range.
    """

    @pytest.fixture(params=["nan_at_top"])
    def scenario(self, request):
        return request.param

    def test_query_in_masked_altitude_returns_finite(self, reader):
        # Mask is on canonical level indices N_LEVEL-3 and N_LEVEL-2, which
        # in ALT_M_GFS_ORDER are 11000 m and 13800 m. Query at 12500 m sits
        # inside the masked altitude band.
        coord = {
            "lat": float(LAT_VALUES_ASC[4]),
            "lon": float(LON_VALUES[4]),
            "alt": 12500.0,
            "timestamp": TIMES_DT[1],
        }
        u, v, _, _ = _wind_interp(reader, coord)
        assert math.isfinite(u), f"u is not finite in masked altitude: {u}"
        assert math.isfinite(v), f"v is not finite in masked altitude: {v}"


class TestAllThreeInterpMethods:
    """For each scenario, all three wind_interpolation methods
    ('linear_neighbors', 'linear_full', 'spline_full') should return finite
    values within wide physical bounds. For all_static specifically, the
    three methods should AGREE to ON_GRID_EPS since the field is constant.
    """

    @pytest.fixture(params=["all_static", "static_by_level"])
    def scenario(self, request):
        return request.param

    def test_three_methods_finite_at_on_grid(
        self, scenario, all_dummies, monkeypatch,
    ):
        from EarthSHAB.Forecast import Forecast
        import EarthSHAB.config_earth as config_earth

        path = all_dummies[scenario]
        coord = _on_grid_query()
        results = {}
        for method in ("linear_neighbors", "linear_full", "spline_full"):
            _patch_config_for_forecast(monkeypatch, path, TIMES_DT[1])
            config_earth.forecast["wind_interpolation"] = method

            reader = Forecast(config_earth.simulation["start_coord"])
            u, v, _, _ = _wind_interp(reader, coord)
            assert math.isfinite(u) and math.isfinite(v), (
                f"non-finite for method={method}, scenario={scenario}: "
                f"u={u}, v={v}"
            )
            results[method] = (u, v)

        if scenario == "all_static":
            for method in ("linear_neighbors", "linear_full", "spline_full"):
                u, v = results[method]
                assert u == pytest.approx(ALL_STATIC_U, abs=ON_GRID_EPS), (
                    f"{method} u mismatch on all_static: {u}"
                )
                assert v == pytest.approx(ALL_STATIC_V, abs=ON_GRID_EPS), (
                    f"{method} v mismatch on all_static: {v}"
                )
