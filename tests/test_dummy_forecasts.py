"""Composite + smoke tests against synthetic forecasts (Q9(c) "composite tier").

The composite under test is ``wind_alt_Interpolate2`` — exercised through
the reader_adapter so the same test body runs against both GFS and ERA5.
For every scenario × reader combination, we assert:

  (1) on-grid query returns the analytical truth within ON_GRID_EPS,
  (2) off-grid query returns the analytical truth within OFF_GRID_EPS,
  (3) GFS and ERA5 produce the same wind at the same physical query
      within CROSS_READER_EPS — this is the consolidation safety net.

The smoke tier asserts ``getNewCoord`` returns finite values in physical
bounds for every scenario × reader, without exact-truth assertions
(getNewCoord adds a geodesic step, which is a separate concern from the
interpolation we're testing here).

A small number of scenario-specific behavioral tests live at the end:
  - bearing_wrap: assert ``linear_neighbors`` gives the wrap-aware answer
    and not the naive midpoint-of-degrees.
  - geopotential_vs_height: assert the ERA5 /g conversion produces the
    same physical altitude as GFS at the same query.
  - nan_at_top: assert fill_missing_data was actually applied (no NaN
    in the returned u, v even though the dummy has masked top levels).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tests.conftest import (
    CROSS_READER_EPS,
    OFF_GRID_EPS,
    ON_GRID_EPS,
    _patch_config_for_era5,
    _patch_config_for_gfs,
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
from tests.tools.reader_adapter import (
    get_new_coord,
    wind_interp,
)
from tests.tools.scenario_truth import truth


# ---------------------------------------------------------------------------
# Helpers: build query coords.
# ---------------------------------------------------------------------------

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
    equidistant-snap behavior that diverges between GFS-asc and ERA5-desc
    lat orderings (a real production cross-reader divergence at midpoints).

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
    u, v, _u_diag, _v_diag = wind_interp(reader, coord)
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
    """For every (scenario, reader_type), the on-grid query should return
    the exact analytical-truth (u, v) within float32 round-trip precision."""

    def test_on_grid_query_matches_truth(self, reader, scenario):
        if scenario == "nan_at_top":
            # nan_at_top reuses static_by_level truth in the valid column;
            # we pick a coord inside the valid (non-masked) altitude range
            # to avoid depending on fill_missing_data's extrapolation choice.
            coord = _on_grid_query()
            # ALT_M_GFS_ORDER[4] = 7400 m, well inside the valid column
            # (masked entries are indices 0-1 = 16500, 13800 m).
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
# Cross-reader equality — the consolidation safety net
# ---------------------------------------------------------------------------

class TestCrossReaderAgreement:
    """For each scenario, GFS and ERA5 readers must return identical wind
    (within CROSS_READER_EPS) for the same physical query. This is the
    invariant that the planned consolidated forecast format must preserve.

    Builds both readers in the same test function to share the same
    monkeypatch lifecycle — the second patch overwrites the first's
    config_earth state, but each reader has captured what it needs at
    construction time, so post-construction calls are stable.
    """

    @pytest.fixture(params=["on_grid", "off_grid"])
    def query_coord(self, request):
        return _on_grid_query() if request.param == "on_grid" else _off_grid_query()

    def test_gfs_and_era5_agree(
        self, scenario, query_coord, all_dummies, monkeypatch,
    ):
        from EarthSHAB.ERA5 import ERA5
        from EarthSHAB.GFS import GFS
        import EarthSHAB.config_earth as config_earth

        path_gfs = all_dummies[(scenario, "gfs")]
        path_era5 = all_dummies[(scenario, "era5")]

        _patch_config_for_gfs(monkeypatch, path_gfs, TIMES_DT[1])
        gfs_reader = GFS(config_earth.simulation["start_coord"])

        _patch_config_for_era5(monkeypatch, path_era5, TIMES_DT[1])
        era5_reader = ERA5(config_earth.simulation["start_coord"])

        u_gfs, v_gfs, *_ = wind_interp(gfs_reader, query_coord)
        u_era5, v_era5, *_ = wind_interp(era5_reader, query_coord)

        assert u_gfs == pytest.approx(u_era5, abs=CROSS_READER_EPS), (
            f"u differs between GFS and ERA5: gfs={u_gfs}, era5={u_era5} "
            f"(scenario={scenario})"
        )
        assert v_gfs == pytest.approx(v_era5, abs=CROSS_READER_EPS), (
            f"v differs between GFS and ERA5: gfs={v_gfs}, era5={v_era5} "
            f"(scenario={scenario})"
        )


# ---------------------------------------------------------------------------
# Smoke: getNewCoord doesn't crash, returns finite + in-bounds values
# ---------------------------------------------------------------------------

# Physical sanity bounds for a balloon trajectory step.
_MAX_WIND_SPEED_MPS = 200.0
_MAX_ALT_M = 80000.0


class TestGetNewCoordSmoke:
    """Every (scenario, reader_type), call getNewCoord once. Assert:
      - it returns 10 values (the documented return tuple)
      - lat, lon, u, v, altitude are finite
      - wind speeds are within physical bounds
      - resolved altitude is within physical bounds
    """

    def test_get_new_coord_returns_finite_in_bounds(self, reader, scenario):
        coord = _on_grid_query()
        out = get_new_coord(reader, coord, dt=1.0)
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
        # ALT_M_GFS_ORDER[4] = 7400, [5] = 5800. Midpoint = 6600 m.
        coord = {
            "lat": float(LAT_VALUES_ASC[4]),
            "lon": float(LON_VALUES[4]),
            "alt": 0.5 * (ALT_M_GFS_ORDER[4] + ALT_M_GFS_ORDER[5]),
            "timestamp": TIMES_DT[1],
        }
        u, v, _, _ = wind_interp(reader, coord)
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
        u, v, _, _ = wind_interp(reader, coord)
        speed = math.hypot(u, v)
        # The scenario uses constant speed = 10 m/s; the interpolated speed
        # should be approximately 10 m/s.
        assert speed == pytest.approx(10.0, abs=OFF_GRID_EPS * 10.0)


class TestGeopotentialVsHeight:
    """For geopotential_vs_height, the readers must produce identical
    physical altitudes — even though one stores h and the other stores h*g."""

    @pytest.fixture(params=["geopotential_vs_height"])
    def scenario(self, request):
        return request.param

    def test_altitude_column_matches_after_g_conversion(self, reader, reader_type):
        """Reader's internal hgtprs/z (after /g) column at any (t, lat, lon)
        cell should match the synthetic altitude column, accounting for the
        ERA5 writer's level-axis reversal.
        """
        # GFS stores ascending altitude with index; ERA5 stores descending.
        expected_column = (
            ALT_M_GFS_ORDER if reader_type == "gfs" else ALT_M_GFS_ORDER[::-1]
        )
        for i, expected_alt in enumerate(expected_column):
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
        # Mask is on level indices 0-1 = altitudes 16500 m, 13800 m.
        # Query at 15000 m — within the masked zone.
        coord = {
            "lat": float(LAT_VALUES_ASC[4]),
            "lon": float(LON_VALUES[4]),
            "alt": 15000.0,
            "timestamp": TIMES_DT[1],
        }
        u, v, _, _ = wind_interp(reader, coord)
        assert math.isfinite(u), f"u is not finite in masked altitude: {u}"
        assert math.isfinite(v), f"v is not finite in masked altitude: {v}"


class TestAllThreeInterpMethods:
    """For every scenario, all three wind_interpolation methods
    ('linear_neighbors', 'linear_full', 'spline_full') should return finite
    values within wide physical bounds. For all_static specifically, the
    three methods should AGREE to ON_GRID_EPS since the field is constant.
    """

    @pytest.fixture(params=["all_static", "static_by_level"])
    def scenario(self, request):
        return request.param

    def test_three_methods_finite_at_on_grid(
        self, scenario, reader_type, all_dummies, monkeypatch,
    ):
        from EarthSHAB.ERA5 import ERA5
        from EarthSHAB.GFS import GFS
        import EarthSHAB.config_earth as config_earth

        path = all_dummies[(scenario, reader_type)]
        coord = _on_grid_query()
        results = {}
        for method in ("linear_neighbors", "linear_full", "spline_full"):
            if reader_type == "gfs":
                _patch_config_for_gfs(monkeypatch, path, TIMES_DT[1])
            else:
                _patch_config_for_era5(monkeypatch, path, TIMES_DT[1])
            # Update wind_interpolation in the patched forecast dict.
            config_earth.forecast["wind_interpolation"] = method

            reader = (
                GFS(config_earth.simulation["start_coord"])
                if reader_type == "gfs"
                else ERA5(config_earth.simulation["start_coord"])
            )
            u, v, _, _ = wind_interp(reader, coord)
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
