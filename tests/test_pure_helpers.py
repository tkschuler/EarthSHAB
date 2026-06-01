"""Tests for pure helper functions that don't need a netCDF file.

Surface tested (per Q9(c) "pure tier"):
    Forecast.closestIdx (alias: closest)
    Forecast.windVectorToBearing
    Forecast.interpolateBearing
    Forecast.interpolateBearingTime
    Forecast.get2NearestAltIdxs
    Forecast.fill_missing_data

These tests instantiate nothing — they call the methods on the class itself
(unbound) where signature allows, or on a stub instance built with
``object.__new__`` to skip __init__. This way we test the pure logic without
needing any forecast file or config_earth patching.

The wind interp composite (``wind_alt_Interpolate2``) and the lookup helpers
that read from instance attributes (``getNearestLatIdx`` etc.) are tested in
test_dummy_forecasts.py / test_dummy_lookups.py against real instances.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from EarthSHAB.Forecast import Forecast


# ---------------------------------------------------------------------------
# Helpers: build a "bare" reader instance with no file open.
# ---------------------------------------------------------------------------

def _bare():
    """Allocate an instance without running __init__. Lets us call pure
    methods that don't touch instance attributes."""
    return object.__new__(Forecast)


@pytest.fixture
def bare():
    return _bare()


# ---------------------------------------------------------------------------
# closestIdx
# ---------------------------------------------------------------------------

class TestClosest:
    def test_exact_match(self, bare):
        arr = [0.0, 1.0, 2.0, 3.0, 4.0]
        assert bare.closestIdx(arr, 2.0) == 2

    def test_below_min_returns_first(self, bare):
        arr = [10.0, 20.0, 30.0]
        assert bare.closestIdx(arr, 0.0) == 0

    def test_above_max_returns_last(self, bare):
        arr = [10.0, 20.0, 30.0]
        assert bare.closestIdx(arr, 1000.0) == 2

    def test_interior_query_picks_nearer(self, bare):
        arr = [0.0, 1.0, 5.0]
        # 0.6 is closer to 1.0 than to 0.0
        assert bare.closestIdx(arr, 0.6) == 1
        # 0.4 is closer to 0.0
        assert bare.closestIdx(arr, 0.4) == 0

    def test_equidistant_returns_first(self, bare):
        # min(...) is stable on first-equal: 1.5 is equidistant from 1.0 and 2.0
        arr = [1.0, 2.0]
        assert bare.closestIdx(arr, 1.5) == 0

    def test_negative_values(self, bare):
        arr = [-10.0, -5.0, 0.0, 5.0]
        assert bare.closestIdx(arr, -7.0) == 1

    def test_closest_alias_matches(self, bare):
        # `closest` is a back-compat alias to closestIdx (GFS used the former,
        # ERA5 the latter — Phase 5 unified them). Verify they're the same.
        arr = [0.0, 1.0, 2.0]
        assert bare.closest(arr, 1.6) == bare.closestIdx(arr, 1.6)


# ---------------------------------------------------------------------------
# windVectorToBearing
# ---------------------------------------------------------------------------

class TestWindVectorToBearing:
    def test_east_wind(self, bare):
        # u=1, v=0 → bearing = arctan2(0, 1) = 0 rad, speed = 1
        bearing, speed = bare.windVectorToBearing(1.0, 0.0)
        assert bearing == pytest.approx(0.0, abs=1e-12)
        assert speed == pytest.approx(1.0, abs=1e-12)

    def test_north_wind(self, bare):
        # u=0, v=1 → bearing = pi/2, speed = 1
        bearing, speed = bare.windVectorToBearing(0.0, 1.0)
        assert bearing == pytest.approx(math.pi / 2, abs=1e-12)
        assert speed == pytest.approx(1.0, abs=1e-12)

    def test_speed_pythagoras(self, bare):
        b, s = bare.windVectorToBearing(3.0, 4.0)
        assert s == pytest.approx(5.0, abs=1e-12)

    def test_zero_vector(self, bare):
        b, s = bare.windVectorToBearing(0.0, 0.0)
        # arctan2(0, 0) is well-defined as 0 in numpy.
        assert s == pytest.approx(0.0, abs=1e-12)

    def test_southwest(self, bare):
        # u=-1, v=-1 → bearing = -3pi/4 (or equivalently 5pi/4)
        b, s = bare.windVectorToBearing(-1.0, -1.0)
        assert s == pytest.approx(math.sqrt(2), abs=1e-12)
        assert b == pytest.approx(-3 * math.pi / 4, abs=1e-12)


# ---------------------------------------------------------------------------
# get2NearestAltIdxs
# ---------------------------------------------------------------------------

class TestGet2NearestAltIdxs:
    def test_query_exactly_on_a_point(self, bare):
        h = np.array([0.0, 1000.0, 2000.0, 3000.0])
        # alt_m == h[1]. Implementation: alt_m > h[h_nearest] is False, so
        # h_idx0 = h_nearest - 1, h_idx1 = h_nearest.
        i0, i1 = bare.get2NearestAltIdxs(h, 1000.0)
        assert (i0, i1) == (0, 1)

    def test_query_between_two_points(self, bare):
        h = np.array([0.0, 1000.0, 2000.0, 3000.0])
        i0, i1 = bare.get2NearestAltIdxs(h, 1500.0)
        # 1500 is closer to either; closest returns 1 (first equidistant) so
        # alt_m > h[1] is True → (1, 2).
        assert (i0, i1) == (1, 2)

    def test_query_below_min(self, bare):
        h = np.array([1000.0, 2000.0, 3000.0])
        # closest returns 0 (h[0]=1000); alt_m=500 < h[0], so h_idx0=-1, h_idx1=0.
        i0, i1 = bare.get2NearestAltIdxs(h, 500.0)
        assert (i0, i1) == (-1, 0)

    def test_query_above_max(self, bare):
        h = np.array([1000.0, 2000.0, 3000.0])
        # closest returns 2; alt_m > h[2], so h_idx0=2, h_idx1=3 (out of bounds!).
        # interpolateBearing fixes the boundary; this just documents the raw return.
        i0, i1 = bare.get2NearestAltIdxs(h, 5000.0)
        assert (i0, i1) == (2, 3)


# ---------------------------------------------------------------------------
# interpolateBearing (raw arrays)
# ---------------------------------------------------------------------------

class TestInterpolateBearing:
    def test_at_exact_lower_altitude(self, bare):
        h = np.array([0.0, 1000.0, 2000.0])
        u = np.array([1.0, 2.0, 3.0])
        v = np.array([0.0, 0.0, 0.0])
        dir_deg, speed = bare.interpolateBearing(h, u, v, 0.0)
        # Pure east at alt 0 → bearing 0°, speed 1.
        assert dir_deg == pytest.approx(0.0, abs=1e-9)
        assert speed == pytest.approx(1.0, abs=1e-9)

    def test_midway_no_wrap(self, bare):
        h = np.array([0.0, 1000.0, 2000.0])
        u = np.array([0.0, 1.0, 0.0])
        v = np.array([1.0, 1.0, 1.0])
        # At alt 500, get2NearestAltIdxs picks (0, 1) → bearing0=90°
        # (north), bearing1=45°. Speeds: 1, sqrt(2). Linear in degrees:
        # midway = (90+45)/2 = 67.5°. Speeds: midway = (1+sqrt(2))/2.
        dir_deg, speed = bare.interpolateBearing(h, u, v, 500.0)
        assert dir_deg == pytest.approx(67.5, abs=1e-6)
        assert speed == pytest.approx((1.0 + math.sqrt(2)) / 2, abs=1e-6)

    def test_bearing_wrap_crosses_zero(self, bare):
        # bearing0 = 350°, bearing1 = 10°. Midpoint should be 0° (or 360°),
        # not 180° (which is what naive linear interp would give without
        # the wrap correction).
        h = np.array([0.0, 1000.0])
        u = np.array([
            math.cos(math.radians(350.0)),
            math.cos(math.radians(10.0)),
        ])
        v = np.array([
            math.sin(math.radians(350.0)),
            math.sin(math.radians(10.0)),
        ])
        dir_deg, _ = bare.interpolateBearing(h, u, v, 500.0)
        # Result modulo 360.
        diff = min(abs(dir_deg), abs(360.0 - dir_deg))
        assert diff == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------------------
# interpolateBearingTime
# ---------------------------------------------------------------------------

class TestInterpolateBearingTime:
    def test_midway_no_wrap(self, bare):
        # bearing0 = 90, bearing1 = 0. hour_index = 0.5 (midway between 0 and 1).
        dir_deg, speed = bare.interpolateBearingTime(
            90.0, 10.0, 0.0, 20.0, 0.5
        )
        assert dir_deg == pytest.approx(45.0, abs=1e-9)
        assert speed == pytest.approx(15.0, abs=1e-9)

    def test_wrap_correction(self, bare):
        dir_deg, _ = bare.interpolateBearingTime(
            350.0, 10.0, 10.0, 10.0, 0.5
        )
        diff = min(abs(dir_deg), abs(360.0 - dir_deg))
        assert diff == pytest.approx(0.0, abs=1e-6)

    def test_at_hour_index_zero(self, bare):
        # Reduces to bearing0.
        dir_deg, speed = bare.interpolateBearingTime(
            42.0, 7.5, 200.0, 1.0, 0.0
        )
        assert dir_deg == pytest.approx(42.0, abs=1e-9)
        assert speed == pytest.approx(7.5, abs=1e-9)


# ---------------------------------------------------------------------------
# fill_missing_data
# ---------------------------------------------------------------------------

class TestFillMissingData:
    def test_no_mask_passthrough(self, bare):
        arr = np.ma.array([1.0, 2.0, 3.0, 4.0], mask=False)
        out = bare.fill_missing_data(arr)
        np.testing.assert_array_equal(out, [1.0, 2.0, 3.0, 4.0])

    def test_interior_nan_linearly_filled(self, bare):
        arr = np.ma.array(
            [1.0, 0.0, 3.0],
            mask=[False, True, False],
        )
        out = bare.fill_missing_data(arr)
        # np.interp on the masked index gives the linear midpoint.
        assert out[1] == pytest.approx(2.0, abs=1e-12)

    def test_trailing_nans_clamped_to_last_valid(self, bare):
        # np.interp with no further valid points extends with the last valid
        # value (clamp behavior). This is the "fill_missing_data at top of
        # GFS profile" path that motivates the spline duplicate-h handling.
        arr = np.ma.array(
            [1.0, 2.0, 0.0, 0.0],
            mask=[False, False, True, True],
        )
        out = bare.fill_missing_data(arr)
        # The two trailing masked entries should both equal the last valid (=2.0).
        assert out[2] == pytest.approx(2.0, abs=1e-12)
        assert out[3] == pytest.approx(2.0, abs=1e-12)

    def test_leading_nans_clamped_to_first_valid(self, bare):
        arr = np.ma.array(
            [0.0, 1.0, 5.0],
            mask=[True, False, False],
        )
        out = bare.fill_missing_data(arr)
        assert out[0] == pytest.approx(1.0, abs=1e-12)

    def test_returns_finite_array(self, bare):
        arr = np.ma.array(
            np.arange(10, dtype=float),
            mask=[False, True, False, True, False, False, True, False, True, False],
        )
        out = bare.fill_missing_data(arr)
        assert np.all(np.isfinite(out))

    def test_plain_ndarray_with_nan(self, bare):
        # Phase 5 fill_missing_data also accepts plain ndarrays with literal
        # NaN values (not just MaskedArrays), because the v2 reader's
        # whole-array reads sometimes yield ndarrays.
        arr = np.array([1.0, np.nan, 3.0])
        out = bare.fill_missing_data(arr)
        assert out[1] == pytest.approx(2.0, abs=1e-12)
