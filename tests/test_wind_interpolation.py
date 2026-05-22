"""test_wind_interpolation.py - Unit tests for the spline_uv helper.

Covers:
  - Spline reproduces linear within sampled-point precision on smooth data
  - Out-of-bounds queries fall back to linear (clamped to endpoints),
    not extrapolated cubic (avoids the spline-overshoot failure mode
    that motivated extrapolate=False).
  - GFS and ERA5 modules expose identical _spline_uv semantics.

Run:
    PYTHONPATH=src pytest tests/test_wind_interpolation.py -v
"""

import math

import numpy as np
import pytest

from EarthSHAB.GFS import _spline_uv as _spline_uv_gfs
from EarthSHAB.ERA5 import _spline_uv as _spline_uv_era5


@pytest.fixture(params=[_spline_uv_gfs, _spline_uv_era5], ids=["GFS", "ERA5"])
def spline(request):
    return request.param


class TestSplineUV:
    def test_passes_through_sample_points_exactly(self, spline):
        h = np.array([0.0, 5000.0, 10000.0, 15000.0, 20000.0, 25000.0])
        u = np.array([1.0, 3.0, 5.0, 2.0, -1.0, 0.5])
        v = np.array([0.0, 1.0, 2.0, 1.5, -0.5, 0.0])
        for i, alt in enumerate(h):
            uu, vv = spline(alt, h, u, v)
            assert uu == pytest.approx(u[i], abs=1e-6)
            assert vv == pytest.approx(v[i], abs=1e-6)

    def test_interior_query_returns_finite(self, spline):
        h = np.array([0.0, 5000.0, 10000.0, 15000.0, 20000.0])
        u = np.array([1.0, 3.0, 5.0, 2.0, -1.0])
        v = np.array([0.0, 1.0, 2.0, 1.5, -0.5])
        uu, vv = spline(12500.0, h, u, v)
        assert math.isfinite(uu) and math.isfinite(vv)

    def test_above_max_altitude_falls_back_to_linear_clamp(self, spline):
        # alt above h_max → np.interp clamps to the highest value (no overshoot).
        h = np.array([0.0, 10000.0, 20000.0])
        u = np.array([10.0, 5.0, 50.0])  # contrived: spline would overshoot wildly above 20km
        v = np.array([0.0, 0.0,  0.0])
        uu, vv = spline(35000.0, h, u, v)
        assert uu == pytest.approx(50.0)  # clamped to last value, NOT splined past
        assert vv == pytest.approx(0.0)

    def test_below_min_altitude_falls_back_to_linear_clamp(self, spline):
        h = np.array([1000.0, 10000.0, 20000.0])
        u = np.array([10.0, 5.0, 1.0])
        v = np.array([1.0,  2.0, 3.0])
        uu, vv = spline(0.0, h, u, v)
        assert uu == pytest.approx(10.0)
        assert vv == pytest.approx(1.0)

    def test_smooth_linear_input_recovers_input_linearly(self, spline):
        # If the input itself is linear in h, spline should reproduce it (within
        # CubicSpline's natural-boundary convention's small wiggle).
        h = np.array([0.0, 1000.0, 2000.0, 3000.0, 4000.0])
        u = 2.0 * h + 1.0  # linear
        v = -0.5 * h + 7.0
        uu, vv = spline(2500.0, h, u, v)
        assert uu == pytest.approx(2.0 * 2500.0 + 1.0, abs=1e-3)
        assert vv == pytest.approx(-0.5 * 2500.0 + 7.0, abs=1e-3)


class TestDuplicateAltitudes:
    """fill_missing_data clamps NaN values at the top of the GFS profile to
    the last valid altitude, producing duplicate h entries. CubicSpline would
    error on this; _spline_uv must dedupe."""

    def test_duplicate_h_at_top_does_not_crash(self, spline):
        h = np.array([0.0, 5000.0, 10000.0, 15000.0, 20000.0, 20000.0, 20000.0])
        u = np.array([1.0, 2.0, 3.0, 2.5, 1.0, 1.0, 1.0])
        v = np.array([0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 0.0])
        uu, vv = spline(12500.0, h, u, v)
        assert math.isfinite(uu) and math.isfinite(vv)

    def test_below_min_strict_increasing_count_under_4_falls_back(self, spline):
        # If dedup leaves <4 points, CubicSpline can't fit — fall back to linear.
        h = np.array([0.0, 10000.0, 10000.0, 10000.0])
        u = np.array([5.0, 1.0, 1.0, 1.0])
        v = np.array([0.0, 0.0, 0.0, 0.0])
        uu, vv = spline(5000.0, h, u, v)
        assert math.isfinite(uu) and math.isfinite(vv)


class TestModuleParity:
    """GFS and ERA5 _spline_uv must produce identical results — same impl."""

    def test_identical_outputs_on_random_data(self):
        rng = np.random.default_rng(42)
        for _ in range(20):
            h = np.sort(rng.uniform(0, 30000, 12))
            u = rng.uniform(-30, 30, 12)
            v = rng.uniform(-30, 30, 12)
            alt = rng.uniform(0, 30000)
            u_gfs, v_gfs   = _spline_uv_gfs(alt, h, u, v)
            u_era, v_era   = _spline_uv_era5(alt, h, u, v)
            assert u_gfs == pytest.approx(u_era, abs=1e-12)
            assert v_gfs == pytest.approx(v_era, abs=1e-12)
