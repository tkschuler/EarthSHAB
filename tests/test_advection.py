"""Unit tests for EarthSHAB.utils.advection (the shared trajectory-advection core).

Covers the two models (geodesic, tangent_plane) and, in particular, guards the
111000/111320 meters-per-degree bug: the tangent-plane step re-derives meters
from the current position against the FIXED launch anchor every step, so if the
forward and inverse conversions use different constants the per-step round-trip
stops being an identity and the trajectory drifts toward the anchor under zero
wind. These are pure tests — no netCDF / config needed.

Run:
    PYTHONPATH=src pytest tests/test_advection.py -v
"""
import math

import numpy as np
import pytest
from geographiclib.geodesic import Geodesic

from EarthSHAB.utils import advection

GEOD = Geodesic.WGS84
ANCHOR_LAT, ANCHOR_LON = 35.0, -100.0
COSC = math.cos(math.radians(ANCHOR_LAT))


def _tan(lats, lons, u, v, dt):
    return advection.advect(
        "tangent_plane", lats, lons, u, v, dt,
        central_lat=ANCHOR_LAT, central_lon=ANCHOR_LON, cos_central=COSC)


def _geo(lats, lons, u, v, dt):
    return advection.advect("geodesic", lats, lons, u, v, dt, geod=GEOD)


# ---------------------------------------------------------------------------
# Zero-wind invariant: a balloon in still air must not move (either model).
# This is the direct regression guard for the 111000/111320 anchor-drift bug.
# ---------------------------------------------------------------------------

class TestZeroWindNoDrift:
    def test_tangent_plane_zero_wind_900_steps(self):
        # A point well away from the anchor, advected 900 times under zero wind.
        lat = np.array([37.0]); lon = np.array([-98.0])
        z = np.array([0.0])
        for _ in range(900):
            lat, lon, _x, _y = _tan(lat, lon, z, z, 60.0)
        assert lat[0] == pytest.approx(37.0, abs=1e-9)
        assert lon[0] == pytest.approx(-98.0, abs=1e-9)

    def test_geodesic_zero_wind_900_steps(self):
        lat = np.array([37.0]); lon = np.array([-98.0])
        z = np.array([0.0])
        for _ in range(900):
            lat, lon, _x, _y = _geo(lat, lon, z, z, 60.0)
        assert lat[0] == pytest.approx(37.0, abs=1e-9)
        assert lon[0] == pytest.approx(-98.0, abs=1e-9)

    def test_tangent_plane_roundtrip_is_identity(self):
        # Forward (lat/lon -> meters) and inverse (meters -> lat/lon) must use the
        # SAME constant, else zero wind maps a point off itself. Check many points.
        rng = np.random.default_rng(0)
        lats = rng.uniform(20, 50, 200)
        lons = rng.uniform(-120, -80, 200)
        z = np.zeros(200)
        lat_new, lon_new, _x, _y = _tan(lats, lons, z, z, 60.0)
        assert np.max(np.abs(lat_new - lats)) < 1e-9
        assert np.max(np.abs(lon_new - lons)) < 1e-9


# ---------------------------------------------------------------------------
# Direction / magnitude sanity
# ---------------------------------------------------------------------------

class TestStepDirection:
    @pytest.mark.parametrize("advect_fn", [_tan, _geo])
    def test_pure_east_wind_moves_east(self, advect_fn):
        lat = np.array([ANCHOR_LAT]); lon = np.array([ANCHOR_LON])
        u = np.array([10.0]); v = np.array([0.0])  # 10 m/s due east
        lat_new, lon_new, _x, _y = advect_fn(lat, lon, u, v, 100.0)
        assert lon_new[0] > lon[0]                    # moved east
        assert lat_new[0] == pytest.approx(lat[0], abs=1e-4)  # ~no N/S motion

    @pytest.mark.parametrize("advect_fn", [_tan, _geo])
    def test_pure_north_wind_moves_north(self, advect_fn):
        lat = np.array([ANCHOR_LAT]); lon = np.array([ANCHOR_LON])
        u = np.array([0.0]); v = np.array([10.0])  # 10 m/s due north
        lat_new, lon_new, _x, _y = advect_fn(lat, lon, u, v, 100.0)
        assert lat_new[0] > lat[0]                     # moved north
        assert lon_new[0] == pytest.approx(lon[0], abs=1e-4)


class TestModelsAgreeForSmallSteps:
    def test_tangent_tracks_geodesic_one_step(self):
        # At the anchor, one moderate step: the two models should agree to well
        # under a meter (their divergence accumulates only over long distances).
        lat = np.array([ANCHOR_LAT]); lon = np.array([ANCHOR_LON])
        u = np.array([15.0]); v = np.array([-8.0])
        tl, to, _x, _y = _tan(lat, lon, u, v, 60.0)
        gl, go, _gx, _gy = _geo(lat, lon, u, v, 60.0)
        dlat_m = (tl[0] - gl[0]) * 111320.0
        dlon_m = (to[0] - go[0]) * 111320.0 * COSC
        # The genuine flat-earth-vs-WGS84 gap is a couple of meters for a ~1 km
        # step (111320 m/deg vs the true ~110941 m/deg meridian at 35 deg); it
        # accumulates over a flight but stays tiny per step. Far below the
        # 100s-of-km the 111000/111320 bug produced.
        assert math.hypot(dlat_m, dlon_m) < 5.0


# ---------------------------------------------------------------------------
# Displacement-return semantics + dispatch
# ---------------------------------------------------------------------------

class TestDisplacementAndDispatch:
    def test_geodesic_displacement_is_incremental(self):
        u = np.array([12.0, -3.0]); v = np.array([-5.0, 7.0])
        lats = np.array([ANCHOR_LAT, 40.0]); lons = np.array([ANCHOR_LON, -95.0])
        _l, _o, x, y = _geo(lats, lons, u, v, 60.0)
        assert np.allclose(x, u * 60.0)
        assert np.allclose(y, v * 60.0)

    def test_tangent_displacement_is_cumulative_from_anchor(self):
        # x_disp/y_disp = current-from-anchor displacement + this step's wind.
        lats = np.array([37.0]); lons = np.array([-98.0])
        u = np.array([10.0]); v = np.array([4.0]); dt = 60.0
        _l, _o, x, y = _tan(lats, lons, u, v, dt)
        cx = (lons[0] - ANCHOR_LON) * advection.M_PER_DEG * COSC + u[0] * dt
        cy = (lats[0] - ANCHOR_LAT) * advection.M_PER_DEG + v[0] * dt
        assert x[0] == pytest.approx(cx)
        assert y[0] == pytest.approx(cy)

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError):
            advection.advect("bogus", np.array([0.0]), np.array([0.0]),
                             np.array([0.0]), np.array([0.0]), 60.0)
