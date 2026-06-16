"""Backend-equivalence tests: numpy / numba / xarray must produce the same wind
(hence the same trajectory step) for a given wind_interpolation method.

The backend is purely an implementation tier (data access + horizontal stencil);
it must not change the physics. numba stores the cube as float32, so agreement is
asserted to a float32-cache tolerance rather than bit-for-bit. numba/xarray are
skipped automatically if not installed.

Also exercises the ``bilinear`` method on every backend.

Run:
    PYTHONPATH=src pytest tests/test_backends.py -v
"""
import math
from datetime import timedelta

import numpy as np
import pytest

import EarthSHAB.config as config
from EarthSHAB.Forecast import Forecast
from tests.tools.make_dummy_forecasts import (
    TIMES_DT, SYNTHETIC_START_COORD, LAT_VALUES_ASC, LON_VALUES, ALT_M_GFS_ORDER,
)

# Available backends (numpy always; numba/xarray only if importable).
BACKENDS = ["numpy"]
for _m in ("numba", "xarray"):
    try:
        __import__(_m)
        BACKENDS.append(_m)
    except Exception:  # pragma: no cover - depends on the install
        pass

METHODS = ["linear_neighbors", "linear_full", "spline_full", "bilinear"]

# Scenarios with real spatial/temporal/vertical variation, so the stencils and
# the alt/time interpolation are actually exercised (a constant field would make
# every method/backend trivially agree).
SCENARIOS = ["linear_ramp_all_dims", "altitude_ramp_only"]

# float32-cache agreement floor (winds here are O(1-10 m/s)).
BACKEND_EPS = 0.05  # m/s


@pytest.fixture
def _no_guard(monkeypatch):
    """Neutralize the config/file res-cadence guard so a dummy file (whatever
    grid it declares) constructs regardless of the live config's netcdf_gfs."""
    monkeypatch.setitem(config.netcdf_gfs, "res", None)
    monkeypatch.setitem(config.netcdf_gfs, "step_hours", None)


def _build(path, backend, wind):
    sc = dict(SYNTHETIC_START_COORD)
    sc["timestamp"] = TIMES_DT[1]
    return Forecast(
        sc, forecast_file=str(path), start_time=TIMES_DT[1], dt=1.0, sim_time=6,
        wind_interpolation=wind, advection="geodesic", backend=backend,
    )


def _query():
    # Off-grid interior point at a time strictly inside a forecast-hour bracket
    # so time interpolation runs. Deliberately NOT the cell midpoint: a 0.3
    # offset keeps the nearest lat/lon cell unambiguous, so the nearest-cell
    # methods can't pick different cells across backends on a tie (that exact-
    # midpoint tie-break ambiguity is real but is not a backend disagreement).
    return {
        "lat": LAT_VALUES_ASC[3] + 0.3 * (LAT_VALUES_ASC[4] - LAT_VALUES_ASC[3]),
        "lon": LON_VALUES[3] + 0.3 * (LON_VALUES[4] - LON_VALUES[3]),
        "alt": 0.5 * (ALT_M_GFS_ORDER[4] + ALT_M_GFS_ORDER[5]),
        "timestamp": TIMES_DT[1] + timedelta(hours=1.0),
    }


@pytest.mark.skipif(len(BACKENDS) < 2, reason="only one backend installed")
@pytest.mark.parametrize("scenario", SCENARIOS)
@pytest.mark.parametrize("method", METHODS)
def test_backends_agree(scenario, method, all_dummies, _no_guard):
    path = all_dummies[scenario]
    coord = _query()
    uv = {}
    for be in BACKENDS:
        f = _build(path, be, method)
        try:
            out = f.getNewCoord(coord, 1.0)
        finally:
            f.close()
        uv[be] = (float(out[2]), float(out[3]))  # x_wind, y_wind

    us = [u for u, v in uv.values()]
    vs = [v for u, v in uv.values()]
    assert all(math.isfinite(x) for x in us + vs), f"non-finite winds: {uv}"
    assert max(us) - min(us) <= BACKEND_EPS, f"u spread too large ({method}/{scenario}): {uv}"
    assert max(vs) - min(vs) <= BACKEND_EPS, f"v spread too large ({method}/{scenario}): {uv}"


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_bilinear_runs_on_default_backend(scenario, all_dummies, _no_guard):
    # bilinear must at least produce finite winds on the default (numpy) backend.
    f = _build(all_dummies[scenario], "numpy", "bilinear")
    try:
        out = f.getNewCoord(_query(), 1.0)
    finally:
        f.close()
    assert math.isfinite(float(out[2])) and math.isfinite(float(out[3]))
