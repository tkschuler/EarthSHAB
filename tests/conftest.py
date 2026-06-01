"""Pytest fixtures shared across test files.

Provides:

* Tolerance constants.
* Session-scoped ``all_dummies`` fixture — writes every scenario synthetic
  .nc once per pytest run into tmp_path_factory.
* ``scenario`` parametrization fixture.
* ``reader`` fixture — monkey-patches ``config_earth`` to point at the
  scenario's dummy file with synthetic start_time / start_coord, then
  instantiates the (single, post-Phase-5) Forecast class.
* ``real_reader`` fixture — leaves ``config_earth`` unpatched, instantiates
  the reader the config currently points at. Used by test_real_forecasts.py.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

import EarthSHAB.config_earth as config_earth
from EarthSHAB.Forecast import Forecast

from tests.tools.make_dummy_forecasts import (
    SCENARIOS,
    SYNTHETIC_START_COORD,
    TIMES_DT,
    build_dummy,
)


# ---------------------------------------------------------------------------
# Tolerance constants
# ---------------------------------------------------------------------------

ON_GRID_EPS = 1e-5     # float32 round-trip precision floor
OFF_GRID_EPS = 1e-4    # one extra interp arithmetic step


# ---------------------------------------------------------------------------
# Synthetic dummy file fixture (session-scoped — built once per pytest run)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def all_dummies(tmp_path_factory):
    """Build one canonical-schema .nc per scenario into a session tmp dir.

    Returns ``dict[scenario_name] -> Path``.
    """
    base = tmp_path_factory.mktemp("dummy_forecasts")
    files = {}
    for scenario in SCENARIOS:
        out = base / f"{scenario}.nc"
        build_dummy(scenario, out)
        files[scenario] = out
    return files


# ---------------------------------------------------------------------------
# Parametrization fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(params=list(SCENARIOS.keys()))
def scenario(request):
    """Parametrize tests over every scenario name."""
    return request.param


# ---------------------------------------------------------------------------
# config_earth patcher
# ---------------------------------------------------------------------------

# Sim runtime kept small but large enough that the synthetic time window
# (4 × 3 hr = 12 hours total) covers it. Tests use TIMES_DT[1] as start.
_SYNTHETIC_SIM_TIME_HOURS = 6


def _synthetic_simulation_dict(start_dt: datetime) -> dict:
    """The dict that would normally live at ``config_earth.simulation``,
    populated with synthetic-frame values."""
    start_coord = dict(SYNTHETIC_START_COORD)
    start_coord["timestamp"] = start_dt
    return dict(
        start_time=start_dt,
        sim_time=_SYNTHETIC_SIM_TIME_HOURS,
        vent=0.0,
        alt_sp=15000.0,
        v_sp=0.0,
        start_coord=start_coord,
        min_alt=start_coord["alt"],
        float=23000,
        dt=1.0,
        balloon_trajectory=None,
    )


def _patch_config_for_forecast(monkeypatch, nc_path: Path, start_dt: datetime) -> None:
    """Point config_earth.forecast['file'] at nc_path and prime simulation state."""
    forecast_dict = dict(
        file=str(nc_path),
        forecast_start_time=start_dt.isoformat(sep=" "),
        forecast_update_interval=60,
        wind_interpolation="linear_neighbors",
    )
    sim_dict = _synthetic_simulation_dict(start_dt)
    monkeypatch.setattr(config_earth, "forecast", forecast_dict, raising=True)
    monkeypatch.setattr(config_earth, "simulation", sim_dict, raising=True)
    monkeypatch.setattr(config_earth, "start_time", start_dt, raising=True)


@pytest.fixture
def reader(scenario, all_dummies, monkeypatch):
    """Construct a Forecast reader bound to the scenario's dummy file.

    Picks ``start_time = TIMES_DT[1]`` (the second time step) so:
      * it's strictly inside the time range (not at t=0 or t=end),
      * leaves at least 2 time steps after it for the simulator's
        ``sim_time`` window check,
      * and lets `wind_alt_Interpolate2` resolve a [t1, t2] bracket
        without hitting the upper boundary.
    """
    nc_path = all_dummies[scenario]
    start_dt = TIMES_DT[1]
    _patch_config_for_forecast(monkeypatch, nc_path, start_dt)
    return Forecast(config_earth.simulation["start_coord"])


@pytest.fixture
def synthetic_start_dt():
    """Convenience accessor for the synthetic start datetime used by `reader`."""
    return TIMES_DT[1]


# ---------------------------------------------------------------------------
# Real-forecast fixture
# ---------------------------------------------------------------------------

@pytest.fixture
def real_reader():
    """Instantiate the Forecast reader against the file the unmodified
    config_earth points at (typically src/EarthSHAB/forecasts/gfs_0p25_*.nc).
    """
    return Forecast(config_earth.simulation["start_coord"])
