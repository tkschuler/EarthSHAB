"""Pytest fixtures shared across test files.

Provides:

* Tolerance constants — single source of truth, per Q11(b).
* Session-scoped ``all_dummies`` fixture — writes every (scenario, schema)
  synthetic .nc once per pytest run into tmp_path_factory.
* ``scenario`` and ``reader_type`` parametrization fixtures.
* ``reader`` fixture — monkey-patches ``config_earth`` to point at the right
  dummy file with synthetic start_time / start_coord / nc_start / etc.,
  then instantiates the corresponding GFS or ERA5 class.
* ``real_reader`` fixture — leaves ``config_earth`` unpatched, instantiates
  whichever class the config currently points at. Used by
  test_real_forecasts.py.

Important: the ``reader`` fixture patches MODULE-LEVEL state in
``EarthSHAB.config_earth``. monkeypatch.setattr restores it after each test.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

import EarthSHAB.config_earth as config_earth
from EarthSHAB.ERA5 import ERA5
from EarthSHAB.GFS import GFS

from tests.tools.make_dummy_forecasts import (
    SCENARIOS,
    SCHEMAS,
    SYNTHETIC_START_COORD,
    TIMES_DT,
    build_dummy,
)


# ---------------------------------------------------------------------------
# Tolerance constants (Q11(b))
# ---------------------------------------------------------------------------

ON_GRID_EPS = 1e-5     # float32 round-trip precision floor
OFF_GRID_EPS = 1e-4    # one extra interp arithmetic step
CROSS_READER_EPS = 1e-5  # both readers from same synthetic source


# ---------------------------------------------------------------------------
# Synthetic dummy file fixture (session-scoped — built once per pytest run)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def all_dummies(tmp_path_factory):
    """Build one .nc per (scenario, schema) into a session tmp dir.

    Returns ``dict[(scenario_name, schema_name)] -> Path``.
    """
    base = tmp_path_factory.mktemp("dummy_forecasts")
    files = {}
    for scenario in SCENARIOS:
        for schema in SCHEMAS:
            out = base / f"{scenario}_{schema}.nc"
            build_dummy(scenario, schema, out)
            files[(scenario, schema)] = out
    return files


# ---------------------------------------------------------------------------
# Parametrization fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(params=list(SCENARIOS.keys()))
def scenario(request):
    """Parametrize tests over every scenario name."""
    return request.param


@pytest.fixture(params=list(SCHEMAS))
def reader_type(request):
    """Parametrize tests over each reader schema."""
    return request.param


# ---------------------------------------------------------------------------
# config_earth patcher (Q5(a))
# ---------------------------------------------------------------------------

# Sim runtime kept small but large enough that the synthetic time window
# (4 × 3hr = 12 hours total) covers it. Q12(a) start time is TIMES_DT[1].
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


def _patch_config_for_gfs(monkeypatch, nc_path: Path, start_dt: datetime) -> None:
    """Set every config_earth global that GFS.__init__ reads."""
    forecast_dict = dict(
        forecast_type="GFS",
        forecast_start_time=start_dt.isoformat(sep=" "),
        GFSrate=60,
        wind_interpolation="linear_neighbors",
    )
    # nc_start must be the FILE'S t=0, not the simulation start. GFS computes
    # hour_index = (coord.timestamp - nc_start) / 3hr. If nc_start != file[0],
    # all queries are off by the offset.
    netcdf_gfs_dict = dict(
        nc_file=str(nc_path),
        nc_start=TIMES_DT[0],
        hourstamp=f"{TIMES_DT[0].hour:02d}",
        res=0.25,
        lat_range=40,
        lon_range=60,
        download_days=1,
    )
    sim_dict = _synthetic_simulation_dict(start_dt)
    monkeypatch.setattr(config_earth, "forecast", forecast_dict, raising=True)
    monkeypatch.setattr(config_earth, "netcdf_gfs", netcdf_gfs_dict, raising=True)
    monkeypatch.setattr(config_earth, "simulation", sim_dict, raising=True)
    monkeypatch.setattr(config_earth, "start_time", start_dt, raising=True)


def _patch_config_for_era5(monkeypatch, nc_path: Path, start_dt: datetime) -> None:
    """Set every config_earth global that ERA5.__init__ reads.

    Important: ERA5.__init__ opens ``"src/EarthSHAB/forecasts/" + filename``
    rather than treating ``filename`` as a full path. To point at a
    tmp_path file, we patch the netCDF4.Dataset call indirectly by setting
    ``netcdf_era5["filename"]`` to a relative path that, combined with the
    hardcoded prefix, resolves to the tmp file. The cleanest way is to
    monkeypatch ERA5's file-open to bypass the prefix.

    Concretely: we replace ``EarthSHAB.ERA5.netCDF4.Dataset`` with a thin
    wrapper that opens ``nc_path`` whenever the filename starts with the
    hardcoded prefix the production code prepends. That keeps production
    code unchanged.
    """
    import EarthSHAB.ERA5 as era5_module

    forecast_dict = dict(
        forecast_type="ERA5",
        forecast_start_time=start_dt.isoformat(sep=" "),
        GFSrate=60,
        wind_interpolation="linear_neighbors",
    )
    netcdf_era5_dict = dict(
        filename=nc_path.name,
        resolution_hr=3,  # synthetic dummy steps every 3 hours
    )
    sim_dict = _synthetic_simulation_dict(start_dt)
    monkeypatch.setattr(config_earth, "forecast", forecast_dict, raising=True)
    monkeypatch.setattr(config_earth, "netcdf_era5", netcdf_era5_dict, raising=True)
    monkeypatch.setattr(config_earth, "simulation", sim_dict, raising=True)
    monkeypatch.setattr(config_earth, "start_time", start_dt, raising=True)

    # Redirect ERA5's hardcoded "src/EarthSHAB/forecasts/" + filename open
    # to the tmp_path file. We patch netCDF4.Dataset only inside the ERA5
    # module so other tests aren't affected.
    original_dataset = era5_module.netCDF4.Dataset

    def _redirected_dataset(path, *args, **kwargs):
        # Any path whose basename matches our synthetic file's basename is
        # redirected to the actual tmp_path file.
        if Path(path).name == nc_path.name:
            return original_dataset(str(nc_path), *args, **kwargs)
        return original_dataset(path, *args, **kwargs)

    monkeypatch.setattr(era5_module.netCDF4, "Dataset", _redirected_dataset)


@pytest.fixture
def reader(scenario, reader_type, all_dummies, monkeypatch):
    """Construct a GFS or ERA5 reader bound to the (scenario, schema) dummy.

    Picks ``start_time = TIMES_DT[1]`` (the second time step) so:
      * it's strictly inside the time range (not at t=0 or t=end),
      * leaves at least 2 time steps after it for the simulator's
        ``sim_time`` window check,
      * and lets `wind_alt_Interpolate2` resolve a [t1, t2] bracket
        without hitting the upper boundary.
    """
    nc_path = all_dummies[(scenario, reader_type)]
    start_dt = TIMES_DT[1]

    if reader_type == "gfs":
        _patch_config_for_gfs(monkeypatch, nc_path, start_dt)
        return GFS(config_earth.simulation["start_coord"])

    _patch_config_for_era5(monkeypatch, nc_path, start_dt)
    return ERA5(config_earth.simulation["start_coord"])


@pytest.fixture
def synthetic_start_dt():
    """Convenience accessor for the synthetic start datetime used by `reader`."""
    return TIMES_DT[1]


# ---------------------------------------------------------------------------
# Real-forecast fixtures (Q6(b))
# ---------------------------------------------------------------------------

@pytest.fixture(params=list(SCHEMAS))
def real_reader_type(request):
    """Parametrize the real-forecast smoke tests over GFS and ERA5."""
    return request.param


@pytest.fixture
def real_reader(real_reader_type):
    """Instantiate whichever real reader the unmodified config_earth points at.

    No patching — exercises the actual config-pointed files in
    src/EarthSHAB/forecasts/. The two files are:
        - gfs_0p25_20220822_12.nc  (GFS, derived from forecast_start_time)
        - SHAB14V_ERA5_20220822_20220823.nc  (ERA5, from netcdf_era5["filename"])
    """
    if real_reader_type == "gfs":
        return GFS(config_earth.simulation["start_coord"])
    return ERA5(config_earth.simulation["start_coord"])
