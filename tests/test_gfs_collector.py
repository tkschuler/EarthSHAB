"""test_gfs_collector.py - Unit tests for the background GFS collector.

The collector defaults to a light 1.0°, 3-hourly grid and is overridden only via
CLI args (--res / --step-hours); it does NOT read config (unlike saveNETCDF.py).
These tests cover that pure configuration/derivation logic plus the grid-aware
helpers (step availability, GRIB/NetCDF filenames, URL building). Network
downloads and cfgrib conversion are not exercised here.

Each test calls ``configure()`` to set/restore the module globals, since import
order leaves them at the defaults.

Run:
    PYTHONPATH=src pytest tests/test_gfs_collector.py -v
"""

import pytest

import EarthSHAB.gfs_collector as gc


@pytest.fixture(autouse=True)
def _reset_defaults():
    """Restore the 1.0°/3h defaults after each test that reconfigures."""
    gc.configure()
    yield
    gc.configure()


# ── Defaults ────────────────────────────────────────────────────────────────────

def test_default_is_1deg_3hr():
    gc.configure()
    assert gc.DEFAULT_RES == 1.0 and gc.DEFAULT_STEP_HOURS == 3
    assert gc.RESOLUTION == "1p00"
    assert gc.STEP_HOURS == 3
    assert gc.FORECAST_HOURS == list(range(0, 384, 3))


def test_does_not_read_config():
    # The collector must be self-contained: importing it should not require the
    # config module (it was removed). Guard against a regression re-adding it.
    assert not hasattr(gc, "config_earth")


# ── configure() override + derivation ─────────────────────────────────────────────

@pytest.mark.parametrize("res,token", [(0.25, "0p25"), (0.5, "0p50"), (1.0, "1p00")])
def test_configure_resolution_token(res, token):
    gc.configure(res=res)
    assert gc.RESOLUTION == token
    assert token in gc.BASE_URL


@pytest.mark.parametrize("step", [1, 2, 3, 6])
def test_configure_step_hours_sets_forecast_hours(step):
    gc.configure(res=0.25, step_hours=step)
    assert gc.STEP_HOURS == step
    assert gc.FORECAST_HOURS == list(range(0, 384, step))


def test_configure_rejects_bad_res():
    with pytest.raises(ValueError):
        gc.configure(res=0.1)


def test_configure_rejects_nonpositive_step():
    with pytest.raises(ValueError):
        gc.configure(step_hours=0)


# ── CLI parsing ───────────────────────────────────────────────────────────────────

def test_parse_args_defaults():
    args = gc.parse_args([])
    assert args.res == 1.0 and args.step_hours == 3


def test_parse_args_overrides():
    args = gc.parse_args(["--res", "0.25", "--step-hours", "1"])
    assert args.res == 0.25 and args.step_hours == 1


def test_parse_args_rejects_invalid_res():
    with pytest.raises(SystemExit):
        gc.parse_args(["--res", "0.1"])


# ── Step availability (grid-aware) ────────────────────────────────────────────────

def test_is_step_available_3hourly_all_grids():
    gc.configure(res=1.0)
    assert gc.is_step_available(0) and gc.is_step_available(384)
    assert gc.is_step_available(3) and not gc.is_step_available(2)


def test_is_step_available_hourly_only_on_0p25_to_f120():
    gc.configure(res=0.25)
    assert gc.is_step_available(1) and gc.is_step_available(120)
    assert not gc.is_step_available(121)  # hourly stops at f120
    assert gc.is_step_available(123)      # but 3-hourly continues


def test_is_step_available_hourly_absent_on_coarse_grid():
    gc.configure(res=1.0)
    assert not gc.is_step_available(1)
    assert not gc.is_step_available(2)


def test_is_step_available_out_of_range():
    assert not gc.is_step_available(-1)
    assert not gc.is_step_available(385)


def test_expected_forecasts_subset_of_forecast_hours():
    gc.configure(res=0.5, step_hours=3)
    expected = gc.expected_forecasts()
    assert expected == [f for f in gc.FORECAST_HOURS if gc.is_step_available(f)]
    assert all(f % 3 == 0 for f in expected)


# ── Filename / path helpers encode res + step ─────────────────────────────────────

def test_grib_filename_encodes_resolution():
    gc.configure(res=1.0)
    assert gc.grib_filename("20260605", "06", 3) == "gfs_20260605_06z_1p00_f003.grib2"
    gc.configure(res=0.25)
    assert gc.grib_filename("20260605", "06", 12) == "gfs_20260605_06z_0p25_f012.grib2"


def test_output_nc_path_encodes_res_and_step():
    gc.configure(res=0.25, step_hours=3)
    assert gc.output_nc_path("20260605", "06").endswith("gfs_20260605_06z_0p25_3h.nc")
    gc.configure(res=1.0, step_hours=6)
    assert gc.output_nc_path("20260605", "06").endswith("gfs_20260605_06z_1p00_6h.nc")


# ── URL building ──────────────────────────────────────────────────────────────────

def test_base_url_format():
    gc.configure(res=0.25)
    url = gc.BASE_URL.format(date="20260605", hour="06", forecast="012")
    assert url == (
        "https://noaa-gfs-bdp-pds.s3.amazonaws.com/"
        "gfs.20260605/06/atmos/gfs.t06z.pgrb2.0p25.f012"
    )
