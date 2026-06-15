"""test_savenetcdf_archive.py - Unit tests for the GFS archive downloader and
the live->archive auto-switch in saveNETCDF.

Covers pure logic only: .idx parsing, byte-range selection, URL building, cycle
validation, the retention check + dispatch routing, the past-date prompt, and the
``step_hours`` temporal-resolution wiring in both downloaders. Network downloads
and cfgrib conversion are not exercised here.

Run:
    PYTHONPATH=src pytest tests/test_savenetcdf_archive.py -v
"""

import builtins
import datetime
import types

import pytest
import xarray as xr

import EarthSHAB.config as config_earth
import EarthSHAB.saveNETCDF as live
import EarthSHAB.saveNETCDF_archive as arch


# Sample .idx body: 4 vars at 500 mb, one UGRD at 0.4 mb, one surface field.
SAMPLE_IDX = (
    "1:0:d=2022082212:UGRD:500 mb:anl:\n"
    "2:100:d=2022082212:VGRD:500 mb:anl:\n"
    "3:250:d=2022082212:HGT:500 mb:anl:\n"
    "4:500:d=2022082212:TMP:500 mb:anl:\n"
    "5:900:d=2022082212:UGRD:0.4 mb:anl:\n"
    "6:1000:d=2022082212:PRMSL:mean sea level:anl:\n"
)


# ── .idx parsing / level / range selection ──────────────────────────────────────

def test_parse_idx_basic():
    rows = arch.parse_idx(SAMPLE_IDX)
    assert len(rows) == 6
    assert rows[0] == (1, 0, "UGRD", "500 mb")
    assert rows[4] == (5, 900, "UGRD", "0.4 mb")


def test_parse_idx_skips_garbage():
    assert arch.parse_idx("not:a:number\n\n:::\n7:x:d=y:UGRD:500 mb:anl:") == []


def test_level_mb():
    assert arch._level_mb("500 mb") == 500.0
    assert arch._level_mb("0.4 mb") == 0.4
    assert arch._level_mb("surface") is None
    assert arch._level_mb("2 m above ground") is None


def test_select_ranges_filters_vars_and_levels():
    ranges = arch.select_ranges(arch.parse_idx(SAMPLE_IDX))
    assert {r[2] for r in ranges} == {"UGRD", "VGRD", "HGT", "TMP"}  # PRMSL excluded
    assert len(ranges) == 5  # surface level excluded


def test_select_ranges_end_offsets_and_open_last():
    ranges = arch.select_ranges(arch.parse_idx(SAMPLE_IDX))
    assert ranges[0][:2] == (0, 99)
    tmp = [r for r in ranges if r[2] == "TMP"][0]
    assert tmp[:2] == (500, 899)


def test_select_ranges_open_ended_final_message():
    idx = "1:0:d=x:TMP:850 mb:anl:\n2:300:d=x:UGRD:850 mb:anl:\n"
    ranges = arch.select_ranges(arch.parse_idx(idx))
    assert ranges[-1][0] == 300 and ranges[-1][1] is None


def test_select_ranges_level_subset():
    ranges = arch.select_ranges(arch.parse_idx(SAMPLE_IDX), want_levels_mb=[0.4])
    assert len(ranges) == 1 and ranges[0][2] == "UGRD" and ranges[0][3] == 0.4


# ── URL construction ────────────────────────────────────────────────────────────

def test_cycle_file_url():
    assert arch.cycle_file_url("20220822", 12, 0) == (
        "https://noaa-gfs-bdp-pds.s3.amazonaws.com/"
        "gfs.20220822/12/atmos/gfs.t12z.pgrb2.0p25.f000"
    )
    assert arch.cycle_file_url("20220822", 6, 27).endswith("gfs.t06z.pgrb2.0p25.f027")


# ── Output filename ───────────────────────────────────────────────────────────

def test_archive_output_path_inserts_marker():
    assert arch._archive_output_path("f/gfs_0p25_3h_20220822_12.nc") == \
        "f/gfs_0p25_3h_20220822_12_archive.nc"


def test_archive_output_path_idempotent():
    p = "f/gfs_0p25_3h_20220822_12_archive.nc"
    assert arch._archive_output_path(p) == p


# ── Cycle validation ────────────────────────────────────────────────────────────

def test_validate_accepts_valid_cycle():
    d = datetime.datetime(2022, 8, 22, 12)
    assert arch._validate_archive_nc_start(d) == d


def test_validate_rejects_pre_archive():
    with pytest.raises(ValueError):
        arch._validate_archive_nc_start(datetime.datetime(2020, 1, 1, 12))


def test_validate_rejects_noncycle_hour():
    with pytest.raises(ValueError):
        arch._validate_archive_nc_start(datetime.datetime(2022, 8, 22, 13))


def test_validate_rejects_future():
    future = (datetime.datetime.utcnow() + datetime.timedelta(days=2)).replace(
        hour=0, minute=0, second=0, microsecond=0
    )
    with pytest.raises(ValueError):
        arch._validate_archive_nc_start(future)


# ── Retention check + dispatch routing ──────────────────────────────────────────

def test_is_past_retention():
    assert live._is_past_retention(datetime.datetime(2022, 8, 22, 12)) is True
    assert live._is_past_retention(datetime.datetime.utcnow() - datetime.timedelta(days=2)) is False


def test_main_routes_to_archive_when_past(monkeypatch):
    monkeypatch.setitem(config_earth.netcdf_gfs, "nc_start", datetime.datetime(2022, 8, 22, 12))
    calls = []
    monkeypatch.setattr(live, "download_gfs_grib_to_netcdf", lambda: calls.append("live"))
    monkeypatch.setattr(arch, "download_gfs_archive_to_netcdf",
                        lambda: calls.append("archive") or "f.nc")
    monkeypatch.setattr(live, "_prompt_proceed_archive", lambda nc: True)
    live.main()
    assert calls == ["archive"]


def test_main_routes_to_live_when_recent(monkeypatch):
    recent = (datetime.datetime.utcnow() - datetime.timedelta(days=2)).replace(
        hour=12, minute=0, second=0, microsecond=0
    )
    monkeypatch.setitem(config_earth.netcdf_gfs, "nc_start", recent)
    calls = []
    monkeypatch.setattr(live, "download_gfs_grib_to_netcdf", lambda: calls.append("live"))
    monkeypatch.setattr(arch, "download_gfs_archive_to_netcdf", lambda: calls.append("archive"))
    live.main()
    assert calls == ["live"]


def test_main_aborts_when_declined(monkeypatch):
    monkeypatch.setitem(config_earth.netcdf_gfs, "nc_start", datetime.datetime(2022, 8, 22, 12))
    monkeypatch.setattr(live, "_prompt_proceed_archive", lambda nc: False)
    monkeypatch.setattr(arch, "download_gfs_archive_to_netcdf",
                        lambda: pytest.fail("must not download when declined"))
    with pytest.raises(SystemExit):
        live.main()


# ── Past-date prompt ────────────────────────────────────────────────────────────

def _patch_tty(monkeypatch, is_tty):
    monkeypatch.setattr(live.sys, "stdin", types.SimpleNamespace(isatty=lambda: is_tty))


def test_prompt_non_interactive_proceeds(monkeypatch):
    _patch_tty(monkeypatch, False)
    assert live._prompt_proceed_archive(datetime.datetime(2022, 8, 22, 12)) is True


@pytest.mark.parametrize("answer,expected", [("y", True), ("yes", True), ("", True),
                                             ("n", False), ("no", False), ("x", False)])
def test_prompt_interactive(monkeypatch, answer, expected):
    _patch_tty(monkeypatch, True)
    monkeypatch.setattr(builtins, "input", lambda *a, **k: answer)
    assert live._prompt_proceed_archive(datetime.datetime(2022, 8, 22, 12)) is expected


# ── step_hours temporal-resolution wiring ───────────────────────────────────────

def test_archive_step_hours_controls_forecast_hours(monkeypatch, tmp_path):
    captured = {}

    def fake_fetch(date_str, cycle_hour, forecast_hours, **kw):
        captured.update(fhrs=list(forecast_hours), date=date_str, cyc=cycle_hour)
        return xr.Dataset({"u": (("valid_time",), [0.0])}, coords={"valid_time": [0]})

    monkeypatch.setattr(arch, "fetch_cycle", fake_fetch)
    monkeypatch.setitem(config_earth.netcdf_gfs, "nc_start", datetime.datetime(2022, 8, 22, 12))
    monkeypatch.setitem(config_earth.netcdf_gfs, "download_days", 1)
    monkeypatch.setitem(config_earth.netcdf_gfs, "nc_file", str(tmp_path / "out.nc"))

    monkeypatch.setitem(config_earth.netcdf_gfs, "step_hours", 1)
    arch.download_gfs_archive_to_netcdf()
    assert captured["fhrs"] == list(range(0, 25, 1))
    assert captured["date"] == "20220822" and captured["cyc"] == 12

    monkeypatch.setitem(config_earth.netcdf_gfs, "step_hours", 3)
    arch.download_gfs_archive_to_netcdf()
    assert captured["fhrs"] == list(range(0, 25, 3))


def test_live_step_hours_controls_forecast_hours(monkeypatch):
    # A valid recent cycle (within retention, past the upload lag).
    base = datetime.datetime.utcnow() - datetime.timedelta(hours=18)
    nc_start = base.replace(hour=(base.hour // 6) * 6, minute=0, second=0, microsecond=0)
    monkeypatch.setitem(config_earth.netcdf_gfs, "nc_start", nc_start)
    monkeypatch.setitem(config_earth.netcdf_gfs, "download_days", 1)
    monkeypatch.setitem(config_earth.netcdf_gfs, "step_hours", 1)

    seen = []
    monkeypatch.setattr(live, "_build_url",
                        lambda date_str, cyc, fhr, *a, **k: seen.append(fhr) or "http://x")
    monkeypatch.setattr(live, "_download_grib", lambda url, dest, retries=3: False)
    monkeypatch.setattr(live.time, "sleep", lambda *a, **k: None)

    # All downloads "fail" -> nothing assembled -> RuntimeError, but URLs were built.
    with pytest.raises(RuntimeError):
        live.download_gfs_grib_to_netcdf()
    assert seen == list(range(0, 25, 1))
