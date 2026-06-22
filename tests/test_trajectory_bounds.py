"""test_trajectory_bounds.py - Tests for the balloon_trajectory time-bounds guard.

A reforecast (config_earth.balloon_trajectory) re-anchors the APRS track to
simulation['start_time'] and steps by each inter-sample dt, so a trajectory
recorded on a *different* day/year than the downloaded forecast does NOT crash —
it silently applies the wrong winds to the truth altitude profile (e.g. a 2022
flight against a 2026 forecast).

simulate._assert_trajectory_within_forecast guards against this by comparing the
trajectory's own recorded timestamps (df['time'], local → UTC) against the
forecast's [model_start, model_end] coverage, and fatally exiting on a mismatch.

Run:
    pytest tests/test_trajectory_bounds.py -v
"""

from datetime import datetime
import textwrap

import pandas as pd
import pytest

from EarthSHAB.simulate import _assert_trajectory_within_forecast, _load_aprs


# ── Helpers ──────────────────────────────────────────────────────────────────

GMT = 7  # config_earth / simulate uses UTC-7 (MST); the loader stores local time.


def _df_local(first_utc, last_utc):
    """Build a trajectory DataFrame whose 'time' column is LOCAL (MST), given the
    desired first/last timestamps in UTC — mirroring _load_aprs output."""
    offset = pd.Timedelta(hours=GMT)
    times_local = [
        pd.Timestamp(first_utc) - offset,
        pd.Timestamp(last_utc) - offset,
    ]
    return pd.DataFrame({"time": pd.to_datetime(times_local)})


# Forecast time window used across tests: 2022-08-22 00:00 → 2022-08-23 23:00 UTC
MODEL_START = datetime(2022, 8, 22, 0, 0, 0)
MODEL_END = datetime(2022, 8, 23, 23, 0, 0)


# ── Tests: in-bounds (no exit) ────────────────────────────────────────────────

class TestTrajectoryWithinForecast:
    """The guard must stay silent when the recorded flight fits inside the forecast."""

    def test_well_within_bounds_passes(self):
        df = _df_local("2022-08-22 12:00:00", "2022-08-22 18:00:00")
        _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=GMT)

    def test_exact_end_boundary_passes(self):
        """traj_end == model_end is allowed (strict > comparison)."""
        df = _df_local("2022-08-22 12:00:00", "2022-08-23 23:00:00")
        _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=GMT)

    def test_exact_start_boundary_passes(self):
        """traj_start == model_start is allowed (strict < comparison)."""
        df = _df_local("2022-08-22 00:00:00", "2022-08-22 06:00:00")
        _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=GMT)


# ── Tests: out-of-bounds (fatal exit) ─────────────────────────────────────────

class TestTrajectoryOutOfForecast:
    """The guard must fatally exit when the recorded flight leaves the forecast."""

    def test_different_year_exits(self):
        """The reported scenario: a 2022 flight against a 2026 forecast."""
        df = _df_local("2022-08-22 14:00:00", "2022-08-23 03:00:00")
        model_start = datetime(2026, 6, 10, 12, 0, 0)
        model_end = datetime(2026, 6, 11, 12, 0, 0)
        with pytest.raises(SystemExit) as exc:
            _assert_trajectory_within_forecast(df, model_start, model_end, gmt_offset_hours=GMT)
        assert exc.value.code != 0

    def test_end_past_forecast_exits(self):
        df = _df_local("2022-08-23 12:00:00", "2022-08-24 06:00:00")
        with pytest.raises(SystemExit) as exc:
            _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=GMT)
        assert exc.value.code != 0

    def test_start_before_forecast_exits(self):
        df = _df_local("2022-08-21 22:00:00", "2022-08-22 06:00:00")
        with pytest.raises(SystemExit) as exc:
            _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=GMT)
        assert exc.value.code != 0

    def test_end_one_hour_past_exits(self):
        """Just over the end boundary must still trip the guard."""
        df = _df_local("2022-08-23 20:00:00", "2022-08-24 00:00:00")  # 00:00 > 23:00
        with pytest.raises(SystemExit):
            _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=GMT)

    def test_message_is_emitted(self, capsys):
        df = _df_local("2022-08-21 00:00:00", "2022-08-21 06:00:00")
        with pytest.raises(SystemExit):
            _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=GMT)
        out = capsys.readouterr().out
        assert "FATAL" in out
        assert "balloon_trajectory" in out


# ── Tests: timezone handling ──────────────────────────────────────────────────

class TestTimezoneOffset:
    """The local→UTC shift (gmt_offset_hours) must be applied before comparing."""

    def test_offset_pushes_inbounds_track_out(self):
        """A track whose LOCAL times sit inside the window is out once shifted to UTC.

        Local span 19:00 → 22:00 (MST) → UTC 02:00 → 05:00 next day, past a
        forecast that ends at 23:00 the same day.
        """
        # Build local times directly (no UTC pre-shift) to isolate the offset.
        df = pd.DataFrame({"time": pd.to_datetime([
            "2022-08-23 19:00:00", "2022-08-23 22:00:00",
        ])})
        with pytest.raises(SystemExit):
            _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=GMT)

    def test_zero_offset_uses_times_as_utc(self):
        """With gmt_offset_hours=0, df['time'] is treated as UTC directly."""
        df = pd.DataFrame({"time": pd.to_datetime([
            "2022-08-22 12:00:00", "2022-08-22 18:00:00",
        ])})
        _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=0)


# ── Tests: integration with _load_aprs (real parse + tz shift) ────────────────

class TestWithLoadedAprs:
    """End-to-end with a real APRS CSV parsed by _load_aprs (which shifts UTC-7)."""

    def _write_csv(self, tmp_path):
        # Flight on 2022-08-22, 17:46 → 19:46 UTC (stored as 10:46 → 12:46 MST).
        csv = textwrap.dedent("""\
            time,lat,lng,altitude,comment
            2022-08-22 17:46:00,34.640,-106.830,1500.0,
            2022-08-22 18:46:00,34.700,-106.780,16000.0,
            2022-08-22 19:46:00,34.750,-106.750,15000.0,
        """)
        p = tmp_path / "flight.csv"
        p.write_text(csv)
        return str(p)

    def test_loaded_track_within_forecast_passes(self, tmp_path):
        df, _ = _load_aprs(self._write_csv(tmp_path))
        _assert_trajectory_within_forecast(df, MODEL_START, MODEL_END, gmt_offset_hours=GMT)

    def test_loaded_track_wrong_year_exits(self, tmp_path):
        df, _ = _load_aprs(self._write_csv(tmp_path))
        model_start = datetime(2026, 6, 10, 0, 0, 0)
        model_end = datetime(2026, 6, 11, 0, 0, 0)
        with pytest.raises(SystemExit):
            _assert_trajectory_within_forecast(df, model_start, model_end, gmt_offset_hours=GMT)
