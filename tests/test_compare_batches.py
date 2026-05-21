"""test_compare_batches.py - Unit tests for compare_batches pure logic.

Covers threshold logic, percentage delta with near-zero baselines, asymmetry
classification, win/loss counting, derived metric computation, and batch_id
resolution. Rendering and plot generation are not covered here (visual review).

Run:
    PYTHONPATH=src pytest tests/test_compare_batches.py -v
"""

import math
import os
from pathlib import Path

import pandas as pd
import pytest

from evaluation.compare_batches import (
    PCT_FLOOR,
    METRICS,
    add_derived_metrics,
    classify_change,
    delta_and_pct,
    find_asymmetries,
    find_missing_metric_columns,
    nan_mean,
    resolve_batch_id,
    win_loss_counts,
)


# ── delta_and_pct ─────────────────────────────────────────────────────────────


class TestDeltaAndPct:
    def test_basic(self):
        d, p = delta_and_pct(100.0, 90.0, abs_floor=1.0)
        assert d == pytest.approx(-10.0)
        assert p == pytest.approx(-10.0)

    def test_positive_delta(self):
        d, p = delta_and_pct(100.0, 120.0, abs_floor=1.0)
        assert d == pytest.approx(20.0)
        assert p == pytest.approx(20.0)

    def test_pct_is_nan_when_a_below_floor(self):
        # |A| < abs_floor → pct undefined but delta still computed
        d, p = delta_and_pct(0.5, 5.0, abs_floor=1.0)
        assert d == pytest.approx(4.5)
        assert math.isnan(p)

    def test_pct_uses_abs_of_a(self):
        # A=-100, B=-50: delta=+50, pct=+50%
        d, p = delta_and_pct(-100.0, -50.0, abs_floor=1.0)
        assert d == pytest.approx(50.0)
        assert p == pytest.approx(50.0)

    def test_nan_inputs_return_nan(self):
        d, p = delta_and_pct(math.nan, 5.0, abs_floor=1.0)
        assert math.isnan(d) and math.isnan(p)
        d, p = delta_and_pct(5.0, math.nan, abs_floor=1.0)
        assert math.isnan(d) and math.isnan(p)

    def test_none_inputs_return_nan(self):
        d, p = delta_and_pct(None, 1.0, abs_floor=1.0)
        assert math.isnan(d) and math.isnan(p)

    def test_non_numeric_string_returns_nan(self):
        d, p = delta_and_pct("bogus", 5.0, abs_floor=1.0)
        assert math.isnan(d) and math.isnan(p)


# ── classify_change ───────────────────────────────────────────────────────────


class TestClassifyChange:
    def test_improvement_passes_both_floors(self):
        # delta=-50 (passes abs floor 1.0), pct=-25% (passes pct floor 5%)
        assert classify_change(-50.0, -25.0, abs_floor=1.0) == "improved"

    def test_regression_passes_both_floors(self):
        assert classify_change(30.0, 15.0, abs_floor=1.0) == "worsened"

    def test_unchanged_when_below_abs_floor(self):
        # delta=0.5 fails abs_floor=1.0 even if pct is huge
        assert classify_change(0.5, 50.0, abs_floor=1.0) == "unchanged"

    def test_unchanged_when_below_pct_floor(self):
        # delta=20 passes abs_floor=1.0 but pct=2% fails PCT_FLOOR=5
        assert classify_change(20.0, 2.0, abs_floor=1.0) == "unchanged"

    def test_significant_when_pct_is_nan_and_abs_passes(self):
        # When baseline was too small for pct, abs alone decides
        assert classify_change(-10.0, math.nan, abs_floor=1.0) == "improved"
        assert classify_change(10.0, math.nan, abs_floor=1.0) == "worsened"

    def test_nan_delta_is_unchanged(self):
        assert classify_change(math.nan, math.nan, abs_floor=1.0) == "unchanged"

    def test_exactly_at_floor_is_unchanged(self):
        # |delta| == abs_floor → unchanged (strict greater-than)
        assert classify_change(1.0, 10.0, abs_floor=1.0) == "unchanged"


# ── win_loss_counts ───────────────────────────────────────────────────────────


class TestWinLossCounts:
    def test_mixed(self):
        # (delta, pct) pairs: 2 improved, 1 unchanged, 1 worsened
        pairs = [(-50.0, -25.0), (30.0, 15.0), (0.5, 50.0), (-100.0, -50.0)]
        c = win_loss_counts(pairs, abs_floor=1.0)
        assert c == {"improved": 2, "unchanged": 1, "worsened": 1}

    def test_empty(self):
        c = win_loss_counts([], abs_floor=1.0)
        assert c == {"improved": 0, "unchanged": 0, "worsened": 0}


# ── nan_mean ──────────────────────────────────────────────────────────────────


class TestNanMean:
    def test_skips_nan(self):
        assert nan_mean([1.0, 2.0, math.nan, 3.0]) == pytest.approx(2.0)

    def test_all_nan_returns_nan(self):
        assert math.isnan(nan_mean([math.nan, math.nan]))

    def test_empty_returns_nan(self):
        assert math.isnan(nan_mean([]))

    def test_skips_none_and_non_numeric(self):
        assert nan_mean([1.0, None, "bad", 3.0]) == pytest.approx(2.0)


# ── find_asymmetries ──────────────────────────────────────────────────────────


def _df(rows):
    return pd.DataFrame(rows)


class TestFindAsymmetries:
    def test_intersection_only(self):
        a = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "ok"}])
        b = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "ok"}])
        out = find_asymmetries(a, b)
        assert out["intersection_ok"] == {("L1", "GFS")}
        assert out["only_in_a"] == set()
        assert out["only_in_b"] == set()
        assert out["failed_in_a"] == set()
        assert out["failed_in_b"] == set()

    def test_only_in_a(self):
        a = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "ok"},
                 {"launch_id": "L2", "forecast_type": "GFS", "status": "ok"}])
        b = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "ok"}])
        out = find_asymmetries(a, b)
        assert out["intersection_ok"] == {("L1", "GFS")}
        assert out["only_in_a"] == {("L2", "GFS")}

    def test_failed_in_b(self):
        a = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "ok"}])
        b = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "error: x"}])
        out = find_asymmetries(a, b)
        assert out["intersection_ok"] == set()
        assert out["failed_in_b"] == {("L1", "GFS")}

    def test_failed_in_both(self):
        a = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "fail"}])
        b = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "fail"}])
        out = find_asymmetries(a, b)
        assert out["failed_in_both"] == {("L1", "GFS")}
        assert out["intersection_ok"] == set()

    def test_forecast_type_asymmetry_within_launch(self):
        # A has GFS+ERA5 for L1, B has only GFS — ERA5 goes to only_in_a
        a = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "ok"},
                 {"launch_id": "L1", "forecast_type": "ERA5", "status": "ok"}])
        b = _df([{"launch_id": "L1", "forecast_type": "GFS", "status": "ok"}])
        out = find_asymmetries(a, b)
        assert out["intersection_ok"] == {("L1", "GFS")}
        assert out["only_in_a"] == {("L1", "ERA5")}


# ── find_missing_metric_columns ───────────────────────────────────────────────


class TestFindMissingMetricColumns:
    def test_missing_in_a(self):
        a = _df([{"launch_id": "L1"}])
        b = _df([{"launch_id": "L1", "landing_distance_km": 100.0}])
        out = find_missing_metric_columns(a, b, ["landing_distance_km"])
        assert out == {"landing_distance_km": "A"}

    def test_missing_in_b(self):
        a = _df([{"launch_id": "L1", "landing_distance_km": 100.0}])
        b = _df([{"launch_id": "L1"}])
        out = find_missing_metric_columns(a, b, ["landing_distance_km"])
        assert out == {"landing_distance_km": "B"}

    def test_missing_in_both(self):
        a = _df([{"launch_id": "L1"}])
        b = _df([{"launch_id": "L1"}])
        out = find_missing_metric_columns(a, b, ["landing_distance_km"])
        assert out == {"landing_distance_km": "both"}

    def test_present_in_both_omitted_from_result(self):
        a = _df([{"launch_id": "L1", "landing_distance_km": 100.0}])
        b = _df([{"launch_id": "L1", "landing_distance_km": 90.0}])
        out = find_missing_metric_columns(a, b, ["landing_distance_km"])
        assert out == {}


# ── add_derived_metrics ───────────────────────────────────────────────────────


class TestAddDerivedMetrics:
    def test_adds_abs_diff_columns(self):
        df = _df([{"sim_float_alt_mean_m": 19000, "truth_float_alt_mean_m": 20000,
                   "sim_ascent_rate_mean_ms": 2.0, "truth_ascent_rate_mean_ms": 2.5,
                   "sim_descent_rate_mean_ms": -2.5, "truth_descent_rate_mean_ms": -3.0,
                   "landing_time_diff_min": -7.5}])
        out = add_derived_metrics(df.copy())
        assert out["float_alt_abs_err_m"].iloc[0] == pytest.approx(1000.0)
        assert out["ascent_rate_abs_err_ms"].iloc[0] == pytest.approx(0.5)
        assert out["descent_rate_abs_err_ms"].iloc[0] == pytest.approx(0.5)
        # use_abs metrics get abs() applied in-place
        assert out["landing_time_diff_min"].iloc[0] == pytest.approx(7.5)

    def test_missing_sim_or_truth_yields_nan(self):
        df = _df([{"sim_float_alt_mean_m": 19000}])  # no truth column
        out = add_derived_metrics(df.copy())
        assert math.isnan(out["float_alt_abs_err_m"].iloc[0])


# ── resolve_batch_id ──────────────────────────────────────────────────────────


class TestResolveBatchId:
    def test_exact_match(self, tmp_path):
        (tmp_path / "2026-01-01T0000_abc1234").mkdir()
        assert resolve_batch_id("2026-01-01T0000_abc1234", str(tmp_path)) \
            == "2026-01-01T0000_abc1234"

    def test_short_hash_match(self, tmp_path):
        (tmp_path / "2026-01-01T0000_abc1234").mkdir()
        assert resolve_batch_id("abc1234", str(tmp_path)) == "2026-01-01T0000_abc1234"

    def test_short_hash_ambiguous_raises(self, tmp_path):
        (tmp_path / "2026-01-01T0000_abc1234").mkdir()
        (tmp_path / "2026-02-02T0000_abc1234").mkdir()
        with pytest.raises(ValueError, match="ambiguous"):
            resolve_batch_id("abc1234", str(tmp_path))

    def test_no_match_raises(self, tmp_path):
        (tmp_path / "2026-01-01T0000_abc1234").mkdir()
        with pytest.raises(FileNotFoundError):
            resolve_batch_id("nonexistent", str(tmp_path))

    def test_missing_dir_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            resolve_batch_id("foo", str(tmp_path / "does_not_exist"))


# ── Sanity check: PCT_FLOOR and METRICS module-level constants ───────────────


class TestModuleConstants:
    def test_pct_floor_is_positive(self):
        assert PCT_FLOOR > 0

    def test_all_metrics_have_required_keys(self):
        for m in METRICS:
            assert "key" in m and "label" in m and "abs_floor" in m and "fmt" in m
            assert m["abs_floor"] > 0
