"""test_evaluate_phases.py - Unit tests for BalloonEvaluator._detect_phases.

Covers two bug fixes made to the phase-detection logic:

  1. NaN altitude handling (nanmax/nanargmax)
     A single corrupted APRS packet produces a NaN altitude.  Before the fix,
     el.max() returned NaN, collapsing the entire float-detection pipeline.

  2. Float-block trimming (front/back)
     The largest "candidate" block was being reported as float even when it
     included the ascending approach into float (rolling mean > 0.3 m/s) or
     the rapid descent exit from float (rolling mean < -0.5 m/s).  The fix
     trims those shoulders off before setting float_mask.

Run:
    pytest tests/test_evaluate_phases.py -v
"""

import numpy as np
import pytest

# ── Import target ────────────────────────────────────────────────────────────

from evaluation.evaluate import BalloonEvaluator

_phases = BalloonEvaluator._detect_phases   # shorthand for brevity


# ── Synthetic flight-profile helpers ─────────────────────────────────────────

def _clean_profile(n_asc=1000, n_flt=1000, n_des=1000, min_alt=500.0,
                   peak_alt=20000.0, v_asc=2.0, v_flt=0.05, v_des=-3.0):
    """Return (el, v, min_alt) for a clean ascent → float → descent profile."""
    el = np.concatenate([
        np.linspace(min_alt, peak_alt, n_asc),
        np.full(n_flt, peak_alt),
        np.linspace(peak_alt, min_alt, n_des),
    ])
    v = np.concatenate([
        np.full(n_asc, v_asc),
        np.full(n_flt, v_flt),
        np.full(n_des, v_des),
    ])
    return el, v, min_alt


# ── Tests: normal detection ───────────────────────────────────────────────────

class TestDetectPhasesNormal:
    def test_float_detected(self):
        el, v, min_alt = _clean_profile()
        _, flt_m, _, _, _ = _phases(el, v, min_alt)
        assert flt_m.any(), "Float should be detected in a clean profile"

    def test_float_altitude_accurate(self):
        el, v, min_alt = _clean_profile(peak_alt=20000.0)
        _, flt_m, _, _, _ = _phases(el, v, min_alt)
        assert el[flt_m].mean() > 19000.0, "Float mean altitude should be near peak"

    def test_ascent_detected(self):
        el, v, min_alt = _clean_profile()
        asc_m, _, _, _, _ = _phases(el, v, min_alt)
        assert asc_m.any(), "Ascent phase should be detected"

    def test_descent_detected(self):
        el, v, min_alt = _clean_profile()
        _, _, des_m, _, _ = _phases(el, v, min_alt)
        assert des_m.any(), "Descent phase should be detected"

    def test_phases_non_overlapping(self):
        el, v, min_alt = _clean_profile()
        asc_m, flt_m, des_m, _, _ = _phases(el, v, min_alt)
        assert not (asc_m & flt_m).any(), "Ascent and float must not overlap"
        assert not (asc_m & des_m).any(), "Ascent and descent must not overlap"
        assert not (flt_m & des_m).any(), "Float and descent must not overlap"

    def test_no_float_for_sharp_turnaround(self):
        """Balloon that ascends then immediately descends should not falsely report float."""
        n = 2000
        el = np.concatenate([np.linspace(500, 20000, 1000), np.linspace(20000, 500, 1000)])
        v  = np.concatenate([np.full(1000, 2.5), np.full(1000, -5.0)])
        _, flt_m, _, _, _ = _phases(el, v, 500.0)
        assert not flt_m.any(), "No float should be detected for an immediate turnaround"


# ── Tests: NaN altitude handling ──────────────────────────────────────────────

class TestDetectPhasesNanAltitude:
    """
    Bug: el.max() / np.argmax(el) returned NaN when any element was NaN,
    cascading to NaN alt_threshold → high_indices empty → broken detection.
    Fix: np.nanmax / np.nanargmax, plus NaN values naturally excluded from
    `el >= alt_threshold` comparisons.
    """

    def test_nan_during_descent_still_detects_float(self):
        """Single NaN altitude during descent must not break float detection."""
        el, v, min_alt = _clean_profile()
        el[2500] = np.nan          # corrupt one reading during descent
        _, flt_m, _, _, _ = _phases(el, v, min_alt)
        assert flt_m.any(), "Float should be detected even with a NaN during descent"

    def test_nan_during_descent_float_altitude_accurate(self):
        el, v, min_alt = _clean_profile(peak_alt=20000.0)
        el[2500] = np.nan
        _, flt_m, _, _, _ = _phases(el, v, min_alt)
        assert np.nanmean(el[flt_m]) > 19000.0, (
            "Float mean altitude should still be near 20000 m with one NaN during descent"
        )

    def test_nan_in_float_region_still_detects_float(self):
        """NaN in the float region itself (corrupted packet at peak) must not crash."""
        el, v, min_alt = _clean_profile()
        el[1500] = np.nan          # NaN right in the middle of float
        _, flt_m, _, _, _ = _phases(el, v, min_alt)
        assert flt_m.any(), "Float should still be detected with a NaN inside the float region"

    def test_nan_excluded_from_high_indices(self):
        """NaN altitudes should not appear in the detected float region."""
        el, v, min_alt = _clean_profile()
        nan_idx = 2500
        el[nan_idx] = np.nan
        _, flt_m, _, _, _ = _phases(el, v, min_alt)
        assert not flt_m[nan_idx], "The NaN-altitude step must not be part of the float mask"

    def test_multiple_nans_in_descent(self):
        """Several NaN altitudes scattered through descent should not break detection."""
        el, v, min_alt = _clean_profile()
        for i in [2100, 2300, 2500, 2700]:
            el[i] = np.nan
        _, flt_m, _, _, _ = _phases(el, v, min_alt)
        assert flt_m.any(), "Float detected with multiple NaNs during descent"
        assert np.nanmean(el[flt_m]) > 19000.0


# ── Tests: float-block trimming ───────────────────────────────────────────────

class TestDetectPhasesFloatTrimming:
    """
    Bug: the largest candidate block in the high-altitude zone was accepted
    as-is, so the 'float' period included:
      • the ascending approach into float (still climbing at ~0.5-1 m/s)
      • the rapid descent exit from float (falling at >0.5 m/s)
    Fix: trim the front of the block until rolling-mean velocity < v_float*0.3,
    and trim the back until rolling-mean velocity > -v_float*0.5.

    Profile
    -------
    ascent        1000 steps   v=2.0 m/s
    approach       200 steps   v=0.9 m/s   ← above front-trim threshold (0.3)
    settled float 2000 steps   v=0.02 m/s  ← true float
    rapid exit     200 steps   v=-2.0 m/s  ← below back-trim threshold (-0.5)
    descent       1000 steps   v=-3.0 m/s
    """

    @pytest.fixture(autouse=True)
    def build_trim_profile(self):
        n_asc, n_app, n_flt, n_exit, n_des = 1000, 200, 2000, 200, 1000
        el = np.concatenate([
            np.linspace(500, 18000, n_asc),
            np.linspace(18000, 20000, n_app),
            np.full(n_flt, 20000.0),
            np.linspace(20000, 18000, n_exit),
            np.linspace(18000, 500, n_des),
        ])
        v = np.concatenate([
            np.full(n_asc, 2.0),
            np.full(n_app, 0.9),    # ascending approach — above trim threshold
            np.full(n_flt, 0.02),   # settled float
            np.full(n_exit, -2.0),  # rapid exit — below trim threshold
            np.full(n_des, -3.0),
        ])
        self.el = el
        self.v  = v
        self.min_alt = 500.0
        self.approach_mid   = n_asc + n_app // 2          # midpoint of approach
        self.settled_float_mid = n_asc + n_app + n_flt // 2  # midpoint of true float
        self.exit_mid       = n_asc + n_app + n_flt + n_exit // 2  # midpoint of exit

    def test_float_is_detected(self):
        _, flt_m, _, _, _ = _phases(self.el, self.v, self.min_alt)
        assert flt_m.any()

    def test_front_trim_excludes_ascending_approach(self):
        """The ascending-approach region (v=0.9, above 0.3 threshold) must not be float."""
        _, flt_m, _, _, _ = _phases(self.el, self.v, self.min_alt)
        assert not flt_m[self.approach_mid], (
            f"Ascending approach at step {self.approach_mid} should not be in float_mask"
        )

    def test_settled_float_center_is_included(self):
        """The settled float center (v=0.02) must be inside the float_mask."""
        _, flt_m, _, _, _ = _phases(self.el, self.v, self.min_alt)
        assert flt_m[self.settled_float_mid], (
            f"Settled float at step {self.settled_float_mid} should be in float_mask"
        )

    def test_back_trim_excludes_rapid_descent_exit(self):
        """The rapid-exit region (v=-2.0, below -0.5 threshold) must not be float."""
        _, flt_m, _, _, _ = _phases(self.el, self.v, self.min_alt)
        assert not flt_m[self.exit_mid], (
            f"Rapid descent exit at step {self.exit_mid} should not be in float_mask"
        )

    def test_float_altitude_not_inflated_by_approach(self):
        """Float mean altitude should be accurate even when an ascending approach is present."""
        _, flt_m, _, _, _ = _phases(self.el, self.v, self.min_alt)
        assert self.el[flt_m].mean() > 19000.0, (
            "Float mean altitude should reflect settled float near 20000 m, "
            "not be dragged down by the ascending approach"
        )

    def test_trim_with_nan_altitude_in_descent(self):
        """Trimming should still work correctly when descent contains a NaN."""
        el_nan = self.el.copy()
        el_nan[-200] = np.nan      # NaN deep in descent
        _, flt_m, _, _, _ = _phases(el_nan, self.v, self.min_alt)
        assert not flt_m[self.approach_mid], "Front trim still works with NaN in descent"
        assert flt_m[self.settled_float_mid], "Settled float still detected with NaN in descent"
