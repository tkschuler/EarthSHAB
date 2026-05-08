"""test_nan_handling.py - Tests for NaN-related bug fixes across the codebase.

Covers two bugs that shared the same root cause: a single corrupted APRS packet
produced a NaN altitude, which then cascaded in two separate ways.

  1. Reforecast NaN cascade  (simulate.py – run_simulation reforecast loop)
     NaN altitude in coord_new["alt"] was passed to the next getNewCoord call,
     which returned NaN lat/lon.  Those NaN coords were then used as the next
     starting position, propagating NaN through all remaining reforecast steps.
     Fix: forward-fill NaN altitude with the previous valid value before storing
     it in coord_new.

  2. HTML trajectory broken by NaN LatLng  (plot_trajectory_map.py – plot_map)
     The reforecast trajectory (orange) contained NaN lat/lon values.  When
     gmplot wrote these as `new google.maps.LatLng(NaN, NaN)`, the browser's
     JavaScript threw a runtime error that halted all subsequent polyline
     rendering — so the primary simulation trajectory (red/blue) never appeared.
     Fix: _drop_nan_coords() strips NaN pairs before any gmap1.plot() call.

Run:
    pytest tests/test_nan_handling.py -v
"""

import math
import textwrap

import numpy as np
import pytest

# ── Import targets ───────────────────────────────────────────────────────────

from EarthSHAB.Plotting.plot_trajectory_map import _drop_nan_coords
from EarthSHAB.simulate import _load_aprs


# ── Helpers ──────────────────────────────────────────────────────────────────

def _list(seq):
    """Materialise any iterable (including zip objects) to a plain list."""
    return list(seq)


# ── Tests: _drop_nan_coords ───────────────────────────────────────────────────

class TestDropNanCoords:
    """Unit tests for EarthSHAB.Plotting.plot_trajectory_map._drop_nan_coords."""

    def test_all_valid_passes_through(self):
        lats = [34.6, 35.0, 35.5]
        lons = [-106.8, -106.5, -106.0]
        out_lats, out_lons = _drop_nan_coords(lats, lons)
        assert _list(out_lats) == lats
        assert _list(out_lons) == lons

    def test_trailing_nan_removed(self):
        """The actual bug: 76 trailing NaN lat/lon from the NaN altitude cascade."""
        lats = [34.6, 35.0, float("nan"), float("nan")]
        lons = [-106.8, -106.5, float("nan"), float("nan")]
        out_lats, out_lons = _drop_nan_coords(lats, lons)
        assert _list(out_lats) == [34.6, 35.0]
        assert _list(out_lons) == [-106.8, -106.5]

    def test_nan_in_middle_removed(self):
        lats = [34.6, float("nan"), 35.5]
        lons = [-106.8, float("nan"), -106.0]
        out_lats, out_lons = _drop_nan_coords(lats, lons)
        assert _list(out_lats) == [34.6, 35.5]
        assert _list(out_lons) == [-106.8, -106.0]

    def test_all_nan_returns_empty(self):
        lats = [float("nan"), float("nan")]
        lons = [float("nan"), float("nan")]
        out_lats, out_lons = _drop_nan_coords(lats, lons)
        assert _list(out_lats) == []
        assert _list(out_lons) == []

    def test_numpy_nan_removed(self):
        """np.nan (numpy.float64 NaN) must be treated the same as float('nan')."""
        lats = np.array([34.6, np.nan, 35.5])
        lons = np.array([-106.8, np.nan, -106.0])
        out_lats, out_lons = _drop_nan_coords(lats, lons)
        result_lats = _list(out_lats)
        assert len(result_lats) == 2
        assert not any(math.isnan(x) for x in result_lats)

    def test_mixed_nan_only_bad_pairs_removed(self):
        """Only pairs where either component is NaN are dropped."""
        lats = [34.6, float("nan"), 35.5, 36.0]
        lons = [-106.8, -106.5, float("nan"), -105.0]
        out_lats, out_lons = _drop_nan_coords(lats, lons)
        # pairs (34.6,-106.8) and (36.0,-105.0) are clean
        assert _list(out_lats) == [34.6, 36.0]
        assert _list(out_lons) == [-106.8, -105.0]

    def test_empty_input_returns_empty(self):
        out_lats, out_lons = _drop_nan_coords([], [])
        assert _list(out_lats) == []
        assert _list(out_lons) == []

    def test_no_nan_count_unchanged(self):
        """Length is preserved when no NaN values are present."""
        lats = list(range(100))
        lons = list(range(100))
        out_lats, out_lons = _drop_nan_coords(lats, lons)
        assert len(_list(out_lats)) == 100
        assert len(_list(out_lons)) == 100


# ── Tests: _load_aprs NaN altitude preservation ───────────────────────────────

class TestLoadAprsNanAltitude:
    """_load_aprs must preserve NaN altitudes rather than silently dropping rows.

    The evaluator's _detect_phases now uses np.nanmax/np.nanargmax to handle
    the NaN gracefully — but only if the NaN actually reaches the function.
    """

    def test_aprsfi_format_nan_altitude_preserved(self, tmp_path):
        csv = textwrap.dedent("""\
            time,lat,lng,altitude,comment
            2024-06-13 17:46:00,34.640,-106.830,1500.0,
            2024-06-13 18:00:00,34.650,-106.820,16000.0,
            2024-06-13 19:22:00,34.700,-106.780,,corrupted packet
            2024-06-13 20:00:00,34.750,-106.750,15000.0,
        """)
        p = tmp_path / "test_aprs.csv"
        p.write_text(csv)

        df, fmt = _load_aprs(str(p))

        assert fmt == "aprs_fi"
        assert len(df) == 4, "All rows including the NaN-altitude row must be loaded"
        assert df["altitude"].isna().sum() == 1, "Exactly one NaN altitude expected"
        assert df["altitude"].isna().iloc[2], "NaN must be at the third row (index 2)"

    def test_aprsfi_format_valid_altitudes_unchanged(self, tmp_path):
        csv = textwrap.dedent("""\
            time,lat,lng,altitude,comment
            2024-06-13 17:46:00,34.640,-106.830,1500.0,
            2024-06-13 20:41:00,35.740,-106.683,1520.0,
        """)
        p = tmp_path / "test_valid.csv"
        p.write_text(csv)

        df, fmt = _load_aprs(str(p))

        assert fmt == "aprs_fi"
        assert df["altitude"].isna().sum() == 0
        assert list(df["altitude"]) == [1500.0, 1520.0]

    def test_aprsfi_time_converted_to_mst(self, tmp_path):
        """Timestamps are shifted from UTC to MST (UTC-7) by _load_aprs."""
        csv = textwrap.dedent("""\
            time,lat,lng,altitude,comment
            2024-06-13 17:46:00,34.640,-106.830,1500.0,
        """)
        p = tmp_path / "test_tz.csv"
        p.write_text(csv)

        df, _ = _load_aprs(str(p))

        # UTC 17:46 → MST 10:46
        assert df["time"].iloc[0].hour == 10
        assert df["time"].iloc[0].minute == 46


# ── Tests: reforecast NaN forward-fill ───────────────────────────────────────

class TestReforecastNanForwardFill:
    """The NaN altitude forward-fill in simulate.py's reforecast loop.

    When alt_aprs[i] is NaN, the fix substitutes coords_aprs[-1]["alt"]
    (the previous valid altitude) before creating coord_new.  Without this,
    the NaN altitude cascades: getNewCoord sees NaN alt → returns NaN lat/lon
    → stored as coord_new lat/lon → used as next starting position → all
    remaining coordinates are NaN.
    """

    # Minimal stand-in for ERA5/GFS.getNewCoord — returns NaN lat/lon when
    # the coordinate's altitude is NaN, otherwise advances by a small step.
    @staticmethod
    def _mock_get_new_coord(coord, dt):
        if np.isnan(coord["alt"]):
            return [np.nan, np.nan] + [np.nan] * 8
        return [coord["lat"] + 0.01, coord["lon"] + 0.01] + [0.0] * 8

    @staticmethod
    def _run_reforecast_loop(alt_aprs, apply_fix):
        """Simulate the reforecast loop from simulate.py for a small array."""
        coords_aprs  = [{"lat": 34.6, "lon": -106.8, "alt": 1500.0, "timestamp": None}]
        lat_aprs_gps = [34.6]
        lon_aprs_gps = [-106.8]

        mock = TestReforecastNanForwardFill._mock_get_new_coord

        for i in range(len(alt_aprs) - 1):
            result   = mock(coords_aprs[i], 60)
            lat_new, lon_new = result[0], result[1]

            alt_i = alt_aprs[i]
            if apply_fix and np.isnan(alt_i):
                alt_i = coords_aprs[-1]["alt"]   # ← the fix

            coord_new = {"lat": lat_new, "lon": lon_new, "alt": alt_i, "timestamp": None}
            coords_aprs.append(coord_new)
            lat_aprs_gps.append(lat_new)
            lon_aprs_gps.append(lon_new)

        return lat_aprs_gps, lon_aprs_gps, coords_aprs

    def test_without_fix_nan_cascades(self):
        """Demonstrates the original bug: one NaN altitude → many NaN lat/lon."""
        alt_aprs = np.array([1500.0, 16000.0, np.nan, 15000.0, 14000.0])
        lats, lons, coords = self._run_reforecast_loop(alt_aprs, apply_fix=False)
        nan_lats = sum(1 for x in lats if isinstance(x, float) and math.isnan(x))
        # Without fix: alt_aprs[2]=NaN → coord[3]["alt"]=NaN → lat[3], lat[4] are NaN
        assert nan_lats > 0, "Without the fix, NaN should cascade into lat_aprs_gps"

    def test_with_fix_no_nan_in_trajectory(self):
        """After the fix, no NaN lat/lon should appear in the reforecast trajectory."""
        alt_aprs = np.array([1500.0, 16000.0, np.nan, 15000.0, 14000.0])
        lats, lons, coords = self._run_reforecast_loop(alt_aprs, apply_fix=True)
        assert not any(isinstance(x, float) and math.isnan(x) for x in lats), (
            "No NaN should appear in lat_aprs_gps after forward-fill fix"
        )
        assert not any(isinstance(x, float) and math.isnan(x) for x in lons), (
            "No NaN should appear in lon_aprs_gps after forward-fill fix"
        )

    def test_with_fix_no_nan_in_coords(self):
        """After the fix, coords_aprs should never carry a NaN altitude."""
        alt_aprs = np.array([1500.0, 16000.0, np.nan, 15000.0, 14000.0])
        _, _, coords = self._run_reforecast_loop(alt_aprs, apply_fix=True)
        nan_alts = [c["alt"] for c in coords if np.isnan(c["alt"])]
        assert not nan_alts, "No NaN altitude should remain in coords_aprs after fix"

    def test_forward_fill_uses_previous_valid_altitude(self):
        """The fill value must be the last valid altitude, not zero or some default."""
        alt_aprs = np.array([5000.0, 18000.0, np.nan, 17000.0])
        _, _, coords = self._run_reforecast_loop(alt_aprs, apply_fix=True)
        # coords[0] alt=5000.0 (initial), coords[1] alt=5000.0 (from alt_aprs[0]),
        # coords[2] alt=18000.0 (from alt_aprs[1]), coords[3] alt=18000.0 (filled).
        assert coords[3]["alt"] == 18000.0, (
            "NaN altitude at index 2 should be replaced with 18000.0 (previous valid)"
        )

    def test_multiple_consecutive_nan_altitudes(self):
        """Multiple consecutive NaN altitudes must all be forward-filled."""
        alt_aprs = np.array([1500.0, 16000.0, np.nan, np.nan, np.nan, 15000.0])
        lats, lons, coords = self._run_reforecast_loop(alt_aprs, apply_fix=True)
        assert not any(isinstance(x, float) and math.isnan(x) for x in lats)
        assert not any(np.isnan(c["alt"]) for c in coords)
