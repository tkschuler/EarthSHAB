"""test_gfs_archive.py - Unit tests for gfs_archive pure logic.

Covers .idx parsing, byte-range selection (variable/level filtering, last-message
open range, requested-level subsetting), and the field diff-stat helper. Network
downloads and cfgrib/xarray conversion are not covered here (integration only).

Run:
    PYTHONPATH=src pytest tests/test_gfs_archive.py -v
"""

import numpy as np

from evaluation.gfs_archive import _level_mb, parse_idx, select_ranges
from evaluation.reproduce_gfs import diff_stats


# Sample .idx body: 6 messages, 4 vars across 2 isobaric levels + a surface field.
SAMPLE_IDX = (
    "1:0:d=2022082212:UGRD:500 mb:anl:\n"
    "2:100:d=2022082212:VGRD:500 mb:anl:\n"
    "3:250:d=2022082212:HGT:500 mb:anl:\n"
    "4:500:d=2022082212:TMP:500 mb:anl:\n"
    "5:900:d=2022082212:UGRD:0.4 mb:anl:\n"
    "6:1000:d=2022082212:PRMSL:mean sea level:anl:\n"
)


def test_parse_idx_basic():
    rows = parse_idx(SAMPLE_IDX)
    assert len(rows) == 6
    assert rows[0] == (1, 0, "UGRD", "500 mb")
    assert rows[4] == (5, 900, "UGRD", "0.4 mb")


def test_parse_idx_skips_garbage():
    rows = parse_idx("not:a:number\n\n:::\n7:abc:d=x:UGRD:500 mb:anl:")
    assert rows == []


def test_level_mb():
    assert _level_mb("500 mb") == 500.0
    assert _level_mb("0.4 mb") == 0.4
    assert _level_mb("surface") is None
    assert _level_mb("2 m above ground") is None


def test_select_ranges_filters_vars_and_levels():
    rows = parse_idx(SAMPLE_IDX)
    # PRMSL excluded (not in WANT_VARS); surface level excluded (not ' mb').
    ranges = select_ranges(rows)
    vars_kept = {r[2] for r in ranges}
    assert vars_kept == {"UGRD", "VGRD", "HGT", "TMP"}
    assert len(ranges) == 5  # 4 at 500mb + 1 UGRD at 0.4mb


def test_select_ranges_end_offsets_and_open_last():
    rows = parse_idx(SAMPLE_IDX)
    ranges = select_ranges(rows)
    # message 1: [0, 99]; message 4 (TMP 500): end = 900-1 = 899
    assert ranges[0][:2] == (0, 99)
    tmp = [r for r in ranges if r[2] == "TMP"][0]
    assert tmp[:2] == (500, 899)
    # last selected message (UGRD 0.4 mb, offset 900) is followed by PRMSL@1000,
    # so its end is 999 (not open) — it is not the final row in the file.
    last_iso = [r for r in ranges if r[3] == 0.4][0]
    assert last_iso[:2] == (900, 999)


def test_select_ranges_open_ended_final_message():
    # When the wanted message IS the final row, end is None (open-ended range).
    idx = "1:0:d=x:TMP:850 mb:anl:\n2:300:d=x:UGRD:850 mb:anl:\n"
    ranges = select_ranges(parse_idx(idx))
    assert ranges[-1][1] is None
    assert ranges[-1][0] == 300


def test_select_ranges_level_subset():
    rows = parse_idx(SAMPLE_IDX)
    ranges = select_ranges(rows, want_levels_mb=[0.4])
    assert len(ranges) == 1
    assert ranges[0][2] == "UGRD" and ranges[0][3] == 0.4


def test_diff_stats_identical_is_zero():
    a = np.array([[1.0, 2.0], [3.0, 4.0]])
    s = diff_stats(a, a.copy())
    assert s["max_abs"] == 0.0
    assert s["mean_abs"] == 0.0
    assert s["rmse"] == 0.0


def test_diff_stats_known_values():
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([3.0, 0.0, 4.0])  # abs diffs 3,0,4
    s = diff_stats(a, b)
    assert s["max_abs"] == 4.0
    assert abs(s["mean_abs"] - (7.0 / 3.0)) < 1e-12
    assert abs(s["rmse"] - np.sqrt((9 + 0 + 16) / 3.0)) < 1e-12


def test_diff_stats_ignores_nan():
    a = np.array([1.0, np.nan, 2.0])
    b = np.array([1.0, 5.0, 2.0])  # only the non-nan pairs (both finite) count
    s = diff_stats(a, b)
    assert s["max_abs"] == 0.0
    assert s["n"] == 2
