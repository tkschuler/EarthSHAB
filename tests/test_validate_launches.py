"""test_validate_launches.py - Validate launches.json and referenced data files.

Tests are organized as one class per launch so pytest groups and indents them
by balloon, making it easy to see which launch has issues at a glance.

Run:
    pytest tests/test_validate_launches.py -v

Extending:
  - Add new per-launch checks by adding a method to the class body inside
    _make_launch_class() below.
  - Add new required fields to REQUIRED_FIELDS.
  - Add new valid balloon shapes to VALID_SHAPES.
  - Add new APRS CSV formats to APRS_COLUMN_SETS.
"""

import json
from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

# ── Paths ─────────────────────────────────────────────────────────────────────

_REPO_ROOT        = Path(__file__).parents[1]
_EVAL_DIR         = _REPO_ROOT / "evaluation"
_BALLOON_DATA_DIR = _REPO_ROOT / "src" / "EarthSHAB" / "balloon_data"
_FORECASTS_DIR    = _REPO_ROOT / "src" / "EarthSHAB" / "forecasts"
_LAUNCHES_JSON    = _EVAL_DIR / "launches.json"
_LAUNCHES_EXAMPLE = _EVAL_DIR / "launches.example.json"

if not _LAUNCHES_JSON.exists():
    pytest.skip(
        f"launches.json not found — copy the example to get started:\n"
        f"  cp evaluation/launches.example.json evaluation/launches.json",
        allow_module_level=True,
    )

# ── Configurable constants ────────────────────────────────────────────────────

REQUIRED_FIELDS = [
    "shab_name", "organization", "launch_time", "sim_time_hr",
    "aprs_file", "launch_lat", "launch_lon", "launch_alt_m",
    "payload_weight_kg", "envelope_weight_kg", "balloon_shape", "balloon_size",
]

FORECAST_FIELDS = ["gfs_file", "era5_file"]

VALID_SHAPES = {"sphere", "trapezoid"}

APRS_COLUMN_SETS = [
    {"time", "lat", "lng", "altitude"},                              # APRS.fi standard
    {"Date", "Time(UTC)", "Latitude", "Longitude", "Altitude(m)"},   # Custom onboard
]

GFS_REQUIRED_VARS  = {"ugrdprs", "vgrdprs", "hgtprs", "lat", "lon", "lev", "time"}
ERA5_REQUIRED_VARS = {"u", "v", "z", "latitude", "longitude", "level", "time"}

# ── Helpers ───────────────────────────────────────────────────────────────────

def _load_launches() -> list[dict]:
    if not _LAUNCHES_JSON.exists():
        return []
    with open(_LAUNCHES_JSON) as f:
        return json.load(f).get("launches", [])


def _lid(launch: dict) -> str:
    return (
        f"{launch.get('organization', '?')}"
        f"_{launch.get('shab_name', '?')}"
        f"_{launch.get('launch_time', '')[:10]}"
    )


def _class_name(launch: dict) -> str:
    return "Test_" + _lid(launch).replace("-", "_")

# ── Global launches.json tests ────────────────────────────────────────────────

def test_launches_json_exists():
    assert _LAUNCHES_JSON.exists(), f"launches.json not found: {_LAUNCHES_JSON}"


def test_launches_json_valid_json():
    assert _LAUNCHES_JSON.exists(), "launches.json missing"
    data = json.loads(_LAUNCHES_JSON.read_text())
    assert isinstance(data, dict)
    assert "launches" in data
    assert isinstance(data["launches"], list)


def test_launches_json_non_empty():
    assert len(_load_launches()) > 0, "launches.json contains no launch entries"


def test_no_duplicate_launches():
    ids = [_lid(l) for l in _load_launches()]
    seen, dupes = set(), []
    for lid in ids:
        if lid in seen:
            dupes.append(lid)
        seen.add(lid)
    assert not dupes, f"Duplicate launch entries found: {dupes}"

# ── Per-launch test class factory ─────────────────────────────────────────────

def _make_launch_class(launch: dict):
    """Return a pytest-discoverable test class for one launch entry.

    Adding a new check: define a new method (name must start with 'test_')
    inside this function and add it to the `methods` dict at the bottom.
    """

    # ── Metadata checks ───────────────────────────────────────────────────────

    def test_required_fields(_):
        missing = [f for f in REQUIRED_FIELDS if f not in launch]
        assert not missing, f"Missing required fields: {missing}"

    def test_has_at_least_one_forecast(_):
        assert any(launch.get(f) for f in FORECAST_FIELDS), (
            f"Must have at least one non-null value among: {FORECAST_FIELDS}"
        )

    def test_launch_time_valid(_):
        lt = launch.get("launch_time", "")
        try:
            datetime.fromisoformat(lt)
        except (ValueError, TypeError):
            pytest.fail(f"launch_time '{lt}' is not a valid ISO datetime string")

    def test_sim_time_positive(_):
        v = launch.get("sim_time_hr")
        assert isinstance(v, (int, float)) and v > 0, (
            f"sim_time_hr must be > 0, got {v}"
        )

    def test_balloon_shape(_):
        shape = launch.get("balloon_shape", "")
        assert shape in VALID_SHAPES, (
            f"balloon_shape '{shape}' not in {VALID_SHAPES}. "
            f"Add it to VALID_SHAPES if intentional."
        )

    def test_balloon_size(_):
        v = launch.get("balloon_size")
        assert isinstance(v, (int, float)) and v > 0, (
            f"balloon_size must be > 0 (diameter in meters), got {v}"
        )

    def test_weights(_):
        for field in ["payload_weight_kg", "envelope_weight_kg"]:
            v = launch.get(field)
            if v is not None:
                assert isinstance(v, (int, float)) and v > 0, (
                    f"{field} must be a positive number, got {v}"
                )

    def test_coordinates(_):
        lat = launch.get("launch_lat")
        lon = launch.get("launch_lon")
        alt = launch.get("launch_alt_m")
        if lat is not None:
            assert -90 <= lat <= 90, f"launch_lat {lat} out of range [-90, 90]"
        if lon is not None:
            assert -180 <= lon <= 180, f"launch_lon {lon} out of range [-180, 180]"
        if alt is not None:
            assert alt >= 0, f"launch_alt_m {alt} must be >= 0"

    # ── APRS file checks ──────────────────────────────────────────────────────

    def test_aprs_exists(_):
        path = _BALLOON_DATA_DIR / launch.get("aprs_file", "")
        assert path.exists(), f"APRS file not found: {path}"

    def test_aprs_columns(_):
        path = _BALLOON_DATA_DIR / launch.get("aprs_file", "")
        if not path.exists():
            pytest.skip("APRS file missing — covered by test_aprs_exists")
        found = set(pd.read_csv(path, nrows=3).columns)
        assert any(req <= found for req in APRS_COLUMN_SETS), (
            f"Does not match any known column format.\n"
            f"  Found:    {found}\n"
            f"  Expected one of: {APRS_COLUMN_SETS}"
        )

    def test_aprs_non_empty(_):
        path = _BALLOON_DATA_DIR / launch.get("aprs_file", "")
        if not path.exists():
            pytest.skip("APRS file missing — covered by test_aprs_exists")
        assert len(pd.read_csv(path)) > 0, f"APRS file has no data rows"

    methods = {
        "test_required_fields":        test_required_fields,
        "test_has_at_least_one_forecast": test_has_at_least_one_forecast,
        "test_launch_time_valid":      test_launch_time_valid,
        "test_sim_time_positive":      test_sim_time_positive,
        "test_balloon_shape":          test_balloon_shape,
        "test_balloon_size":           test_balloon_size,
        "test_weights":                test_weights,
        "test_coordinates":            test_coordinates,
        "test_aprs_exists":            test_aprs_exists,
        "test_aprs_columns":           test_aprs_columns,
        "test_aprs_non_empty":         test_aprs_non_empty,
    }

    # ── GFS checks (only if gfs_file is set) ─────────────────────────────────

    if launch.get("gfs_file"):

        def test_gfs_exists(_):
            path = _FORECASTS_DIR / launch["gfs_file"]
            assert path.exists(), f"GFS forecast not found: {path}"

        def test_gfs_filename_parseable(_):
            parts = launch["gfs_file"].replace(".nc", "").split("_")
            assert len(parts) >= 4, (
                f"GFS filename must match gfs_0p25_YYYYMMDD_HH[_suffix].nc"
            )
            assert len(parts[2]) == 8 and parts[2].isdigit(), (
                f"Date part '{parts[2]}' must be 8 digits (YYYYMMDD)"
            )
            assert parts[3].isdigit() and 0 <= int(parts[3]) <= 23, (
                f"Hour part '{parts[3]}' must be 0–23"
            )

        def test_gfs_variables(_):
            path = _FORECASTS_DIR / launch["gfs_file"]
            if not path.exists():
                pytest.skip("GFS file missing — covered by test_gfs_exists")
            try:
                import netCDF4 as nc
                ds = nc.Dataset(path)
                found = set(ds.variables.keys())
                ds.close()
                missing = GFS_REQUIRED_VARS - found
                assert not missing, f"Missing variables: {missing}. Found: {found}"
            except ImportError:
                pytest.skip("netCDF4 not installed")

        methods["test_gfs_exists"]              = test_gfs_exists
        methods["test_gfs_filename_parseable"]  = test_gfs_filename_parseable
        methods["test_gfs_variables"]           = test_gfs_variables

    # ── ERA5 checks (only if era5_file is set) ────────────────────────────────

    if launch.get("era5_file"):

        def test_era5_exists(_):
            path = _FORECASTS_DIR / launch["era5_file"]
            assert path.exists(), f"ERA5 forecast not found: {path}"

        def test_era5_variables(_):
            path = _FORECASTS_DIR / launch["era5_file"]
            if not path.exists():
                pytest.skip("ERA5 file missing — covered by test_era5_exists")
            try:
                import netCDF4 as nc
                ds = nc.Dataset(path)
                found = set(ds.variables.keys())
                ds.close()
                missing = ERA5_REQUIRED_VARS - found
                assert not missing, f"Missing variables: {missing}. Found: {found}"
            except ImportError:
                pytest.skip("netCDF4 not installed")

        methods["test_era5_exists"]    = test_era5_exists
        methods["test_era5_variables"] = test_era5_variables

    cls = type(_class_name(launch), (), methods)
    cls.__qualname__ = cls.__name__
    return cls

# ── Inject one class per launch into module scope so pytest discovers them ────

for _launch in _load_launches():
    _cls = _make_launch_class(_launch)
    globals()[_cls.__name__] = _cls
