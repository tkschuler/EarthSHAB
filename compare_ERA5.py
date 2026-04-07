# compare_ERA5.py
# pytest-based comparison tests for two forecast NetCDF files that should contain
# the same meteorological data but use different formatting / naming conventions.
#
# Run:
#   pytest -q compare_ERA5.py
# or for print output:
#   pytest -q -s compare_ERA5.py
#
# Notes:
# - This script normalizes:
#     lev   -> level
#     lat   -> latitude
#     lon   -> longitude
#     tmpprs  -> t
#     ugrdprs -> u
#     vgrdprs -> v
#     hgtprs  -> z
# - It converts longitudes from 0..360 to -180..180 when needed.
# - It sorts coordinates into a consistent ascending order.
# - It compares the exact shared grid between the two files.
# - Tolerances are chosen from your debug output, where the shared-grid differences
#   were very small but not exactly zero.


from __future__ import annotations

import numpy as np
import pytest
import xarray as xr


# --------------------------------------------------------------------------
# EDIT THESE PATHS
# --------------------------------------------------------------------------
FILE_A = "../HAB-COM/GFS_DATA/COMBINED/gfs_20260407_06z_combined_processed.nc"
FILE_B = "/home/schuler/EarthSHAB/forecasts/gfs_0p25_20260407_06.nc"


# --------------------------------------------------------------------------
# Full-field tolerances on exact overlap
# --------------------------------------------------------------------------
FIELD_TOLS = {
    "u": {"atol": 0.11, "rtol": 1e-4},
    "v": {"atol": 0.11, "rtol": 1e-4},
    "t": {"atol": 0.02, "rtol": 1e-5},
    "z": {"atol": 0.05, "rtol": 1e-5},
}

# --------------------------------------------------------------------------
# Summary-stat tolerances
#
# These need to be looser than field tolerances for extrema, because
# np.nanmin / np.nanmax compare only one element and can differ by a bit more.
# --------------------------------------------------------------------------
STAT_TOLS = {
    "u": {"atol": 0.03, "rtol": 1e-4},
    "v": {"atol": 0.06, "rtol": 1e-4},
    "t": {"atol": 0.005, "rtol": 5e-5},
    "z": {"atol": 0.05, "rtol": 5e-5},
}


def open_dataset(path: str) -> xr.Dataset:
    return xr.open_dataset(path, engine="netcdf4")


def normalize_dataset(ds: xr.Dataset) -> xr.Dataset:
    rename_map: dict[str, str] = {}

    if "lev" in ds.coords or "lev" in ds.dims:
        rename_map["lev"] = "level"
    if "lat" in ds.coords or "lat" in ds.dims:
        rename_map["lat"] = "latitude"
    if "lon" in ds.coords or "lon" in ds.dims:
        rename_map["lon"] = "longitude"

    if "tmpprs" in ds.data_vars:
        rename_map["tmpprs"] = "t"
    if "ugrdprs" in ds.data_vars:
        rename_map["ugrdprs"] = "u"
    if "vgrdprs" in ds.data_vars:
        rename_map["vgrdprs"] = "v"
    if "hgtprs" in ds.data_vars:
        rename_map["hgtprs"] = "z"

    if rename_map:
        ds = ds.rename(rename_map)

    if "valid_time" in ds.coords:
        ds = ds.assign_coords(time=ds["valid_time"].values)

    if "longitude" in ds.coords:
        lon = ds["longitude"].values
        if np.nanmax(lon) > 180:
            new_lon = ((lon + 180) % 360) - 180
            ds = ds.assign_coords(longitude=new_lon)

    for coord in ("time", "level", "latitude", "longitude"):
        if coord in ds.coords:
            ds = ds.sortby(coord)

    target_order = ("time", "level", "latitude", "longitude")
    for var in ("u", "v", "t", "z"):
        if var in ds.data_vars:
            ds[var] = ds[var].transpose(*target_order)

    return ds


def assert_strictly_increasing(values: np.ndarray, coord_name: str) -> None:
    diffs = np.diff(values)
    assert np.all(diffs > 0), (
        f"{coord_name} is not strictly increasing. "
        f"First values: {values[:5]} Last values: {values[-5:]}"
    )


def overlapping_coords(ds_a: xr.Dataset, ds_b: xr.Dataset):
    times = np.intersect1d(ds_a.time.values, ds_b.time.values)
    levels = np.intersect1d(ds_a.level.values, ds_b.level.values)
    lats = np.intersect1d(ds_a.latitude.values, ds_b.latitude.values)
    lons = np.intersect1d(ds_a.longitude.values, ds_b.longitude.values)
    return times, levels, lats, lons


def restrict_to_overlap(
    ds: xr.Dataset,
    times: np.ndarray,
    levels: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
) -> xr.Dataset:
    return ds.sel(
        time=times,
        level=levels,
        latitude=lats,
        longitude=lons,
    )


def print_dataset_summary(name: str, ds: xr.Dataset) -> None:
    print("\n" + "=" * 80)
    print(name)
    print("dims:", ds.dims)

    for coord in ("time", "level", "latitude", "longitude"):
        vals = ds[coord].values
        print(f"\n{coord}:")
        print("  first 5:", vals[:5])
        print("  last  5:", vals[-5:])
        print("  len:", len(vals))

    for var in ("u", "v", "t", "z"):
        arr = ds[var]
        print(f"\n{var}: dims={arr.dims}, shape={arr.shape}, dtype={arr.dtype}")
        if hasattr(arr, "attrs") and arr.attrs:
            print(f"  attrs: {arr.attrs}")


def print_overlap_summary(ds_a: xr.Dataset, ds_b: xr.Dataset) -> None:
    times, levels, lats, lons = overlapping_coords(ds_a, ds_b)

    print("\n" + "=" * 80)
    print("OVERLAP SUMMARY")
    print("times :", len(times), times[:5], "...", times[-5:] if len(times) else [])
    print("levels:", len(levels), levels[:10], "...", levels[-10:] if len(levels) else [])
    print("lats  :", len(lats), lats[:10], "...", lats[-10:] if len(lats) else [])
    print("lons  :", len(lons), lons[:10], "...", lons[-10:] if len(lons) else [])


@pytest.fixture(scope="module")
def ds_a() -> xr.Dataset:
    return normalize_dataset(open_dataset(FILE_A))


@pytest.fixture(scope="module")
def ds_b() -> xr.Dataset:
    return normalize_dataset(open_dataset(FILE_B))


@pytest.fixture(scope="module")
def ds_overlap(ds_a: xr.Dataset, ds_b: xr.Dataset):
    times, levels, lats, lons = overlapping_coords(ds_a, ds_b)

    assert len(times) > 0, "No overlapping times"
    assert len(levels) > 0, "No overlapping levels"
    assert len(lats) > 0, "No overlapping latitudes"
    assert len(lons) > 0, "No overlapping longitudes"

    a = restrict_to_overlap(ds_a, times, levels, lats, lons)
    b = restrict_to_overlap(ds_b, times, levels, lats, lons)
    return a, b


def test_debug_dataset_summaries(ds_a: xr.Dataset, ds_b: xr.Dataset):
    print_dataset_summary("DATASET A", ds_a)
    print_dataset_summary("DATASET B", ds_b)
    print_overlap_summary(ds_a, ds_b)
    assert True


def test_required_variables_exist(ds_a: xr.Dataset, ds_b: xr.Dataset):
    for ds in (ds_a, ds_b):
        for var in ("u", "v", "t", "z"):
            assert var in ds.data_vars, f"Missing variable {var}"


def test_standard_dims(ds_a: xr.Dataset, ds_b: xr.Dataset):
    for ds in (ds_a, ds_b):
        for var in ("u", "v", "t", "z"):
            assert ds[var].dims == ("time", "level", "latitude", "longitude")


def test_coords_are_strictly_increasing(ds_a: xr.Dataset, ds_b: xr.Dataset):
    for ds in (ds_a, ds_b):
        assert_strictly_increasing(ds.time.values, "time")
        assert_strictly_increasing(ds.level.values, "level")
        assert_strictly_increasing(ds.latitude.values, "latitude")
        assert_strictly_increasing(ds.longitude.values, "longitude")


def test_longitudes_are_normalized(ds_a: xr.Dataset, ds_b: xr.Dataset):
    for ds in (ds_a, ds_b):
        lon = ds.longitude.values
        assert np.nanmin(lon) >= -180
        assert np.nanmax(lon) <= 180


def test_overlap_exists(ds_a: xr.Dataset, ds_b: xr.Dataset):
    times, levels, lats, lons = overlapping_coords(ds_a, ds_b)
    assert len(times) > 0
    assert len(levels) > 0
    assert len(lats) > 0
    assert len(lons) > 0


def test_overlap_shapes_match(ds_overlap):
    a, b = ds_overlap
    for var in ("u", "v", "t", "z"):
        assert a[var].shape == b[var].shape


def test_overlap_coords_match_exactly(ds_overlap):
    a, b = ds_overlap
    np.testing.assert_array_equal(a.time.values, b.time.values)
    np.testing.assert_array_equal(a.level.values, b.level.values)
    np.testing.assert_array_equal(a.latitude.values, b.latitude.values)
    np.testing.assert_array_equal(a.longitude.values, b.longitude.values)


@pytest.mark.parametrize("var", ["u", "v", "t", "z"])
def test_values_close_on_exact_overlap(ds_overlap, var: str):
    a, b = ds_overlap

    np.testing.assert_allclose(
        a[var].values,
        b[var].values,
        rtol=FIELD_TOLS[var]["rtol"],
        atol=FIELD_TOLS[var]["atol"],
        err_msg=(
            f"{var} differs more than allowed on exact overlap "
            f"(rtol={FIELD_TOLS[var]['rtol']}, atol={FIELD_TOLS[var]['atol']})"
        ),
    )


@pytest.mark.parametrize("var", ["u", "v", "t", "z"])
def test_nan_masks_match(ds_overlap, var: str):
    a, b = ds_overlap
    np.testing.assert_array_equal(
        np.isnan(a[var].values),
        np.isnan(b[var].values),
        err_msg=f"NaN masks differ for {var}",
    )


@pytest.mark.parametrize("var", ["u", "v", "t", "z"])
def test_summary_stats_close(ds_overlap, var: str):
    a, b = ds_overlap

    av = a[var].values
    bv = b[var].values

    rtol = STAT_TOLS[var]["rtol"]
    atol = STAT_TOLS[var]["atol"]

    np.testing.assert_allclose(
        np.nanmean(av), np.nanmean(bv), rtol=rtol, atol=atol,
        err_msg=f"{var} mean differs too much"
    )
    np.testing.assert_allclose(
        np.nanstd(av), np.nanstd(bv), rtol=rtol, atol=atol,
        err_msg=f"{var} std differs too much"
    )
    np.testing.assert_allclose(
        np.nanmin(av), np.nanmin(bv), rtol=rtol, atol=atol,
        err_msg=f"{var} min differs too much"
    )
    np.testing.assert_allclose(
        np.nanmax(av), np.nanmax(bv), rtol=rtol, atol=atol,
        err_msg=f"{var} max differs too much"
    )


@pytest.mark.parametrize("var", ["u", "v", "t", "z"])
def test_max_abs_difference_within_expected(ds_overlap, var: str):
    a, b = ds_overlap
    max_abs_diff = np.nanmax(np.abs(a[var].values - b[var].values))

    expected_max = {
        "u": 0.11,
        "v": 0.11,
        "t": 0.02,
        "z": 0.05,
    }[var]

    assert max_abs_diff <= expected_max, (
        f"{var} max abs diff {max_abs_diff} exceeds expected {expected_max}"
    )


@pytest.mark.parametrize("var", ["u", "v", "t", "z"])
def test_print_comparison_stats(ds_overlap, var: str):
    a, b = ds_overlap
    av = a[var].values
    bv = b[var].values
    diff = av - bv

    print("\n" + "-" * 80)
    print(f"{var} exact-overlap comparison")
    print(
        "A mean/std/min/max:",
        float(np.nanmean(av)),
        float(np.nanstd(av)),
        float(np.nanmin(av)),
        float(np.nanmax(av)),
    )
    print(
        "B mean/std/min/max:",
        float(np.nanmean(bv)),
        float(np.nanstd(bv)),
        float(np.nanmin(bv)),
        float(np.nanmax(bv)),
    )
    print(
        "diff mean/std/min/max:",
        float(np.nanmean(diff)),
        float(np.nanstd(diff)),
        float(np.nanmin(diff)),
        float(np.nanmax(diff)),
    )
    print("max abs diff:", float(np.nanmax(np.abs(diff))))

    assert True


def nearest_index(values: np.ndarray, x: float) -> int:
    return int(np.abs(values - x).argmin())


@pytest.mark.parametrize(
    "sample",
    [
        {"time_idx": 0, "level": 100.0, "lat": 30.0, "lon": -120.0},
        {"time_idx": 2, "level": 500.0, "lat": 35.0, "lon": -110.0},
        {"time_idx": 4, "level": 850.0, "lat": 40.0, "lon": -100.0},
    ],
)
def test_point_samples_are_close(ds_overlap, sample):
    a, b = ds_overlap

    time_val = a.time.values[sample["time_idx"]]
    level_val = a.level.values[nearest_index(a.level.values, sample["level"])]
    lat_val = a.latitude.values[nearest_index(a.latitude.values, sample["lat"])]
    lon_val = a.longitude.values[nearest_index(a.longitude.values, sample["lon"])]

    a_point = a.sel(time=time_val, level=level_val, latitude=lat_val, longitude=lon_val)
    b_point = b.sel(time=time_val, level=level_val, latitude=lat_val, longitude=lon_val)

    point_atol = {
        "u": 0.11,
        "v": 0.11,
        "t": 0.02,
        "z": 0.05,
    }

    for var in ("u", "v", "t", "z"):
        np.testing.assert_allclose(
            a_point[var].values,
            b_point[var].values,
            rtol=0.0,
            atol=point_atol[var],
            err_msg=(
                f"Point sample mismatch for {var} at "
                f"time={time_val}, level={level_val}, lat={lat_val}, lon={lon_val}"
            ),
        )


def pytest_sessionfinish(session, exitstatus):
    return