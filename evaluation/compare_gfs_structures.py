"""
Compare the NetCDF variable/dimension structure of two GFS forecast files
to diagnose KeyError: 'time' failures.

Usage:
    python -m evaluation.compare_gfs_structures
"""

import netCDF4 as nc
import numpy as np
import os

FORECASTS_DIR = os.path.join(os.path.dirname(__file__), "..", "src", "EarthSHAB", "forecasts")

FILES = {
    "SHAB14V (working)":  "gfs_0p25_20220822_12.nc",
    "SHAB1   (failing)":  "gfs_0p25_20201001_12_SHAB1.nc",
}


def describe_var(v):
    shape = tuple(v.shape)
    dtype = str(v.dtype)
    dims  = tuple(v.dimensions)
    attrs = {k: getattr(v, k) for k in v.ncattrs()}
    return {"shape": shape, "dtype": dtype, "dims": dims, "attrs": attrs}


def print_separator(char="─", width=70):
    print(char * width)


def compare():
    datasets = {}
    for label, fname in FILES.items():
        path = os.path.join(FORECASTS_DIR, fname)
        ds = nc.Dataset(path, "r")
        datasets[label] = ds
        print(f"\n{'='*70}")
        print(f"  {label}")
        print(f"  {path}")
        print(f"{'='*70}")

        print("\n  Dimensions:")
        for name, dim in ds.dimensions.items():
            print(f"    {name:30s}  size={len(dim)}  unlimited={dim.isunlimited()}")

        print("\n  Variables:")
        for name, var in ds.variables.items():
            info = describe_var(var)
            print(f"    {name:30s}  shape={str(info['shape']):20s}  dtype={info['dtype']:10s}  dims={info['dims']}")
            if name.lower() in ("time", "reftime", "time1", "forecast_reference_time"):
                units = info["attrs"].get("units", "<no units>")
                calendar = info["attrs"].get("calendar", "<no calendar>")
                try:
                    vals = var[:]
                    sample = vals[:5].tolist() if len(vals) >= 5 else vals[:].tolist()
                except Exception:
                    sample = "<unreadable>"
                print(f"      └ units={units!r}  calendar={calendar!r}  sample={sample}")

        print("\n  Global attributes:")
        for attr in ds.ncattrs():
            val = getattr(ds, attr)
            print(f"    {attr:30s}  {str(val)[:60]}")

    print(f"\n\n{'='*70}")
    print("  STRUCTURAL DIFF")
    print(f"{'='*70}")

    labels = list(datasets.keys())
    a_vars = set(datasets[labels[0]].variables.keys())
    b_vars = set(datasets[labels[1]].variables.keys())

    only_a = a_vars - b_vars
    only_b = b_vars - a_vars
    common = a_vars & b_vars

    if only_a:
        print(f"\n  Only in {labels[0]}:")
        for v in sorted(only_a):
            print(f"    + {v}")

    if only_b:
        print(f"\n  Only in {labels[1]}:")
        for v in sorted(only_b):
            print(f"    + {v}")

    if not only_a and not only_b:
        print("\n  Variable names match exactly.")

    print(f"\n  Shape/dim differences in common variables:")
    found_diff = False
    for v in sorted(common):
        da = describe_var(datasets[labels[0]].variables[v])
        db = describe_var(datasets[labels[1]].variables[v])
        diffs = []
        if da["shape"] != db["shape"]:
            diffs.append(f"shape {da['shape']} vs {db['shape']}")
        if da["dims"] != db["dims"]:
            diffs.append(f"dims {da['dims']} vs {db['dims']}")
        if da["dtype"] != db["dtype"]:
            diffs.append(f"dtype {da['dtype']} vs {db['dtype']}")
        if diffs:
            print(f"    {v}: {', '.join(diffs)}")
            found_diff = True

    if not found_diff:
        print("    (none)")

    for ds in datasets.values():
        ds.close()


if __name__ == "__main__":
    compare()
