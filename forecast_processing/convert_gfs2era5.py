# convert_processed_to_gfs_style.py

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc

INPUT_FILE = "../HAB-COM/GFS_DATA/COMBINED/gfs_20260407_06z_combined_processed.nc"
OUTPUT_FILE = "forecasts/converted_gfs_style.nc"


ds = xr.open_dataset(INPUT_FILE, engine="netcdf4")
print(ds)

ds = xr.open_dataset(OUTPUT_FILE, engine="netcdf4")
print(ds)
dfgdfg

def convert_to_gfs_style(input_file: str, output_file: str) -> None:
    ds = xr.open_dataset(input_file, engine="netcdf4")

    required_coords = {"time", "level", "latitude", "longitude"}
    required_vars = {"u", "v", "t", "z"}

    missing_coords = required_coords - set(ds.coords)
    missing_vars = required_vars - set(ds.data_vars)

    if missing_coords:
        raise ValueError(f"Missing coords: {sorted(missing_coords)}")
    if missing_vars:
        raise ValueError(f"Missing vars: {sorted(missing_vars)}")

    # Rename to GFS-style names
    ds = ds.rename(
        {
            "level": "lev",
            "latitude": "lat",
            "longitude": "lon",
            "u": "ugrdprs",
            "v": "vgrdprs",
            "t": "tmpprs",
            "z": "hgtprs",
        }
    )

    # Convert lon from [-180, 180) to [0, 360)
    lon = ds["lon"].values.astype(float)
    lon = lon % 360.0
    ds = ds.assign_coords(lon=lon)
    ds = ds.sortby("lon")

    # Match GFS-style ordering
    ds = ds.sortby("lat", ascending=True)     # File 2 has lat ascending
    ds = ds.sortby("lev", ascending=False)    # File 2 has pressure descending

    # Standard dim order
    ds = ds.transpose("time", "lev", "lat", "lon")

    # Variable order
    var_order = ["hgtprs", "tmpprs", "ugrdprs", "vgrdprs"]
    ds = xr.Dataset({v: ds[v] for v in var_order}, coords=ds.coords)

    # Coordinate dtypes to look more like File 2
    ds["lev"] = ds["lev"].astype(np.float64)
    ds["lat"] = ds["lat"].astype(np.float64)
    ds["lon"] = ds["lon"].astype(np.float64)

    # Convert time to numeric "days since 0001-01-01"
    time_values = ds["time"].values
    datetime_objects = [pd.Timestamp(t).to_pydatetime() for t in time_values]

    julian_dates = nc.date2num(
        datetime_objects,
        units="days since 0001-01-01",
        calendar="standard",
        has_year_zero=True,
    )

    ds = ds.assign_coords(time=julian_dates.astype(np.float64))
    ds["time"].attrs = {
        "units": "days since 0001-01-01",
        "calendar": "standard",
    }

    ds["lev"].attrs.update({"units": "hPa", "long_name": "pressure level"})
    ds["lat"].attrs.update({"units": "degrees_north"})
    ds["lon"].attrs.update({"units": "degrees_east"})

    ds.attrs = dict(ds.attrs)
    ds.attrs["Conventions"] = "CF-1.7"

    encoding = {var: {"zlib": True, "complevel": 4} for var in ds.data_vars}
    encoding["time"] = {"_FillValue": None}

    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out, engine="netcdf4", encoding=encoding)

    print(f"Saved: {out}")
    print(ds)
    print("Dimensions:", dict(ds.sizes))


if __name__ == "__main__":
    convert_to_gfs_style(INPUT_FILE, OUTPUT_FILE)