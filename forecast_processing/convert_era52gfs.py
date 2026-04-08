# convert_gfs_style_to_processed.py

from __future__ import annotations
from pathlib import Path

import numpy as np
import xarray as xr
import netCDF4 as nc

INPUT_FILE = "forecasts/gfs_0p25_20260407_06.nc"
OUTPUT_FILE = "forecasts/converted_processed_style.nc"


def decode_time_if_needed(ds: xr.Dataset) -> xr.Dataset:
    """
    If time is numeric and has CF-style units/calendar attrs, decode to datetime64.
    Otherwise leave it alone.
    """
    if "time" not in ds.coords:
        raise ValueError("Dataset is missing time coordinate.")

    if np.issubdtype(ds["time"].dtype, np.datetime64):
        return ds

    time_attrs = ds["time"].attrs
    units = time_attrs.get("units")
    calendar = time_attrs.get("calendar", "standard")

    if units is None:
        return ds

    time_vals = ds["time"].values
    decoded = nc.num2date(time_vals, units=units, calendar=calendar, only_use_cftime_datetimes=False)
    ds = ds.assign_coords(time=np.array(decoded, dtype="datetime64[ns]"))
    return ds


def convert_gfs_style_to_processed(input_file: str, output_file: str) -> None:
    ds = xr.open_dataset(input_file, engine="netcdf4", decode_times=False)

    required_coords = {"time", "lev", "lat", "lon"}
    required_vars = {"hgtprs", "tmpprs", "ugrdprs", "vgrdprs"}

    missing_coords = required_coords - set(ds.coords)
    missing_vars = required_vars - set(ds.data_vars)

    if missing_coords:
        raise ValueError(f"Missing coords: {sorted(missing_coords)}")
    if missing_vars:
        raise ValueError(f"Missing vars: {sorted(missing_vars)}")

    # Decode numeric CF-style time to datetime64 if needed
    ds = decode_time_if_needed(ds)

    # Rename to processed-style names
    ds = ds.rename(
        {
            "lev": "level",
            "lat": "latitude",
            "lon": "longitude",
            "hgtprs": "z",
            "tmpprs": "t",
            "ugrdprs": "u",
            "vgrdprs": "v",
        }
    )

    # Convert lon from [0, 360) to [-180, 180)
    lon = ds["longitude"].values.astype(float)
    lon = ((lon + 180.0) % 360.0) - 180.0
    ds = ds.assign_coords(longitude=lon)
    ds = ds.sortby("longitude")

    # Match processed-file style:
    # - latitude descending
    # - level ascending
    ds = ds.sortby("latitude", ascending=False)
    ds = ds.sortby("level", ascending=True)

    # Standard dim order
    ds = ds.transpose("time", "level", "latitude", "longitude")

    # Variable order
    var_order = ["u", "v", "t", "z"]
    ds = xr.Dataset({v: ds[v] for v in var_order}, coords=ds.coords)

    # Cast to processed-style dtypes
    for v in var_order:
        ds[v] = ds[v].astype(np.float32)

    ds["level"] = ds["level"].astype(np.float64)
    ds["latitude"] = ds["latitude"].astype(np.float64)
    ds["longitude"] = ds["longitude"].astype(np.float64)

    # Minimal attrs
    ds.attrs = dict(ds.attrs)
    ds.attrs["Conventions"] = "CF-1.7"

    ds["level"].attrs.update({"units": "hPa", "long_name": "pressure level"})
    ds["latitude"].attrs.update({"units": "degrees_north"})
    ds["longitude"].attrs.update({"units": "degrees_east"})

    encoding = {var: {"zlib": True, "complevel": 4} for var in ds.data_vars}
    encoding["time"] = {"_FillValue": None}

    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out, engine="netcdf4", encoding=encoding)

    print(f"Saved: {out}")
    print(ds)
    print("Dimensions:", dict(ds.sizes))


if __name__ == "__main__":
    convert_gfs_style_to_processed(INPUT_FILE, OUTPUT_FILE)