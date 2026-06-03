# author: TKS
# date: Mar 02 2025
#
# GFS GRIB downloader → v2 canonical NetCDF converter.
#
# Downloads raw GRIB2 files from NOAA's AWS S3 bucket, then combines them into
# a single v2-schema NetCDF file. The v2 schema is the canonical EarthSHAB format:
# CF-1.7 compliant, with dimensions (valid_time, pressure_level, latitude, longitude)
# and variables (u, v, z, t).
#
# GFS Update Cycle:
#   Frequency: 4 times/day, every 6 hours starting at 00:00 UTC
#   Upload delay: typically 3-4 hours after cycle time
#   Forecast range: up to 8 days ahead
#   Temporal resolution:
#     - First 12 hours: 2-hour frequency
#     - 12 hours–2 days: 3-hour frequency
#     - Days 2–8: 6-hour frequency
#
# Script workflow:
#   1. Check for latest available GFS cycle on NOAA's server
#   2. Download missing GRIB2 files for that cycle
#   3. Combine downloaded GRIBs into a single v2 NetCDF file
#   4. Clean up old GRIB2 files and previous-cycle NetCDF files

import os
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from requests.exceptions import RequestException
import time
import datetime
from tqdm import tqdm
from termcolor import colored
import functools

from pathlib import Path
import xarray as xr

# Configuration
BASE_URL = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.{date}/{hour}/atmos/gfs.t{hour}z.pgrb2.1p00.f{forecast}"
SAVE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "GFS_DATA")  # Project root / GFS_DATA
GRIBS_DIR = os.path.join(SAVE_DIR, "GRIBS")
OUTPUT_DIR = os.path.join(SAVE_DIR, "NETCDF")
FORECAST_HOURS = list(range(0, 8*24, 3))  # Up to 8 days ahead, 3-hour resolution
MAX_RETRIES = 3
TIMEOUT = 10
CHECK_INTERVAL = 5 * 60  # Check for new cycles every 5 minutes

# Which fields to extract from the GRIBs → tweak as needed
# Each dict is passed to cfgrib's filter_by_keys
GRIB_FILTERS = [
    # Pressure-level winds & temperature (example)
    {"typeOfLevel": "isobaricInhPa", "shortName": "u"},
    {"typeOfLevel": "isobaricInhPa", "shortName": "v"},
    {"typeOfLevel": "isobaricInhPa", "shortName": "t"},
    # Geopotential height on pressure levels (sometimes 'gh' or 'z', GFS uses 'gh')
    {"typeOfLevel": "isobaricInhPa", "shortName": "gh"},

    # Surface 10m winds (uncomment if you want them)
    # {"typeOfLevel": "surface", "shortName": "10u"},
    # {"typeOfLevel": "surface", "shortName": "10v"},
]


def make_session() -> requests.Session:
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=0.8,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("HEAD", "GET"),
        raise_on_status=False,
    )
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.mount("http://", HTTPAdapter(max_retries=retry))
    return s

SESSION = make_session()

# Create output directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(GRIBS_DIR, exist_ok=True)

# --- v2 NetCDF output helpers ---
def output_nc_path(date: str, hour: str) -> str:
    return os.path.join(OUTPUT_DIR, f"gfs_{date}_{hour}z.nc")

def has_output(date: str, hour: str) -> bool:
    p = output_nc_path(date, hour)
    return os.path.exists(p) and os.path.getsize(p) > 0

def cleanup_old_output(current_date: str, current_hour: str):
    if not os.path.isdir(OUTPUT_DIR):
        return
    prefix = f"gfs_{current_date}_{current_hour}z"
    for fn in os.listdir(OUTPUT_DIR):
        if fn.endswith(".nc") and not fn.startswith(prefix):
            try:
                os.remove(os.path.join(OUTPUT_DIR, fn))
            except OSError:
                pass
            


def retry(max_retries=3, delay=5, exceptions=(requests.exceptions.RequestException,)):
    """Decorator to automatically retry failed network operations."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    print(colored(f"WARNING: Attempt {attempt + 1}/{max_retries} failed: {e}. Retrying in {delay} sec...","yellow"))
                    time.sleep(delay)
            print(colored(f"ERROR: {func.__name__} failed after {max_retries} retries.","red"))
            return None
        return wrapper
    return decorator

def is_gfs_cycle_available(date, hour):
    """Check if the GFS data for a given cycle exists on NOAA's server yet."""
    url = BASE_URL.format(date=date, hour=hour, forecast="000")  # dummy forecast hour
    try:
        response = SESSION.head(url, timeout=TIMEOUT, allow_redirects=True)
        return response.status_code == 200
    except RequestException as e:
        print(colored(f"WARNING: HEAD check failed for {url}: {e}", "yellow"))
        return False


def get_latest_cycle(now):
    """Get the latest available GFS cycle based on current UTC time."""
    print(f"Current UTC Time: {now}")

    # Round down current time to the nearest GFS cycle
    gfs_hours = [18, 12, 6, 0]
    nearest_hour = max([h for h in gfs_hours if h <= now.hour])

    # Create a timestamp for the nearest cycle today
    cycle_time = now.replace(hour=nearest_hour, minute=0, second=0, microsecond=0)

    # Format and return cycle time
    date = cycle_time.strftime("%Y%m%d")
    hour = f"{cycle_time.hour:02d}"

    # Try to find the latest available cycle
    # Make sure latest cycle is available (3-4 hour delay).  Pretty sure this should never be more than 6 hours. So just 1 check
    if is_gfs_cycle_available(date, hour):
        return date, hour

    # If not found, go back one cycle (6 hours)
    print(colored(f"WARNING: Cycle time {cycle_time} not uploaded yet. Going back 6 hours", "yellow"))
    cycle_time -= datetime.timedelta(hours=6)
    date = cycle_time.strftime("%Y%m%d")
    hour = f"{cycle_time.hour:02d}"

    return date, hour

@retry(max_retries=3, delay=5)
def download_gfs_file(date, hour, forecast):
    """Download a GFS GRIB2 file with resume support."""
    url = BASE_URL.format(date=date, hour=hour, forecast=f"{forecast:03d}")
    filename = f"gfs_{date}_{hour}z_f{forecast:03d}.grib2"
    file_path = os.path.join(GRIBS_DIR, filename)

    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        print(f"{filename} already exists. Skipping download.")
        return file_path

    temp_file = file_path + ".part"
    resume_byte_pos = os.path.getsize(temp_file) if os.path.exists(temp_file) else 0
    headers = {"Range": f"bytes={resume_byte_pos}-"} if resume_byte_pos > 0 else {}

    print(f"Downloading {filename} from {url}... (resuming at {resume_byte_pos} bytes)")

    response = requests.get(url, stream=True, headers=headers, timeout=TIMEOUT)

    if response.status_code == 404:
        print(colored(f"File not found: {url}","yellow"))
        return None

    if response.status_code not in [200, 206]:
        print(colored(f"Unexpected status code {response.status_code} for {filename}","red"))
        return None

    total_size = int(response.headers.get("content-length", 0)) + resume_byte_pos
    block_size = 8192
    progress_bar = tqdm(total=total_size, unit="B", unit_scale=True, desc=filename, initial=resume_byte_pos)

    with open(temp_file, "ab") as file:
        for chunk in response.iter_content(chunk_size=block_size):
            if chunk:
                file.write(chunk)
                progress_bar.update(len(chunk))

    progress_bar.close()
    os.rename(temp_file, file_path)
    print(colored(f"Download complete: {file_path}","green"))
    return file_path


def get_missing_files(date, hour, forecasts):
    """Check which GRIB2 files are missing."""
    missing_files = []
    for forecast in forecasts:
        filename = f"gfs_{date}_{hour}z_f{forecast:03d}.grib2"
        file_path = os.path.join(GRIBS_DIR, filename)
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            missing_files.append(forecast)
    return missing_files


def cleanup_old_gribs(current_date, current_hour):
    """Remove old GRIB2 files, keeping only the latest cycle."""
    if not os.path.isdir(GRIBS_DIR):
        return
    prefix = f"gfs_{current_date}_{current_hour}z"
    for filename in os.listdir(GRIBS_DIR):
        if not filename.startswith(prefix):
            try:
                os.remove(os.path.join(GRIBS_DIR, filename))
            except OSError:
                pass


def main():
    """Continuously check for new GFS cycles and download → convert to v2 NetCDF."""
    last_checked_cycle = None

    while True:
        now = datetime.datetime.now(datetime.timezone.utc)
        date, hour = get_latest_cycle(now)

        if last_checked_cycle != (date, hour):
            print(f"New GFS cycle detected: {date} {hour}Z. Downloading...")

            # Download all forecast hours for this cycle
            for forecast in FORECAST_HOURS:
                download_gfs_file(date, hour, forecast)

            # Verify all files downloaded
            missing_files = get_missing_files(date, hour, FORECAST_HOURS)

            if not missing_files:
                print(f"All GFS files downloaded for {date} {hour}Z!")

                # Build v2 NetCDF only if it doesn't already exist
                if has_output(date, hour):
                    print(colored("[v2 NetCDF] Output file already exists — skipping rebuild.", "cyan"))
                else:
                    try:
                        combine_gfs_gribs_to_netcdf(date, hour)
                    except Exception as e:
                        print(colored(f"[v2 NetCDF] Failed to build output file: {e}", "red"))
                    else:
                        print(colored("[v2 NetCDF] Output file ready.", "green"))

                # Update cycle marker and clean old GRIB2 files
                last_checked_cycle = (date, hour)
                cleanup_old_gribs(date, hour)
                cleanup_old_output(date, hour)

                print(f"Checking again in 5 minutes...")

        else:
            print(f"No new data available. Checking again in 5 minutes...")

        time.sleep(CHECK_INTERVAL)


def _open_grib_filtered(grib_path: str, filter_keys: dict):
    """
    Open one GRIB file as an xarray.Dataset filtered to a single 'group' of fields
    (e.g., one shortName + level type). Returns None if that combination isn't present.
    """
    try:
        ds = xr.open_dataset(
            grib_path,
            engine="cfgrib",
            backend_kwargs={
                "filter_by_keys": filter_keys,
                # Avoid creating .idx files next to each GRIB
                "indexpath": "",
            },
        )
        # Add/normalize a unified valid_time coordinate for easy concat
        # GFS files typically have the analysis 'time' and a forecast 'step'
        if "time" in ds and "step" in ds:
            valid_time = (ds["time"] + ds["step"]).astype("datetime64[ns]")
            # Make valid_time a dimension for concatenation
            ds = ds.assign_coords(valid_time=valid_time)
            # Some fields keep 'step' as a length-1 dimension; squeeze it out
            for dim in list(ds.dims):
                if ds.sizes[dim] == 1 and dim != "isobaricInhPa":
                    ds = ds.squeeze(dim, drop=True)
            # Ensure valid_time is a dimension we can concat over
            if "valid_time" not in ds.dims:
                ds = ds.expand_dims("valid_time")
        elif "time" in ds:
            # Fallback if step is missing for any reason
            ds = ds.assign_coords(valid_time=ds["time"])
            if "valid_time" not in ds.dims:
                ds = ds.expand_dims("valid_time")

        # Drop now-redundant coords if present (optional)
        for maybe in ("time", "step"):
            if maybe in ds:
                try:
                    ds = ds.drop_vars(maybe)
                except Exception:
                    pass

        return ds
    except Exception as e:
        # Not all filter combos exist in every GRIB; silence those
        # but keep truly unexpected errors visible during debugging if you want.
        # print(f"[DEBUG] {grib_path} missing {filter_keys}: {e}")
        return None


def combine_gfs_gribs_to_netcdf(date: str, hour: str, filters: list[dict] = None) -> str:
    """
    Combine all downloaded GRIBs into a v2-schema NetCDF file.

    v2 schema: dimensions (valid_time, pressure_level, latitude, longitude);
    variables (u, v, z, t); CF-1.7 compliant with proper metadata.

    Returns path to the produced NetCDF.
    """
    filters = filters or GRIB_FILTERS
    pattern = f"gfs_{date}_{hour}z_f*.grib2"
    files = sorted(Path(GRIBS_DIR).glob(pattern))
    if not files:
        raise FileNotFoundError(f"No GRIB files found for cycle {date} {hour}Z")

    # Collect datasets per filter (one per variable + level type)
    per_filter_datasets = []
    for filt in filters:
        ds_list = []
        for grib in files:
            ds = _open_grib_filtered(str(grib), filt)
            if ds is not None:
                ds_list.append(ds)

        if ds_list:
            ds_cat = xr.concat(ds_list, dim="valid_time")
            ds_cat = ds_cat.sortby("valid_time")
            per_filter_datasets.append(ds_cat)

    if not per_filter_datasets:
        raise RuntimeError("No datasets opened from GRIB filters.")

    # Merge all variables into one dataset
    ds = xr.merge(per_filter_datasets, compat="no_conflicts")

    # Rename dimensions to v2 schema
    rename_map = {"isobaricInhPa": "pressure_level", "latitude": "latitude", "longitude": "longitude"}
    ds = ds.rename(rename_map)

    # Convert longitude from [0, 360) to [-180, 180)
    if "longitude" in ds.coords:
        lon = ds.coords["longitude"].values
        lon = ((lon + 180) % 360) - 180
        ds = ds.assign_coords(longitude=("longitude", lon))
        ds = ds.sortby("longitude")

    # Sort latitude descending (north to south) for v2 schema
    if "latitude" in ds.coords:
        ds = ds.sortby("latitude", ascending=False)

    # Ensure pressure_level is descending (high hPa → low hPa)
    if "pressure_level" in ds.coords:
        ds = ds.sortby("pressure_level", ascending=False)

    # Standardize variable names to v2: u, v, z, t
    var_map = {}
    for var in ds.data_vars:
        if var in ("u", "v", "z", "t"):
            continue  # Already correct
        # Map common GFS/GRIB names
        if "ugrd" in var.lower() or var.lower() == "10u":
            var_map[var] = "u"
        elif "vgrd" in var.lower() or var.lower() == "10v":
            var_map[var] = "v"
        elif any(x in var.lower() for x in ("hgt", "gh", "z")):
            var_map[var] = "z"
        elif "tmp" in var.lower() or "t" == var.lower():
            var_map[var] = "t"
    if var_map:
        ds = ds.rename(var_map)

    # Add/fix v2 coordinate attributes
    if "valid_time" in ds.coords:
        ds["valid_time"].attrs.update({
            "units": "seconds since 1970-01-01",
            "calendar": "proleptic_gregorian",
            "standard_name": "time",
            "long_name": "time",
        })

    if "pressure_level" in ds.coords:
        ds["pressure_level"].attrs.update({
            "units": "hPa",
            "positive": "down",
            "stored_direction": "decreasing",
            "standard_name": "air_pressure",
            "long_name": "pressure",
        })

    if "latitude" in ds.coords:
        ds["latitude"].attrs.update({
            "units": "degrees_north",
            "standard_name": "latitude",
            "long_name": "latitude",
            "stored_direction": "decreasing",
        })

    if "longitude" in ds.coords:
        ds["longitude"].attrs.update({
            "units": "degrees_east",
            "standard_name": "longitude",
            "long_name": "longitude",
        })

    # Add v2 variable metadata
    if "u" in ds:
        ds["u"].attrs.update({
            "standard_name": "eastward_wind",
            "long_name": "U component of wind",
            "units": "m s**-1",
        })
    if "v" in ds:
        ds["v"].attrs.update({
            "standard_name": "northward_wind",
            "long_name": "V component of wind",
            "units": "m s**-1",
        })
    if "z" in ds:
        ds["z"].attrs.update({
            "standard_name": "geopotential",
            "long_name": "Geopotential",
            "units": "m**2 s**-2",
        })
    if "t" in ds:
        ds["t"].attrs.update({
            "standard_name": "air_temperature",
            "long_name": "Temperature",
            "units": "K",
        })

    # Add global attributes (v2 schema requirement)
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    ds.attrs.update({
        "Conventions": "CF-1.7",
        "institution": "NOAA/NCEP (GFS)",
        "history": f"{now} gfs_collector.py",
    })

    # Compression
    encoding = {var: {"zlib": True, "complevel": 4, "shuffle": True} for var in ds.data_vars}

    # Write v2 NetCDF
    out_path = output_nc_path(date, hour)
    ds.to_netcdf(out_path, engine="netcdf4", encoding=encoding)
    ds.close()

    print(colored(f"[v2 NetCDF] Wrote v2-schema dataset → {out_path}", "green"))
    return out_path



if __name__ == "__main__":
    main()