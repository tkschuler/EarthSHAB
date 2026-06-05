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
import numpy as np
import warnings

from pathlib import Path
import xarray as xr
try:
    import cfgrib
except ImportError:
    cfgrib = None

# Suppress xarray FutureWarnings about compat defaults and other warnings
xr.set_options(use_new_combine_kwarg_defaults=True)
warnings.filterwarnings("ignore", category=FutureWarning)

# Configuration

# GFS horizontal resolution. Valid AWS products: "1p00" (1°), "0p50" (0.5°),
# "0p25" (0.25°). The global collector defaults to 1° (~2.6 GB / cycle) to keep
# whole-globe files manageable. "0p25" matches saveNETCDF.py fidelity but is
# ~16× larger and far slower to download — only practical for short horizons.
# Note: trajectories computed from a 1° file will differ from 0.25° runs because
# the reader (Forecast.py) samples the single nearest grid cell with no spatial
# interpolation, so a coarser grid means a more distant sample point.
RESOLUTION = "0p25"

BASE_URL = (
    "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.{date}/{hour}/atmos/"
    "gfs.t{hour}z.pgrb2." + RESOLUTION + ".f{forecast}"
)
SAVE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "GFS_DATA")  # Project root / GFS_DATA
GRIBS_DIR = os.path.join(SAVE_DIR, "GRIBS")
OUTPUT_DIR = os.path.join(SAVE_DIR, "NETCDF")

# Download entire 8-day GFS forecast available on AWS.
# AWS GFS: 3-hourly throughout [0, 3, 6, 9, ... 384h].
FORECAST_HOURS = list(range(0, 384, 3))

MAX_RETRIES = 3
TIMEOUT = 10
CHECK_INTERVAL = 5 * 60  # Check for new cycles every 5 minutes


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

# --- GRIB filename helpers ---
# Resolution is baked into the GRIB filename so that changing RESOLUTION does not
# silently reuse stale files of a different grid (download_gfs_file skips any file
# that already exists on disk).
def grib_prefix(date: str, hour: str) -> str:
    return f"gfs_{date}_{hour}z_{RESOLUTION}"

def grib_filename(date: str, hour: str, forecast: int) -> str:
    return f"{grib_prefix(date, hour)}_f{forecast:03d}.grib2"

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

def download_gfs_file(date, hour, forecast):
    """Download a GFS GRIB2 file with resume support. Skip missing files immediately."""
    url = BASE_URL.format(date=date, hour=hour, forecast=f"{forecast:03d}")
    filename = grib_filename(date, hour, forecast)
    file_path = os.path.join(GRIBS_DIR, filename)

    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        return file_path

    # Check if file exists on server (head request, no retry for 404)
    try:
        head_resp = SESSION.head(url, timeout=TIMEOUT)
        if head_resp.status_code == 404:
            return None  # File doesn't exist; skip silently
        if head_resp.status_code != 200:
            print(colored(f"[{filename}] HEAD returned {head_resp.status_code}; skipping","yellow"))
            return None
    except RequestException as e:
        print(colored(f"[{filename}] HEAD check failed: {e}; will retry on download","yellow"))

    # Download with retry for transient errors (not 404)
    for attempt in range(MAX_RETRIES):
        try:
            temp_file = file_path + ".part"
            resume_byte_pos = os.path.getsize(temp_file) if os.path.exists(temp_file) else 0
            headers = {"Range": f"bytes={resume_byte_pos}-"} if resume_byte_pos > 0 else {}

            response = requests.get(url, stream=True, headers=headers, timeout=TIMEOUT)

            if response.status_code == 404:
                return None  # Already checked above; shouldn't happen
            if response.status_code not in [200, 206]:
                print(colored(f"[{filename}] GET returned {response.status_code}","yellow"))
                if attempt < MAX_RETRIES - 1:
                    time.sleep(5)
                    continue
                else:
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

        except RequestException as e:
            print(colored(f"[{filename}] Download failed (attempt {attempt + 1}/{MAX_RETRIES}): {e}","yellow"))
            if attempt < MAX_RETRIES - 1:
                time.sleep(5)
            else:
                return None

    return None


def get_missing_files(date, hour, forecasts):
    """Check which GRIB2 files are missing."""
    missing_files = []
    for forecast in forecasts:
        filename = grib_filename(date, hour, forecast)
        file_path = os.path.join(GRIBS_DIR, filename)
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            missing_files.append(forecast)
    return missing_files


def cleanup_old_gribs(current_date, current_hour):
    """Remove old GRIB2 files, keeping only the latest cycle."""
    if not os.path.isdir(GRIBS_DIR):
        return
    prefix = grib_prefix(current_date, current_hour)
    for filename in os.listdir(GRIBS_DIR):
        if not filename.startswith(prefix):
            try:
                os.remove(os.path.join(GRIBS_DIR, filename))
            except OSError:
                pass


def main():
    """Continuously check for new GFS cycles and download → convert to v2 NetCDF.

    Resilient to crashes/stops: if all files are present but conversion incomplete,
    will retry conversion on next run.
    """
    last_downloaded_cycle = None

    while True:
        now = datetime.datetime.now(datetime.timezone.utc)
        date, hour = get_latest_cycle(now)

        # Download new cycle if detected
        if last_downloaded_cycle != (date, hour):
            print(f"New GFS cycle detected: {date} {hour}Z. Downloading...")

            # Download all forecast hours for this cycle
            for forecast in FORECAST_HOURS:
                download_gfs_file(date, hour, forecast)

            last_downloaded_cycle = (date, hour)

        # Always check if output exists; if not and we have files, try to convert
        # (handles resuming after crashes/stops)
        conversion_successful = False

        if has_output(date, hour):
            print(colored(f"[v2 NetCDF] Output {date} {hour}Z already exists — skipping rebuild.", "cyan"))
            conversion_successful = True  # Mark as successful since output exists
        else:
            # Check what files we have for this cycle
            missing_files = get_missing_files(date, hour, FORECAST_HOURS)
            downloaded_count = len(FORECAST_HOURS) - len(missing_files)

            if downloaded_count > 0:
                print(f"Converting {downloaded_count}/{len(FORECAST_HOURS)} forecast files...")
                if missing_files and missing_files[:5]:
                    print(colored(f"  (missing {len(missing_files)} files beyond AWS availability)", "yellow"))

                try:
                    combine_gfs_gribs_to_netcdf(date, hour)
                    print(colored("[v2 NetCDF] Output file ready.", "green"))
                    conversion_successful = True

                except Exception as e:
                    print(colored(f"[v2 NetCDF] Failed to build output file: {e}", "red"))
                    print(colored("  Will retry on next cycle check...", "yellow"))
            else:
                print(colored(f"[v2 NetCDF] No files for {date} {hour}Z yet; will retry...", "yellow"))

        # Only clean up old data if conversion succeeded
        if conversion_successful:
            cleanup_old_gribs(date, hour)
            cleanup_old_output(date, hour)

        print(f"Checking again in 5 minutes...")
        time.sleep(CHECK_INTERVAL)


def _grib_to_dataset(grib_path: str):
    """
    Open a GRIB2 file with cfgrib using filter_by_keys to extract only u, v, t, gh.
    Much faster than reading the whole file.
    """
    filters = [
        {"typeOfLevel": "isobaricInhPa", "shortName": "u"},
        {"typeOfLevel": "isobaricInhPa", "shortName": "v"},
        {"typeOfLevel": "isobaricInhPa", "shortName": "t"},
        {"typeOfLevel": "isobaricInhPa", "shortName": "gh"},
    ]

    ds_list = []
    for filt in filters:
        try:
            ds = xr.open_dataset(
                str(grib_path),
                engine="cfgrib",
                backend_kwargs={
                    "filter_by_keys": filt,
                    "indexpath": "",
                },
            )
            # Clear encoding IMMEDIATELY after opening to avoid conflicts later
            for var in list(ds.coords) + list(ds.data_vars):
                if var in ds.variables:
                    ds[var].encoding.clear()
            ds_list.append(ds)
        except Exception:
            pass  # Variable may not exist in this file

    if not ds_list:
        return None

    # Merge filtered datasets
    try:
        merged = xr.merge(ds_list, compat="override", join="override")
    except Exception:
        merged = ds_list[0]

    # Clear encoding again after merge (it can be reattached during merge)
    for var in list(merged.coords) + list(merged.data_vars):
        if var in merged.variables:
            merged[var].encoding.clear()

    # Add unified valid_time coordinate
    if "time" in merged and "step" in merged:
        valid_time = (merged["time"] + merged["step"]).astype("datetime64[ns]")
        merged = merged.assign_coords(valid_time=valid_time)
        merged = merged.drop_vars(["time", "step"], errors="ignore")
    elif "time" in merged:
        merged = merged.assign_coords(valid_time=merged["time"])
        merged = merged.drop_vars("time", errors="ignore")

    # Ensure valid_time is a dimension
    if "valid_time" in merged.coords and "valid_time" not in merged.dims:
        merged = merged.expand_dims("valid_time")

    return merged


def combine_gfs_gribs_to_netcdf(date: str, hour: str) -> str:
    """
    Combine all downloaded GRIBs into a v2-schema NetCDF file.

    v2 schema: dimensions (valid_time, pressure_level, latitude, longitude);
    variables (u, v, z, t); CF-1.7 compliant with proper metadata.

    Returns path to the produced NetCDF.
    """
    if cfgrib is None:
        raise RuntimeError("cfgrib is required but not installed. Install via conda: conda install -c conda-forge cfgrib eccodes")

    pattern = f"{grib_prefix(date, hour)}_f*.grib2"
    files = sorted(Path(GRIBS_DIR).glob(pattern))
    if not files:
        raise FileNotFoundError(f"No GRIB files found for cycle {date} {hour}Z")

    # Open and merge all GRIB files (fast with filter_by_keys)
    hourly_datasets = []
    for grib_path in tqdm(files, desc="Opening GRIBs", position=0, leave=True):
        ds = _grib_to_dataset(str(grib_path))
        if ds is not None:
            hourly_datasets.append(ds)

    if not hourly_datasets:
        raise RuntimeError("No datasets successfully opened from GRIB files.")

    # Concatenate all forecast hours
    ds = xr.concat(hourly_datasets, dim="valid_time")
    ds = ds.sortby("valid_time")

    # Clear all encoding metadata to avoid conflicts when setting attributes later
    for var in list(ds.coords) + list(ds.data_vars):
        if var in ds.variables:
            ds[var].encoding.clear()

    # Rename dimensions to v2 schema
    rename_map = {"isobaricInhPa": "pressure_level"}
    ds = ds.rename({k: v for k, v in rename_map.items() if k in ds})

    # Convert longitude from [0, 360) to [-180, 180)
    if "longitude" in ds.coords:
        lon_vals = ds.coords["longitude"].values
        lon_vals = np.where(lon_vals >= 180, lon_vals - 360, lon_vals)
        ds = ds.assign_coords(longitude=lon_vals)
        # Re-sort ascending: the conversion leaves longitude non-monotonic
        # ([0..179, -180..-1]). The v2 schema (and Forecast.py's documented
        # contract) require ascending longitude. Mirrors saveNETCDF.py.
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

    # Convert geopotential height (meters) to geopotential (m²/s²) for v2 schema
    g = 9.80665
    if "z" in ds:
        ds["z"] = ds["z"] * g

    # Set v2 schema attributes (clear encoding for each var right before setting attrs)
    for var in ["valid_time", "pressure_level", "latitude", "longitude"]:
        if var in ds.coords:
            ds[var].encoding.clear()

    if "valid_time" in ds.coords:
        # valid_time is a datetime64 coordinate. CF time units/calendar must go in
        # encoding (xarray writes them during serialization); putting them in attrs
        # collides with that and raises "Key 'units' already exists in attrs".
        # The actual units/calendar encoding is applied in the encoding dict at
        # write time (below), since encoding is cleared again before writing.
        ds["valid_time"].attrs = {
            "standard_name": "time",
            "long_name": "time",
        }
    if "pressure_level" in ds.coords:
        ds["pressure_level"].attrs = {
            "units": "hPa",
            "positive": "down",
            "stored_direction": "decreasing",
            "standard_name": "air_pressure",
            "long_name": "pressure",
        }
    if "latitude" in ds.coords:
        ds["latitude"].attrs = {
            "units": "degrees_north",
            "standard_name": "latitude",
            "long_name": "latitude",
            "stored_direction": "decreasing",
        }
    if "longitude" in ds.coords:
        ds["longitude"].attrs = {
            "units": "degrees_east",
            "standard_name": "longitude",
            "long_name": "longitude",
        }

    # Set variable attributes (clear encoding first)
    for var in ["u", "v", "z", "t"]:
        if var in ds:
            ds[var].encoding.clear()

    if "u" in ds:
        ds["u"].attrs = {
            "standard_name": "eastward_wind",
            "long_name": "U component of wind",
            "units": "m s**-1",
        }
    if "v" in ds:
        ds["v"].attrs = {
            "standard_name": "northward_wind",
            "long_name": "V component of wind",
            "units": "m s**-1",
        }
    if "z" in ds:
        ds["z"].attrs = {
            "standard_name": "geopotential",
            "long_name": "Geopotential",
            "units": "m**2 s**-2",
        }
    if "t" in ds:
        ds["t"].attrs = {
            "standard_name": "air_temperature",
            "long_name": "Temperature",
            "units": "K",
        }

    # Add global attributes (v2 schema requirement)
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    ds.attrs.update({
        "Conventions": "CF-1.7",
        "institution": "NOAA/NCEP (GFS)",
        "history": f"{now} gfs_collector.py",
    })

    # Clear all encoding to avoid conflicts with attributes
    for var in list(ds.coords) + list(ds.data_vars):
        if var in ds.variables:
            ds[var].encoding.clear()

    # Set compression encoding only (no other encoding that could conflict)
    encoding = {var: {"zlib": True, "complevel": 4, "shuffle": True} for var in ds.data_vars}

    # CF time encoding for the datetime coordinate. Lives in encoding (not attrs)
    # so xarray serializes valid_time as numeric seconds-since-epoch without colliding.
    if "valid_time" in ds.coords:
        encoding["valid_time"] = {
            "units": "seconds since 1970-01-01",
            "calendar": "proleptic_gregorian",
            "dtype": "int64",
        }

    # Write v2 NetCDF
    out_path = output_nc_path(date, hour)
    ds.to_netcdf(out_path, engine="netcdf4", encoding=encoding)
    ds.close()

    print(colored(f"[v2 NetCDF] Wrote v2-schema dataset → {out_path}", "green"))
    return out_path



if __name__ == "__main__":
    main()