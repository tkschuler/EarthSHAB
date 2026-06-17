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
import argparse
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
import netCDF4 as nc
try:
    import cfgrib
except ImportError:
    cfgrib = None

# Suppress xarray FutureWarnings about compat defaults and other warnings
xr.set_options(use_new_combine_kwarg_defaults=True)
warnings.filterwarnings("ignore", category=FutureWarning)

# Configuration

# Unlike saveNETCDF.py (which honours the netcdf_gfs config knobs for a single
# requested cycle), the collector is a long-running background downloader, so it
# defaults to the cheapest sensible grid/cadence — 1.0° at 3-hourly steps — and
# only deviates when explicitly asked via CLI args (--res / --step-hours). This
# keeps an always-on station from silently pulling ~16× more data than intended.
# Valid AWS products: "1p00" (1°), "0p50" (0.5°), "0p25" (0.25°). Note:
# trajectories from a 1° file differ from 0.25° runs because the reader
# (Forecast.py) samples the single nearest grid cell with no spatial
# interpolation, so a coarser grid means a more distant sample point.
DEFAULT_RES = 1.0          # degrees: 0.25, 0.5, or 1.0
DEFAULT_STEP_HOURS = 3     # hours between forecast steps

SAVE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "GFS_DATA")  # Project root / GFS_DATA
GRIBS_DIR = os.path.join(SAVE_DIR, "GRIBS")
OUTPUT_DIR = os.path.join(SAVE_DIR, "NETCDF")

# These are (re)assigned by configure(); module-level defaults let the module be
# imported without first parsing CLI args.
RESOLUTION = None    # e.g. "1p00" — used in URLs and filenames
BASE_URL = None
STEP_HOURS = None    # (h) temporal step
FORECAST_HOURS = None


def configure(res: float = DEFAULT_RES, step_hours: int = DEFAULT_STEP_HOURS):
    """Set the grid/cadence-derived globals from a resolution and temporal step.

    Recomputes RESOLUTION, BASE_URL, STEP_HOURS, and FORECAST_HOURS. Called once
    at import with the defaults and again from main() with any CLI overrides;
    every downstream function reads these globals at call time, so reassigning
    them here is sufficient.
    """
    global RESOLUTION, BASE_URL, STEP_HOURS, FORECAST_HOURS

    if res not in (0.25, 0.5, 1.0):
        raise ValueError(f"res must be one of 0.25, 0.5, 1.0 (got {res})")
    if step_hours < 1:
        raise ValueError(f"step_hours must be a positive integer (got {step_hours})")

    RESOLUTION = ("%.2f" % res).replace(".", "p")  # 0.25->0p25, 0.5->0p50, 1.0->1p00
    BASE_URL = (
        "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.{date}/{hour}/atmos/"
        "gfs.t{hour}z.pgrb2." + RESOLUTION + ".f{forecast}"
    )
    STEP_HOURS = int(step_hours)
    # Forecast horizon (full 8-day GFS available on AWS) at the chosen step.
    # AWS GFS is 3-hourly throughout [0..384h]; 0.5°/1.0° grids have no hourly
    # steps, and 0.25° hourly only exists to f120 — missing steps are skipped on
    # download (see is_step_available()).
    FORECAST_HOURS = list(range(0, 384, STEP_HOURS))


# Initialise globals with the defaults so plain `import gfs_collector` works.
configure()


def is_step_available(forecast: int) -> bool:
    """Whether NOAA's AWS bucket ever publishes this forecast hour for the grid.

    All grids carry 3-hourly steps over [0, 384]. The 0.25° grid additionally
    carries hourly steps over [0, 120]. Any other step 404s forever, so it must
    be excluded from the completeness check — otherwise the NetCDF would never be
    considered "complete" and would never build.
    """
    if forecast < 0 or forecast > 384:
        return False
    if forecast % 3 == 0:
        return True
    if RESOLUTION == "0p25" and forecast <= 120:
        return True
    return False


def expected_forecasts():
    """Forecast hours in FORECAST_HOURS that actually exist on the server.

    These are the GRIBs that must ALL be on disk before we build the NetCDF and
    replace the previous forecast. A step that is merely not-yet-uploaded (the
    cycle publishes f000→f384 progressively) is still "expected" and keeps the
    cycle incomplete until it lands.
    """
    return [f for f in FORECAST_HOURS if is_step_available(f)]

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
# Filename encodes the grid resolution and temporal step, e.g.
# gfs_20260605_06z_0p25_3h.nc, so different grids/cadences don't collide.
def output_nc_path(date: str, hour: str) -> str:
    return os.path.join(OUTPUT_DIR, f"gfs_{date}_{hour}z_{RESOLUTION}_{STEP_HOURS}h.nc")

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

    Resilient to crashes/stops: each loop re-attempts any missing GRIBs (already
    downloaded files are skipped cheaply) and only builds the NetCDF once every
    server-available forecast hour for the cycle is on disk.

    Critically, the previous forecast (old GRIBs + old NetCDF) is only deleted
    AFTER a complete cycle has been successfully converted. A partially-uploaded
    or partially-downloaded cycle never triggers a build or a cleanup, so we
    never replace a good forecast with a truncated one.
    """
    while True:
        now = datetime.datetime.now(datetime.timezone.utc)
        date, hour = get_latest_cycle(now)

        # If we already built this cycle's output, nothing to do but housekeeping.
        if has_output(date, hour):
            print(colored(f"[v2 NetCDF] Output {date} {hour}Z already exists — skipping rebuild.", "cyan"))
            cleanup_old_gribs(date, hour)
            cleanup_old_output(date, hour)
            print(f"Checking again in 5 minutes...")
            time.sleep(CHECK_INTERVAL)
            continue

        # Attempt to download every requested forecast hour. Existing files are
        # returned immediately; off-grid steps 404 and are skipped; not-yet-
        # uploaded steps will simply be retried on the next loop.
        print(f"Cycle {date} {hour}Z: downloading any missing GRIBs...")
        for forecast in FORECAST_HOURS:
            download_gfs_file(date, hour, forecast)

        # Completeness is measured against the steps the server actually carries,
        # not all of FORECAST_HOURS (many hourly steps 404 by design).
        expected = expected_forecasts()
        missing = get_missing_files(date, hour, expected)
        have = len(expected) - len(missing)

        if missing:
            # Cycle not fully available yet (still uploading) or some downloads
            # failed. Do NOT build and do NOT delete the previous forecast —
            # wait and retry. This is the guard against converting after only a
            # handful of files have arrived.
            print(colored(
                f"[v2 NetCDF] Cycle {date} {hour}Z incomplete: {have}/{len(expected)} "
                f"expected GRIBs present, {len(missing)} still missing. "
                f"Holding off on build (keeping previous forecast).", "yellow"))
            print(f"Checking again in 5 minutes...")
            time.sleep(CHECK_INTERVAL)
            continue

        # All expected GRIBs are present → safe to build and replace the old one.
        print(f"All {len(expected)} expected GRIBs present. Building v2 NetCDF...")
        try:
            combine_gfs_gribs_to_netcdf(date, hour)
            print(colored("[v2 NetCDF] Output file ready.", "green"))
            cleanup_old_gribs(date, hour)
            cleanup_old_output(date, hour)
        except Exception as e:
            print(colored(f"[v2 NetCDF] Failed to build output file: {e}", "red"))
            print(colored("  Will retry on next cycle check...", "yellow"))

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


# v2 data-variable metadata: name → (standard_name, long_name, units).
_VAR_META = {
    "u": ("eastward_wind", "U component of wind", "m s**-1"),
    "v": ("northward_wind", "V component of wind", "m s**-1"),
    "z": ("geopotential", "Geopotential", "m**2 s**-2"),
    "t": ("air_temperature", "Temperature", "K"),
}


def _epoch_seconds(dt64) -> int:
    """Convert a numpy datetime64 to integer seconds since the Unix epoch."""
    return int(np.datetime64(dt64, "s").astype("int64"))


def _normalize_step(ds):
    """Apply the v2 coordinate/variable conventions to one forecast-hour dataset.

    Renames the level dim, converts longitude to [-180, 180) (restoring ascending
    order), sorts latitude/pressure descending, maps GRIB names to (u, v, z, t),
    and converts geopotential height to geopotential. Returns the dataset with
    dims (valid_time, pressure_level, latitude, longitude).
    """
    g = 9.80665

    if "isobaricInhPa" in ds:
        ds = ds.rename({"isobaricInhPa": "pressure_level"})

    # Map GRIB variable names to the v2 set.
    var_map = {}
    for var in list(ds.data_vars):
        if var in ("u", "v", "z", "t"):
            continue
        low = var.lower()
        if "ugrd" in low or low == "10u":
            var_map[var] = "u"
        elif "vgrd" in low or low == "10v":
            var_map[var] = "v"
        elif "hgt" in low or "gh" in low or low == "z":
            var_map[var] = "z"
        elif "tmp" in low or low == "t":
            var_map[var] = "t"
    if var_map:
        ds = ds.rename(var_map)

    # Convert longitude [0, 360) → [-180, 180) and restore ascending order.
    if "longitude" in ds.coords:
        lon_vals = ds.coords["longitude"].values
        lon_vals = np.where(lon_vals >= 180, lon_vals - 360, lon_vals)
        ds = ds.assign_coords(longitude=lon_vals).sortby("longitude")
    if "latitude" in ds.coords:
        ds = ds.sortby("latitude", ascending=False)
    if "pressure_level" in ds.coords:
        ds = ds.sortby("pressure_level", ascending=False)

    # Geopotential height (m) → geopotential (m² s⁻²); keep everything float32.
    if "z" in ds:
        ds["z"] = (ds["z"] * g).astype("float32")
    for name in ("u", "v", "t"):
        if name in ds:
            ds[name] = ds[name].astype("float32")

    return ds


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

    # ── Stream each forecast hour straight to disk ──────────────────────────
    # Earlier versions opened every forecast hour, xr.concat'd them into one
    # array, and wrote at the end. For a global 0.25° cycle that array is ~66 GB
    # (and decoding/concat copies push the peak past that), which OOM-kills the
    # process. Instead we create the NetCDF from the first GRIB and append one
    # valid_time slice per file, so peak memory is a single step (~0.5 GB at
    # 0.25°) regardless of how many forecast hours or what resolution.
    out_path = output_nc_path(date, hour)
    tmp_path = out_path + ".tmp"
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()

    ncf = None
    data_vars = {}
    vt = None
    ref_shape = None
    written = 0
    try:
        for grib_path in tqdm(files, desc="Opening GRIBs", position=0, leave=True):
            ds = _grib_to_dataset(str(grib_path))
            if ds is None:
                continue
            ds = _normalize_step(ds)

            lev = ds["pressure_level"].values
            lat = ds["latitude"].values
            lon = ds["longitude"].values

            if ncf is None:
                # Build the file skeleton from the first step's grid.
                ncf = nc.Dataset(tmp_path, "w", format="NETCDF4")
                ncf.createDimension("valid_time", None)   # unlimited: append per step
                ncf.createDimension("pressure_level", lev.size)
                ncf.createDimension("latitude", lat.size)
                ncf.createDimension("longitude", lon.size)

                vt = ncf.createVariable("valid_time", "i8", ("valid_time",))
                vt.units = "seconds since 1970-01-01"
                vt.calendar = "proleptic_gregorian"
                vt.standard_name = "time"
                vt.long_name = "time"

                pl = ncf.createVariable("pressure_level", "f8", ("pressure_level",))
                pl.units = "hPa"
                pl.positive = "down"
                pl.stored_direction = "decreasing"
                pl.standard_name = "air_pressure"
                pl.long_name = "pressure"
                pl[:] = lev

                la = ncf.createVariable("latitude", "f8", ("latitude",))
                la.units = "degrees_north"
                la.standard_name = "latitude"
                la.long_name = "latitude"
                la.stored_direction = "decreasing"
                la[:] = lat

                lo = ncf.createVariable("longitude", "f8", ("longitude",))
                lo.units = "degrees_east"
                lo.standard_name = "longitude"
                lo.long_name = "longitude"
                lo[:] = lon

                chunks = (1, lev.size, lat.size, lon.size)
                for name, (sn, ln, units) in _VAR_META.items():
                    v = ncf.createVariable(
                        name, "f4",
                        ("valid_time", "pressure_level", "latitude", "longitude"),
                        zlib=True, complevel=4, shuffle=True, chunksizes=chunks,
                    )
                    v.standard_name = sn
                    v.long_name = ln
                    v.units = units
                    data_vars[name] = v

                ncf.Conventions = "CF-1.7"
                ncf.institution = "NOAA/NCEP (GFS)"
                ncf.history = f"{now} gfs_collector.py"

                ref_shape = (lev.size, lat.size, lon.size)
            elif (lev.size, lat.size, lon.size) != ref_shape:
                # Grid changed mid-cycle (should never happen) — skip the oddball.
                print(colored(f"  [!] {grib_path.name}: grid {(lev.size, lat.size, lon.size)} "
                              f"!= {ref_shape}; skipping", "yellow"))
                ds.close()
                continue

            # Append this step's slice. .values[0] drops the size-1 valid_time axis.
            vt[written] = _epoch_seconds(ds["valid_time"].values[0])
            for name, v in data_vars.items():
                if name in ds:
                    v[written] = ds[name].values[0]
            written += 1
            ds.close()

        if ncf is None or written == 0:
            raise RuntimeError("No datasets successfully opened from GRIB files.")
    finally:
        if ncf is not None:
            ncf.close()

    os.replace(tmp_path, out_path)
    print(colored(f"[v2 NetCDF] Wrote v2-schema dataset ({written} steps) → {out_path}", "green"))
    return out_path



def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="EarthSHAB.gfs_collector",
        description="Continuously download the latest GFS cycle and build a v2 "
                    "NetCDF. Defaults to a 1.0° grid at 3-hourly steps.",
    )
    p.add_argument(
        "--res", type=float, default=DEFAULT_RES, choices=(0.25, 0.5, 1.0),
        help="Grid resolution in degrees (default: %(default)s).",
    )
    p.add_argument(
        "--step-hours", type=int, default=DEFAULT_STEP_HOURS, metavar="H",
        help="Hours between forecast steps (default: %(default)s). Note: only "
             "the 0.25° grid carries sub-3-hourly steps, and only to f120.",
    )
    return p.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    configure(res=args.res, step_hours=args.step_hours)
    print(f"GFS collector: {RESOLUTION} grid, {STEP_HOURS}-hourly steps")
    main()