"""
saveNETCDF.py  –  EarthSHAB GFS forecast downloader (GRIB-filter edition)
==========================================================================
NOAA retired the OpenDAP/DODS interface that EarthSHAB originally used.
This replacement fetches the same data through NOAA's GRIB-filter service
(https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl), saves each
forecast hour as a temporary GRIB2 file, and merges everything into a
single NetCDF-4 file whose structure matches what GFS.py already expects.

Variables downloaded (all on isobaric pressure levels):
  UGRD  – U-component of wind  (m/s)
  VGRD  – V-component of wind  (m/s)
  TMP   – Temperature          (K)
  HGT   – Geopotential height  (gpm)

Requirements (add to requirements.txt):
  cfgrib        – pip install cfgrib
  eccodes       – conda install -c conda-forge eccodes   (or apt/brew)
  xarray        – already in EarthSHAB requirements
  netCDF4       – already in EarthSHAB requirements
  requests      – usually available; pip install requests

Usage:
  python saveNETCDF.py
  (uses config_earth.py for lat/lon centre, range, date, and download_days)
"""

import sys
import time
import datetime
import tempfile
import requests
import numpy as np
import xarray as xr
try:
    import cfgrib
except ImportError:  # cfgrib needs eccodes (conda); keep the module importable without it
    cfgrib = None
import netCDF4 as nc
from pathlib import Path
from termcolor import colored
import pandas as pd

# ── EarthSHAB config ──────────────────────────────────────────────────────────
from EarthSHAB.config_earth import netcdf_gfs, simulation

xr.set_options(use_new_combine_kwarg_defaults=True) # cfgrib may raise FutureWarnings about combine_attrs; this silences them for now...

# ── Pressure levels available in GFS 0.25° pgrb2 isobaric data ───────────────
PRESSURE_LEVELS_MB = [
    1000, 975, 950, 925, 900, 850, 800, 750, 700, 650,
     600, 550, 500, 450, 400, 350, 300, 250, 200, 150,
     100,  70,  50,  40,  30,  20,  15,  10,
]

# ── GRIB filter base URL ───────────────────────────────────────────────────────
FILTER_BASE = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl"

# Seconds to wait between fetches (NOAA asks for ≥10 s between requests)
FETCH_DELAY = 2


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _level_params(levels_mb):
    """Build the query-string fragments that turn on each pressure level."""
    return "&".join(f"lev_{p}_mb=on" for p in levels_mb)


def _build_url(date_str, cycle_hour, forecast_hour, lat_range, lon_range,
               center_lat, center_lon):
    """
    Construct a GRIB-filter URL for one GFS forecast step.

    date_str      : 'YYYYMMDD'
    cycle_hour    : 0, 6, 12, or 18  (model initialisation hour)
    forecast_hour : 0, 3, 6, … 240+  (hours into forecast)
    lat/lon_range : index counts (each = res degrees) from config
    center_lat/lon: centre of the bounding box; lon in [-180, 180]
    NOTE: GRIB filter requires longitudes in [-180, 180], NOT [0, 360].
    """
    fhh = f"{int(forecast_hour):03d}"
    cyc = f"{int(cycle_hour):02d}"

    # Bounding box – clamp to valid globe extents and round to 4 dp
    top    = round(min( 90.0, center_lat + lat_range / 2.0), 4)
    bottom = round(max(-90.0, center_lat - lat_range / 2.0), 4)
    left   = round(max(-180.0, center_lon - lon_range / 2.0), 4)
    right  = round(min( 180.0, center_lon + lon_range / 2.0), 4)

    level_str = _level_params(PRESSURE_LEVELS_MB)

    params = (
        f"file=gfs.t{cyc}z.pgrb2.0p25.f{fhh}"
        f"&{level_str}"
        f"&var_UGRD=on&var_VGRD=on&var_TMP=on&var_HGT=on"
        f"&subregion="
        f"&toplat={top}&leftlon={left}&rightlon={right}&bottomlat={bottom}"
        f"&dir=%2Fgfs.{date_str}%2F{cyc}%2Fatmos"
    )
    return f"{FILTER_BASE}?{params}"


def _adjust_run_date(start_time):
    """
    Adjust the run_date to match an available forecast cycle on NOAA's server.
    """
    # GFS model runs at 00, 06, 12, and 18 UTC
    cycle_hour = (start_time.hour // 6) * 6  # Nearest GFS cycle (00, 06, 12, 18)
    run_date = start_time.replace(hour=cycle_hour, minute=0, second=0, microsecond=0)

    # Check if the calculated run_date is within the NOAA retention window
    now_utc = datetime.datetime.utcnow()
    GFS_RETENTION_DAYS = 9
    oldest_available = (now_utc - datetime.timedelta(days=GFS_RETENTION_DAYS)).replace(
        hour=0, minute=0, second=0, microsecond=0
    )

    if run_date < oldest_available:
        print(colored("\n[WARNING] Calculated run_date is too far in the past for available forecasts.", "yellow"))
        print(colored(f"  Adjusting run_date to the oldest available forecast: {oldest_available}", "yellow"))
        run_date = oldest_available

    print(f"Calculated run_date for GFS download: {run_date} (cycle hour: {cycle_hour:02d}Z)")

    return run_date


def _download_grib(url, dest_path, retries=3):
    """Download a GRIB2 file from *url* to *dest_path*. Returns True on success."""
    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(url, timeout=120, stream=True)
            resp.raise_for_status()
            content_type = resp.headers.get("Content-Type", "")
            if "html" in content_type.lower():
                print(f"  [!] Server returned HTML (likely bad request). Skipping.")
                return False
            with open(dest_path, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1 << 16):
                    fh.write(chunk)
            return True
        except requests.RequestException as exc:
            print(f"  [!] Attempt {attempt}/{retries} failed: {exc}")
            if attempt < retries:
                time.sleep(FETCH_DELAY)
    return False


def _grib_to_dataset(grib_path):
    """
    Open a GRIB2 file with cfgrib and return a merged xarray Dataset
    containing ugrd, vgrd, t, and gh on isobaric levels.
    cfgrib may split a GRIB2 into multiple datasets; we merge them.
    """
    if cfgrib is None:
        raise RuntimeError(
            "cfgrib is required to download GFS forecasts but is not installed. "
            "Install it via conda: `conda install -c conda-forge cfgrib eccodes`."
        )
    datasets = []
    try:
        datasets = cfgrib.open_datasets(
            str(grib_path),
            backend_kwargs={"indexpath": ""},   # don't write .idx files
        )
    except Exception as exc:
        print(f"  [!] cfgrib error opening {grib_path}: {exc}")
        return None

    if not datasets:
        return None

    isobaric_ds = [ds for ds in datasets if "isobaricInhPa" in ds.coords]
    if not isobaric_ds:
        print(f"  [!] No isobaric data found in {grib_path}")
        return None

    merged = xr.merge(isobaric_ds, compat="override")
    return merged

def _validate_requested_nc_start(nc_start, availability_lag_hours=4):
    """
    Validate the user-requested forecast cycle from config.

    Rules:
    - nc_start must be exactly on a GFS cycle hour: 00, 06, 12, or 18 UTC
    - it must not be older than NOAA's retention window
    - it must not be too recent / in the future relative to expected upload lag

    availability_lag_hours:
        Approximate delay after cycle time before files are likely available.
        3-4 hours is common; defaulting to 4 is safer.
    """
    run_date = nc_start.replace(minute=0, second=0, microsecond=0)

    if run_date.hour not in (0, 6, 12, 18):
        raise ValueError(
            f"netcdf_gfs['nc_start'] must be a GFS cycle hour "
            f"(00, 06, 12, 18 UTC), got {run_date}"
        )

    now_utc = datetime.datetime.utcnow().replace(tzinfo=None)

    # Too old
    GFS_RETENTION_DAYS = 9
    oldest_available = (now_utc - datetime.timedelta(days=GFS_RETENTION_DAYS)).replace(
        hour=0, minute=0, second=0, microsecond=0
    )
    if run_date < oldest_available:
        raise ValueError(
            f"netcdf_gfs['nc_start']={run_date} is older than NOAA retention window. "
            f"Oldest likely available is about {oldest_available}."
        )

    # Too new / likely not uploaded yet
    latest_likely_available = now_utc - datetime.timedelta(hours=availability_lag_hours)
    latest_cycle_hour = (latest_likely_available.hour // 6) * 6
    latest_cycle_dt = latest_likely_available.replace(
        hour=latest_cycle_hour, minute=0, second=0, microsecond=0
    )

    if run_date > latest_cycle_dt:
        raise ValueError(
            f"netcdf_gfs['nc_start']={run_date} is too recent to reliably exist yet. "
            f"With an assumed ~{availability_lag_hours} hour upload lag, the latest likely "
            f"available cycle is about {latest_cycle_dt}."
        )

    return run_date


# ─────────────────────────────────────────────────────────────────────────────
# Main download routine
# ─────────────────────────────────────────────────────────────────────────────

def download_gfs_grib_to_netcdf():
    """
    Download GFS GRIB2 subsets via the NOAA GRIB filter and save a single
    NetCDF file compatible with EarthSHAB's GFS.py reader.
    """
    start_dt = simulation["start_time"]
    lat_center = simulation["start_coord"]["lat"]
    lon_center = simulation["start_coord"]["lon"]
    lat_range = netcdf_gfs["lat_range"]
    lon_range = netcdf_gfs["lon_range"]
    download_days = netcdf_gfs["download_days"]
    step_hours = netcdf_gfs.get("step_hours", 3)
    out_filename = netcdf_gfs["nc_file"]
    forecast_start_time = netcdf_gfs["nc_start"]

    out_path = Path(out_filename)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Use nc_start as the requested forecast cycle
    run_date = _validate_requested_nc_start(forecast_start_time, availability_lag_hours=4)
    date_str = run_date.strftime("%Y%m%d")
    cycle_hour = run_date.hour

    total_hours = download_days * 24
    forecast_hours = list(range(0, total_hours + 1, step_hours))

    print(f"\nEarthSHAB GFS downloader (GRIB-filter edition)")
    print(f"  Model run   : {date_str} {cycle_hour:02d}Z")
    print(f"  Forecast hrs: {forecast_hours[0]}–{forecast_hours[-1]} (every {step_hours} h)")
    print(f"  Bounding box: lat {lat_center}±{lat_range/2}°, "
          f"lon {lon_center}±{lon_range/2}°")
    print(f"  Output      : {out_path}\n")
    #print(f"Start Time: {start_time} UTC")

    hourly_datasets = []

    with tempfile.TemporaryDirectory(prefix="earthshab_grib_") as tmpdir:
        for fhr in forecast_hours:
            url = _build_url(date_str, cycle_hour, fhr, lat_range, lon_range, lat_center, lon_center)
            grb_name = f"gfs_{date_str}_{cycle_hour:02d}z_f{fhr:03d}.grb2"
            grb_path = Path(tmpdir) / grb_name

            print(f"  Downloading f{fhr:03d} … ", end="", flush=True)
            ok = _download_grib(url, grb_path)
            if not ok:
                print("FAILED – skipping")
                time.sleep(FETCH_DELAY)
                continue
            print(f"OK ({grb_path.stat().st_size / 1024:.1f} kB)")

            ds = _grib_to_dataset(grb_path)
            if ds is None:
                time.sleep(FETCH_DELAY)
                continue

            # cfgrib creates datasets with 'time', 'step', and 'valid_time' coords.
            # Drop 'step' and keep 'valid_time' for the canonical format.
            if "step" in ds.coords:
                ds = ds.drop_vars("step")
            if "time" in ds.coords and "valid_time" in ds.coords:
                # Drop the forecast reference time, keep valid_time
                ds = ds.drop_vars("time")

            # Ensure valid_time is a dimension
            if "valid_time" in ds.coords and "valid_time" not in ds.dims:
                ds = ds.expand_dims("valid_time")

            hourly_datasets.append(ds)

            time.sleep(FETCH_DELAY)

    if not hourly_datasets:
        raise RuntimeError("No data was successfully downloaded. Check your config and network connectivity.")

    print("\nMerging all forecast hours …")
    combined = xr.concat(hourly_datasets, dim="valid_time")

    # Convert to canonical v2 format
    # Rename variables: gh → z, u/v/t stay the same
    rename_map = {
        "gh": "z",
        "isobaricInhPa": "pressure_level",
        "latitude": "latitude",  # keep name
        "longitude": "longitude",  # keep name
    }
    combined = combined.rename({k: v for k, v in rename_map.items() if k in combined})

    # Convert geopotential height (m) to geopotential (m²/s²) by multiplying by g
    g = 9.80665
    combined["z"] = combined["z"] * g
    combined["z"].attrs["units"] = "m**2 s**-2"
    combined["z"].attrs["long_name"] = "Geopotential"
    combined["z"].attrs["standard_name"] = "geopotential"

    # Convert longitude from [0, 360) to [-180, 180)
    lon_vals = combined["longitude"].values
    lon_vals = np.where(lon_vals >= 180, lon_vals - 360, lon_vals)
    combined = combined.assign_coords(longitude=lon_vals)

    # Sort by longitude to maintain ascending order after conversion
    combined = combined.sortby("longitude")

    # Flip latitude to descending order (north to south)
    combined = combined.sortby("latitude", ascending=False)

    # Set coordinate attributes
    combined["pressure_level"].attrs["units"] = "hPa"
    combined["pressure_level"].attrs["long_name"] = "pressure"
    combined["pressure_level"].attrs["standard_name"] = "air_pressure"
    combined["pressure_level"].attrs["positive"] = "down"
    combined["pressure_level"].attrs["stored_direction"] = "decreasing"

    combined["latitude"].attrs["units"] = "degrees_north"
    combined["latitude"].attrs["long_name"] = "latitude"
    combined["latitude"].attrs["standard_name"] = "latitude"
    combined["latitude"].attrs["stored_direction"] = "decreasing"

    combined["longitude"].attrs["units"] = "degrees_east"
    combined["longitude"].attrs["long_name"] = "longitude"
    combined["longitude"].attrs["standard_name"] = "longitude"

    # Set u/v/t attributes
    combined["u"].attrs["units"] = "m s**-1"
    combined["u"].attrs["long_name"] = "U component of wind"
    combined["u"].attrs["standard_name"] = "eastward_wind"

    combined["v"].attrs["units"] = "m s**-1"
    combined["v"].attrs["long_name"] = "V component of wind"
    combined["v"].attrs["standard_name"] = "northward_wind"

    combined["t"].attrs["units"] = "K"
    combined["t"].attrs["long_name"] = "Temperature"
    combined["t"].attrs["standard_name"] = "air_temperature"

    # Convert time to seconds since 1970-01-01 (Unix epoch)
    time_values = combined["valid_time"].values
    datetime_objects = [pd.Timestamp(t).to_pydatetime() for t in time_values]

    print(f"time values (datetime): {datetime_objects}")

    epoch_seconds = nc.date2num(
        datetime_objects,
        units="seconds since 1970-01-01",
        calendar="proleptic_gregorian",
    )

    combined = combined.assign_coords(valid_time=epoch_seconds.astype(np.int64))
    combined["valid_time"].attrs = {
        "units": "seconds since 1970-01-01",
        "calendar": "proleptic_gregorian",
        "standard_name": "time",
        "long_name": "time",
    }

    # Add global CF-1.7 convention and provenance
    combined.attrs["Conventions"] = "CF-1.7"
    combined.attrs["institution"] = "NOAA/NCEP (GFS)"
    combined.attrs["history"] = f"{pd.Timestamp.now().isoformat()} - Downloaded and converted by EarthSHAB saveNETCDF.py"

    # Reorder dimensions to canonical order: valid_time, pressure_level, latitude, longitude
    combined = combined.transpose("valid_time", "pressure_level", "latitude", "longitude")

    # Order variables: u, v, z, t
    variable_order = ["u", "v", "z", "t"]
    combined = xr.Dataset(
        {var: combined[var] for var in variable_order if var in combined},
        coords=combined.coords
    )

    # Remove unnecessary coordinates (already handled above but double-check)
    for coord in ["step", "time"]:
        if coord in combined.coords and coord not in combined.dims:
            combined = combined.drop_vars(coord)


    print(f"Writing {out_path} …")
    encoding = {var: {"zlib": True, "complevel": 4} for var in combined.data_vars}
    encoding["valid_time"] = {"_FillValue": None}  # Remove _FillValue from the time variable
    combined.to_netcdf(str(out_path), encoding=encoding)
    print(f"\nDone!  Saved to: {out_path}")
    print(f"  Dimensions : {dict(combined.sizes)}")
    print(f"  Variables  : {list(combined.data_vars)}")
    print(f"  Time values (epoch seconds): {epoch_seconds}")
    return str(out_path)


# ─────────────────────────────────────────────────────────────────────────────
# Live-vs-archive dispatch
# ─────────────────────────────────────────────────────────────────────────────

def _is_past_retention(nc_start, retention_days=9):
    """True if *nc_start* is older than NOAA's live (NOMADS) retention window."""
    now_utc = datetime.datetime.utcnow()
    oldest_available = (now_utc - datetime.timedelta(days=retention_days)).replace(
        hour=0, minute=0, second=0, microsecond=0
    )
    return nc_start.replace(minute=0, second=0, microsecond=0) < oldest_available


def _prompt_proceed_archive(nc_start):
    """Warn (yellow) that the cycle is past the live window and ask to use the archive.

    Returns True to proceed with the AWS archive download. In a non-interactive
    session (no TTY) it defaults to proceeding, since the archive is the only way
    to obtain a past cycle — it prints a notice so the choice is visible in logs.
    """
    print(colored(
        f"\n[WARNING] Requested cycle {nc_start} is older than NOAA's ~9-day live "
        f"retention window, so it is not available from NOMADS.", "yellow"))
    print(colored(
        "          A historical copy is available from the AWS GFS archive "
        "(noaa-gfs-bdp-pds).", "yellow"))
    if not sys.stdin.isatty():
        print(colored("          Non-interactive session — proceeding with the "
                      "archived forecast.", "yellow"))
        return True
    answer = input(colored(
        "          Download the archived forecast instead? [Y/n]: ", "yellow")).strip().lower()
    return answer in ("", "y", "yes")


def main():
    """Download the configured GFS cycle, auto-switching to the AWS archive
    (saveNETCDF_archive.py) when the cycle predates NOAA's live retention window."""
    nc_start = netcdf_gfs["nc_start"]
    if _is_past_retention(nc_start):
        if not _prompt_proceed_archive(nc_start):
            print(colored("Aborted. Choose a cycle within the last ~9 days for a "
                          "live download, or confirm the archive next time.", "red"))
            sys.exit(1)
        # Auto-switch to the archive downloader.
        from EarthSHAB.saveNETCDF_archive import download_gfs_archive_to_netcdf
        return download_gfs_archive_to_netcdf()
    return download_gfs_grib_to_netcdf()


if __name__ == "__main__":
    main()