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

import os
import time
import datetime
import tempfile
import requests
import numpy as np
import xarray as xr
import cfgrib
import netCDF4 as nc
from pathlib import Path
from termcolor import colored

# ── EarthSHAB config ──────────────────────────────────────────────────────────
from config_earth import netcdf_gfs, simulation

# ── Pressure levels available in GFS 0.25° pgrb2 isobaric data ───────────────
# This is the exact set accepted by the GRIB filter. Levels not in this list
# (e.g. 875, 825, 775 mb) will cause an "invalid parameter" error.
# Source: https://nomads.ncep.noaa.gov/gribfilter.php?ds=gfs_0p25
# EarthSHAB SHAB flights typically reach 15,000–25,000 m (~55–25 mb),
# so all levels from 1000 mb down to 10 mb are included here.
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
    # to prevent floating point garbage appearing in the URL.
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


def _download_grib(url, dest_path, retries=3):
    """Download a GRIB2 file from *url* to *dest_path*. Returns True on success."""
    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(url, timeout=120, stream=True)
            resp.raise_for_status()
            # NOAA returns HTML error pages (not HTTP error codes) on bad requests
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

    # Keep only datasets that have an isobaricInhPa vertical coordinate
    isobaric_ds = [ds for ds in datasets if "isobaricInhPa" in ds.coords]
    if not isobaric_ds:
        print(f"  [!] No isobaric data found in {grib_path}")
        return None

    merged = xr.merge(isobaric_ds, compat="override")
    return merged


# ─────────────────────────────────────────────────────────────────────────────
# Main download routine
# ─────────────────────────────────────────────────────────────────────────────

def download_gfs_grib_to_netcdf():
    """
    Download GFS GRIB2 subsets via the NOAA GRIB filter and save a single
    NetCDF file compatible with EarthSHAB's GFS.py reader.
    """
    # ── Pull settings from config_earth.py ────────────────────────────────────
    start_dt      = simulation["start_time"]          # datetime object
    lat_center    = simulation["start_coord"]["lat"]  # centre latitude
    lon_center    = simulation["start_coord"]["lon"]  # centre longitude, [-180,180]
    lat_range     = netcdf_gfs["lat_range"]           # index counts (× res = degrees)
    lon_range     = netcdf_gfs["lon_range"]           # index counts (× res = degrees)
    download_days = netcdf_gfs["download_days"]       # number of days
    out_filename  = netcdf_gfs["nc_file"]             # full path e.g. "forecasts/gfs_....nc"

    out_path = Path(out_filename)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # GFS runs at 00, 06, 12, 18 UTC – pick the most recent one before start
    cycle_hour = (start_dt.hour // 6) * 6
    run_date   = start_dt.replace(hour=cycle_hour, minute=0, second=0,
                                  microsecond=0)
    date_str   = run_date.strftime("%Y%m%d")

    # Forecast hours to cover: every 3 h for download_days days
    total_hours    = download_days * 24
    forecast_hours = list(range(0, total_hours + 1, 3))

    print(f"\nEarthSHAB GFS downloader (GRIB-filter edition)")
    print(f"  Model run   : {date_str} {cycle_hour:02d}Z")
    print(f"  Forecast hrs: {forecast_hours[0]}–{forecast_hours[-1]} (every 3 h)")
    print(f"  Bounding box: lat {lat_center}±{lat_range/2}°, "
          f"lon {lon_center}±{lon_range/2}°")
    print(f"  Output      : {out_path}\n")

    hourly_datasets = []

    # Safety check: ensure requested GFS forecast still exists on NOAA's servers
    now_utc = datetime.datetime.utcnow()

    # NOAA GFS retention window (conservative)
    GFS_RETENTION_DAYS = 9
    oldest_available = (now_utc - datetime.timedelta(days=GFS_RETENTION_DAYS)).replace(
    hour=0, minute=0, second=0, microsecond=0
)

    if start_dt < oldest_available:
        print(colored("\n[ERROR] Requested simulation start time is too far in the past for downloadable forecast.", "red"))
        print(colored(f"  Simulation start : {start_dt}", "yellow"))
        print(colored(f"  Oldest available : {oldest_available}", "yellow"))
        raise ValueError()

    with tempfile.TemporaryDirectory(prefix="earthshab_grib_") as tmpdir:
        for fhr in forecast_hours:
            url       = _build_url(date_str, cycle_hour, fhr,
                                   lat_range, lon_range,
                                   lat_center, lon_center)
            grb_name  = f"gfs_{date_str}_{cycle_hour:02d}z_f{fhr:03d}.grb2"
            grb_path  = Path(tmpdir) / grb_name

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

            # Stamp with valid time so we can concatenate along time later
            valid_time = run_date + datetime.timedelta(hours=fhr)
            ds = ds.expand_dims("time").assign_coords(
                time=[np.datetime64(valid_time, "ns")]
            )
            hourly_datasets.append(ds)

            time.sleep(FETCH_DELAY)   # be polite to NOAA servers

    if not hourly_datasets:
        raise RuntimeError("No data was successfully downloaded. "
                           "Check your config and network connectivity.")

    print("\nMerging all forecast hours …")
    combined = xr.concat(hourly_datasets, dim="time")

    # ── Rename variables to EarthSHAB conventions ─────────────────────────────
    # cfgrib uses ECMWF short names; map to the names GFS.py looks for.
    rename_map = {}
    name_mapping = {
        "u":              "ugrdprs",
        "v":              "vgrdprs",
        "t":              "tmpprs",
        "gh":             "hgtprs",
        "isobaricInhPa":  "lev",
        "latitude":       "lat",
        "longitude":      "lon",
    }
    for cfgrib_name, earthshab_name in name_mapping.items():
        if cfgrib_name in combined:
            rename_map[cfgrib_name] = earthshab_name
        if cfgrib_name in combined.coords or cfgrib_name in combined.dims:
            rename_map[cfgrib_name] = earthshab_name

    combined = combined.rename(rename_map)

    # Ensure pressure levels are in hPa (cfgrib already stores them that way)
    if "lev" in combined.coords:
        combined["lev"].attrs["units"] = "hPa"
        combined["lev"].attrs["long_name"] = "pressure level"

    # latitude / longitude dimension names
    for old, new in [("latitude", "lat"), ("longitude", "lon")]:
        if old in combined.dims:
            combined = combined.rename({old: new})

    # ── Write NetCDF ──────────────────────────────────────────────────────────
    print(f"Writing {out_path} …")
    encoding = {
        var: {"zlib": True, "complevel": 4}
        for var in combined.data_vars
    }
    combined.to_netcdf(str(out_path), encoding=encoding)
    print(f"\nDone!  Saved to: {out_path}")
    print(f"  Dimensions : {dict(combined.dims)}")
    print(f"  Variables  : {list(combined.data_vars)}")
    return str(out_path)


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    download_gfs_grib_to_netcdf()