"""
saveNETCDF_archive.py  –  EarthSHAB GFS forecast downloader (AWS archive edition)
=================================================================================
NOAA's NOMADS GRIB-filter service (used by ``saveNETCDF.py``) only retains the
last ~9 days of GFS cycles, so historical forecasts can't be fetched live. This
module is the historical counterpart: it pulls the *identical* ``pgrb2.0p25``
product from the AWS Open Data mirror (``noaa-gfs-bdp-pds``), which goes back to
~early 2021, using the per-file ``.idx`` sidecar to grab only the GRIB2 messages
we need via HTTP byte-range requests.

It is config-compatible with ``saveNETCDF.py`` — it reads the same
``config_earth.netcdf_gfs`` / ``simulation`` settings, downloads the same
variables (UGRD/VGRD/TMP/HGT) on the same 28 isobaric levels every 3 hours over
``download_days``, and writes a NetCDF-4 file in the same canonical v2 schema —
so a file produced here is a drop-in replacement for a live download.

``saveNETCDF.py`` auto-switches to this module when the requested cycle predates
NOAA's live retention window (see ``saveNETCDF.main``).

Variables downloaded (all on isobaric pressure levels):
  UGRD – U-component of wind  (m/s)    VGRD – V-component of wind (m/s)
  TMP  – Temperature          (K)      HGT  – Geopotential height (gpm)

Requirements: cfgrib, eccodes, xarray, netCDF4, requests (same as saveNETCDF.py).

Usage:
  python -m EarthSHAB.saveNETCDF_archive
  (uses config_earth.py for lat/lon centre, range, date, and download_days)
"""

import datetime
import os
import sys
import tempfile
import time

import netCDF4 as nc
import numpy as np
import pandas as pd
import requests
import xarray as xr
from termcolor import colored

try:
    import cfgrib
    # cfgrib emits FutureWarnings about combine_attrs during merge; silence them.
    xr.set_options(use_new_combine_kwarg_defaults=True)
except ImportError:  # pragma: no cover - cfgrib required at runtime
    cfgrib = None

from EarthSHAB.config_earth import netcdf_gfs, simulation

# ── Constants ──────────────────────────────────────────────────────────────────

AWS_BASE = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"

WANT_VARS = ("UGRD", "VGRD", "HGT", "TMP")

G = 9.80665  # geopotential height (gpm) -> geopotential (m^2/s^2)

# Earliest cycle for which pgrb2.0p25 exists in the AWS bucket (approx).
ARCHIVE_START = datetime.datetime(2021, 3, 1)

# Match saveNETCDF.py exactly so an archive file is a drop-in for a live one.
PRESSURE_LEVELS_MB = [
    1000, 975, 950, 925, 900, 850, 800, 750, 700, 650,
     600, 550, 500, 450, 400, 350, 300, 250, 200, 150,
     100,  70,  50,  40,  30,  20,  15,  10,
]

_session = requests.Session()
_session.headers.update({"User-Agent": "EarthSHAB-saveNETCDF-archive/1.0"})

# Seconds between fetches (be polite; AWS has no published rate limit).
FETCH_DELAY = 0.0


# ── URL construction ────────────────────────────────────────────────────────────

def cycle_file_url(date_str, cycle_hour, forecast_hour, base=AWS_BASE):
    """S3 object URL for one GFS forecast step (date_str='YYYYMMDD')."""
    cyc = f"{int(cycle_hour):02d}"
    fhh = f"{int(forecast_hour):03d}"
    return f"{base}/gfs.{date_str}/{cyc}/atmos/gfs.t{cyc}z.pgrb2.0p25.f{fhh}"


# ── .idx parsing and byte-range selection ───────────────────────────────────────

def parse_idx(text):
    """Parse a GFS ``.idx`` body into ``(msgnum, offset, var, level)`` tuples."""
    rows = []
    for line in text.splitlines():
        parts = line.split(":")
        if len(parts) < 5:
            continue
        try:
            msgnum = int(parts[0])
            offset = int(parts[1])
        except ValueError:
            continue
        rows.append((msgnum, offset, parts[3], parts[4]))
    return rows


def _level_mb(level_str):
    """Float pressure (mb) for an isobaric level string, else None."""
    if not level_str.endswith(" mb"):
        return None
    try:
        return float(level_str[:-3])
    except ValueError:
        return None


def select_ranges(rows, want_vars=WANT_VARS, want_levels_mb=None):
    """Select ``(start, end, var, level_mb)`` byte ranges for the wanted fields.

    ``end`` is the inclusive last byte (next message offset - 1), or None for the
    final message in the file (open-ended range).
    """
    want_vars = set(want_vars)
    want_set = [float(x) for x in want_levels_mb] if want_levels_mb is not None else None
    offsets = [r[1] for r in rows]
    out = []
    for i, (_, offset, var, level) in enumerate(rows):
        if var not in want_vars:
            continue
        mb = _level_mb(level)
        if mb is None:
            continue
        if want_set is not None and not any(abs(mb - w) < 1e-6 for w in want_set):
            continue
        end = offsets[i + 1] - 1 if i + 1 < len(rows) else None
        out.append((offset, end, var, mb))
    return out


# ── Download ────────────────────────────────────────────────────────────────────

def _get(url, headers=None, timeout=300, retries=3, delay=2.0):
    last = None
    for attempt in range(1, retries + 1):
        try:
            resp = _session.get(url, headers=headers, timeout=timeout)
            resp.raise_for_status()
            return resp
        except requests.RequestException as exc:
            last = exc
            if attempt < retries:
                time.sleep(delay)
    raise RuntimeError(f"GET failed after {retries} attempts: {url} ({last})")


def fetch_grib_bytes(base_url, want_levels_mb=None):
    """Fetch the selected GRIB2 messages for one forecast step; return (bytes, n)."""
    idx_resp = _get(base_url + ".idx", timeout=120)
    ranges = select_ranges(parse_idx(idx_resp.text), WANT_VARS, want_levels_mb)
    if not ranges:
        raise RuntimeError(f"No matching GRIB messages in {base_url}.idx")
    buf = bytearray()
    for start, end, _var, _mb in ranges:
        rng = f"bytes={start}-{end}" if end is not None else f"bytes={start}-"
        buf += _get(base_url, headers={"Range": rng}).content
    return bytes(buf), len(ranges)


# ── Canonical v2 conversion (mirrors saveNETCDF.py / migrate_v1.py) ──────────────

def _grib_to_isobaric(grib_path):
    """Open a GRIB2 file with cfgrib and merge isobaric datasets (u/v/t/gh).

    GFS stores levels >= 1 hPa under ``isobaricInhPa`` (hPa) and levels < 1 hPa
    under ``isobaricInPa`` (Pa), so cfgrib returns two datasets. Normalise both
    to ``isobaricInhPa`` (hPa) and concatenate along the level dimension.
    """
    datasets = cfgrib.open_datasets(str(grib_path), backend_kwargs={"indexpath": ""})
    parts = []
    for d in datasets:
        if "isobaricInhPa" in d.coords:
            parts.append(d)
        elif "isobaricInPa" in d.coords:
            d = d.rename({"isobaricInPa": "isobaricInhPa"})
            d = d.assign_coords(isobaricInhPa=d["isobaricInhPa"].values / 100.0)
            parts.append(d)
    if not parts:
        raise RuntimeError(f"No isobaric data in {grib_path}")
    merged = parts[0] if len(parts) == 1 else xr.concat(parts, dim="isobaricInhPa")
    merged = merged.load()
    if "step" in merged.coords:
        merged = merged.drop_vars("step")
    if "time" in merged.coords and "valid_time" in merged.coords:
        merged = merged.drop_vars("time")
    if "valid_time" in merged.coords and "valid_time" not in merged.dims:
        merged = merged.expand_dims("valid_time")
    return merged


def to_canonical_v2(combined):
    """Convert a merged cfgrib Dataset to EarthSHAB's canonical v2 schema."""
    rename_map = {"gh": "z", "isobaricInhPa": "pressure_level"}
    combined = combined.rename({k: v for k, v in rename_map.items() if k in combined})

    combined["z"] = combined["z"] * G
    combined["z"].attrs.update({
        "units": "m**2 s**-2", "long_name": "Geopotential", "standard_name": "geopotential",
    })

    lon_vals = combined["longitude"].values
    lon_vals = np.where(lon_vals >= 180, lon_vals - 360, lon_vals)
    combined = combined.assign_coords(longitude=lon_vals).sortby("longitude")
    combined = combined.sortby("latitude", ascending=False)
    combined = combined.sortby("pressure_level", ascending=False)

    combined["pressure_level"].attrs.update({
        "units": "hPa", "long_name": "pressure", "standard_name": "air_pressure",
        "positive": "down", "stored_direction": "decreasing",
    })
    combined["latitude"].attrs.update({
        "units": "degrees_north", "long_name": "latitude", "standard_name": "latitude",
        "stored_direction": "decreasing",
    })
    combined["longitude"].attrs.update({
        "units": "degrees_east", "long_name": "longitude", "standard_name": "longitude",
    })
    combined["u"].attrs.update({"units": "m s**-1", "long_name": "U component of wind", "standard_name": "eastward_wind"})
    combined["v"].attrs.update({"units": "m s**-1", "long_name": "V component of wind", "standard_name": "northward_wind"})
    combined["t"].attrs.update({"units": "K", "long_name": "Temperature", "standard_name": "air_temperature"})

    time_values = combined["valid_time"].values
    datetime_objects = [pd.Timestamp(t).to_pydatetime() for t in time_values]
    epoch_seconds = nc.date2num(
        datetime_objects, units="seconds since 1970-01-01", calendar="proleptic_gregorian",
    )
    combined = combined.assign_coords(valid_time=np.atleast_1d(epoch_seconds).astype(np.int64))
    combined["valid_time"].attrs = {
        "units": "seconds since 1970-01-01", "calendar": "proleptic_gregorian",
        "standard_name": "time", "long_name": "time",
    }

    combined.attrs["Conventions"] = "CF-1.7"
    combined.attrs["institution"] = "NOAA/NCEP (GFS)"
    combined.attrs["history"] = (
        f"{pd.Timestamp.now().isoformat()} - Downloaded from AWS noaa-gfs-bdp-pds "
        f"and converted by EarthSHAB saveNETCDF_archive.py"
    )

    combined = combined.transpose("valid_time", "pressure_level", "latitude", "longitude")
    variable_order = [v for v in ("u", "v", "z", "t") if v in combined]
    combined = xr.Dataset({v: combined[v] for v in variable_order}, coords=combined.coords)
    return combined


def _crop(ds, bbox):
    eps = 1e-3
    return ds.sel(
        latitude=slice(bbox["lat_max"] + eps, bbox["lat_min"] - eps),
        longitude=slice(bbox["lon_min"] - eps, bbox["lon_max"] + eps),
    )


def fetch_cycle(date_str, cycle_hour, forecast_hours, bbox=None,
                levels_mb=None, base=AWS_BASE, delay=FETCH_DELAY, verbose=True):
    """Fetch a forecast cycle from the archive; return a canonical-v2 Dataset.

    Each step is converted and cropped immediately so memory stays bounded
    regardless of how many forecast hours are requested.
    """
    if cfgrib is None:
        raise RuntimeError("cfgrib is required for fetch_cycle but is not installed")

    steps = []
    total_bytes = 0
    for fhr in forecast_hours:
        url = cycle_file_url(date_str, cycle_hour, fhr, base)
        if verbose:
            print(f"  Downloading f{int(fhr):03d} … ", end="", flush=True)
        try:
            data, n_msgs = fetch_grib_bytes(url, levels_mb)
        except RuntimeError as exc:
            print(colored(f"FAILED – {exc}", "red"))
            time.sleep(delay)
            continue
        total_bytes += len(data)
        if verbose:
            print(f"OK ({n_msgs} msgs, {len(data)/1024:.0f} kB)")
        with tempfile.NamedTemporaryFile(suffix=".grb2", delete=False) as tmp:
            tmp.write(data)
            tmp_path = tmp.name
        try:
            step = to_canonical_v2(_grib_to_isobaric(tmp_path))
            if bbox is not None:
                step = _crop(step, bbox)
            steps.append(step.load())
        finally:
            os.unlink(tmp_path)
        if delay:
            time.sleep(delay)

    if not steps:
        raise RuntimeError(
            "No forecast hours downloaded. Check connectivity and that the cycle "
            "exists in the AWS archive."
        )

    ds = xr.concat(steps, dim="valid_time") if len(steps) > 1 else steps[0]
    if verbose:
        print(f"  Total downloaded: {total_bytes/1e6:.1f} MB")
    return ds


# ── Validation ───────────────────────────────────────────────────────────────────

def _validate_archive_nc_start(nc_start):
    """Validate the requested cycle for an archive download.

    Must be exactly on a GFS cycle hour (00/06/12/18 UTC), not in the future,
    and not before the AWS pgrb2.0p25 archive start.
    """
    run_date = nc_start.replace(minute=0, second=0, microsecond=0)
    if run_date.hour not in (0, 6, 12, 18):
        raise ValueError(
            f"netcdf_gfs['nc_start'] must be a GFS cycle hour (00, 06, 12, 18 UTC), "
            f"got {run_date}"
        )
    if run_date < ARCHIVE_START:
        raise ValueError(
            f"netcdf_gfs['nc_start']={run_date} predates the AWS GFS archive "
            f"(~{ARCHIVE_START.date()}). Pre-2021 cycles need NCAR RDA ds084.1."
        )
    now_utc = datetime.datetime.utcnow().replace(tzinfo=None)
    if run_date > now_utc:
        raise ValueError(f"netcdf_gfs['nc_start']={run_date} is in the future.")
    return run_date


# ── Main download routine ─────────────────────────────────────────────────────────

def download_gfs_archive_to_netcdf():
    """Download a historical GFS cycle from the AWS archive into a v2 NetCDF.

    Reads the same config as saveNETCDF.py and produces a structurally identical
    file (same 28 levels, every 3 h over download_days, same bounding box).
    """
    lat_center = simulation["start_coord"]["lat"]
    lon_center = simulation["start_coord"]["lon"]
    lat_range = netcdf_gfs["lat_range"]
    lon_range = netcdf_gfs["lon_range"]
    download_days = netcdf_gfs["download_days"]
    step_hours = netcdf_gfs.get("step_hours", 3)
    out_filename = netcdf_gfs["nc_file"]
    forecast_start_time = netcdf_gfs["nc_start"]

    run_date = _validate_archive_nc_start(forecast_start_time)
    date_str = run_date.strftime("%Y%m%d")
    cycle_hour = run_date.hour

    total_hours = download_days * 24
    forecast_hours = list(range(0, total_hours + 1, step_hours))

    bbox = {
        "lat_max": round(min(90.0, lat_center + lat_range / 2.0), 4),
        "lat_min": round(max(-90.0, lat_center - lat_range / 2.0), 4),
        "lon_min": round(max(-180.0, lon_center - lon_range / 2.0), 4),
        "lon_max": round(min(180.0, lon_center + lon_range / 2.0), 4),
    }

    out_dir = os.path.dirname(out_filename)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(colored("\nEarthSHAB GFS downloader (AWS archive edition)", "blue", attrs=["bold"]))
    print(f"  Source      : {AWS_BASE} (noaa-gfs-bdp-pds, pgrb2.0p25)")
    print(f"  Model run   : {date_str} {cycle_hour:02d}Z")
    print(f"  Forecast hrs: {forecast_hours[0]}–{forecast_hours[-1]} (every {step_hours} h)")
    print(f"  Bounding box: lat[{bbox['lat_min']},{bbox['lat_max']}] "
          f"lon[{bbox['lon_min']},{bbox['lon_max']}]")
    print(f"  Output      : {out_filename}\n")

    ds = fetch_cycle(date_str, cycle_hour, forecast_hours, bbox=bbox,
                     levels_mb=PRESSURE_LEVELS_MB)

    encoding = {var: {"zlib": True, "complevel": 4} for var in ds.data_vars}
    encoding["valid_time"] = {"_FillValue": None}
    ds.to_netcdf(str(out_filename), encoding=encoding)

    print(colored(f"\nDone!  Saved to: {out_filename}", "green"))
    print(f"  Dimensions : {dict(ds.sizes)}")
    print(f"  Variables  : {list(ds.data_vars)}")
    return str(out_filename)


if __name__ == "__main__":
    download_gfs_archive_to_netcdf()
