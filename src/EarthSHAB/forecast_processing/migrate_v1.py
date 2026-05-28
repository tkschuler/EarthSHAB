"""migrate_v1: convert pre-v2 EarthSHAB forecast files to the canonical v2 schema.

Phase 4 of the v2.0 strangler-fig refactor. The v2 reader rejects every legacy
layout it encounters (see ForecastFormatError in ERA5.py); this script is the
single supported path to upgrade an archived file. There is intentionally no
silent auto-conversion inside the reader — see docs/migration-v2.md.

Two legacy formats are detected and converted:

  1. v1 GFS
     - signature: data variables include ugrdprs / vgrdprs / hgtprs / tmpprs
     - dims  : time, lev, lat, lon
     - time  : 'days since 0001-01-01' (sometimes implicit; GFS.py hardcoded it)
     - lat   : ASCENDING
     - lon   : [0, 360)
     - hgtprs: geopotential height in METERS (gpm)

  2. v1 ERA5
     - signature: vars u/v/z/t with dims time / level / latitude / longitude
     - level : ASCENDING hPa
     - z     : geopotential (m**2 s**-2)
     - Covers both the convert_era52gfs.py 'processed' output (CF-1.7,
       'seconds since 1970-01-01') and pre-Sep-2024 raw CDS downloads
       (CF-1.6, 'hours since 1900-01-01'). xarray decodes either form on
       open; the converter handles them identically.

For each input <name>.nc:
  - the original is renamed to <name>.v1.nc as a safety backup
  - the canonical file is written in its place at <name>.nc

The operation is idempotent. If <name>.v1.nc already exists, or the file is
already in v2 canonical layout, it is skipped with a clear message.

Other layouts are not recognized and are skipped with ERROR.

CLI:
    python -m EarthSHAB.forecast_processing.migrate_v1 <file_or_dir> [...]
    python -m EarthSHAB.forecast_processing.migrate_v1 --dry-run <target>

When a directory is given, every *.nc file at the top level is processed
(non-recursive). Existing *.v1.nc backups are skipped automatically.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm


G = 9.80665  # standard gravity, used to convert hgtprs (m) -> z (m^2 s^-2)


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

FORMAT_V2 = "v2"
FORMAT_V1_GFS = "v1_gfs"
FORMAT_V1_ERA5 = "v1_era5"
FORMAT_UNKNOWN = "unknown"

V2_REQUIRED_DIMS = {"valid_time", "pressure_level", "latitude", "longitude"}
# `t` is intentionally optional: _finalize() drops it for archived GFS files
# whose source layout had no `tmpprs` variable.
V2_REQUIRED_VARS = {"u", "v", "z"}
V1_GFS_VARS = {"ugrdprs", "vgrdprs", "hgtprs"}
V1_ERA5_DIMS = {"time", "level", "latitude", "longitude"}


def detect_format(path: Path) -> str:
    """Return one of the FORMAT_* tags by sniffing variables and dimensions."""
    with netCDF4.Dataset(str(path)) as ds:
        vars_ = set(ds.variables.keys())
        dims_ = set(ds.dimensions.keys())

    if V2_REQUIRED_DIMS.issubset(dims_) and V2_REQUIRED_VARS.issubset(vars_):
        return FORMAT_V2
    if V1_GFS_VARS.issubset(vars_):
        return FORMAT_V1_GFS
    if V1_ERA5_DIMS.issubset(dims_) and V2_REQUIRED_VARS.issubset(vars_):
        return FORMAT_V1_ERA5
    return FORMAT_UNKNOWN


# ---------------------------------------------------------------------------
# Canonical attribute / encoding helpers
# ---------------------------------------------------------------------------

def _encode_valid_time(times) -> tuple[np.ndarray, dict]:
    """Encode a sequence of datetime-like values as 'seconds since 1970-01-01'."""
    dts = [pd.Timestamp(t).to_pydatetime() for t in times]
    epoch_seconds = netCDF4.date2num(
        dts,
        units="seconds since 1970-01-01",
        calendar="proleptic_gregorian",
    )
    attrs = {
        "units": "seconds since 1970-01-01",
        "calendar": "proleptic_gregorian",
        "standard_name": "time",
        "long_name": "time",
    }
    return np.asarray(epoch_seconds, dtype=np.int64), attrs


def _apply_canonical_attrs(ds: xr.Dataset) -> xr.Dataset:
    """Set canonical CF attrs on whichever data vars are present in `ds`."""
    ds["pressure_level"].attrs.update({
        "units": "hPa",
        "long_name": "pressure",
        "standard_name": "air_pressure",
        "positive": "down",
        "stored_direction": "decreasing",
    })
    ds["latitude"].attrs.update({
        "units": "degrees_north",
        "long_name": "latitude",
        "standard_name": "latitude",
        "stored_direction": "decreasing",
    })
    ds["longitude"].attrs.update({
        "units": "degrees_east",
        "long_name": "longitude",
        "standard_name": "longitude",
    })
    data_attrs = {
        "u": {"units": "m s**-1", "long_name": "U component of wind", "standard_name": "eastward_wind"},
        "v": {"units": "m s**-1", "long_name": "V component of wind", "standard_name": "northward_wind"},
        "z": {"units": "m**2 s**-2", "long_name": "Geopotential", "standard_name": "geopotential"},
        "t": {"units": "K", "long_name": "Temperature", "standard_name": "air_temperature"},
    }
    for name, attrs in data_attrs.items():
        if name in ds.data_vars:
            ds[name].attrs.update(attrs)
    return ds


def _finalize(ds: xr.Dataset, source_format: str) -> xr.Dataset:
    """Common end-of-pipeline: order dims, cast to float32, set attrs."""
    ds = ds.transpose("valid_time", "pressure_level", "latitude", "longitude")
    # `t` is required by the canonical schema but neither reader actually
    # uses temperature today; older archived GFS files (pre-tmpprs) round-
    # trip cleanly if we drop it from the canonical write.
    present_data_vars = [v for v in ("u", "v", "z", "t") if v in ds.data_vars]
    canon = xr.Dataset(
        {v: ds[v].astype(np.float32) for v in present_data_vars},
        coords=ds.coords,
    )
    canon = _apply_canonical_attrs(canon)
    canon.attrs["Conventions"] = "CF-1.7"
    canon.attrs["institution"] = ds.attrs.get("institution", "")
    history_prev = ds.attrs.get("history", "").strip()
    stamp = datetime.now(timezone.utc).isoformat()
    new_line = (
        f"{stamp} - Migrated from {source_format} to v2 canonical by "
        "EarthSHAB.forecast_processing.migrate_v1"
    )
    canon.attrs["history"] = f"{history_prev}\n{new_line}".strip()
    return canon


# ---------------------------------------------------------------------------
# Per-format converters
# ---------------------------------------------------------------------------

def _decode_gfs_time(t_var: xr.DataArray) -> np.ndarray:
    """Reproduce GFS.py's hardcoded decode so timestamps match what the v1
    reader would have surfaced — even when the file omits the time attrs.

    netCDF4.num2date with has_year_zero=True returns cftime objects regardless
    of only_use_cftime_datetimes; convert each to a stdlib datetime by reading
    its calendar fields directly so downstream encoding can treat them
    uniformly as proleptic_gregorian.
    """
    units = t_var.attrs.get("units") or "days since 0001-01-01"
    calendar = t_var.attrs.get("calendar") or "standard"
    decoded = netCDF4.num2date(
        np.asarray(t_var.values),
        units=units,
        calendar=calendar,
        has_year_zero=True,
    )
    return np.asarray([
        datetime(d.year, d.month, d.day, d.hour, d.minute, d.second, d.microsecond)
        for d in decoded
    ])


def _valid_subset_indices(src: Path) -> dict:
    """Find the unmasked rectangular bounding box (time, lat, lon) in a v1 GFS file.

    Archived global GFS downloads typically carry data for only a small
    region and a short time window — the rest of the grid is masked. The
    readers do exactly this at load time (GFS.determineRanges); doing it
    here before xarray sees the file lets us subset away >99% of empty
    cells before any rearrange or rewrite work.
    """
    with netCDF4.Dataset(str(src)) as ds:
        # u[:, 0, :, :] mirrors the slice GFS.determineRanges samples.
        sample = ds.variables["ugrdprs"][:, 0, :, :]
        if not np.ma.is_masked(sample):
            n_time = ds.dimensions["time"].size
            n_lat = ds.dimensions["lat"].size
            n_lon = ds.dimensions["lon"].size
            return {"time": slice(0, n_time), "lat": slice(0, n_lat), "lon": slice(0, n_lon)}

        valid = ~np.asarray(sample.mask)
        if not valid.any():
            raise ValueError(f"{src}: ugrdprs is fully masked; nothing to migrate")
        t_any = valid.any(axis=(1, 2))
        lat_any = valid.any(axis=(0, 2))
        lon_any = valid.any(axis=(0, 1))
        t_idx = np.where(t_any)[0]
        lat_idx = np.where(lat_any)[0]
        lon_idx = np.where(lon_any)[0]
        return {
            "time": slice(int(t_idx.min()), int(t_idx.max()) + 1),
            "lat":  slice(int(lat_idx.min()),  int(lat_idx.max())  + 1),
            "lon":  slice(int(lon_idx.min()),  int(lon_idx.max())  + 1),
        }


def _convert_v1_gfs(src: Path) -> xr.Dataset:
    """xarray-based path with mask-aware subsetting.

    Global archived GFS files store >99% masked cells outside the original
    download region. Subsetting to the valid bounding box first keeps the
    in-memory footprint of the xarray pipeline proportional to the actual
    data, not the global grid.
    """
    subset = _valid_subset_indices(src)
    raw = xr.open_dataset(src, engine="netcdf4", decode_times=False).isel(
        time=subset["time"], lat=subset["lat"], lon=subset["lon"],
    )
    raw = raw.assign_coords(time=_decode_gfs_time(raw["time"]))

    # v1 GFS files omit the _FillValue attribute but still store the CF-1.7
    # default float32 fill (9.9692e+36) at unset cells (e.g., above model top
    # or outside the original download bbox). The v1 reader picked these up
    # via netCDF4's auto-mask using that default; xarray here read them as
    # raw data. Convert to NaN explicitly so downstream interpolation does
    # the right thing and the new file's _FillValue=NaN actually applies.
    CF_DEFAULT_FILL_F32 = 9.9692099683868690e36
    for v in ("ugrdprs", "vgrdprs", "hgtprs", "tmpprs"):
        if v in raw.data_vars:
            raw[v] = raw[v].where(raw[v] < CF_DEFAULT_FILL_F32 * 0.5)

    rename_map = {
        "time": "valid_time",
        "lev": "pressure_level",
        "lat": "latitude",
        "lon": "longitude",
        "ugrdprs": "u",
        "vgrdprs": "v",
        "hgtprs": "z",
    }
    if "tmpprs" in raw.data_vars:
        rename_map["tmpprs"] = "t"
    ds = raw.rename(rename_map)

    # hgtprs (geopotential meters) -> z (geopotential m^2/s^2)
    ds["z"] = ds["z"] * G

    # lon [0, 360) -> [-180, 180), then sort ascending
    lon = ds["longitude"].values.astype(np.float64)
    lon = ((lon + 180.0) % 360.0) - 180.0
    ds = ds.assign_coords(longitude=lon).sortby("longitude")

    # lat ascending -> descending
    ds = ds.sortby("latitude", ascending=False)
    # pressure descending (no-op for native GFS but guarantees the invariant)
    ds = ds.sortby("pressure_level", ascending=False)

    epoch, attrs = _encode_valid_time(ds["valid_time"].values)
    ds = ds.assign_coords(valid_time=epoch)
    ds["valid_time"].attrs = attrs
    return _finalize(ds, FORMAT_V1_GFS)


def _convert_v1_era5(src: Path) -> xr.Dataset:
    raw = xr.open_dataset(src, engine="netcdf4")
    ds = raw.rename({"time": "valid_time", "level": "pressure_level"})
    # ascending hPa -> descending hPa (data reverses along the axis)
    ds = ds.sortby("pressure_level", ascending=False)
    ds = ds.sortby("longitude")
    ds = ds.sortby("latitude", ascending=False)

    epoch, attrs = _encode_valid_time(ds["valid_time"].values)
    ds = ds.assign_coords(valid_time=epoch)
    ds["valid_time"].attrs = attrs
    return _finalize(ds, FORMAT_V1_ERA5)


CONVERTERS: dict[str, Callable[[Path], xr.Dataset]] = {
    FORMAT_V1_GFS: _convert_v1_gfs,
    FORMAT_V1_ERA5: _convert_v1_era5,
}


# ---------------------------------------------------------------------------
# Write + filesystem orchestration
# ---------------------------------------------------------------------------

def _write_canonical(ds: xr.Dataset, dst: Path) -> None:
    encoding = {var: {"zlib": True, "complevel": 4} for var in ds.data_vars}
    encoding["valid_time"] = {"_FillValue": None, "dtype": "int64"}
    ds.to_netcdf(str(dst), engine="netcdf4", encoding=encoding)


def _backup_path(path: Path) -> Path:
    """Map foo.nc -> foo.v1.nc (alongside the original)."""
    return path.with_name(path.stem + ".v1.nc")


def migrate_file(path: Path, dry_run: bool = False) -> str:
    """Migrate one file in place. Returns a single-line status string.

    Status prefix vocabulary (stable; parseable by tooling):
      OK      - converted
      SKIP    - skipped (already v2, backup exists, or not a target)
      ERROR   - inspection or conversion failed
      DRY-RUN - would have done the action described
    """
    if not path.is_file():
        return f"SKIP {path}: not a file"
    if path.suffix != ".nc":
        return f"SKIP {path}: not a .nc file"
    if path.name.endswith(".v1.nc"):
        return f"SKIP {path}: looks like a .v1.nc backup"

    backup = _backup_path(path)
    if backup.exists():
        return f"SKIP {path}: backup {backup.name} already exists (previously migrated)"

    try:
        fmt = detect_format(path)
    except Exception as e:
        return f"ERROR {path}: cannot inspect ({e})"

    if fmt == FORMAT_V2:
        return f"SKIP {path}: already v2 canonical"
    if fmt == FORMAT_UNKNOWN:
        return f"ERROR {path}: format not recognized as v1 GFS or v1 ERA5"

    if dry_run:
        return f"DRY-RUN {path}: would convert from {fmt} (backup -> {backup.name})"

    tmp = path.with_name(path.stem + ".v2.tmp.nc")
    try:
        ds = CONVERTERS[fmt](path)
        try:
            _write_canonical(ds, tmp)
        finally:
            ds.close()
        path.rename(backup)
        tmp.rename(path)
        return f"OK {path}: converted from {fmt} (backup: {backup.name})"
    except Exception as e:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        return f"ERROR {path}: conversion failed ({e})"


def _iter_targets(paths: list[Path]):
    for p in paths:
        if p.is_dir():
            yield from sorted(p.glob("*.nc"))
        else:
            yield p


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m EarthSHAB.forecast_processing.migrate_v1",
        description=(
            "Migrate legacy EarthSHAB v1 forecast files to the v2 canonical "
            "netCDF schema (see docs/forecast-schema-v2.md). Backs the "
            "original up as <name>.v1.nc and writes the canonical file in "
            "place. Idempotent."
        ),
    )
    parser.add_argument(
        "targets", nargs="+", type=Path,
        help="One or more .nc files, or directories containing .nc files "
             "(non-recursive).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Report what would happen, but make no filesystem changes.",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print per-file status. Default mode shows a progress bar and "
             "prints only ERROR lines plus an end-of-run summary.",
    )
    args = parser.parse_args(argv)

    targets = list(_iter_targets(args.targets))
    counts = {"OK": 0, "SKIP": 0, "ERROR": 0, "DRY-RUN": 0}

    bar = tqdm(targets, unit="file", disable=args.verbose, leave=True)
    for tgt in bar:
        bar.set_postfix_str(tgt.name[:50], refresh=False)
        msg = migrate_file(tgt, dry_run=args.dry_run)
        prefix = msg.split(" ", 1)[0]
        counts[prefix] = counts.get(prefix, 0) + 1
        if args.verbose:
            print(msg)
        elif prefix == "ERROR":
            bar.write(msg)
    bar.close()

    print(
        f"\nMigration summary: {counts['OK']} migrated, "
        f"{counts['SKIP']} skipped, {counts['ERROR']} errors, "
        f"{counts['DRY-RUN']} would-migrate (dry-run)"
    )
    return 1 if counts["ERROR"] else 0


if __name__ == "__main__":
    sys.exit(main())
