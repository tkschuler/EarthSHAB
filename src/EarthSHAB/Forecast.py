"""Forecast reader for the v2 canonical netcdf forecast schema.

Single class for both GFS- and ERA5-sourced netcdf files.
This module replaces both ERA5.py and GFS.py.

Standard v2 file format::

    dims:       valid_time, pressure_level, latitude, longitude
    data vars:  u, v, z, t  (z is geopotential m**2 s**-2; t may be absent)
    coords:     latitude descending, longitude in [-180, 180) ascending,
                pressure_level descending hPa,
                valid_time 'seconds since 1970-01-01' / proleptic_gregorian
    global:     Conventions = 'CF-1.7'
                institution = 'NOAA/NCEP (GFS)' or 'ECMWF (ERA5)'

Every netcdf file is a bounded forecast subset (no full-world arrays with
mask-based subsetting), although full world forecasts can still be downloaded
by applying full world lat/lon bounds. The reader slices the entire stored array.

Hot-swappable wind-lookup design
--------------------------------
Three orthogonal axes select the wind field, each a constructor kwarg with a
``config.forecast[...]`` fallback:

  * ``backend``            -- 'xarray' (slow/original) | 'numpy' | 'numba' (default)
  * ``wind_interpolation`` -- 'bilinear' | 'linear_neighbors' | 'linear_full' | 'spline_full'
  * ``advection``          -- 'geodesic' (WGS84) | 'tangent_plane' (fixed-anchor spherical)

The skeleton (this class) owns the public API -- ``getNewCoord`` (scalar) and
``getNewCoords`` (batch) -- plus advection. The **tier** swaps only the wind
lookup, implemented by the backend classes in ``utils.forecast_backends``; the
vertical altitude/time math is shared (``utils.wind_interp``). ``getNewCoord``
is a thin length-1 wrapper over ``getNewCoords`` so the two share one code path.

The shape mirrors EarthSHAB's ``Forecast.py`` (same ``Forecast(start_coord)``
construction and ``getNewCoord`` return) so this is upstreamable.
"""
import math
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4
from geographiclib.geodesic import Geodesic
from pytz import timezone
from termcolor import colored

import EarthSHAB.config as config
import EarthSHAB.utils.transform as transform
import EarthSHAB.utils.forecast_backends as forecast_backends


class ForecastFormatError(ValueError):
    """Raised when a forecast file is not in the v2 canonical schema.

    The message includes the exact migrate_v1 CLI invocation a user should
    run to upgrade the file in place. There is no auto-conversion in the
    reader by design (explicit > silent data mutation).
    """


# Signatures used to recognise legacy formats and emit a targeted migration
# message.
_V1_GFS_VARS = {"ugrdprs", "vgrdprs", "hgtprs"}
_V1_ERA5_DIMS = {"time", "level", "latitude", "longitude"}
_V2_REQUIRED_DIMS = {"valid_time", "pressure_level", "latitude", "longitude"}


def _raise_if_v1(file: netCDF4.Dataset, path: str) -> None:
    vars_ = set(file.variables.keys())
    dims_ = set(file.dimensions.keys())

    if _V2_REQUIRED_DIMS.issubset(dims_):
        return  # canonical v2 — proceed

    if _V1_GFS_VARS.issubset(vars_):
        legacy = "v1 GFS (ugrdprs/vgrdprs/hgtprs variables)"
    elif _V1_ERA5_DIMS.issubset(dims_):
        legacy = "v1 ERA5 (time/level dim naming)"
    else:
        legacy = "unrecognized layout (neither v2 canonical nor a known v1 format)"

    raise ForecastFormatError(
        f"{path}: detected {legacy}; HAB-COM requires the canonical v2 schema "
        "documented in forecast-schema-v2.md. Re-download the forecast with "
        "collectors/gfs_collector.py (which writes v2 directly), or convert the "
        "file to the v2 layout (dims: valid_time, pressure_level, latitude, "
        "longitude; vars: u, v, z, t with z as geopotential m**2 s**-2)."
    )


def _detect_source(file: netCDF4.Dataset, path: str) -> str:
    """Resolve self.source from `institution` attr, falling back to filename.

    Returns 'GFS', 'ERA5', or 'unknown'.
    """
    inst = ""
    if "institution" in file.ncattrs():
        inst = (file.getncattr("institution") or "").upper()
    if "NOAA" in inst or "GFS" in inst or "NCEP" in inst:
        return "GFS"
    if "ECMWF" in inst or "ERA5" in inst:
        return "ERA5"

    name = Path(path).name.lower()
    if "gfs" in name:
        return "GFS"
    if "era5" in name:
        return "ERA5"
    return "unknown"


def _epoch_seconds(times):
    """Convert num2date datetimes (cftime or stdlib) to int64 unix seconds."""
    epoch = datetime(1970, 1, 1)
    out = np.empty(len(times), dtype=np.int64)
    for i, t in enumerate(times):
        out[i] = int((datetime(t.year, t.month, t.day, t.hour, t.minute,
                               int(t.second)) - epoch).total_seconds())
    return out


class Forecast:
    """Reader for v2 canonical forecast files (GFS- or ERA5-sourced)."""

    def __init__(self, start_coord=None, *, forecast_file=None, start_time=None,
                 dt=None, sim_time=None, wind_interpolation=None, advection=None,
                 backend=None):
        # Each parameter falls back to config.py when not explicitly provided.
        # The Monte Carlo runners pass overrides so they can drive Forecast.py
        # without mutating global config.
        self.dt = dt if dt is not None else config.simulation['dt']
        self.sim_time = sim_time if sim_time is not None else config.simulation['sim_time']

        # --- the three orthogonal wind-field selectors ---
        self.wind_interpolation = (
            wind_interpolation if wind_interpolation is not None
            else config.forecast.get('wind_interpolation', 'linear_neighbors')
        )
        # Advection model used by getNewCoord(s):
        #   'geodesic'      -> WGS84 ellipsoid step (geographiclib), re-anchored
        #                      at the balloon each step.
        #   'tangent_plane' -> spherical tangent-plane step about a FIXED anchor
        #                      (the launch coordinate).
        self.advection = (
            advection if advection is not None
            else config.forecast.get('advection', 'geodesic')
        )
        # Wind-lookup tier: 'xarray' (slow/original) | 'numpy' | 'numba'.
        self.backend_name = (
            backend if backend is not None
            else config.forecast.get('backend', 'numba')
        )

        forecast_path = forecast_file if forecast_file is not None else config.forecast['file']
        try:
            self.file = netCDF4.Dataset(forecast_path)
        except OSError:
            print(colored(f"Unable to locate netcdf file: {forecast_path}", "red"))
            sys.exit(1)

        _raise_if_v1(self.file, forecast_path)
        self.source = _detect_source(self.file, forecast_path)

        self.geod = Geodesic.WGS84
        self.start_coord = start_coord if start_coord is not None else config.simulation["start_coord"]
        self.start_time = start_time if start_time is not None else config.simulation['start_time']
        self.min_alt_m = self.start_coord['alt']

        # Fixed tangent-plane anchor (the launch coordinate). Only used when
        # advection == 'tangent_plane'.
        self.central_lat = float(self.start_coord['lat'])
        self.central_lon = transform.wrap_lon(float(self.start_coord['lon']))
        self._cos_central = math.cos(math.radians(self.central_lat))

        time_arr = self.file.variables['valid_time']
        time_units = time_arr.units if hasattr(time_arr, 'units') else "seconds since 1970-01-01"
        time_calendar = time_arr.calendar if hasattr(time_arr, 'calendar') else "proleptic_gregorian"
        self.time_convert = netCDF4.num2date(time_arr[:], units=time_units, calendar=time_calendar)
        self.model_start_datetime = self.time_convert[0]
        self.model_end_datetime = self.time_convert[-1]
        self._time_unix_s = _epoch_seconds(self.time_convert)

        # Cadence (hours per time step) derived from the file itself rather
        # than a config key — ERA5 is typically 1 hr, GFS typically 3 hr.
        if len(self.time_convert) >= 2:
            self.resolution_hr = (
                self.time_convert[1] - self.time_convert[0]
            ).total_seconds() / 3600.0
        else:
            self.resolution_hr = 1.0

        # v2 standard format: file IS a tight bounding-box subset. No mask-
        # based subsetting, no index bookkeeping. Whole-array reads.
        g = 9.80665  # convert geopotential (m^2/s^2) to geopotential height (m)
        self.lat = np.asarray(self.file.variables['latitude'][:], dtype=float)
        self.lon = np.asarray(self.file.variables['longitude'][:], dtype=float)
        self.levels = self.file.variables['pressure_level'][:]
        self.ugdrps0 = self.file.variables['u'][:]
        self.vgdrps0 = self.file.variables['v'][:]
        self.hgtprs = self.file.variables['z'][:] / g

        # Bounds for trajectory containment checks / region polygons.
        self.LAT_LOW  = float(np.min(self.lat))
        self.LAT_HIGH = float(np.max(self.lat))
        self.LON_LOW  = float(np.min(self.lon))
        self.LON_HIGH = float(np.max(self.lon))

        # Actual horizontal grid spacing, derived from the file itself. This is
        # the resolution the simulation truly runs at, independent of any config.
        self.resolution_deg = (
            abs(float(self.lat[1] - self.lat[0])) if len(self.lat) >= 2 else None
        )

        print(colored("Forecast Information (Parsed from netcdf file):", "blue", attrs=['bold']))
        print(f"  source: {self.source}")
        print(f"  backend: {self.backend_name}  wind: {self.wind_interpolation}  advection: {self.advection}")
        print(f"  LAT RANGE: min: {self.LAT_LOW}  max: {self.LAT_HIGH}  size: {len(self.lat)}")
        print(f"  LON RANGE: min: {self.LON_LOW}  max: {self.LON_HIGH}  size: {len(self.lon)}")
        if self.resolution_deg is not None:
            print(f"  resolution: {self.resolution_deg:g}° grid")
        print(f"  cadence:   {self.resolution_hr:.2f} hr")
        print()

        # Guard against a silent config/file mismatch in grid spacing or cadence.
        # config.forecast['file'] is normally built from netcdf_gfs['res'] and
        # ['step_hours'], so a stale knob can quietly load a file with a different
        # grid or cadence than declared — the sim still runs, but not on the data
        # the config describes. Compare the file's actual spacing/cadence against
        # the declared config values and flag any drift.
        declared_res  = config.netcdf_gfs.get("res")
        declared_step = config.netcdf_gfs.get("step_hours")

        mismatches = []
        if (declared_res is not None and self.resolution_deg is not None
                and abs(self.resolution_deg - float(declared_res)) > 1e-3):
            mismatches.append(
                f"grid spacing: config netcdf_gfs['res'] = {declared_res}° but the "
                f"file is a {self.resolution_deg:g}° grid"
            )
        if (declared_step is not None
                and abs(self.resolution_hr - float(declared_step)) > 1e-3):
            mismatches.append(
                f"cadence: config netcdf_gfs['step_hours'] = {declared_step} h but "
                f"the file is {self.resolution_hr:g} h"
            )

        if mismatches:
            detail = "\n  - ".join(mismatches)
            if self.source == "GFS":
                # GFS filenames are built from these knobs, so a mismatch means
                # config.forecast['file'] points at the wrong forecast — fail hard
                # rather than silently run on the wrong data.
                raise ValueError(
                    f"Forecast config/file mismatch:\n  - {detail}\n"
                    f"  file: {forecast_path}\n"
                    f"These config keys build the GFS forecast filename, so this "
                    f"usually means config.forecast['file'] points at the wrong "
                    f"file. Update netcdf_gfs['res']/['step_hours'] to match the "
                    f"file, or point forecast['file'] at the matching forecast."
                )
            elif self.source == "ERA5":
                # ERA5 files are supplied manually, not produced from the GFS
                # download knobs, so a mismatch is informational rather than fatal
                # — but still surfaced so a stale config doesn't silently
                # misdescribe the run. The sim uses the file's actual values.
                print(colored(
                    f"[WARNING] The loaded ERA5 forecast does not match the "
                    f"resolution/cadence declared in config.netcdf_gfs:\n"
                    f"  - {detail}\n"
                    f"  file: {forecast_path}\n"
                    f"  The simulation will use the file's actual resolution and "
                    f"cadence shown above. Update netcdf_gfs if this is unexpected.",
                    "yellow"))

        # Check requested simulation window fits inside the forecast.
        desired_simulation_end_time = self.start_time + timedelta(hours=self.sim_time)
        diff_time = (self.model_end_datetime - self.start_time).total_seconds()
        print("Sim start time: ", self.start_time)
        print("NetCDF end time:", self.model_end_datetime)
        print("Max sim runtime:", diff_time // 3600, "hours")
        print("Des sim runtime:", self.sim_time, "hours")
        print()

        if not desired_simulation_end_time <= self.model_end_datetime:
            print(colored(
                f"Desired simulation run time of {self.sim_time} hours is out of "
                "bounds of downloaded forecast. Check simulation start time "
                "and/or download a new forecast.", "red"))
            sys.exit()

        # Build the selected wind-lookup backend. All backends take the same
        # kwargs and ignore what they don't need (see utils.forecast_backends).
        self.backend = forecast_backends.make_backend(
            self.backend_name,
            path=forecast_path,
            lat=self.lat, lon=self.lon, hgt_m=self.hgtprs,
            u=self.ugdrps0, v=self.vgdrps0, time_unix_s=self._time_unix_s,
        )

    # ------------------------------------------------------------------
    # Nearest-index helpers (used for the diagnostic return fields)
    # ------------------------------------------------------------------

    def closestIdx(self, arr, k):
        """Given an ordered array and a value, return the index of the closest item."""
        return min(range(len(arr)), key=lambda i: abs(arr[i] - k))

    # Back-compat alias (GFS used `closest`, ERA5 used `closestIdx`).
    closest = closestIdx

    def getNearestLatIdx(self, lat, *_unused):
        return self.closestIdx(self.lat, float(lat))

    def getNearestLonIdx(self, lon, *_unused):
        return self.closestIdx(self.lon, float(lon))

    # Back-compat aliases (GFS naming).
    getNearestLat = getNearestLatIdx
    getNearestLon = getNearestLonIdx

    def getNearestAltbyIndex(self, int_hr_idx, lat_i, lon_i, alt_m):
        """Nearest altitude index at a specific (time, lat, lon) cell."""
        return self.closestIdx(self.hgtprs[int_hr_idx, :, lat_i, lon_i], alt_m)

    def getNearestAlt(self, hour_index, lat, lon, alt):
        """GFS-style signature: (hour_index, lat_deg, lon_deg, alt_m)."""
        lat_i = self.getNearestLatIdx(lat)
        lon_i = self.getNearestLonIdx(lon)
        return self.getNearestAltbyIndex(int(hour_index), lat_i, lon_i, alt)

    def fill_missing_data(self, data):
        """Linearly interpolate over missing entries in a 1-D profile column."""
        data = data.filled(np.nan) if np.ma.isMaskedArray(data) else np.asarray(data, dtype=float)
        nans, x = np.isnan(data), lambda z: z.nonzero()[0]
        if nans.any() and (~nans).any():
            data[nans] = np.interp(x(nans), x(~nans), data[~nans])
        return data

    # ------------------------------------------------------------------
    # Time / advection helpers
    # ------------------------------------------------------------------

    def _hour_index(self, timestamp):
        """Fractional forecast-hour index for a timestamp."""
        diff = timestamp - self.model_start_datetime
        return (diff.days * 24 + diff.seconds / 3600.0) / self.resolution_hr

    def _advect(self, lats, lons, u, v, dt):
        """Advance (lats, lons) under wind (u, v) for dt seconds → (lat_new, lon_new).

        Both branches are vectorized over members; geodesic loops member-wise
        because geographiclib isn't array-aware.
        """
        if self.advection == 'tangent_plane':
            # Spherical tangent-plane step about the fixed launch anchor (exact
            # same arithmetic as utils.transform.{latlon_to_meters,
            # meters_to_latlon}_spherical, minus their NaN guards).
            cosc = self._cos_central
            cx = (lons - self.central_lon) * 111000.0 * cosc
            cy = (lats - self.central_lat) * 111000.0
            x_new = cx + u * dt
            y_new = cy + v * dt
            lat_new = np.clip(self.central_lat + y_new / 111320.0, -89.999999, 89.999999)
            lon_new = ((self.central_lon + x_new / (111320.0 * cosc) + 180.0) % 360.0) - 180.0
            return lat_new, lon_new

        # Geodesic: WGS84 ellipsoid step, re-anchored at the balloon each step.
        M = len(lats)
        lat_new = np.empty(M, dtype=np.float64)
        lon_new = np.empty(M, dtype=np.float64)
        for m in range(M):
            bearing = 90.0 - math.degrees(math.atan2(v[m], u[m]))
            d = math.hypot(u[m], v[m]) * dt
            g = self.geod.Direct(float(lats[m]), float(lons[m]), bearing, d)
            lat_new[m] = g['lat2']
            lon_new[m] = g['lon2']
        return lat_new, lon_new

    # ------------------------------------------------------------------
    # Top-level coordinate propagation (batch + scalar)
    # ------------------------------------------------------------------

    def getNewCoords(self, lats, lons, alts, timestamp, dt, *, diagnostic=True):
        """Advance a *batch* of coords by `dt` seconds under the local wind field.

        :param lats/lons/alts: 1-D arrays (length M) of current positions
        :param timestamp: shared datetime for all members
        :param dt: integration interval in seconds
        :param diagnostic: also compute the linear_full diagnostic winds and the
            nearest-cell diagnostic fields (needed for the full scalar return /
            EarthSHAB contract). Set False in hot batch loops to skip the extra
            wind sample + per-member nearest-index lookups.
        :returns: list of 10 length-M arrays mirroring getNewCoord's 10-tuple::
            [lat_new, lon_new, x_wind, y_wind, x_wind_diag, y_wind_diag,
             bearing, nearest_lat, nearest_lon, nearest_alt]
        """
        lats = np.asarray(lats, dtype=np.float64)
        lons = np.asarray(lons, dtype=np.float64)
        alts = np.asarray(alts, dtype=np.float64)
        M = lats.size

        hour_index = self._hour_index(timestamp)

        # (1) wind lookup — the swappable tier.
        u, v = self.backend.sample_winds(self.wind_interpolation, lats, lons, alts, hour_index)

        # diagnostic baseline: always the linear_full path (matches the
        # original wind_alt_Interpolate2 [u, v, u_lf, v_lf] contract).
        if diagnostic:
            if self.wind_interpolation == 'linear_full':
                u_lf, v_lf = u, v
            else:
                u_lf, v_lf = self.backend.sample_winds('linear_full', lats, lons, alts, hour_index)
        else:
            u_lf = np.full(M, np.nan)
            v_lf = np.full(M, np.nan)

        bearing = 90.0 - np.degrees(np.arctan2(v, u))  # wind components → compass bearing

        # (2) advection.
        lat_new, lon_new = self._advect(lats, lons, u, v, dt)

        # Out-of-bounds warning (once, on the advected position).
        oob = ((lat_new < self.LAT_LOW) | (lat_new > self.LAT_HIGH)
               | (lon_new < self.LON_LOW) | (lon_new > self.LON_HIGH))
        if oob.any():
            print(colored("WARNING: Trajectory is out of bounds of downloaded netcdf forecast", "yellow"))

        # Launch pin: members at/below the launch altitude stay put.
        pinned = alts <= self.min_alt_m
        lat_out = np.where(pinned, lats, lat_new)
        lon_out = np.where(pinned, lons, lon_new)

        # (3) nearest-cell diagnostic fields.
        if diagnostic:
            int_hr = int(hour_index)
            nlat = np.empty(M, dtype=np.float64)
            nlon = np.empty(M, dtype=np.float64)
            nalt = np.empty(M, dtype=np.float64)
            for m in range(M):
                li = self.getNearestLatIdx(lats[m])
                oi = self.getNearestLonIdx(lons[m])
                zi = self.getNearestAltbyIndex(int_hr, li, oi, alts[m])
                nlat[m] = self.lat[li]
                nlon[m] = self.lon[oi]
                nalt[m] = self.hgtprs[0, zi, li, oi]
        else:
            nlat = np.full(M, np.nan)
            nlon = np.full(M, np.nan)
            nalt = np.full(M, np.nan)

        return [lat_out, lon_out, u, v, u_lf, v_lf, bearing, nlat, nlon, nalt]

    def getNewCoord(self, coord, dt):
        """Advance a single balloon coord by `dt` seconds (EarthSHAB-shaped).

        Thin length-1 wrapper over getNewCoords so scalar and batch share one
        code path (and are bit-identical at M=1).

        :param coord: dict with keys lat, lon, alt, timestamp
        :param dt:    integration interval in seconds
        :returns:     [lat_new, lon_new, x_wind_vel, y_wind_vel,
                       x_wind_vel_diag, y_wind_vel_diag, bearing,
                       nearest_lat, nearest_lon, nearest_alt]
        """
        try:
            res = self.getNewCoords(
                np.array([coord['lat']], dtype=np.float64),
                np.array([coord['lon']], dtype=np.float64),
                np.array([coord['alt']], dtype=np.float64),
                coord['timestamp'], dt)
        except Exception:
            print(colored(
                "Mismatch with simulation and forecast timstamps. Check simulation "
                "start time and/or download a new forecast.", "red"))
            sys.exit(1)
        return [field[0] for field in res]

    # ------------------------------------------------------------------
    # Back-compat scalar interpolation entry (parity with EarthSHAB's API)
    # ------------------------------------------------------------------

    def wind_alt_Interpolate2(self, alt_m, diff_time, lat_idx, lon_idx, lat=None, lon=None):
        """Scalar wind interpolation -> [u, v, u_diag, v_diag].

        Kept for API parity with EarthSHAB's Forecast.py. Routes through the
        selected backend; `(u_diag, v_diag)` is always the linear_full path.
        For 'bilinear' the actual `lat`/`lon` are required (nearest-cell methods
        fall back to the grid coordinates of the supplied indices).
        """
        hour_index = (diff_time.days * 24 + diff_time.seconds / 3600.0) / self.resolution_hr
        if lat is None:
            lat = float(self.lat[lat_idx])
        if lon is None:
            lon = float(self.lon[lon_idx])
        lats = np.array([lat], dtype=np.float64)
        lons = np.array([lon], dtype=np.float64)
        alts = np.array([alt_m], dtype=np.float64)
        u, v = self.backend.sample_winds(self.wind_interpolation, lats, lons, alts, hour_index)
        if self.wind_interpolation == 'linear_full':
            u_lf, v_lf = u, v
        else:
            u_lf, v_lf = self.backend.sample_winds('linear_full', lats, lons, alts, hour_index)
        return [float(u[0]), float(v[0]), float(u_lf[0]), float(v_lf[0])]

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def close(self):
        """Close the underlying netCDF dataset (and any backend dataset).

        CPython's __del__ on netCDF4.Dataset is unreliable enough that
        long-running batch evaluators accumulate HDF5 chunk-cache state across
        launches and segfault. Callers should `close()` explicitly. Idempotent.
        """
        b = getattr(self, "backend", None)
        if b is not None:
            try:
                b.close()
            except Exception:
                pass
        f = getattr(self, "file", None)
        if f is not None:
            try:
                f.close()
            except Exception:
                pass
            self.file = None

    def __del__(self):
        # Best-effort cleanup; exceptions during interpreter shutdown swallowed.
        try:
            self.close()
        except Exception:
            pass

    def datetime_epoch(self, timedate_obj):
        """Convert a datetime (or ISO string) to GMT epoch seconds."""
        gmt = timezone('GMT')
        try:
            dt = timedate_obj.strftime("%Y-%m-%d %H:%M:%S")
            naive_ts = datetime.strptime(dt, "%Y-%m-%d %H:%M:%S")
        except AttributeError:
            naive_ts = datetime.strptime(timedate_obj, "%Y-%m-%d %H:%M:%S")
        local_ts = gmt.localize(naive_ts)
        return int(local_ts.timestamp())
