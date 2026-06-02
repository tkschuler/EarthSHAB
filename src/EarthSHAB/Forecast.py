"""Forecast reader for the v2 canonical schema.

Single class for both GFS- and ERA5-sourced files. Source-specific quirks
were eliminated in Phases 1–4; this module replaces both ERA5.py and GFS.py.

Standard v2 file format (locked in by Phase 5)::

    dims:       valid_time, pressure_level, latitude, longitude
    data vars:  u, v, z, t  (z is geopotential m**2 s**-2; t may be absent)
    coords:     latitude descending, longitude in [-180, 180) ascending,
                pressure_level descending hPa,
                valid_time 'seconds since 1970-01-01' / proleptic_gregorian
    global:     Conventions = 'CF-1.7'
                institution = 'NOAA/NCEP (GFS)' or 'ECMWF (ERA5)'

Every file is a tight bounding-box subset — no full-world arrays with
mask-based subsetting. The reader slices the entire stored array.

Source provenance (self.source) is read from the ``institution`` global
attribute. When that attribute is missing or empty (the case for already-
migrated files predating Phase 5), the reader falls back to filename pattern::

    - basename containing 'gfs' (case-insensitive) -> 'GFS'
    - basename containing 'era5'                   -> 'ERA5'
    - otherwise                                    -> 'unknown'

Contributing authors: Craig Motell and Michael Rodriguez of NIWC Pacific
Edited and integrated: Tristan Schuler
"""
import math
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4
from geographiclib.geodesic import Geodesic
from pytz import timezone
from scipy import interpolate
from scipy.interpolate import CubicSpline
from termcolor import colored

import EarthSHAB.config_earth as config_earth


class ForecastFormatError(ValueError):
    """Raised when a forecast file is not in the v2 canonical schema.

    The message includes the exact migrate_v1 CLI invocation a user should
    run to upgrade the file in place. There is no auto-conversion in the
    reader by design (explicit > silent data mutation).
    """


# Signatures used to recognise legacy formats and emit a targeted migration
# message. Mirrors detection in EarthSHAB.forecast_processing.migrate_v1.
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
        f"{path}: detected {legacy}; EarthSHAB v2.0+ requires the canonical "
        "schema documented in docs/forecast-schema-v2.md. Convert this file "
        "in place with:\n"
        f"    python -m EarthSHAB.forecast_processing.migrate_v1 {path}\n"
        "The original will be backed up alongside it as <name>.v1.nc. "
        "See docs/migration-v2.md for details."
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


def _spline_uv(alt_m, h_ascending, u_ascending, v_ascending):
    """CubicSpline on u and v independently across the altitude profile.

    extrapolate=False; when alt_m is outside [h_min, h_max] falls back to
    np.interp (which clamps to endpoints), avoiding spline overshoot above
    the highest pressure level where the balloon often actually flies.

    h_ascending, u_ascending, v_ascending must be in ascending order of h.
    Duplicate h values (which can occur after fill_missing_data clamps
    NaN-extrapolation at the top of the profile to the last valid altitude)
    are filtered out — CubicSpline requires strictly increasing x.
    """
    h_arr = np.asarray(h_ascending)
    u_arr = np.asarray(u_ascending)
    v_arr = np.asarray(v_ascending)
    keep = np.concatenate(([True], np.diff(h_arr) > 0))
    h_s, u_s, v_s = h_arr[keep], u_arr[keep], v_arr[keep]

    h_lo, h_hi = h_s[0], h_s[-1]
    if alt_m < h_lo or alt_m > h_hi or len(h_s) < 4:
        return (
            float(np.interp(alt_m, h_s, u_s)),
            float(np.interp(alt_m, h_s, v_s)),
        )
    cs_u = CubicSpline(h_s, u_s, extrapolate=False)
    cs_v = CubicSpline(h_s, v_s, extrapolate=False)
    return float(cs_u(alt_m)), float(cs_v(alt_m))


class Forecast:
    """Reader for v2 canonical forecast files (GFS- or ERA5-sourced)."""

    def __init__(self, start_coord):
        self.dt = config_earth.simulation['dt']
        self.sim_time = config_earth.simulation['sim_time']

        forecast_path = config_earth.forecast['file']
        try:
            self.file = netCDF4.Dataset(forecast_path)
        except OSError:
            print(colored(f"Unable to locate netcdf file: {forecast_path}", "red"))
            sys.exit(1)

        _raise_if_v1(self.file, forecast_path)
        self.source = _detect_source(self.file, forecast_path)

        self.geod = Geodesic.WGS84
        self.start_coord = config_earth.simulation["start_coord"]
        self.start_time = config_earth.simulation['start_time']
        self.min_alt_m = self.start_coord['alt']

        time_arr = self.file.variables['valid_time']
        time_units = time_arr.units if hasattr(time_arr, 'units') else "seconds since 1970-01-01"
        time_calendar = time_arr.calendar if hasattr(time_arr, 'calendar') else "proleptic_gregorian"
        self.time_convert = netCDF4.num2date(time_arr[:], units=time_units, calendar=time_calendar)
        self.model_start_datetime = self.time_convert[0]
        self.model_end_datetime = self.time_convert[-1]

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

        print(colored("Forecast Information (Parsed from netcdf file):", "blue", attrs=['bold']))
        print(f"  source: {self.source}")
        print(f"  LAT RANGE: min: {self.LAT_LOW}  max: {self.LAT_HIGH}  size: {len(self.lat)}")
        print(f"  LON RANGE: min: {self.LON_LOW}  max: {self.LON_HIGH}  size: {len(self.lon)}")
        print(f"  cadence:   {self.resolution_hr:.2f} hr")
        print()

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

    # ------------------------------------------------------------------
    # Nearest-index helpers
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

    def get2NearestAltIdxs(self, h, alt_m):
        """Two enclosing altitude indices for bearing/speed interpolation.

        Index wrap from 0 to -1 is handled in interpolateBearing().
        """
        h_nearest = self.closestIdx(h, alt_m)
        if alt_m > h[h_nearest]:
            return h_nearest, h_nearest + 1
        return h_nearest - 1, h_nearest

    def getNearestAltbyIndex(self, int_hr_idx, lat_i, lon_i, alt_m):
        """Nearest altitude index at a specific (time, lat, lon) cell."""
        return self.closestIdx(self.hgtprs[int_hr_idx, :, lat_i, lon_i], alt_m)

    def getNearestAlt(self, hour_index, lat, lon, alt):
        """GFS-style signature: (hour_index, lat_deg, lon_deg, alt_m)."""
        lat_i = self.getNearestLatIdx(lat)
        lon_i = self.getNearestLonIdx(lon)
        return self.getNearestAltbyIndex(int(hour_index), lat_i, lon_i, alt)

    # ------------------------------------------------------------------
    # Interpolation
    # ------------------------------------------------------------------

    def fill_missing_data(self, data):
        """Linearly interpolate over missing entries in a 1-D profile column."""
        data = data.filled(np.nan) if np.ma.isMaskedArray(data) else np.asarray(data, dtype=float)
        nans, x = np.isnan(data), lambda z: z.nonzero()[0]
        if nans.any() and (~nans).any():
            data[nans] = np.interp(x(nans), x(~nans), data[~nans])
        return data

    def windVectorToBearing(self, u, v):
        bearing = np.arctan2(v, u)
        speed = np.power((np.power(u, 2) + np.power(v, 2)), .5)
        return [bearing, speed]

    def interpolateBearing(self, h, u, v, alt_m):
        """Linear interpolation between adjacent pressure levels on
        bearing+speed, with 0/360° axis-crossover correction.
        """
        h_idx0, h_idx1 = self.get2NearestAltIdxs(h, alt_m)
        if h_idx0 == -1:
            h_idx0 = 0
        if h_idx1 == 0:
            h_idx1 = -1

        bearing0, speed0 = self.windVectorToBearing(u[h_idx0], v[h_idx0])
        bearing1, speed1 = self.windVectorToBearing(u[h_idx1], v[h_idx1])

        angle1 = np.degrees(bearing0) % 360
        angle2 = np.degrees(bearing1) % 360
        if abs(angle2 - angle1) > 180:
            if angle2 > angle1:
                angle1 += 360
            else:
                angle2 += 360

        interp_speed = np.interp(alt_m, [h[h_idx0], h[h_idx1]], [speed0, speed1])
        interp_dir_deg = np.interp(alt_m, [h[h_idx0], h[h_idx1]], [angle1, angle2]) % 360
        return (interp_dir_deg, interp_speed)

    def interpolateBearingTime(self, bearing0, speed0, bearing1, speed1, hour_index):
        """Interpolate bearing+speed in time between two adjacent forecast hours."""
        angle1 = bearing0 % 360
        angle2 = bearing1 % 360
        if abs(angle2 - angle1) > 180:
            if angle2 > angle1:
                angle1 += 360
            else:
                angle2 += 360
        fp = [int(hour_index), int(hour_index) + 1]
        interp_speed = np.interp(hour_index, fp, [speed0, speed1])
        interp_dir_deg = np.interp(hour_index, fp, [angle1, angle2]) % 360
        return (interp_dir_deg, interp_speed)

    def wind_alt_Interpolate2(self, alt_m, diff_time, lat_idx, lon_idx):
        """Two-step interpolation: altitude across pressure levels at each of
        the two enclosing forecast times, then linear in time.

        Returns [u, v, u_diag, v_diag] where (u_diag, v_diag) is always the
        linear_full path, returned for diagnostic comparison regardless of
        which wind_interpolation method is configured.
        """
        hour_index = (diff_time.days * 24 + diff_time.seconds / 3600.0) / self.resolution_hr

        v_0 = self.fill_missing_data(self.vgdrps0[int(hour_index), :, lat_idx, lon_idx])
        u_0 = self.fill_missing_data(self.ugdrps0[int(hour_index), :, lat_idx, lon_idx])
        h0  = self.fill_missing_data(self.hgtprs[int(hour_index), :, lat_idx, lon_idx])

        v_1 = self.fill_missing_data(self.vgdrps0[int(hour_index) + 1, :, lat_idx, lon_idx])
        u_1 = self.fill_missing_data(self.ugdrps0[int(hour_index) + 1, :, lat_idx, lon_idx])
        h1  = self.fill_missing_data(self.hgtprs[int(hour_index) + 1, :, lat_idx, lon_idx])

        # Diagnostic baseline: linear_full path on u,v across the full
        # altitude profile, then linear in time.
        u0_lf = np.interp(alt_m, h0, u_0)
        v0_lf = np.interp(alt_m, h0, v_0)
        u1_lf = np.interp(alt_m, h1, u_1)
        v1_lf = np.interp(alt_m, h1, v_1)
        fp = [int(hour_index), int(hour_index) + 1]
        u_lf = np.interp(hour_index, fp, [u0_lf, u1_lf])
        v_lf = np.interp(hour_index, fp, [v0_lf, v1_lf])

        method = config_earth.forecast.get('wind_interpolation', 'linear_neighbors')
        if method == 'linear_neighbors':
            bearing_t0, speed_t0 = self.interpolateBearing(h0, u_0, v_0, alt_m)
            bearing_t1, speed_t1 = self.interpolateBearing(h1, u_1, v_1, alt_m)
            bearing_interpolated, speed_interpolated = self.interpolateBearingTime(
                bearing_t0, speed_t0, bearing_t1, speed_t1, hour_index)
            u = speed_interpolated * np.cos(np.radians(bearing_interpolated))
            v = speed_interpolated * np.sin(np.radians(bearing_interpolated))
        elif method == 'linear_full':
            u, v = u_lf, v_lf
        elif method == 'spline_full':
            u0_sp, v0_sp = _spline_uv(alt_m, h0, u_0, v_0)
            u1_sp, v1_sp = _spline_uv(alt_m, h1, u_1, v_1)
            u = np.interp(hour_index, fp, [u0_sp, u1_sp])
            v = np.interp(hour_index, fp, [v0_sp, v1_sp])
        else:
            raise ValueError(
                f"Unknown wind_interpolation: {method!r}. "
                "Expected 'linear_neighbors', 'linear_full', or 'spline_full'."
            )

        return [u, v, u_lf, v_lf]

    # ------------------------------------------------------------------
    # Top-level coordinate propagation
    # ------------------------------------------------------------------

    def getNewCoord(self, coord, dt):
        """Advance a balloon coord by `dt` seconds under the local wind field.

        :param coord: dict with keys lat, lon, alt, timestamp
        :param dt:    integration interval in seconds
        :returns:     [lat_new, lon_new, x_wind_vel, y_wind_vel,
                       x_wind_vel_diag, y_wind_vel_diag, bearing,
                       nearest_lat, nearest_lon, nearest_alt]
        """
        diff_time = coord["timestamp"] - self.model_start_datetime
        hour_index = (diff_time.days * 24 + diff_time.seconds / 3600.0) / self.resolution_hr
        int_hr_idx = int(hour_index)

        lat_idx = self.getNearestLatIdx(coord["lat"])
        lon_idx = self.getNearestLonIdx(coord["lon"])
        z = self.getNearestAltbyIndex(int_hr_idx, lat_idx, lon_idx, coord["alt"])

        try:
            x_wind_vel, y_wind_vel, x_wind_vel_old, y_wind_vel_old = \
                self.wind_alt_Interpolate2(coord['alt'], diff_time, lat_idx, lon_idx)
            if x_wind_vel is None:
                return [None] * 10
        except Exception:
            print(colored(
                "Mismatch with simulation and forecast timstamps. Check simulation "
                "start time and/or download a new forecast.", "red"))
            sys.exit(1)

        bearing = math.degrees(math.atan2(y_wind_vel, x_wind_vel))
        bearing = 90 - bearing  # 90-degree rotation: wind components → compass bearing
        d = math.pow((math.pow(y_wind_vel, 2) + math.pow(x_wind_vel, 2)), .5) * dt
        g = self.geod.Direct(coord["lat"], coord["lon"], bearing, d)

        if (g['lat2'] < self.LAT_LOW or g['lat2'] > self.LAT_HIGH
                or g['lon2'] < self.LON_LOW or g['lon2'] > self.LON_HIGH):
            print(colored("WARNING: Trajectory is out of bounds of downloaded netcdf forecast", "yellow"))

        if coord["alt"] <= self.min_alt_m:
            # Pin to the launch coordinate while the balloon hasn't lifted off.
            return [coord['lat'], coord['lon'], x_wind_vel, y_wind_vel,
                    x_wind_vel_old, y_wind_vel_old, bearing,
                    self.lat[lat_idx], self.lon[lon_idx],
                    self.hgtprs[0, z, lat_idx, lon_idx]]
        return [g['lat2'], g['lon2'], x_wind_vel, y_wind_vel,
                x_wind_vel_old, y_wind_vel_old, bearing,
                self.lat[lat_idx], self.lon[lon_idx],
                self.hgtprs[0, z, lat_idx, lon_idx]]

    def close(self):
        """Close the underlying netCDF dataset.

        CPython's __del__ on netCDF4.Dataset is unreliable enough that
        long-running batch evaluators (evaluation/run_batch.py) accumulate
        HDF5 chunk-cache state across launches and segfault around the
        6th–7th iteration. Callers should `close()` explicitly when done.
        Idempotent.
        """
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
