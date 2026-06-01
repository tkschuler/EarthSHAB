"""
This is file generates a 3d windrose plot for a particular coordinate and timestamp.
The polar plot displays information on wind speed and direction  at
various altitudes in a visual format

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from termcolor import colored
from scipy.interpolate import CubicSpline
import sys

import EarthSHAB.config_earth as config_earth
from EarthSHAB.Forecast import Forecast, _spline_uv

class Windmap:
    def __init__(self):

        self.forecast = Forecast(config_earth.simulation['start_coord'])
        self.file = self.forecast.file
        self.source = self.forecast.source

        self.start_time = config_earth.simulation["start_time"]
        self.coord = config_earth.simulation['start_coord']
        self.LAT = self.coord["lat"]
        self.LON = self.coord["lon"]

        # v2 canonical schema: same variable names regardless of source.
        self.lat = self.file.variables['latitude'][:]
        self.lon = self.file.variables['longitude'][:]
        self.levels = self.file.variables['pressure_level'][:]
        self.vgrdprs = self.file.variables['v']
        self.ugrdprs = self.file.variables['u']
        self.hgtprs = self.file.variables['z']
        try:
            self.tmpprs = self.file.variables['t']
        except KeyError:
            print(colored("Temperature not present in this forecast", "yellow"))
            self.tmpprs = None
        self.hour_index, self.new_timestamp = self.getHourIndex(self.start_time)

    def closest(self,arr, k):
        return min(range(len(arr)), key = lambda i: abs(arr[i]-k))

    def time_in_range(self,start, end, x):
        """Return true if x is in the range [start, end]"""
        if start <= end:
            return start <= x <= end
        else:
            return start <= x or x <= end

    def getHourIndex(self,start_time):
        times = self.forecast.time_convert
        #Check is simulation start time is within netcdf file
        if not self.time_in_range(times[0], times[-1], self.start_time):
            print(colored("Simulation start time " + str(self.start_time) + " is not within netcdf timerange of " + str(times[0]) + " - " + str(times[-1]) , "red"))
            sys.exit()
        else:
            print(colored("Simulation start time " + str(self.start_time) + " is within netcdf timerange of " + str(times[0]) + " - " + str(times[-1]) , "green"))

        #Find closest time using lambda function:
        closest_time = min(times, key=lambda sub: abs(sub - start_time))

        #Need to write some exception handling code
        #Find the corresponding index to a matching key (closest time) in the array (full list of timestamps in netcdf)
        hour_index = [i for i ,e in enumerate(times) if e == closest_time][0]

        return hour_index, closest_time


    def windVectorToBearing(self, u, v, h):
        """ Converts U-V wind data at specific heights to angular and radial
        components for polar plotting.

        :param u: U-Vector Wind Component from Forecast
        :type u: float64 array
        :param v: V-Vector Wind Component from Forecast
        :type v: float64 array
        :param h: Corresponding Converted Altitudes (m) from Forecast
        :type h: float64 array
        :returns: Array of bearings, radius, colors, and color map for plotting
        :rtype: array

        """

        # Calculate altitude
        bearing = np.arctan2(v,u)
        bearing = np.unwrap(bearing)
        r = np.power((np.power(u,2)+np.power(v,2)),.5)

        # Set up Color Bar
        colors = h
        cmap=mpl.colors.ListedColormap(colors)

        return [bearing, r , colors, cmap]


    def plotWind2(self, hour_index, lat, lon, num_interpolations=100, interpolation='linear'):
        """ Plots a 3D windrose by interpolating wind speed and direction between pressure levels.

        :param hour_index: Time index from forecast file
        :type hour_index: int
        :param lat: Latitude coordinate [deg]
        :type lat: float
        :param lon: Longitude coordinate [deg]
        :type lon: float
        :param num_interpolations: Points to interpolate between adjacent pressure levels (linear)
            or total sample points across the full altitude range (spline)
        :type num_interpolations: int
        :param interpolation: ``'linear'`` interpolates speed and direction between adjacent pressure
            levels with angle-wrap correction (default); ``'spline'`` fits a CubicSpline across all
            pressure levels for smoother, continuous curves
        :type interpolation: str
        :returns: 3D windrose plot
        :rtype: matplotlib.plot

        """

        lat_i = self.forecast.getNearestLatIdx(lat)
        lon_i = self.forecast.getNearestLonIdx(lon)
        # self.forecast.hgtprs is already geopotential / g (meters); use it directly.
        ugrd, vgrd, hgt = self.forecast.ugdrps0, self.forecast.vgdrps0, self.forecast.hgtprs

        u = ugrd[hour_index, :, lat_i, lon_i]
        v = vgrd[hour_index, :, lat_i, lon_i]
        h = hgt[hour_index, :, lat_i, lon_i]

        # Remove missing data
        if np.ma.isMaskedArray(u):
            u = u.filled(np.nan)
        if np.ma.isMaskedArray(v):
            v = v.filled(np.nan)
        nans = ~np.isnan(u)
        u = u[nans]
        v = v[nans]
        h = np.asarray(h)[nans]

        bearing, r, _, _ = self.windVectorToBearing(u, v, h)

        if interpolation == 'spline':
            # Fit CubicSpline across all pressure levels on the unwrapped bearing and speed
            h_new = np.linspace(h[0], h[-1], num_interpolations * (len(h) - 1))
            cs_r = CubicSpline(h, r)
            cs_bearing = CubicSpline(h, bearing)  # bearing is already np.unwrap'd
            interpolated_speeds = cs_r(h_new)
            interpolated_directions_deg = np.degrees(cs_bearing(h_new)) % 360
            interpolated_altitudes = h_new
        else:
            # Linear interpolation between adjacent pressure levels with angle-wrap correction
            interpolated_altitudes = []
            interpolated_speeds = []
            interpolated_directions_deg = []

            for i in range(len(h) - 1):
                angle1 = np.degrees(bearing[i]) % 360
                angle2 = np.degrees(bearing[i + 1]) % 360

                if abs(angle2 - angle1) > 180:
                    if angle2 > angle1:
                        angle1 += 360
                    else:
                        angle2 += 360

                for j in range(num_interpolations + 1):
                    alpha = j / num_interpolations
                    interp_alt = h[i] + alpha * (h[i + 1] - h[i])
                    interpolated_altitudes.append(interp_alt)
                    interpolated_speeds.append(np.interp(interp_alt, [h[i], h[i + 1]], [r[i], r[i + 1]]))
                    interpolated_directions_deg.append(
                        np.interp(interp_alt, [h[i], h[i + 1]], [angle1, angle2]) % 360
                    )

        forecast_type = self.forecast.source
        interp_label = "Spline" if interpolation == 'spline' else "Linear"

        fig = plt.figure(figsize=(10, 8))
        ax1 = fig.add_subplot(111, projection='polar')
        sc = ax1.scatter(np.radians(interpolated_directions_deg), interpolated_altitudes, c=interpolated_speeds, cmap='winter', s=2)
        ax1.title.set_text(f"{forecast_type} 3D Windrose for ({self.LAT}, {self.LON}) on {self.new_timestamp}")
        plt.colorbar(sc, label='Wind Speed (m/s)')
        fig.suptitle(f"Wind Interpolation using Speed/Direction {interp_label} Interpolation")

    def _profile_at(self, hour_index, lat, lon):
        """Return (u, v, h) pressure-level wind profile at (time, lat, lon),
        with masks/NaNs removed and geopotential converted to meters.
        """
        g = 9.80665
        lat_i = self.closest(self.lat, float(lat))
        lon_i = self.closest(self.lon, float(lon))

        u_raw = self.ugrdprs[hour_index, :, lat_i, lon_i]
        v_raw = self.vgrdprs[hour_index, :, lat_i, lon_i]
        h_raw = self.hgtprs[hour_index, :, lat_i, lon_i]
        u = u_raw.filled(np.nan) if np.ma.isMaskedArray(u_raw) else np.asarray(u_raw, dtype=float)
        v = v_raw.filled(np.nan) if np.ma.isMaskedArray(v_raw) else np.asarray(v_raw, dtype=float)
        h = np.asarray(h_raw, dtype=float)
        nans = ~np.isnan(u)
        u, v, h = u[nans], v[nans], h[nans]
        h = h / g
        return u, v, h

    @staticmethod
    def _interp_linear_neighbors(h_query, h, u, v):
        """Reproduces GFS.interpolateBearing semantics across a query grid:
        for each h_query, find enclosing pressure levels, linearly interp
        bearing+speed (with 0/360° wrap correction), convert back to u,v.
        """
        u_out = np.empty_like(h_query, dtype=float)
        v_out = np.empty_like(h_query, dtype=float)
        bearing = np.degrees(np.arctan2(v, u)) % 360.0
        speed   = np.hypot(u, v)
        for q_idx, hq in enumerate(h_query):
            # find enclosing indices
            if hq <= h[0]:
                i0, i1 = 0, 1
            elif hq >= h[-1]:
                i0, i1 = len(h) - 2, len(h) - 1
            else:
                i1 = int(np.searchsorted(h, hq))
                i0 = i1 - 1
            a1, a2 = bearing[i0], bearing[i1]
            if abs(a2 - a1) > 180:
                if a2 > a1:
                    a1 += 360
                else:
                    a2 += 360
            sp = np.interp(hq, [h[i0], h[i1]], [speed[i0], speed[i1]])
            dr = np.interp(hq, [h[i0], h[i1]], [a1, a2]) % 360
            u_out[q_idx] = sp * np.cos(np.radians(dr))
            v_out[q_idx] = sp * np.sin(np.radians(dr))
        return u_out, v_out

    def _interp_method(self, method, h_query, h, u, v):
        """Compute (u, v) on h_query using the named production interpolation.

        :param method: one of ``'linear_neighbors'``, ``'linear_full'``, ``'spline_full'``
        """
        if method == 'linear_neighbors':
            return self._interp_linear_neighbors(h_query, h, u, v)
        if method == 'linear_full':
            return np.interp(h_query, h, u), np.interp(h_query, h, v)
        if method == 'spline_full':
            u_out = np.empty_like(h_query, dtype=float)
            v_out = np.empty_like(h_query, dtype=float)
            for i, hq in enumerate(h_query):
                u_out[i], v_out[i] = _spline_uv(hq, h, u, v)
            return u_out, v_out
        raise ValueError(
            f"Unknown wind_interpolation method {method!r}. "
            "Expected 'linear_neighbors', 'linear_full', or 'spline_full'."
        )

    def plotWindMethod(self, hour_index, lat, lon, method='linear_full',
                       alt_step=100.0, ax=None, vmin=None, vmax=None,
                       alt_max=None, show_title=True):
        """3D polar windrose for a single production interpolation method.

        Same visual conventions as :func:`plotWind2` (polar projection,
        radius = altitude, angle = bearing, color = wind speed), but the
        wind profile is interpolated using the named production method
        from ``config_earth.forecast['wind_interpolation']`` instead of
        the ad-hoc visualization-only schemes in ``plotWind2``. This
        makes the panel a faithful picture of what the simulator
        actually sees with that method.

        :param method: ``'linear_neighbors'``, ``'linear_full'``, or ``'spline_full'``
        :type method: str
        :param alt_step: Altitude sample spacing [m] for the dense query grid
        :type alt_step: float
        :param ax: Optional pre-existing polar matplotlib Axes
        :param vmin: Color-scale lower bound for wind speed [m/s] (shared bounds across panels)
        :param vmax: Color-scale upper bound for wind speed [m/s]
        :param alt_max: Outer radial limit [m] (shared across panels)
        :returns: (Figure, PathCollection) — the scatter handle lets the
                  caller attach a colorbar
        """
        u_lvl, v_lvl, h_lvl = self._profile_at(hour_index, lat, lon)
        if len(h_lvl) < 2:
            raise RuntimeError(f"Not enough valid pressure levels at ({lat}, {lon}).")
        h_dense = np.arange(h_lvl[0], h_lvl[-1] + alt_step / 2.0, alt_step)

        u, v = self._interp_method(method, h_dense, h_lvl, u_lvl, v_lvl)
        speed = np.hypot(u, v)
        bearing_deg = np.degrees(np.arctan2(v, u)) % 360.0

        if ax is None:
            fig = plt.figure(figsize=(7, 7))
            ax = fig.add_subplot(111, projection='polar')
        else:
            fig = ax.figure

        sc = ax.scatter(np.radians(bearing_deg), h_dense, c=speed,
                        cmap='winter', s=4, vmin=vmin, vmax=vmax)
        # Overlay the raw pressure-level samples (larger black-edged dots).
        lvl_speed   = np.hypot(u_lvl, v_lvl)
        lvl_bearing = np.degrees(np.arctan2(v_lvl, u_lvl)) % 360.0
        ax.scatter(np.radians(lvl_bearing), h_lvl, c=lvl_speed, cmap='winter',
                   s=45, edgecolor='black', linewidth=0.6,
                   vmin=vmin, vmax=vmax, zorder=10)
        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels(['E', '', 'N', '', 'W', '', 'S', ''])
        if alt_max is not None:
            ax.set_ylim(0, alt_max)
        if show_title:
            ax.set_title(method, fontsize=12, pad=14)
        return fig, sc

    def plotWindMethodsComparison(self, hour_index, lat, lon, alt_step=100.0):
        """Three side-by-side polar windrose panels — one per production
        interpolation method — with shared color and radial scales.
        """
        u_lvl, v_lvl, h_lvl = self._profile_at(hour_index, lat, lon)
        if len(h_lvl) < 2:
            raise RuntimeError(f"Not enough valid pressure levels at ({lat}, {lon}).")
        h_dense = np.arange(h_lvl[0], h_lvl[-1] + alt_step / 2.0, alt_step)

        # Pre-compute speeds across all methods so panels share color bounds.
        all_speeds = [np.hypot(u_lvl, v_lvl)]
        for m in ('linear_neighbors', 'linear_full', 'spline_full'):
            uu, vv = self._interp_method(m, h_dense, h_lvl, u_lvl, v_lvl)
            all_speeds.append(np.hypot(uu, vv))
        vmin = 0.0
        vmax = float(np.max(np.concatenate(all_speeds)))
        alt_max = float(h_dense[-1])

        forecast_type = self.forecast.source
        fig = plt.figure(figsize=(18, 7))
        axes = [fig.add_subplot(1, 3, i + 1, projection='polar') for i in range(3)]
        sc_last = None
        for ax, method in zip(axes,
                              ('linear_neighbors', 'linear_full', 'spline_full')):
            _, sc_last = self.plotWindMethod(
                hour_index, lat, lon, method=method, alt_step=alt_step,
                ax=ax, vmin=vmin, vmax=vmax, alt_max=alt_max,
            )
        fig.suptitle(
            f"{forecast_type} wind interpolation comparison — "
            f"({lat:.3f}, {lon:.3f}) @ {self.new_timestamp}\n"
            "radius = altitude (m), angle = bearing, color = speed (m/s)",
            fontsize=12,
        )
        # Single shared colorbar on the right.
        fig.colorbar(sc_last, ax=axes, orientation='vertical',
                     shrink=0.85, pad=0.05, label='Wind Speed (m/s)')
        return fig
