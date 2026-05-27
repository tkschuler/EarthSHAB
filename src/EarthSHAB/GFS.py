""" GFS extracts meteorological data from a NOAA netcdf (.nc file) data set.  To speed up predicting trajectories,
saveNETCDF.py should be run before main.py. This way, the large data set doesn't need to be redownloaded each time the trajectory is run.
For now, only wind velocity is used for the simulation prediction.  Atmposheric properties such as temperarature and pressure are based off of
the U.S. Standard Atmosphere tables from 1976, and the fluids library is used for these, which can be seen in radiation3.py

"""

import numpy as np
import netCDF4
from scipy.interpolate import CubicSpline
from termcolor import colored
import math
from geographiclib.geodesic import Geodesic
from datetime import timedelta
import sys
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()   #Hacky solution for Python 3.6 to use ISO format Strings

import EarthSHAB.config_earth as config_earth


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
    # Keep only strictly increasing samples.
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


class GFS:
    def __init__(self, centered_coord):
        # import variables from configuration file
        self.centered_coord = centered_coord
        self.min_alt = config_earth.simulation['min_alt']
        self.start_time = config_earth.simulation['start_time']
        self.sim_time = config_earth.simulation['sim_time']
        #self.hours3 = config_earth.netcdf_gfs['hours3']

        self.file = netCDF4.Dataset(config_earth.netcdf_gfs["nc_file"])  # Only accepting manual uploads for now
        self.gfs_time = config_earth.netcdf_gfs['nc_start']
        self.res = config_earth.netcdf_gfs['res']

        self.geod = Geodesic.WGS84

        # Phase 3: Read canonical v2 format
        time_arr = self.file.variables['valid_time']
        time_units = time_arr.units if hasattr(time_arr, 'units') else "seconds since 1970-01-01"
        time_calendar = time_arr.calendar if hasattr(time_arr, 'calendar') else "proleptic_gregorian"

        self.time_convert = netCDF4.num2date(time_arr[:], units=time_units, calendar=time_calendar)

        # Determine Index values from netcdf4 subset (using canonical variable name 'u')
        netcdf_ranges = self.file.variables['u'][:,0,:,:]
        self.determineRanges(netcdf_ranges)

        self.lat = np.asarray(
            self.file.variables['latitude'][self.lat_min_idx:self.lat_max_idx + 1],
            dtype=float
        )
        self.lon = np.asarray(
            self.file.variables['longitude'][self.lon_min_idx:self.lon_max_idx + 1],
            dtype=float
        )

        self.ugdrps0 = self.file.variables['u'][
            self.start_time_idx:self.end_time_idx + 1,
            :,
            self.lat_min_idx:self.lat_max_idx + 1,
            self.lon_min_idx:self.lon_max_idx + 1
        ]
        self.vgdrps0 = self.file.variables['v'][
            self.start_time_idx:self.end_time_idx + 1,
            :,
            self.lat_min_idx:self.lat_max_idx + 1,
            self.lon_min_idx:self.lon_max_idx + 1
        ]

        # Phase 3: Read geopotential (z) and convert to height in meters by dividing by g
        g = 9.80665
        self.hgtprs = self.file.variables['z'][
            self.start_time_idx:self.end_time_idx + 1,
            :,
            self.lat_min_idx:self.lat_max_idx + 1,
            self.lon_min_idx:self.lon_max_idx + 1
        ] / g

        # min/max lat/lon degree values from netcdf4 subset
        # Phase 3: canonical v2 format has descending latitude, so min/max are swapped in the array
        lat_values = self.file.variables['latitude'][self.lat_min_idx:self.lat_max_idx + 1]
        lon_values = self.file.variables['longitude'][self.lon_min_idx:self.lon_max_idx + 1]

        self.LAT_LOW  = float(np.min(lat_values))
        self.LAT_HIGH = float(np.max(lat_values))
        self.LON_LOW  = float(np.min(lon_values))
        self.LON_HIGH = float(np.max(lon_values))

        print("LAT RANGE: min: " + str(self.LAT_LOW), " (deg) max: " + str(self.LAT_HIGH) + " (deg) array size: " + str(self.lat_max_idx-self.lat_min_idx+1))
        print("LON RANGE: min: " + str(self.LON_LOW), " (deg) max: " + str(self.LON_HIGH) + " (deg) array size: " + str(self.lon_max_idx-self.lon_min_idx+1))

        # Import the netcdf4 subset to speed up table lookup in this script
        self.levels = self.file.variables['pressure_level'][:]

        #print("Data downloaded.\n\n")
        print()

        #Check if number of hours will fit in simulation time
        desired_simulation_end_time = self.start_time + timedelta(hours=self.sim_time)
        diff_time = (self.time_convert[self.end_time_idx] - self.start_time).total_seconds() #total number of seconds between 2 timestamps

        print("Sim start time: ", self.start_time)
        print("NetCDF end time:", self.time_convert[self.end_time_idx])
        print("Max sim runtime:", diff_time//3600, "hours")
        print("Des sim runtime:", self.sim_time, "hours")
        print()

        if not desired_simulation_end_time <= self.time_convert[self.end_time_idx]:
            print(colored("Desired simulation run time of " + str(self.sim_time)  +
            " hours is out of bounds of downloaded forecast. " +
            "Check simulation start time and/or download a new forecast.", "red"))
            sys.exit()


    def determineRanges(self,netcdf_ranges):
        """
        Determine the following variable ranges (min/max) within the netcdf file:

        -time
        -lat
        -lon

        .. note:: the levels variable is uncessary for knowing the ranges for, because they vary from
            coordinate to coordinate, and all levels in the array are always be used.

        """


        print(colored("Forecast Information (Parsed from netcdf file):", "blue", attrs=['bold']))

        results = np.all(~netcdf_ranges.mask)
        #Results will almost always be false,  unless an entire netcdf of the world is downloaded. Or if the netcdf is downloaded via another method with lat/lon bounds
        if results == False:
            timerange, latrange, lonrange = np.nonzero(~netcdf_ranges.mask)

            self.start_time_idx = timerange.min()
            self.end_time_idx = timerange.max()
            self.lat_min_idx = latrange.min() #Min/Max are switched compared to with ERA5
            self.lat_max_idx = latrange.max()
            self.lon_min_idx = lonrange.min()
            self.lon_max_idx = lonrange.max()
        else: #This might be broken for time
            #Need to double check this for an entire world netcdf download
            t_size, lat_size, lon_size = netcdf_ranges.shape
            self.start_time_idx = 0
            self.end_time_idx   = t_size   - 1
            self.lat_min_idx    = 0
            self.lat_max_idx    = lat_size - 1
            self.lon_min_idx    = 0
            self.lon_max_idx    = lon_size - 1

    def closest(self, arr, k):
        """ Given an ordered array and a value, determines the index of the closest item contained in the array.

        """
        return min(range(len(arr)), key = lambda i: abs(arr[i]-k))

    def getNearestLat(self, lat, min_unused=None, max_unused=None):
        """
        Determine the nearest latitude index using the actual stored latitude array.
        """
        return self.closest(self.lat, float(lat))

    def getNearestLon(self, lon, min_unused=None, max_unused=None):
        """
        Determine the nearest longitude index using the actual stored longitude array.
        Phase 3: canonical v2 longitudes are in [-180, 180], so normalize the query the same way.
        """
        lon_normalized = float(lon)
        # Normalize to [-180, 180)
        while lon_normalized >= 180:
            lon_normalized -= 360
        while lon_normalized < -180:
            lon_normalized += 360
        return self.closest(self.lon, lon_normalized)

    def getNearestAlt(self,hour_index,lat,lon,alt):
        """ Determines the nearest altitude based off of geo potential height of a .25 degree lat/lon area.

        """

        lat_i = self.getNearestLat(lat,self.LAT_LOW,self.LAT_HIGH)
        lon_i = self.getNearestLon(lon,self.LON_LOW,self.LON_HIGH)
        i = self.closest(self.hgtprs[int(hour_index),:,lat_i,lon_i], alt)
        return i

    def get2NearestAltIdxs(self, h, alt_m):
        """ Determines 2 nearest indexes for altitude for interpolating angles.
        It does index wrap from 0 to -1, which is taken care of in `ERA5.interpolateBearing()`

        """
        h_nearest = self.closest(h, alt_m)

        if alt_m > h[h_nearest]:
            h_idx0 = h_nearest
            h_idx1 = h_nearest +1
        else:
            h_idx0 = h_nearest - 1
            h_idx1 = h_nearest

        return h_idx0, h_idx1

    def wind_alt_Interpolate(self, coord):
        """
        This function performs a 2-step linear interpolation to determine horizontal wind velocity at a
        3d desired coordinate and timestamp.

        The figure below shows a visual representation of how wind data is stored in netcdf forecasts based on
        lat, lon, and geopotential height. The data forms a non-uniform grid, that also changes in time. Therefore
        we performs a 2-step linear interpolation to determine horizontal wind velocity at a desired 3D coordinate
        and particular timestamp.

        To start, the two nearest .25 degree lat/lon areas to the desired coordinate are looked up along with the
        2 closest timestamps t0 and t1. This produces 6 arrays: u-wind, v-wind, and geopotential heights at the lower
        and upper closest timestamps (t0 and t1).

        Next, the geopotential height is converted to altitude (m) for each timestamp. For the first interpolation,
        the u-v wind components at the desired altitude are determined (1a and 1b) using np.interp.

        Then, once the wind speeds at matching altitudes for t0 and t1 are detemined, a second linear interpolation
        is performed with respect to time (t0 and t1).

        .. image:: ../../img/netcdf-2step-interpolation.png

        :param coord: Coordinate of balloon
        :type coord: dict
        :returns: [u_wind_vel, v_wind_vel]
        :rtype: array

        """

        diff = coord["timestamp"] - self.gfs_time
        hour_index = (diff.days*24 + diff.seconds / 3600.)/3

        lat_i = self.getNearestLat(coord["lat"],self.LAT_LOW,self.LAT_HIGH)
        lon_i = self.getNearestLon(coord["lon"],self.LON_LOW,self.LON_HIGH)
        z_low = self.getNearestAlt(hour_index,coord["lat"],coord["lon"],coord["alt"]) #fix this for lower and Upper


        # First interpolate wind speeds between 2 closest time steps to match altitude estimates (hgtprs), which can change with time
        v_0 = self.vgdrps0[int(hour_index),:,lat_i,lon_i] # Round hour index to nearest int
        v_0 = self.fill_missing_data(v_0) # Fill the missing wind data. With netcdf4 there are always 3 missing values at the higher elevations
        u_0 = self.ugdrps0[int(hour_index),:,lat_i,lon_i]
        u_0 = self.fill_missing_data(u_0)

        u0 = np.interp(coord["alt"],self.hgtprs[int(hour_index),:,lat_i,lon_i],u_0)
        v0 = np.interp(coord["alt"],self.hgtprs[int(hour_index),:,lat_i,lon_i],v_0)

        # Next interpolate the wind velocities with respect to time.
        v_1 = self.vgdrps0[int(hour_index)+1,:,lat_i,lon_i]
        v_1 = self.fill_missing_data(v_1)

        u_1 = self.ugdrps0[int(hour_index)+1,:,lat_i,lon_i]
        u_1 = self.fill_missing_data(u_1)

        u1 = np.interp(coord["alt"],self.hgtprs[int(hour_index)+1,:,lat_i,lon_i],u_1) # Round hour index to next timestep
        v1 = np.interp(coord["alt"],self.hgtprs[int(hour_index)+1,:,lat_i,lon_i],v_1)

        u = np.interp(hour_index,[int(hour_index),int(hour_index)+1],[u0,u1])
        v = np.interp(hour_index,[int(hour_index),int(hour_index)+1],[v0,v1])

        return[u,v]

    def wind_alt_Interpolate2(self, coord):

        diff = coord["timestamp"] - self.gfs_time
        hour_index = (diff.days*24 + diff.seconds / 3600.)/3

        lat_i = self.getNearestLat(coord["lat"],self.LAT_LOW,self.LAT_HIGH)
        lon_i = self.getNearestLon(coord["lon"],self.LON_LOW,self.LON_HIGH)


        z_low = self.getNearestAlt(hour_index,coord["lat"],coord["lon"],coord["alt"]) #fix this for lower and Upper


        # First interpolate wind speeds between 2 closest time steps to match altitude estimates (hgtprs), which can change with time
        #t0
        v_0 = self.vgdrps0[int(hour_index),:,lat_i,lon_i] # Round hour index to nearest int
        v_0 = self.fill_missing_data(v_0) # Fill the missing wind data. With netcdf4 there are always 3 missing values at the higher elevations

        u_0 = self.ugdrps0[int(hour_index),:,lat_i,lon_i]
        u_0 = self.fill_missing_data(u_0)

        h0 = self.hgtprs[int(hour_index),:,lat_i,lon_i]
        h0 = self.fill_missing_data(h0)

        #t1

        v_1 = self.vgdrps0[int(hour_index)+1,:,lat_i,lon_i]
        v_1 = self.fill_missing_data(v_1)

        u_1 = self.ugdrps0[int(hour_index)+1,:,lat_i,lon_i]
        u_1 = self.fill_missing_data(u_1)

        h1 = self.hgtprs[int(hour_index)+1,:,lat_i,lon_i]
        h1 = self.fill_missing_data(h1)

        # linear_full path (np.interp on u,v across full altitude profile,
        # then linear in time). Always computed for diagnostic comparison
        # (returned as u_diag, v_diag below).
        u0_lf = np.interp(coord["alt"], h0, u_0)
        v0_lf = np.interp(coord["alt"], h0, v_0)
        u1_lf = np.interp(coord["alt"], h1, u_1)
        v1_lf = np.interp(coord["alt"], h1, v_1)
        u_lf  = np.interp(hour_index, [int(hour_index), int(hour_index)+1], [u0_lf, u1_lf])
        v_lf  = np.interp(hour_index, [int(hour_index), int(hour_index)+1], [v0_lf, v1_lf])

        method = config_earth.forecast.get('wind_interpolation', 'linear_neighbors')

        if method == 'linear_neighbors':
            bearing_t0, speed_t0 = self.interpolateBearing(h0, u_0, v_0, coord["alt"])
            bearing_t1, speed_t1 = self.interpolateBearing(h1, u_1, v_1, coord["alt"])
            bearing_interpolated, speed_interpolated = self.interpolateBearingTime(
                bearing_t0, speed_t0, bearing_t1, speed_t1, hour_index)
            u = speed_interpolated * np.cos(np.radians(bearing_interpolated))
            v = speed_interpolated * np.sin(np.radians(bearing_interpolated))
        elif method == 'linear_full':
            u, v = u_lf, v_lf
        elif method == 'spline_full':
            u0_sp, v0_sp = _spline_uv(coord["alt"], h0, u_0, v_0)
            u1_sp, v1_sp = _spline_uv(coord["alt"], h1, u_1, v_1)
            u = np.interp(hour_index, [int(hour_index), int(hour_index)+1], [u0_sp, u1_sp])
            v = np.interp(hour_index, [int(hour_index), int(hour_index)+1], [v0_sp, v1_sp])
        else:
            raise ValueError(
                f"Unknown wind_interpolation: {method!r}. "
                "Expected 'linear_neighbors', 'linear_full', or 'spline_full'."
            )

        # Return selected (u, v) for sim use, plus linear_full as diagnostic baseline.
        return [u, v, u_lf, v_lf]

    def windVectorToBearing(self, u, v):
        """Helper function to conver u-v wind components to bearing and speed.

        """
        bearing = np.arctan2(v,u)
        speed = np.power((np.power(u,2)+np.power(v,2)),.5)

        return [bearing, speed]

    def interpolateBearing(self, h, u, v, alt_m ): #h, bearing0, speed0, bearing1, speed1):
        """Given altitude, u_velocity and v_velocity arrays as well as a desired altitude, perform a linear interpolation
        between the 2 altitudes (h0 and h1) while accounting for possible 0/360degree axis crossover.

        """

        h_idx0, h_idx1 = self.get2NearestAltIdxs(h, alt_m)

        #Check to make sure altitude isn't outside bounds of altitude array for interpolating
        if h_idx0 == -1:
            h_idx0 = 0

        if h_idx1 == 0: #this one should most likely never trigger because altitude forecasts go so high.
            h_idx1 = -1

        bearing0, speed0 = self.windVectorToBearing(u[h_idx0], v[h_idx0])
        bearing1, speed1 = self.windVectorToBearing(u[h_idx1], v[h_idx1])

        interp_dir_deg = 0
        angle1 = np.degrees(bearing0) %360
        angle2 = np.degrees(bearing1) %360
        angular_difference = abs(angle2-angle1)

        if angular_difference > 180:
            if (angle2 > angle1):
                angle1 += 360
            else:
                angle2 += 360

        interp_speed = np.interp(alt_m, [h[h_idx0], h[h_idx1]], [speed0, speed1])
        interp_dir_deg = np.interp(alt_m, [h[h_idx0], h[h_idx1]], [angle1, angle2]) % 360 #make sure in the range (0, 360)

        return (interp_dir_deg, interp_speed)

    def interpolateBearingTime(self, bearing0, speed0, bearing1, speed1, hour_index ): #h, bearing0, speed0, bearing1, speed1):
        """Similar to `ERA5.interpolateBearing()` however bearings and speeds are already known
        and linearly then interpolated with respect to time (t0 and t1)

        """

        interp_dir_deg = 0
        angle1 = bearing0 %360
        angle2 = bearing1 %360
        angular_difference = abs(angle2-angle1)

        if angular_difference > 180:
            if (angle2 > angle1):
                angle1 += 360
            else:
                angle2 += 360

        fp = [int(hour_index), int(hour_index) + 1]
        interp_speed = np.interp(hour_index, fp, [speed0, speed1])
        interp_dir_deg = np.interp(hour_index, fp, [angle1, angle2]) % 360 #make sure in the range (0, 360)

        return (interp_dir_deg, interp_speed)

    def fill_missing_data(self, data):
        """Helper function to fill in linearly interpolate and fill in missing data

        """

        data = data.filled(np.nan)
        nans, x = np.isnan(data), lambda z: z.nonzero()[0]
        data[nans]= np.interp(x(nans), x(~nans), data[~nans])
        return data

    def getNewCoord(self, coord, dt):
        """ Determines the new coordinates every second due to wind velocity
        :param coord: Coordinate of balloon
        :type coord: dict
        :returns: [lat_new, lon_new, x_wind_vel, y_wind_vel, bearing, closest_lat, closest_lon, closest alt]
        :rtype: array

        """

        diff = coord["timestamp"] - self.gfs_time
        hour_index = (diff.days*24 + diff.seconds / 3600.)/3

        i = self.getNearestLat(coord["lat"],self.LAT_LOW,self.LAT_HIGH)
        j = self.getNearestLon(coord["lon"],self.LON_LOW,self.LON_HIGH)
        z = self.getNearestAlt(int(hour_index),coord["lat"],coord["lon"],coord["alt"])

        #Get wind estimate for current coordiante

        #add some error handling here?
        #old option:
        #x_wind_vel,y_wind_vel = self.wind_alt_Interpolate(coord)
        #x_wind_vel_old, y_wind_vel_old = None, None

        x_wind_vel,y_wind_vel, x_wind_vel_old, y_wind_vel_old = self.wind_alt_Interpolate2(coord)

        bearing = math.degrees(math.atan2(y_wind_vel,x_wind_vel))
        bearing = 90 - bearing # perform 90 degree rotation for bearing from wind data
        d = math.pow((math.pow(y_wind_vel,2)+math.pow(x_wind_vel,2)),.5) * dt #dt multiplier
        g = self.geod.Direct(coord["lat"], coord["lon"], bearing, d)

        # Phase 3: Check bounds with canonical v2 longitude in [-180, 180]
        lon2_normalized = g['lon2']
        while lon2_normalized >= 180:
            lon2_normalized -= 360
        while lon2_normalized < -180:
            lon2_normalized += 360

        if g['lat2'] < self.LAT_LOW or g['lat2'] > self.LAT_HIGH or lon2_normalized < self.LON_LOW or lon2_normalized > self.LON_HIGH:
            print(colored("WARNING: Trajectory is out of bounds of downloaded netcdf forecast", "yellow"))

        if coord["alt"] <= self.min_alt:
            # Balloon should remain stationary if it's reached the minimum altitude
            return [coord['lat'],coord['lon'],x_wind_vel,y_wind_vel, x_wind_vel_old, y_wind_vel_old, bearing, self.lat[i], self.lon[j], self.hgtprs[0,z,i,j]] # hgtprs doesn't matter here so is set to 0
        else:
            return [g['lat2'],g['lon2'],x_wind_vel,y_wind_vel, x_wind_vel_old, y_wind_vel_old, bearing, self.lat[i], self.lon[j], self.hgtprs[0,z,i,j]]
