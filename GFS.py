""" GFS extracts meteorological data from a NOAA netcdf (.nc file) data set.  To speed up predicting trajectories,
saveNETCDF.py should be run before main.py. This way, the large data set doesn't need to be redownloaded each time the trajectory is run.
For now, only wind velocity is used for the simulation prediction.  Atmposheric properties such as temperarature and pressure are based off of
the U.S. Standard Atmosphere tables from 1976, and the fluids library is used for these, which can be seen in radiation3.py

"""

import numpy as np
import netCDF4
from termcolor import colored
import math
from geographiclib.geodesic import Geodesic
import datetime
from datetime import timedelta
import sys
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()   #Hacky solution for Python 3.6 to use ISO format Strings

import config_earth

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

        time_arr = self.file.variables['time']

        #Manually add time units, not imported with units formatted in saveNETCDF.py
        self.time_convert = netCDF4.num2date(time_arr[:], units="days since 0001-01-01", has_year_zero=True)


        #Determine Index values from netcdf4 subset
        netcdf_ranges = self.file.variables['ugrdprs'][:,0,:,:]
        self.determineRanges(netcdf_ranges)

        # smaller array of downloaded forecast subset
        self.lat  = self.file.variables['lat'][self.lat_min_idx:self.lat_max_idx]
        self.lon  = self.file.variables['lon'][self.lon_min_idx:self.lon_max_idx]

        # min/max lat/lon degree values from netcdf4 subset
        self.LAT_LOW  = self.file.variables['lat'][self.lat_min_idx]
        self.LON_LOW  = self.file.variables['lon'][self.lon_min_idx]
        self.LAT_HIGH = self.file.variables['lat'][self.lat_max_idx]
        self.LON_HIGH = self.file.variables['lon'][self.lon_max_idx]

        print("LAT RANGE: min: " + str(self.LAT_LOW), " (deg) max: " + str(self.LAT_HIGH) + " (deg) array size: " + str(self.lat_max_idx-self.lat_min_idx+1))
        print("LON RANGE: min: " + str(self.LON_LOW), " (deg) max: " + str(self.LON_HIGH) + " (deg) array size: " + str(self.lon_max_idx-self.lon_min_idx+1))

        # Import the netcdf4 subset to speed up table lookup in this script
        self.levels = self.file.variables['lev'][:]
        self.ugdrps0 = self.file.variables['ugrdprs'][self.start_time_idx:self.end_time_idx+1,:,self.lat_min_idx:self.lat_max_idx,self.lon_min_idx:self.lon_max_idx]
        self.vgdrps0 = self.file.variables['vgrdprs'][self.start_time_idx:self.end_time_idx+1,:,self.lat_min_idx:self.lat_max_idx,self.lon_min_idx:self.lon_max_idx]
        self.hgtprs  = self.file.variables['hgtprs'][self.start_time_idx:self.end_time_idx+1,:,self.lat_min_idx:self.lat_max_idx,self.lon_min_idx:self.lon_max_idx]

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
            '''
            self.start_time_idx = 0
            self.end_time_idx = len(self.time_convert)-1
            lati, loni = netcdf_ranges.shape
            self.lat_min_idx = lati
            self.lat_max_idx = 0
            self.lon_max_idx = loni
            self.lon_min_idx = 0
            '''

    def closest(self, arr, k):
        """ Given an ordered array and a value, determines the index of the closest item contained in the array.

        """
        return min(range(len(arr)), key = lambda i: abs(arr[i]-k))

    def getNearestLat(self,lat,min,max):
        """ Determines the nearest lattitude (to .25 degrees)

        """
        arr = np.arange(start=min, stop=max, step=self.res)
        i = self.closest(arr, lat)
        return i

    def getNearestLon(self,lon,min,max):
        """ Determines the nearest longitude (to .25 degrees)

        """

        lon = lon % 360 #convert from -180-180 to 0-360
        arr = np.arange(start=min, stop=max, step=self.res)
        i = self.closest(arr, lon)
        return i

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

        u0 = np.interp(coord["alt"],h0,u_0)
        v0 = np.interp(coord["alt"],h0,v_0)

        bearing_t0, speed_t0 = self.interpolateBearing(h0, u_0, v_0, coord["alt"])

        #t1

        v_1 = self.vgdrps0[int(hour_index)+1,:,lat_i,lon_i]
        v_1 = self.fill_missing_data(v_1)

        u_1 = self.ugdrps0[int(hour_index)+1,:,lat_i,lon_i]
        u_1 = self.fill_missing_data(u_1)

        h1 = self.hgtprs[int(hour_index)+1,:,lat_i,lon_i]
        h1 = self.fill_missing_data(h1)

        u1 = np.interp(coord["alt"],h1,u_1) # Round hour index to next timestep
        v1 = np.interp(coord["alt"],h1,v_1)

        bearing_t1, speed_t1 = self.interpolateBearing(h1, u_1, v_1, coord["alt"])


        #Old Method
        u_old = np.interp(hour_index,[int(hour_index),int(hour_index)+1],[u0,u1])
        v_old = np.interp(hour_index,[int(hour_index),int(hour_index)+1],[v0,v1])

        #New Method (Bearing/speed  linear interpolation):
        bearing_interpolated, speed_interpolated =  self.interpolateBearingTime(bearing_t0, speed_t0, bearing_t1, speed_t1, hour_index )

        u_new = speed_interpolated * np.cos(np.radians(bearing_interpolated))
        v_new = speed_interpolated * np.sin(np.radians(bearing_interpolated))

        return [u_new, v_new, u_old, v_old]

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

        if g['lat2'] < self.LAT_LOW or g['lat2']  > self.LAT_HIGH or (g['lon2'] % 360) < self.LON_LOW or (g['lon2'] % 360) > self.LON_HIGH:
            print(colored("WARNING: Trajectory is out of bounds of downloaded netcdf forecast", "yellow"))

        if coord["alt"] <= self.min_alt:
            # Balloon should remain stationary if it's reached the minimum altitude
            return [coord['lat'],coord['lon'],x_wind_vel,y_wind_vel, x_wind_vel_old, y_wind_vel_old, bearing, self.lat[i], self.lon[j], self.hgtprs[0,z,i,j]] # hgtprs doesn't matter here so is set to 0
        else:
            return [g['lat2'],g['lon2'],x_wind_vel,y_wind_vel, x_wind_vel_old, y_wind_vel_old, bearing, self.lat[i], self.lon[j], self.hgtprs[0,z,i,j]]

        if g['lat2'] < self.LAT_LOW or g['lat2'] > self.LAT_HIGH or g['lon2'] < self.LON_LOW or g['lon2'] > self.LON_HIGH:
            print(colored("WARNING: Trajectory is out of bounds of downloaded netcdf forecast", "yellow"))
