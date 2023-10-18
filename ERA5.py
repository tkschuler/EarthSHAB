"""
This version was edited to integrate with EarthSHAB.

Contributing authors: Craig Motell and Michael Rodriguez of NIWC Pacific
Edited and integrated: Tristan Schuler
"""
import logging
import numpy as np
import netCDF4
from termcolor import colored
import math
from geographiclib.geodesic import Geodesic
import sys
import os
from scipy import interpolate
from pytz import timezone
from datetime import datetime, timedelta

import config_earth #integrate with EARTHSHAB

class ERA5:
    #def __init__(self, start_coord, end_coord, input_file, dt_sec, sim_time_hours, use_time):
    def __init__(self, start_coord): #CHANGED INPUT VARIABLES
        """Create a class object containing information about an ERA5. Along the way,
        it checks whether it can take a subset of your model data. For example, if you have
        data before or after your simulation times, you can ignore these data through saving
        indexes.

        :param start_coord: starting coordinate of balloon for simulation
        :type coord: dict

        .. note:: Similar to class GFS which stores NOAA GFS data from NOMAD server

        .. seealso:: GFS

        """

        self.dt = config_earth.simulation['dt']
        self.sim_time = config_earth.simulation['sim_time']

        try:
            self.file = netCDF4.Dataset("forecasts/" + config_earth.netcdf_era5["filename"])  # Only accepting manual uploads for now

        except:
            print(colored((f'Unable to locate netcdf file'),"red"))
            sys.exit(1)

        self.geod = Geodesic.WGS84
        self.resolution_hr = config_earth.netcdf_era5['resolution_hr']
        self.start_coord = config_earth.simulation["start_coord"]

        self.init_with_time()

    #def init_without_time(self, start_coord, end_coord, dt_sec, sim_time_hours): # modded the number of varibales for troubleshootin
    def init_with_time(self):

        ''' adds to class object containing information about an ERA5 but without any time
        variables.

        Returns
        -------
        ERA5 : ERA5 Class Object
            Initialize class object with our NetCDF data.

        Notes
        -----
        Similar to class GFS which stores NOAA GFS data from NOMAD server

        See Also
        --------

        '''

        self.start_time = config_earth.simulation['start_time']
        self.min_alt_m = self.start_coord['alt']

        time_arr = self.file.variables['time']
        #Convert from epoch to human readable time. Different than GFS for now.
        self.time_convert = netCDF4.num2date(time_arr[:], time_arr.units, time_arr.calendar)
        self.model_start_datetime = self.time_convert[0] #reanalysis should be same time as gfs predictions for now
        self.model_end_datetime = self.time_convert[-1]

        # Determine Index values from netcdf4 subset
        netcdf_ranges = self.file.variables['u'][0, 0, :, :]
        self.determineRanges(netcdf_ranges)

        # smaller array of downloaded forecast subset
        self.lat = self.file.variables['latitude'][self.lat_max_idx:self.lat_min_idx]
        self.lon = self.file.variables['longitude'][self.lon_min_idx:self.lon_max_idx]

        # min/max lat/lon degree values from netcdf4 subset
        self.LAT_LOW  = self.file.variables['latitude'][self.lat_min_idx-1]
        self.LON_LOW  = self.file.variables['longitude'][self.lon_min_idx]
        self.LAT_HIGH = self.file.variables['latitude'][self.lat_max_idx]
        self.LON_HIGH = self.file.variables['longitude'][self.lon_max_idx-1]

        print("LAT RANGE: min:" + str(self.file.variables['latitude'][self.lat_min_idx-1]), " max: " + str(self.file.variables['latitude'][self.lat_max_idx]) + " size: " + str(self.lat_min_idx-self.lat_max_idx))
        print("LON RANGE: min:" + str(self.file.variables['longitude'][self.lon_min_idx]), " max: " + str(self.file.variables['longitude'][self.lon_max_idx-1]) + " size: " + str(self.lon_max_idx-self.lon_min_idx))

        # Import the netcdf4 subset to speed up table lookup in this script
        g = 9.80665 # gravitation constant used to convert geopotential height to height

        self.levels = self.file.variables['level'][:]
        #these are typos, fix later.  Should be ugrdprs0
        self.ugdrps0 = self.file.variables['u'][self.start_time_idx:self.end_time_idx+1, :, self.lat_max_idx:self.lat_min_idx, self.lon_min_idx:self.lon_max_idx]
        self.vgdrps0 = self.file.variables['v'][self.start_time_idx:self.end_time_idx+1, :, self.lat_max_idx:self.lat_min_idx, self.lon_min_idx:self.lon_max_idx]
        self.hgtprs = self.file.variables['z'][self.start_time_idx:self.end_time_idx+1, :, self.lat_max_idx:self.lat_min_idx, self.lon_min_idx:self.lon_max_idx] / g #what is this divide by g???

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

    def determineRanges(self, netcdf_ranges):
        """
        Determine the dimensions of actual data. If you have columns or rows with missing data
        as indicated by NaN, then this function will return your actual shape size so you can
        resize your data being used.

        """

        results = np.all(~netcdf_ranges.mask)
        print(results)
        if results == False:
            timerange, latrange, lonrange = np.nonzero(~netcdf_ranges.mask)

            self.start_time_idx = timerange.min()
            self.end_time_idx = timerange.max()
            self.lat_min_idx = latrange.min() #Min/Max are switched compared to with ERA5
            self.lat_max_idx = latrange.max()
            self.lon_min_idx = lonrange.min()
            self.lon_max_idx = lonrange.max()
        else: #This might be broken for time
            self.start_time_idx = 0
            self.end_time_idx = len(self.time_convert)-1
            lati, loni = netcdf_ranges.shape
            #THese are backwards???
            self.lat_min_idx = lati
            self.lat_max_idx = 0
            self.lon_max_idx = loni
            self.lon_min_idx = 0

    def closestIdx(self, arr, k):
        """ Given an ordered array and a value, determines the index of the closest item contained in the array.

        """
        return min(range(len(arr)), key=lambda i: abs(arr[i] - k))

    def getNearestLatIdx(self, lat, min, max):
        """ Determines the nearest latitude index (to .25 degrees),  which will be integer index starting at 0 where you data is starting

        """
        # arr = np.arange(start=min, stop=max, step=self.res)
        # i = self.closest(arr, lat)
        i = self.closestIdx(self.lat, lat)
        return i

    def getNearestLonIdx(self, lon, min, max):
        """ Determines the nearest longitude (to .25 degrees)

        """

        # lon = lon % 360 #convert from -180-180 to 0-360
        # arr = np.arange(start=min, stop=max, step=self.res)
        i = self.closestIdx(self.lon, lon)
        return i

    def getNearestAltbyIndex(self, int_hr_idx, lat_i, lon_i, alt_m):
        """ Determines the nearest altitude based off of geo potential height of a .25 degree lat/lon area.
            at a given latitude, longitude index

        """
        i = self.closestIdx(self.hgtprs[int_hr_idx, :, lat_i, lon_i], alt_m)

        return i

    def wind_alt_Interpolate(self, alt_m, diff_time, lat_idx, lon_idx):
        """
        Performs a 2 step linear interpolation to determine horizontal wind velocity.  First the altitude is interpolated between the two nearest
        .25 degree lat/lon areas.  Once altitude is matched, a second linear interpolation is performed with respect to time.

        .. note:: I performed a lot of clean up of this code from what Craig originally sent me.  I'm not exactly sure what
            about all the difference from GFS to ERA5 for this function.  I will do more cleaning and variable renaming in a later version.
            The results are the same right now thoug.

        Parameters
        ----------
        alt_m : float
            Contains current altitude of balloon object.
        diff_time: ?
            This parameter is calculated in two different places, need to clean up this parameter
        lon_idx : integer
            Closest index into longitude array
        lat_idx : integer
            Closest index into latitude array

        Returns:
            u : float
                u component wind vector
            v : float
                v component wind vector

        """

        #we need to clean this up, ends up we check for the hour_index before we call this function
        #diff = coord["timestamp"] - self.model_start_datetime
        hour_index = (diff_time.days * 24 + diff_time.seconds / 3600.) / self.resolution_hr

        hgt_m = alt_m

        # First interpolate wind speeds between 2 closest time steps to match altitude estimates (hgtprs), which can change with time
        v_0 = self.vgdrps0[int(hour_index), :, lat_idx, lon_idx]  # Round hour index to nearest int
        v_0 = self.fill_missing_data(v_0)  # Fill the missing wind data if occurs at upper heights

        u_0 = self.ugdrps0[int(hour_index), :, lat_idx, lon_idx]
        u_0 = self.fill_missing_data(u_0)

        #t_0 = self.temp0[int_hr_idx, :, lat_idx, lon_idx]
        #t_0 = self.fill_missing_data(t_0)

        xp_0 = self.hgtprs[int(hour_index), :, lat_idx, lon_idx]
        xp_0 = self.fill_missing_data(xp_0)


        #f = interpolate.interp1d(xp_0, t_0, assume_sorted=False, fill_value="extrapolate")
        #t0_pt = f(hgt_m)
        f = interpolate.interp1d(xp_0, u_0, assume_sorted=False, fill_value="extrapolate")
        u0_pt = f(hgt_m)
        f = interpolate.interp1d(xp_0, v_0, assume_sorted=False, fill_value="extrapolate")
        v0_pt = f(hgt_m)

        # Next interpolate the wind velocities with respect to time.
        v_1 = self.vgdrps0[int(hour_index) + 1, :, lat_idx, lon_idx]
        v_1 = self.fill_missing_data(v_1)
        u_1 = self.ugdrps0[int(hour_index) + 1, :, lat_idx, lon_idx]
        u_1 = self.fill_missing_data(u_1)
        #t_1 = self.temp0[int_hr_idx2, :, lat_idx, lon_idx]
        #t_1 = self.fill_missing_data(t_1)
        xp_1 = self.hgtprs[int(hour_index) +1 , :, lat_idx, lon_idx]
        xp_1 = self.fill_missing_data(xp_1)

        #f = interpolate.interp1d(xp_1, t_1, assume_sorted=False, fill_value="extrapolate" )

        #t1_pt = f(hgt_m)
        f = interpolate.interp1d(xp_1, u_1, assume_sorted=False, fill_value="extrapolate")
        u1_pt = f(hgt_m)
        f = interpolate.interp1d(xp_1, v_1, assume_sorted=False, fill_value="extrapolate")
        v1_pt = f(hgt_m)

        fp = [int(hour_index), int(hour_index) + 1]
        u = np.interp(hour_index, fp, [u0_pt, u1_pt])
        v = np.interp(hour_index, fp, [v0_pt, v1_pt])
        #t = np.interp(hour_index, fp, [t0_pt, t1_pt])

        return [u, v]

    def fill_missing_data(self, data):
        """Helper function to fill in linearly interpolate and fill in missing data

        """

        data = data.filled(np.nan)
        nans, x = np.isnan(data), lambda z: z.nonzero()[0]
        data[nans] = np.interp(x(nans), x(~nans), data[~nans])
        return data

    def getNewCoord(self, coord, dt):
        """
        Determines the new coordinates of the balloon based on the effects of U and V wind components.

        Parameters
        ----------
        coord : dict
            Contains current position, altitude and time of balloon object.
        dt : float
            Integration time

        Returns:
            lat_new, lon_new, x_wind_vel, y_wind_vel, bearing, closest_lat, closest_lon, closest alt

        """

        diff_time = coord["timestamp"] - self.model_start_datetime
        hour_index = (diff_time.days * 24 + diff_time.seconds / 3600.) / self.resolution_hr
        int_hr_idx = int(hour_index)

        lat_idx = self.getNearestLatIdx(coord["lat"], self.lat_max_idx, self.lat_min_idx)
        lon_idx = self.getNearestLonIdx(coord["lon"], self.lon_min_idx, self.lon_max_idx)
        #z = self.getNearestAltbyIndex(lat_idx, lon_idx, coord["alt_m"])  # fix this for lower and Upper
        z = self.getNearestAltbyIndex(int_hr_idx, lat_idx, lon_idx, coord["alt"])  # fix this for lower and Upper

        x_wind_vel, y_wind_vel = self.wind_alt_Interpolate(coord['alt'], diff_time, lat_idx, lon_idx)

        try:
            x_wind_vel, y_wind_vel = self.wind_alt_Interpolate(coord['alt'], diff_time, lat_idx, lon_idx)  # for now hour index is 0
            if x_wind_vel == None:
                return [None, None, None, None, None, None, None, None]
                # x_wind_vel, y_wind_vel = self.wind_interpolate_by_index(coord, i, j, z)  # for now hour index is 0
        except:
            print(colored(
                "Mismatch with simulation and forecast timstamps. Check simulation start time and/or download a new forecast.",
                "red"))
            sys.exit(1)

        bearing = math.degrees(math.atan2(y_wind_vel, x_wind_vel))
        bearing = 90 - bearing  # perform 90 degree rotation for bearing from wind data
        d = math.pow((math.pow(y_wind_vel, 2) + math.pow(x_wind_vel, 2)), .5) * dt  # dt multiplier
        g = self.geod.Direct(coord["lat"], coord["lon"], bearing, d)

        #GFS is %360 here.  The files are organized a bit differently
        if g['lat2'] < self.LAT_LOW or g['lat2']  > self.LAT_HIGH or (g['lon2']) < self.LON_LOW or (g['lon2']) > self.LON_HIGH:
            print(colored("WARNING: Trajectory is out of bounds of downloaded netcdf forecast", "yellow"))

        '''Changed to alt instead of alt_m'''
        if coord["alt"] <= self.min_alt_m:
            # Balloon should remain stationary if it's reached the minimum altitude
            return [coord['lat'], coord['lon'], x_wind_vel, y_wind_vel, bearing, self.lat[lat_idx], self.lon[lon_idx],
                    self.hgtprs[0, z, lat_idx, lon_idx]]  # hgtprs doesn't matter here so is set to 0
        else:
            return [g['lat2'], g['lon2'], x_wind_vel, y_wind_vel, bearing, self.lat[lat_idx], self.lon[lon_idx],
                    self.hgtprs[0, z, lat_idx, lon_idx]]

    #This function is unused now after code cleanup, but leaving for now because Craig included it with the orginal code
    '''CHECK THIS'''
    def datetime_epoch(self, timedate_obj):
        gmt = timezone('GMT')
        try:
            dt = timedate_obj.strftime("%Y-%m-%d %H:%M:%S")
            naive_ts = datetime.strptime(dt, "%Y-%m-%d %H:%M:%S")
        except:
            naive_ts = datetime.strptime(timedate_obj, "%Y-%m-%d %H:%M:%S")

        local_ts = gmt.localize(naive_ts)
        epoch_ts = local_ts.timestamp()

        return int(epoch_ts)
