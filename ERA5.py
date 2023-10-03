"""
This version was edited to integrate with EarthSHAB.

Version 1.2 (02-25-2022)
Authors: MR and CEM
Edited and Integrated: Tristan Schuler
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
from datetime import datetime

import config_earth #integrate with EARTHSAHB

logging.basicConfig(filename='era5.log',
                    level=logging.WARNING, # Changed from Info for Demo
                    format='%(asctime)s | %(levelname)s | %(message)s',
                    datefmt="%Y-%m-%d %H:%M")

class ERA5:
    #def __init__(self, start_coord, end_coord, input_file, dt_sec, sim_time_hours, use_time):
    def __init__(self, start_coord): #CHANGED INPUT VARIABLES
        """Create a class object containing information about an ERA5. Along the way,
        it checks whether it can take a subset of your model data. For example, if you have
        data before or after your simulation times, you can ignore these data through saving
        indexes.

        Parameters
        ----------
        start_coord : dict
            Contains information about the starting location.

        Returns
        -------
        ERA5 : Class Object
            Object with our NetCDF data.

        Notes
        -----
        Similar to class GFS which stores NOAA GFS data from NOMAD server

        See Also
        --------
        GFS
        """





        self.dt = config_earth.dt

        '''
        #CHANGES
        file_path = os.path.realpath(__file__)
        logging.debug(f"Initializing ERA5 class in File path = {file_path}")
        full_file_path = netcdf_prop.forecast_dir + input_file
        '''

        file_path = os.path.realpath(__file__)
        logging.debug(f"Initializing ERA5 class in File path = {file_path}")
        full_file_path = "forecasts/" + config_earth.netcdf_era5['filename']

        try:
            logging.debug(f" Extracting netCDF data from {full_file_path}")
            self.file = netCDF4.Dataset(full_file_path, 'r')
        except:
            logging.error(f'Unable to locate netcdf file {full_file_path}, terminating program')
            logging.debug(f'Unable to locate netcdf file {full_file_path}, terminating program')
            print(colored((f'Unable to locate netcdf file {full_file_path}'),"red"))
            sys.exit(1)

        self.geod = Geodesic.WGS84
        self.resolution_hr = config_earth.netcdf_era5['resolution_hr']
        self.num_model_hrs = sim_time_hours = config_earth.simulation["sim_time"]
        self.start_coord = config_earth.simulation["start_coord"]

        '''
        #CHANGES
        if use_time:
            self.init_with_time(start_coord, end_coord, dt_sec, sim_time_hours)
        else:
            self.init_without_time(start_coord, end_coord, dt_sec, sim_time_hours)
        '''

        self.dt = config_earth.dt

        self.init_with_time()


    '''THIS FUNCTION IS UNCESCCARY FOR EarthSHAB?'''
    #def init_without_time(self, start_coord, end_coord, dt_sec, sim_time_hours): # modded the number of varibales for troubleshootin
    '''
    def init_without_time(self, start_coord, end_coord, dt_sec, sim_time_hours): # modded the number of varibales for troubleshootin
        """ adds to class object containing information about an ERA5 but without any time
        variables.

        Parameters
        ----------
        start_coord : dict
            Contains information about the starting location.
        end_coord : dict
            Contains information about the ending location, can contain
            all None or empty values
        dt_sec : float
            Sampling rate of simulation in seconds

        Returns
        -------
        ERA5 : ERA5 Class Object
            Initialize class objectwith our NetCDF data.

        Notes
        -----
        Similar to class GFS which stores NOAA GFS data from NOMAD server

        See Also
        --------

        """
        logging.debug("Now saving model data independent of time. Happens when we use time averaged data.")
        print("Now saving model data independent of time. Happens when we use time averaged data.")


        balloon_launch_time = start_coord['timestamp']
        balloon_last_time = end_coord['end_datetime']

        logging.debug(f"Ballon flight times are from {balloon_launch_time} to {balloon_last_time}")
        logging.info(f"Ballon flight times are from {balloon_launch_time} to {balloon_last_time}")
        print(f"Ballon flight times are from {balloon_launch_time} to {balloon_last_time}")
        print(f"Ballon flight times are from {balloon_launch_time} to {balloon_last_time}")

        self.start_coord = start_coord
        self.min_alt_m = start_coord['alt_m']

        # Initialize min/max lat/lon index values from netcdf4 subset

        # Determine Index values from netcdf4 subset
        lla = self.file.variables['u'][0, :, :]
        self.latlon_index_range(lla)
        logging.debug(f"By default will use the following latitude index boundaries: {self.lat_top_idx}, {self.lat_bot_idx}")
        logging.debug(f"In degrees this covers : {self.file.variables['latitude'][self.lat_top_idx]} to  {self.file.variables['latitude'][self.lat_bot_idx-1]}")
        logging.debug(f"By default will use the following longitude index boundaries: {self.lon_left_idx}, {self.lon_right_idx}")
        logging.debug(f"In degrees this covers : {self.file.variables['longitude'][self.lon_left_idx]} to  {self.file.variables['longitude'][self.lon_right_idx-1]}")
        print(f"By default will use the following latitude index boundaries: {self.lat_top_idx}, {self.lat_bot_idx}")
        print(f"In degrees this covers : {self.file.variables['latitude'][self.lat_top_idx]} to  {self.file.variables['latitude'][self.lat_bot_idx-1]} degrees")
        print(f"By default will use the following longitude index boundaries: {self.lon_left_idx}, {self.lon_right_idx} degrees")
        print(f"In degrees this covers : {self.file.variables['longitude'][self.lon_left_idx]} to  {self.file.variables['longitude'][self.lon_right_idx-1]}")

        # smaller array of downloaded forecast subset
        self.test = self.file.variables['latitude']
        logging.debug(f"Current shape of the latitude variable: {self.file.variables['latitude'].shape}")
        self.lat = self.file.variables['latitude'][self.lat_top_idx:self.lat_bot_idx]
        self.lon = self.file.variables['longitude'][self.lon_left_idx:self.lon_right_idx]
        logging.debug(f"Final shape of the latitude variable: {self.lat.shape}")

        # min/max lat/lon degree values from netcdf4 subset
        self.lat_bot_deg = self.file.variables['latitude'][self.lat_bot_idx-1]
        self.lon_left_deg = self.file.variables['longitude'][self.lon_left_idx]
        self.lat_top_deg = self.file.variables['latitude'][self.lat_top_idx]
        self.lon_right_deg = self.file.variables['longitude'][self.lon_right_idx-1]

        message = "LAT RANGE - min:" + str(self.lat_bot_deg) + " max: " + str(self.lat_top_deg) + " size: " + str(
            abs(self.lat_bot_idx - self.lat_top_idx))
        logging.info(message)
        message = "LON RANGE - min:" + str(self.lon_left_deg) + " max: " + str(self.lon_right_deg) + " size: " + str(
            abs(self.lon_left_idx - self.lon_right_idx))
        logging.info(message)

        # Import the netcdf4 subset to speed up table lookup in this script
        self.levels = self.file.variables['level'][:]
        level_dim, lat_dim, lon_dim = self.file.variables['u'].shape

        # currently, using balloon start and end time to get simulation points.
        # should retest this method below to see if still works
        if self.num_model_hrs != None:
            simulation_points = int(self.num_model_hrs * int(3600 * (1 / dt_sec)))  # Simulation time in seconds
        else:
            simulation_points = None


        g = 9.80665 # gravitation constant used to convert geopotential height to height
        logging.debug(f"Original model wind shape:  {self.file.variables['u'].shape}")
        self.ugdrps0 = self.file.variables['u'][:, self.lat_top_idx:self.lat_bot_idx,
                       self.lon_left_idx:self.lon_right_idx]
        logging.debug(f"Reduced or final model wind shape:  {self.ugdrps0.shape}")
        self.vgdrps0 = self.file.variables['v'][:, self.lat_top_idx:self.lat_bot_idx,
                       self.lon_left_idx:self.lon_right_idx]
        self.hgtprs = self.file.variables['z'][:, self.lat_top_idx:self.lat_bot_idx,
                      self.lon_left_idx:self.lon_right_idx] / g

        #self.temp0 = self.file.variables['t'][start_time_idx:end_time_idx+1, :, self.lat_top_idx:self.lat_bot_idx,
        #             self.lon_left_idx:self.lon_right_idx]

        logging.info("ERA5 Data Finished downloading.\n\n")
    '''
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









        try:
            time_arr = self.file.variables['time']  # do not cast to numpy array yet, time units are different than in GFS
            print(time_arr)
        except:
            time_arr = [0.00, 0.01]
        time_convert = netCDF4.num2date(time_arr[:], time_arr.units, time_arr.calendar)
        time_model = []

        self.model_start_datetime = time_convert[0]
        self.model_end_datetime = time_convert[-1]

        balloon_launch_time = self.start_coord['timestamp']
        #balloon_last_time = end_coord['timestamp']
        start_time_idx = 0
        '''
        #CHANGES
        logging.debug(f"Original model start and times are {self.model_start_datetime} to {self.model_end_datetime}")
        logging.debug(f"Balloon flight times are from {balloon_launch_time} to {balloon_last_time}")
        logging.info(f"Model data spans the time from {time_convert[0]} to {time_convert[-1]}")
        logging.info(f"Balloon flight times are from {balloon_launch_time} to {balloon_last_time}")
        print(f"Model data spans the time from {time_convert[0]} to {time_convert[-1]}")
        print(f"Balloon flight times are from {balloon_launch_time} to {balloon_last_time}")

        # Depending on the dates chosen, we can sometimes remove unused model data the is prior or after
        # our simulation run as performed below

        original = self.model_start_datetime
        for i in range(0, len(time_convert)):
            logging.debug(f"Check if can remove model start times not used #{i}: model={time_convert[i]}, balloon start time: {balloon_launch_time}")
            if time_convert[i] > balloon_launch_time:
                start_time_idx = i-1
                logging.debug(f"Set model time index start to {start_time_idx} with model time of {time_convert[start_time_idx]}")

                if start_time_idx < 0:
                    start_time_idx = 0  # just make sure that your model data start isn't negative
                    logging.debug("Resetting model start time to first instance of model data, did you get the right model data?")

                self.model_start_datetime = time_convert[start_time_idx]
                if original != self.model_start_datetime:
                    logging.debug(f"Changed model start time to {self.model_start_datetime} from {original}")
                break

        if balloon_last_time == None:
            logging.debug("Because you did not specify a balloon last time, run until hits ground or out of model data...")
            end_time_idx = len(time_convert) - 1
        else:
            end_time_idx = len(time_convert) - 1
            original = end_time_idx
            for i in range(start_time_idx, len(time_convert)):
                logging.debug(f"Check if can remove model times not used #{i}: {time_convert[i]}, end time: {balloon_last_time}")
                if time_convert[i] > balloon_last_time:
                    end_time_idx = i
                    self.model_end_datetime = time_convert[end_time_idx]

                    logging.debug(f"Set model time index end to {end_time_idx}, original was set to {original}")
                    logging.debug(f"Model end time set to: {self.model_end_datetime}")
                    break

        ok = self.ok_model_balloon_times(balloon_launch_time, self.model_start_datetime)
        if ok is False:
            sys.exit(2)
        '''
        end_time_idx = len(time_convert) - 1 #Only need this variable from above?

        for i in range(start_time_idx, end_time_idx+1):
            t = time_convert[i]
            dt = t.strftime("%Y-%m-%d %H:%M:%S")
            time_model.append(dt)
            #epoch_time = time.mktime(t.timetuple()) will give an answer dependent on your time zone
            #gmt = timezone('GMT') # already defined above
            #naive_ts = datetime.strptime(dt,"%Y-%m-%d %H:%M:%S")
            #local_ts = gmt.localize(naive_ts)
            #epoch_ts = local_ts.timestamp()
            #epoch_model.append(int(epoch_ts))

        model_hours = (self.datetime_epoch(self.model_end_datetime) - self.datetime_epoch(self.model_start_datetime))/3600

        self.hour_index_check = 0
        self.diff_check = 0.0
        logging.debug(f"First datetime {time_model[0]}, Next datetime {time_model[1]}, last datetime {time_model[-1]}")


        #self.start_coord = start_coord #don't need to define this again
        self.min_alt_m = self.start_coord['alt']

        # Initialize min/max lat/lon index values from netcdf4 subset
        #lat_top_idx = netcdf_prop.netcdf["lat_top"]
        #lat_bot_idx = netcdf_prop.netcdf["lat_bot"]
        #lon_left_idx = netcdf_prop.netcdf["lon_left"]
        #lon_right_idx = netcdf_prop.netcdf["lon_right"]

        # Determine Index values from netcdf4 subset
        lla = self.file.variables['u'][0, 0, :, :]
        self.latlon_index_range(lla)
        logging.debug(f"By default will use the following latitude index boundaries: {self.lat_top_idx}, {self.lat_bot_idx}")
        logging.debug(f"In degrees this covers : {self.file.variables['latitude'][self.lat_top_idx]} to  {self.file.variables['latitude'][self.lat_bot_idx-1]} degrees")
        logging.debug(f"By default will use the following longitude index boundaries: {self.lon_left_idx}, {self.lon_right_idx} degrees")
        logging.debug(f"In degrees this covers : {self.file.variables['longitude'][self.lon_left_idx]} to  {self.file.variables['longitude'][self.lon_right_idx-1]}")

        print(f"By default will use the following latitude index boundaries: {self.lat_top_idx}, {self.lat_bot_idx}")
        print(f"In degrees this covers : {self.file.variables['latitude'][self.lat_top_idx]} to  {self.file.variables['latitude'][self.lat_bot_idx-1]} degrees")
        print(f"By default will use the following longitude index boundaries: {self.lon_left_idx}, {self.lon_right_idx} degrees")
        print(f"In degrees this covers : {self.file.variables['longitude'][self.lon_left_idx]} to  {self.file.variables['longitude'][self.lon_right_idx-1]}")

        # smaller array of downloaded forecast subset
        self.test = self.file.variables['latitude']
        logging.debug(f"Current shape of the latitude variable: {self.file.variables['latitude'].shape}")
        self.lat = self.file.variables['latitude'][self.lat_top_idx:self.lat_bot_idx]
        self.lon = self.file.variables['longitude'][self.lon_left_idx:self.lon_right_idx]
        logging.debug(f"Final shape of the latitude variable: {self.lat.shape}")

        # min/max lat/lon degree values from netcdf4 subset
        self.lat_bot_deg = self.file.variables['latitude'][self.lat_bot_idx-1]
        self.lon_left_deg = self.file.variables['longitude'][self.lon_left_idx]
        self.lat_top_deg = self.file.variables['latitude'][self.lat_top_idx]
        self.lon_right_deg = self.file.variables['longitude'][self.lon_right_idx-1]

        message = "LAT RANGE - min:" + str(self.lat_bot_deg) + " max: " + str(self.lat_top_deg) + " size: " + str(
            abs(self.lat_bot_idx - self.lat_top_idx))
        logging.info(message)
        message = "LON RANGE - min:" + str(self.lon_left_deg) + " max: " + str(self.lon_right_deg) + " size: " + str(
            abs(self.lon_left_idx - self.lon_right_idx))
        logging.info(message)

        # Import the netcdf4 subset to speed up table lookup in this script
        self.levels = self.file.variables['level'][:]
        time_dim, level_dim, lat_dim, lon_dim = self.file.variables['u'].shape
        hour_index = 0  # start simulating at netcdf file start time
        if self.num_model_hrs == None:
            self.simulation_points = int((self.datetime_epoch(self.model_end_datetime) - self.datetime_epoch(self.model_start_datetime))/self.dt)
            end_index = int(hour_index + model_hours+1)
        else:
            end_index = int(hour_index + self.num_model_hrs + 1)
            simulation_points = int(self.num_model_hrs * int(3600 * (1 / self.dt)))  # Simulation time in seconds

        if end_index > time_dim:
            logging.debug(f"Warning: Number of hours to simulate are set at {end_index} but your model data only contains {time_dim} hours. Will reduce to hours to match available model data")
            logging.warning(
                f"Warning: Number of hours to simulate are set at {end_index} but your model data only contains {time_dim} hours. Will reduce to hours to match available model data")
            end_index = time_dim

        g = 9.80665 # gravitation constant used to convert geopotential height to height
        logging.debug(f"Original model wind shape:  {self.file.variables['u'].shape}")
        self.ugdrps0 = self.file.variables['u'][start_time_idx:end_time_idx+1, :, self.lat_top_idx:self.lat_bot_idx,
                       self.lon_left_idx:self.lon_right_idx]
        logging.debug(f"Reduced or final model wind shape:  {self.ugdrps0.shape}")
        self.vgdrps0 = self.file.variables['v'][start_time_idx:end_time_idx+1, :, self.lat_top_idx:self.lat_bot_idx,
                       self.lon_left_idx:self.lon_right_idx]
        self.hgtprs = self.file.variables['z'][start_time_idx:end_time_idx+1, :, self.lat_top_idx:self.lat_bot_idx,
                      self.lon_left_idx:self.lon_right_idx] / g

        #self.temp0 = self.file.variables['t'][start_time_idx:end_time_idx+1, :, self.lat_top_idx:self.lat_bot_idx,
        #             self.lon_left_idx:self.lon_right_idx]

        logging.info("ERA5 Data Finished downloading.\n\n")
    """
    #REMOVED THIS FUNCTION
    def ok_model_balloon_times(self, balloon_launch_time, model_start_datetime):
        '''
        This function checks if we have model data to cover are expect flight times

        :param balloon_launch_time:
        :param model_start_datetime:
        :return:
        '''
        logging.debug("For valid model data, the model start should be before or equal to your simulation start")
        logging.debug("Further, cannot be more than one 59 minutes before")
        logging.debug(f"Model start date is {model_start_datetime}, simulation start is {balloon_launch_time}")
        minutes_diff = (balloon_launch_time - model_start_datetime).total_seconds() / 60.0
        logging.debug(f"Difference in minutes is {minutes_diff}")
        if minutes_diff > 0 and minutes_diff <= 59:
            logging.debug("Model data is in correct time window and passes test.")
            return True
        else:
            logging.error("Model data not cover time of balloon flight and program will exit.")
            logging.error(f"Model start date is {model_start_datetime}, simulation start is {balloon_launch_time}")
            print("Model data not cover time of balloon flight and program will exit.")
            print(f"Model start date is {model_start_datetime}, simulation start is {balloon_launch_time}")
            return False
    """

    def latlon_index_range(self, lla):
        """
        Determine the dimensions of actual data. If you have columns or rows with missing data
        as indicated by NaN, then this function will return your actual shape size so you can
        resize your data being used.
        """
        logging.debug("Scraping given netcdf4 forecast file\n (" + str(self.file) + "\n for subset size")
        logging.info("Scraping given netcdf4 forecast file (" + str(self.file) + " for subset size")
        logging.info(colored('...', 'white', attrs=['blink']))

        logging.debug(f"Shape of u data compent: {lla.shape}")
        results = np.all(~lla.mask)
        if results == False:
            logging.debug("Found missing data inside the latitude, longitude, grid will determine range and set new latitude, longitude boundary indexes.")
            rows, columns = np.nonzero(~lla.mask)
            logging.debug('Row values :', (rows.min(), rows.max()))  # print the min and max rows
            logging.debug('Column values :', (columns.min(), columns.max()))  # print the min and max columns
            self.lat_top_idx = rows.min()
            self.lat_bot_idx = rows.max()
            self.lon_left_idx = columns.min()
            self.lon_right_idx = columns.max()
        else:
            lati, loni = lla.shape
            self.lat_bot_idx = lati
            self.lat_top_idx = 0
            self.lon_right_idx = loni
            self.lon_left_idx = 0

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
    """
    #Removed This Function
    def getNearestAltbyIndexNoTime(self, lat_i, lon_i, alt_m):
        ''' Determines the nearest altitude based off of geopotential height of a .25 degree lat/lon area.
            at a given latitude, longitude index. note, this version has no dimension

            :returns integer index into nearest height
        '''
        i = self.closestIdx(self.hgtprs[:, lat_i, lon_i], alt_m)

        return i
    """

    def wind_alt_Interpolate(self, alt_m, diff_time, lat_idx, lon_idx):
        """
        Performs a 2 step linear interpolation to determine horizontal wind velocity.  First the altitude is interpolated between the two nearest
        .25 degree lat/lon areas.  Once altitude is matched, a second linear interpolation is performed with respect to time.

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
        int_hr_idx = int(hour_index)
        self.diff_check = diff_time
        self.hour_index_check = int_hr_idx
        int_hr_idx2 = int_hr_idx + 1

        hgt_m = alt_m

        #Check to see if we have exceeded available model time
        time_indices, height_indices, lat_indices, lon_indices = self.vgdrps0.shape
        if int_hr_idx2 >= time_indices:
            return [None, None]

        # First interpolate wind speeds between 2 closest time steps to match altitude estimates (hgtprs), which can change with time
        v_0 = self.vgdrps0[int_hr_idx, :, lat_idx, lon_idx]  # Round hour index to nearest int
        v_0 = self.fill_missing_data(v_0)  # Fill the missing wind data if occurs at upper heights

        u_0 = self.ugdrps0[int_hr_idx, :, lat_idx, lon_idx]
        u_0 = self.fill_missing_data(u_0)

        #t_0 = self.temp0[int_hr_idx, :, lat_idx, lon_idx]
        #t_0 = self.fill_missing_data(t_0)

        xp_0 = self.hgtprs[int_hr_idx, :, lat_idx, lon_idx]
        xp_0 = self.fill_missing_data(xp_0)

        #f = interpolate.interp1d(xp_0, t_0, assume_sorted=False, fill_value="extrapolate")
        #t0_pt = f(hgt_m)
        f = interpolate.interp1d(xp_0, u_0, assume_sorted=False, fill_value="extrapolate")
        u0_pt = f(hgt_m)
        f = interpolate.interp1d(xp_0, v_0, assume_sorted=False, fill_value="extrapolate")
        v0_pt = f(hgt_m)

        # Next interpolate the wind velocities with respect to time.
        v_1 = self.vgdrps0[int_hr_idx2, :, lat_idx, lon_idx]
        v_1 = self.fill_missing_data(v_1)
        u_1 = self.ugdrps0[int_hr_idx2, :, lat_idx, lon_idx]
        u_1 = self.fill_missing_data(u_1)
        #t_1 = self.temp0[int_hr_idx2, :, lat_idx, lon_idx]
        #t_1 = self.fill_missing_data(t_1)
        xp_1 = self.hgtprs[int_hr_idx2, :, lat_idx, lon_idx]
        xp_1 = self.fill_missing_data(xp_1)

        #f = interpolate.interp1d(xp_1, t_1, assume_sorted=False, fill_value="extrapolate" )

        #t1_pt = f(hgt_m)
        f = interpolate.interp1d(xp_1, u_1, assume_sorted=False, fill_value="extrapolate")
        u1_pt = f(hgt_m)
        f = interpolate.interp1d(xp_1, v_1, assume_sorted=False, fill_value="extrapolate")
        v1_pt = f(hgt_m)

        fp = [int_hr_idx, int_hr_idx2]
        u = np.interp(hour_index, fp, [u0_pt, u1_pt])
        v = np.interp(hour_index, fp, [v0_pt, v1_pt])
        #t = np.interp(hour_index, fp, [t0_pt, t1_pt])

        return [u, v]

    '''
    #Removed THis Function
    def wind_alt_InterpolateNoTime(self, alt_m, lat_idx, lon_idx):
        """
        Performs a 2 step linear interpolation to determine horizontal wind velocity.  First the altitude is interpolated between the two nearest
        .25 degree lat/lon areas.  Once altitude is matched, a second linear interpolation is performed with respect to time.

        Parameters
        ----------
        alt_m : float
            Contains altitude of balloon object.
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

        hgt_m = alt_m


        # First interpolate wind speeds between 2 closest time steps to match altitude estimates (hgtprs), which can change with time
        v_0 = self.vgdrps0[:, lat_idx, lon_idx]  # Round hour index to nearest int
        v_0 = self.fill_missing_data(v_0)  # Fill the missing wind data if occurs at upper heights

        u_0 = self.ugdrps0[:, lat_idx, lon_idx]
        u_0 = self.fill_missing_data(u_0)

        # t_0 = self.temp0[int_hr_idx, :, lat_idx, lon_idx]
        # t_0 = self.fill_missing_data(t_0)

        xp_0 = self.hgtprs[:, lat_idx, lon_idx]
        xp_0 = self.fill_missing_data(xp_0)

        # f = interpolate.interp1d(xp_0, t_0, assume_sorted=False, fill_value="extrapolate")
        # t0_pt = f(hgt_m)
        f = interpolate.interp1d(xp_0, u_0, assume_sorted=False, fill_value="extrapolate")
        u0_pt = f(hgt_m)
        f = interpolate.interp1d(xp_0, v_0, assume_sorted=False, fill_value="extrapolate")
        v0_pt = f(hgt_m)

        return [u0_pt, v0_pt]
    '''

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
        self.diff_check = diff_time
        self.hour_index_check = int_hr_idx

        lat_idx = self.getNearestLatIdx(coord["lat"], self.lat_top_idx, self.lat_bot_idx)
        lon_idx = self.getNearestLonIdx(coord["lon"], self.lon_left_idx, self.lon_right_idx)
        #z = self.getNearestAltbyIndex(lat_idx, lon_idx, coord["alt_m"])  # fix this for lower and Upper
        z = self.getNearestAltbyIndex(int_hr_idx, lat_idx, lon_idx, coord["alt"])  # fix this for lower and Upper

        try:
            x_wind_vel, y_wind_vel = self.wind_alt_Interpolate(coord['alt'], diff_time, lat_idx, lon_idx)  # for now hour index is 0
            if x_wind_vel == None:
                return [None, None, None, None, None, None, None, None]
                # x_wind_vel, y_wind_vel = self.wind_interpolate_by_index(coord, i, j, z)  # for now hour index is 0
        except:
            logging.error(f'Simulation Timestamp:    , {coord["timestamp"]}')
            logging.error(f"Model Forecast Start time: {self.model_start_datetime}")
            logging.error(f"Last valid difference: {self.diff_check}")
            logging.error(f"Last valid hour index {self.hour_index_check}")
            logging.error(colored(
                "Mismatch with simulation and forecast timstamps. Check simulation start time and/or download a new forecast.",
                "red"))

            sys.exit(1)

        bearing = math.degrees(math.atan2(y_wind_vel, x_wind_vel))
        bearing = 90 - bearing  # perform 90 degree rotation for bearing from wind data
        d = math.pow((math.pow(y_wind_vel, 2) + math.pow(x_wind_vel, 2)), .5) * dt  # dt multiplier
        g = self.geod.Direct(coord["lat"], coord["lon"], bearing, d)

        '''Changed to alt instead of alt_m'''
        if coord["alt"] <= self.min_alt_m:
            # Balloon should remain stationary if it's reached the minimum altitude
            return [coord['lat'], coord['lon'], x_wind_vel, y_wind_vel, bearing, self.lat[lat_idx], self.lon[lon_idx],
                    self.hgtprs[0, z, lat_idx, lon_idx]]  # hgtprs doesn't matter here so is set to 0
        else:
            return [g['lat2'], g['lon2'], x_wind_vel, y_wind_vel, bearing, self.lat[lat_idx], self.lon[lon_idx],
                    self.hgtprs[0, z, lat_idx, lon_idx]]
    '''
    def getNewCoordNoTime(self, coord, dt):
        """
        Determines the new coordinates of the balloon based on the effects of U and V wind components.
        Here the model data does not have a time dimension so model value changes are only a function
        of height (or level) and latitude, longitude.

        Parameters
        ----------
        coord : dict
            Contains current position, altitude and time of balloon object.
        dt : float
            Integration time

        Returns:
            lat_new, lon_new, x_wind_vel, y_wind_vel, bearing, closest_lat, closest_lon, closest alt

        """

        #diff_time = coord["timestamp"] - self.model_start_datetime
        #hour_index = (diff_time.days * 24 + diff_time.seconds / 3600.) / self.resolution_hr
        #int_hr_idx = int(hour_index)
        #self.diff_check = diff_time
        #self.hour_index_check = int_hr_idx

        lat_idx = self.getNearestLatIdx(coord["lat"], self.lat_top_idx, self.lat_bot_idx)
        lon_idx = self.getNearestLonIdx(coord["lon"], self.lon_left_idx, self.lon_right_idx)
        z = self.getNearestAltbyIndexNoTime(lat_idx, lon_idx, coord["alt_m"])  # fix this for lower and Upper

        try:
            x_wind_vel, y_wind_vel = self.wind_alt_InterpolateNoTime(coord['alt_m'], lat_idx, lon_idx)  # for now hour index is 0
            if x_wind_vel == None:
                return [None, None, None, None, None, None, None, None]
                # x_wind_vel, y_wind_vel = self.wind_interpolate_by_index(coord, i, j, z)  # for now hour index is 0
        except:
            logging.error(f'Simulation Timestamp:    , {coord["timestamp"]}')
            logging.error(f"Last valid difference: {self.diff_check}")
            logging.error(f"Last valid hour index {self.hour_index_check}")
            logging.error(colored(
                "Mismatch with simulation and forecast timstamps. Check simulation start time and/or download a new forecast.",
                "red"))

            sys.exit(1)

        bearing = math.degrees(math.atan2(y_wind_vel, x_wind_vel))
        bearing = 90 - bearing  # perform 90 degree rotation for bearing from wind data
        d = math.pow((math.pow(y_wind_vel, 2) + math.pow(x_wind_vel, 2)), .5) * dt  # dt multiplier
        g = self.geod.Direct(coord["lat"], coord["lon"], bearing, d)

        if coord["alt_m"] <= self.min_alt_m:
            # Balloon should remain stationary if it's reached the minimum altitude
            return [coord['lat'], coord['lon'], x_wind_vel, y_wind_vel, bearing, self.lat[lat_idx], self.lon[lon_idx],
                    self.hgtprs[z, lat_idx, lon_idx]]  # hgtprs doesn't matter here so is set to 0
        else:
            return [g['lat2'], g['lon2'], x_wind_vel, y_wind_vel, bearing, self.lat[lat_idx], self.lon[lon_idx],
                    self.hgtprs[z, lat_idx, lon_idx]]
    '''

    #I don't remember why I made this function for Now
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
