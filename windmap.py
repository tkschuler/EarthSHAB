"""
This is file generates a 3d windrose plot for a particular coordinate and timestamp.
The polar plot displays information on wind speed and direction  at
various altitudes in a visual format

"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import ticker
#import netCDF4
from termcolor import colored
from scipy.interpolate import CubicSpline
import fluids
import config_earth
from datetime import datetime, timedelta
import sys

import GFS
import ERA5

class Windmap:
    def __init__(self):

        if config_earth.forecast['forecast_type'] == "GFS":
            self.gfs = GFS.GFS(config_earth.simulation['start_coord'])
            self.file = self.gfs.file
        else:
            self.era5 = ERA5.ERA5(config_earth.simulation['start_coord'])
            self.file = self.era5.file


        #FIX THE DEFINING OF THESE Variables
        #*******************
        #self.res = config_earth.netcdf_gfs['res'] #fix this
        self.nc_start = config_earth.netcdf_gfs["nc_start"] #fix this

        self.start_time = config_earth.simulation["start_time"]
        self.coord = config_earth.simulation['start_coord']
        self.LAT = self.coord["lat"]
        self.LON = self.coord["lon"]

        # VERIFY TIMESTAMP INFO MAKES SENSE
        hours = None
        self.new_timestamp = None

        #EDIT TO REMOVE THESE
        #Should I move these to ERA5 and GFS for organization?
        if config_earth.forecast['forecast_type'] == "GFS":
            self.lat  = self.file.variables['lat'][:]
            self.lon  = self.file.variables['lon'][:]
            self.levels = self.file.variables['lev'][:]

            self.vgrdprs = self.file.variables['vgrdprs']
            self.ugrdprs = self.file.variables['ugrdprs']
            self.hgtprs = self.file.variables['hgtprs']
            try:
                self.tmpprs = self.file.variables['tmpprs']
            except:
                print(colored("Temperature not downloaded in this GFS forecast", "yellow"))
                self.tmpprs = None
            self.hour_index, self.new_timestamp = self.getHourIndex(self.start_time)

        if config_earth.forecast['forecast_type'] == "ERA5":

            self.lat  = self.file.variables['latitude'][:]
            self.lon  = self.file.variables['longitude'][:]
            self.levels = self.file.variables['level'][:]

            self.vgrdprs = self.file.variables['v']
            self.ugrdprs = self.file.variables['u']
            self.hgtprs = self.file.variables['z']
            self.tmpprs = self.file.variables['time'] #are t and time the same? #THIS IS WRONG SHOULD BE TEMPERATURE

            self.hour_index, self.new_timestamp = self.getHourIndex(self.start_time)

    def closest(self,arr, k):
        return min(range(len(arr)), key = lambda i: abs(arr[i]-k))

    def time_in_range(self,start, end, x):
        """Return true if x is in the range [start, end]"""
        if start <= end:
            return start <= x <= end
        else:
            return start <= x or x <= end

    #THIS IS ASSUMING THE CORRECT ERA5 TIME PERIOD HAS BEEN DOWNLOADED
    def getHourIndex(self,start_time):
        if config_earth.forecast['forecast_type'] == "GFS":
            times = self.gfs.time_convert
        else:
            times = self.era5.time_convert
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


    def getWind(self,hour_index,lat_i,lon_i, interpolation_frequency = 1):
        """ Calculates a wind vector estimate at a particular 3D coordinate and timestamp
        using a 2-step linear interpolation approach.

        Currently using scipy.interpolat.CubicSpline instead of np.interp like in GFS and ERA5.

        See also :meth:`GFS.GFS.wind_alt_Interpolate`

        :param hour_index: Time index from forecast file
        :type hour_index: int
        :param lat_i: Array index for corresponding netcdf lattitude array
        :type lat_i: int
        :param lon_i: Array index for corresponding netcdf laongitude array
        :type lon_i: int
        :returns: [U, V]
        :rtype: float64 2d array

        """

        g = 9.80665 # gravitation constant used to convert geopotential height to height

        u = self.ugrdprs[hour_index,:,lat_i,lon_i]
        v = self.vgrdprs[hour_index,:,lat_i,lon_i]
        h = self.hgtprs[hour_index,:,lat_i,lon_i]

        # Remove missing data
        u = u.filled(np.nan)
        v = v.filled(np.nan)
        nans = ~np.isnan(u)
        u= u[nans]
        v= v[nans]
        h = h[nans]

        #for ERA5, need to reverse all array so h is increasing.
        if config_earth.forecast['forecast_type'] == "ERA5":
            u = np.flip(u)
            v = np.flip(v)
            h = np.flip(h)
            h = h / g #have to do this for ERA5

        #Fix this interpolation method later, espcially for ERA5
        cs_u = CubicSpline(h, u)
        cs_v = CubicSpline(h, v)
        if config_earth.forecast['forecast_type'] == "GFS":
            h_new = np.arange(0, h[-1], interpolation_frequency) # New altitude range
        elif config_earth.forecast['forecast_type'] == "ERA5":
            h_new = np.arange(0, h[-1], interpolation_frequency) # New altitude range
        u = cs_u(h_new)
        v = cs_v(h_new)

        return self.windVectorToBearing(u, v, h_new)


    def plotWind2(self,hour_index,lat,lon, num_interpolations = 100):
        """ Calculates a wind vector estimate at a particular 3D coordinate and timestamp
        using a 2-step linear interpolation approach.

        I believe this is more accurate than linear interpolating the U-V wind components. Using arctan2
        creates a non-linear distribution of points, however they should still be accurate. arctan2 causes issues
        when there are 180 degree opposing winds because it is undefined, and the speed goes to 0 and back up to the new speed.

        Therefore instead we interpolate the wind speed and direction, and use an angle wrapping check to see if the shortest
        distance around the circle crosses the 0/360 degree axis. This prevents speed down to 0, as well as undefined regions.

        1. Get closest timesteps (t0 and t1) and lat/lon indexes for desired time
        2. Convert the entire altitude range at t0 and t1 U-V wind vector components to wind speed and direction (rad and m/s)
        3. Perform a linear interpolation of wind speed and direction. Fill in num_iterations

        :param hour_index: Time index from forecast file
        :type hour_index: int
        :param lat: Latitude coordinate [deg]
        :type lat_i: float
        :param lon_i: Longitude coordinate [deg]
        :type lon_i: cloast
        :returns: 3D windrose plot
        :rtype: matplotlib.plot

        """

        if config_earth.forecast['forecast_type'] == "GFS":
            lat_i = self.gfs.getNearestLat(lat, -90, 90.01 ) #I think instead of min max this is because download netcdf downloads the whole world, but many of the spots are empty.
            lon_i = self.gfs.getNearestLon(lon, 0, 360 )

        elif config_earth.forecast['forecast_type'] == "ERA5":
            lat_i = self.era5.getNearestLatIdx(lat, self.era5.lat_min_idx, self.era5.lat_max_idx)
            lon_i = self.era5.getNearestLonIdx(lon, self.era5.lon_min_idx, self.era5.lon_max_idx)


        g = 9.80665 # gravitation constant used to convert geopotential height to height

        u = self.ugrdprs[hour_index,:,lat_i,lon_i]
        v = self.vgrdprs[hour_index,:,lat_i,lon_i]
        h = self.hgtprs[hour_index,:,lat_i,lon_i]

        # Remove missing data
        u = u.filled(np.nan)
        v = v.filled(np.nan)
        nans = ~np.isnan(u)
        u= u[nans]
        v= v[nans]
        h = h[nans]

        #for ERA5, need to reverse all array so h is increasing.
        if config_earth.forecast['forecast_type'] == "ERA5":
            u = np.flip(u)
            v = np.flip(v)
            h = np.flip(h)
            h = h / g #have to do this for ERA5

        bearing, r , colors, cmap = self.windVectorToBearing(u, v, h)

        # Create interpolated altitudes and corresponding wind data
        interpolated_altitudes = []
        interpolated_speeds = []
        interpolated_directions_deg = []


        for i in range(len(h) - 1):

            #Do some angle wrapping checks
            interp_dir_deg = 0
            angle1 = np.degrees(bearing[i]) %360
            angle2 = np.degrees(bearing[i + 1]) %360
            angular_difference = abs(angle2-angle1)

            if angular_difference > 180:
                if (angle2 > angle1):
                    angle1 += 360
                else:
                    angle2 += 360


            for j in range(num_interpolations + 1):
                alpha = j / num_interpolations
                interp_alt = h[i] + alpha * (h[i + 1] - h[i])
                interp_speed = np.interp(interp_alt, [h[i], h[i + 1]], [r[i], r[i + 1]])

                interp_dir_deg = np.interp(interp_alt, [h[i], h[i + 1]], [angle1, angle2]) % 360 #make sure in the range (0, 360)

                interpolated_altitudes.append(interp_alt)
                interpolated_speeds.append(interp_speed)
                interpolated_directions_deg.append(interp_dir_deg)

        fig = plt.figure(figsize=(10, 8))
        ax1 = fig.add_subplot(111, projection='polar')

        if config_earth.forecast['forecast_type'] == "GFS":
            sc = ax1.scatter(np.radians(interpolated_directions_deg), interpolated_altitudes, c=interpolated_speeds, cmap='winter', s=2)
            ax1.title.set_text("GFS 3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))

        elif config_earth.forecast['forecast_type'] == "ERA5":
            sc = ax1.scatter(np.radians(interpolated_directions_deg), interpolated_altitudes, c=interpolated_speeds, cmap='winter', s=2)
            ax1.title.set_text("ERA5 3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))

        cbar = plt.colorbar(sc, label='Wind Speed (m/s)')
        #plt.scatter(np.radians(interpolated_directions_deg), interpolated_altitudes)

        # Set title
        fig.suptitle("Wind Interpolation using Wind Speed and Directio Linear Interpolation")
        #plt.title('Windmap with Wind Angles Interpolated')

    '''
    def plotWindOLD(self,hour_index,lat,lon, num_interpolations = 100):

        if config_earth.forecast['forecast_type'] == "GFS":
            lat_i = self.gfs.getNearestLat(lat, -90, 90.01 ) #I think instead of min max this is because download netcdf downloads the whole world, but many of the spots are empty.
            lon_i = self.gfs.getNearestLon(lon, 0, 360 )

        elif config_earth.forecast['forecast_type'] == "ERA5":
            lat_i = self.era5.getNearestLatIdx(lat, self.era5.lat_min_idx, self.era5.lat_max_idx)
            lon_i = self.era5.getNearestLonIdx(lon, self.era5.lon_min_idx, self.era5.lon_max_idx)


        g = 9.80665 # gravitation constant used to convert geopotential height to height

        u = self.ugrdprs[hour_index,:,lat_i,lon_i]
        v = self.vgrdprs[hour_index,:,lat_i,lon_i]
        h = self.hgtprs[hour_index,:,lat_i,lon_i]

        # Remove missing data
        u = u.filled(np.nan)
        v = v.filled(np.nan)
        nans = ~np.isnan(u)
        u= u[nans]
        v= v[nans]
        h = h[nans]

        #for ERA5, need to reverse all array so h is increasing.
        if config_earth.forecast['forecast_type'] == "ERA5":
            u = np.flip(u)
            v = np.flip(v)
            h = np.flip(h)
            h = h / g #have to do this for ERA5


        # Create interpolated altitudes and corresponding wind data
        interpolated_altitudes = []
        interpolated_u = []
        interpolated_v = []

        for i in range(len(h) - 1):
            for j in range(num_interpolations + 1):
                alpha = j / num_interpolations
                interp_alt = h[i] + alpha * (h[i + 1] - h[i])
                interp_u = np.interp(interp_alt, [h[i], h[i + 1]], [u[i], u[i + 1]])
                interp_v = np.interp(interp_alt, [h[i], h[i + 1]], [v[i], v[i + 1]])

                interpolated_altitudes.append(interp_alt)
                interpolated_u.append(interp_u)
                interpolated_v.append(interp_v)

        bearing, r , colors, cmap = self.windVectorToBearing(interpolated_u, interpolated_v, interpolated_altitudes)
        #bearing, r , colors, cmap = self.windVectorToBearing(np.full(len(interpolated_altitudes), 3), np.full(len(interpolated_altitudes), 2), interpolated_altitudes)


        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='polar')

        # Create a scatter plot where radius is altitude, angle is wind direction (in radians), and color represents wind speed
        sc = ax.scatter(bearing, colors, c=r, cmap='winter', s=2)
        cbar = plt.colorbar(sc, label='Wind Speed (m/s)')

        #plt.scatter(np.radians(interpolated_directions_deg), interpolated_altitudes)

        # Set title
        fig.suptitle("Wind Interpolation using OLDDDDDDDD")
        #plt.title('Windmap with Wind Angles Interpolated')
    '''

    def plotWindVelocity(self,hour_index,lat,lon, interpolation_frequency = 1):
        """ Plots a 3D Windrose for a particular coordinate and timestamp from a downloaded forecast.

        :param hour_index: Time index from forecast file
        :type hour_index: int
        :param lat: Latitude
        :type lat: float
        :param lon: Longitude
        :type lon: float
        :returns:

        """

        #Should I remove arguments for the function and just use the initialized functions from config?


        # Find location in data
        if config_earth.forecast['forecast_type'] == "GFS":
            lat_i = self.gfs.getNearestLat(lat, -90, 90.01 ) #I think instead of min max this is because download netcdf downloads the whole world, but many of the spots are empty.
            lon_i = self.gfs.getNearestLon(lon, 0, 360 )
            bearing1, r1 , colors1, cmap1 = self.getWind(hour_index,lat_i,lon_i, interpolation_frequency)

        elif config_earth.forecast['forecast_type'] == "ERA5":
            lat_i = self.era5.getNearestLatIdx(lat, self.era5.lat_min_idx, self.era5.lat_max_idx)
            lon_i = self.era5.getNearestLonIdx(lon, self.era5.lon_min_idx, self.era5.lon_max_idx)
            bearing1, r1 , colors1, cmap1 = self.getWind(hour_index,lat_i,lon_i, interpolation_frequency)

        # Plot figure and legend
        fig = plt.figure(figsize=(10, 8))
        fig.suptitle("Wind Interpolation using  Spline and U-V wind components")
        ax1 = fig.add_subplot(111, projection='polar')
        if config_earth.forecast['forecast_type'] == "GFS":
            sc2 = ax1.scatter(bearing1, colors1, c=r1, cmap='winter', alpha=0.75, s = 2)
            ax1.title.set_text("GFS 3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))

        elif config_earth.forecast['forecast_type'] == "ERA5":
            sc2 = ax1.scatter(bearing1, colors1, c=r1, cmap='winter', alpha=0.75, s = 2)
            ax1.title.set_text("ERA5 3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))

        #ax1.set_xticks([0,9,90,180, 180,270,270,360]) #Fixes the FixedLocator Warning for the line below
        ax1.set_xticks(ax1.get_xticks())
        ax1.set_xticklabels(['E', '', 'N', '', 'W', '', 'S', ''])

        plt.colorbar(sc2, ax=ax1, label=" Wind Velocity (m/s)")


        plt.figure()
        plt.plot(self.ugrdprs[hour_index,:,lat_i,lon_i], self.hgtprs[hour_index,:,lat_i,lon_i])
        plt.plot(self.vgrdprs[hour_index,:,lat_i,lon_i], self.hgtprs[hour_index,:,lat_i,lon_i])
        plt.title("U V Wind Plot")

    """
        #Update this
    def plotTempAlt(self,hour_index,lat,lon):
        if config_earth.forecast['forecast_type'] == "ERA5":
            hour_index = 0 #not right, hard_coded for now

        # Find nearest lat/lon in ncdf4 resolution
        plt.figure(figsize=(10, 8))
        lat_i = self.getNearestLat(lat)
        lon_i = self.getNearestLon(lon)

        # Extract relevant u/v wind velocity, and altitude
        T = self.tmpprs[hour_index,:,lat_i,lon_i]
        h = self.hgtprs[hour_index,:,lat_i,lon_i]

        '''
        # Forecast data is sparse, so use a cubic spline to add more points
        cs_T = CubicSpline(h, T)
        h_new = np.arange(0, 50000, 10) # New altitude range
        T = cs_T(h_new)
        '''

        # ISA Temperature Model
        el = np.arange(0, 50000, 10)
        T_atm = []

        for e in el:
            atm = fluids.atmosphere.ATMOSPHERE_1976(e)
            T_atm.append(atm.T)


        # Formatting
        plt.xlabel("Temperature (K)")
        plt.ylabel("Altitude (m)")
        plt.title('Atmospheric Temperature Profile for (' + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))



        plt.plot(T,h, label = "GFS Forecast")
        plt.plot(T_atm,el, label = "ISA Model")
        plt.legend(loc='upper right')
    """

    def makePlots(self):
        print(self.hour_index)
        self.plotWindVelocity(self.hour_index,self.LAT,self.LON, interpolation_frequency = 100)
        self.plotWind2(self.hour_index,self.LAT,self.LON, num_interpolations = 10)
        #self.plotWindOLD(self.hour_index,self.LAT,self.LON, num_interpolations = 100)
        wind.file.close()
        plt.show()


#wind = Windmap()
#wind.makePlots()
