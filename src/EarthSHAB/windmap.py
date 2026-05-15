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
import EarthSHAB.GFS as GFS
import EarthSHAB.ERA5 as ERA5

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


    def getWind(self, hour_index, lat_i, lon_i, interpolation_frequency=1, interpolation='spline'):
        """ Calculates a wind vector estimate at a particular 3D coordinate and timestamp.

        :param hour_index: Time index from forecast file
        :type hour_index: int
        :param lat_i: Array index for corresponding netcdf lattitude array
        :type lat_i: int
        :param lon_i: Array index for corresponding netcdf longitude array
        :type lon_i: int
        :param interpolation: ``'spline'`` uses CubicSpline on U-V components (default);
            ``'linear'`` uses piecewise linear interpolation
        :type interpolation: str
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

        h_new = np.arange(0, h[-1], interpolation_frequency)

        if interpolation == 'spline':
            cs_u = CubicSpline(h, u)
            cs_v = CubicSpline(h, v)
            u = cs_u(h_new)
            v = cs_v(h_new)
        else:
            u = np.interp(h_new, h, u)
            v = np.interp(h_new, h, v)

        return self.windVectorToBearing(u, v, h_new)


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

        if config_earth.forecast['forecast_type'] == "GFS":
            lat_i = self.gfs.getNearestLat(lat, -90, 90.01)
            lon_i = self.gfs.getNearestLon(lon, 0, 360)
        elif config_earth.forecast['forecast_type'] == "ERA5":
            lat_i = self.era5.getNearestLatIdx(lat, self.era5.lat_min_idx, self.era5.lat_max_idx)
            lon_i = self.era5.getNearestLonIdx(lon, self.era5.lon_min_idx, self.era5.lon_max_idx)

        g = 9.80665

        u = self.ugrdprs[hour_index,:,lat_i,lon_i]
        v = self.vgrdprs[hour_index,:,lat_i,lon_i]
        h = self.hgtprs[hour_index,:,lat_i,lon_i]

        # Remove missing data
        u = u.filled(np.nan)
        v = v.filled(np.nan)
        nans = ~np.isnan(u)
        u = u[nans]
        v = v[nans]
        h = h[nans]

        if config_earth.forecast['forecast_type'] == "ERA5":
            u = np.flip(u)
            v = np.flip(v)
            h = np.flip(h)
            h = h / g

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

        forecast_type = config_earth.forecast['forecast_type']
        interp_label = "Spline" if interpolation == 'spline' else "Linear"

        fig = plt.figure(figsize=(10, 8))
        ax1 = fig.add_subplot(111, projection='polar')
        sc = ax1.scatter(np.radians(interpolated_directions_deg), interpolated_altitudes, c=interpolated_speeds, cmap='winter', s=2)
        ax1.title.set_text(f"{forecast_type} 3D Windrose for ({self.LAT}, {self.LON}) on {self.new_timestamp}")
        plt.colorbar(sc, label='Wind Speed (m/s)')
        fig.suptitle(f"Wind Interpolation using Speed/Direction {interp_label} Interpolation")

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

    def plotWindVelocity(self, hour_index, lat, lon, interpolation_frequency=1, interpolation='spline'):
        """ Plots a 3D Windrose for a particular coordinate and timestamp from a downloaded forecast.

        :param hour_index: Time index from forecast file
        :type hour_index: int
        :param lat: Latitude
        :type lat: float
        :param lon: Longitude
        :type lon: float
        :param interpolation: ``'spline'`` uses CubicSpline on U-V components (default);
            ``'linear'`` uses piecewise linear interpolation
        :type interpolation: str

        """

        if config_earth.forecast['forecast_type'] == "GFS":
            lat_i = self.gfs.getNearestLat(lat, -90, 90.01)
            lon_i = self.gfs.getNearestLon(lon, 0, 360)
        elif config_earth.forecast['forecast_type'] == "ERA5":
            lat_i = self.era5.getNearestLatIdx(lat, self.era5.lat_min_idx, self.era5.lat_max_idx)
            lon_i = self.era5.getNearestLonIdx(lon, self.era5.lon_min_idx, self.era5.lon_max_idx)

        bearing1, r1, colors1, _ = self.getWind(hour_index, lat_i, lon_i, interpolation_frequency, interpolation)

        forecast_type = config_earth.forecast['forecast_type']
        interp_label = "Spline" if interpolation == 'spline' else "Linear"

        fig = plt.figure(figsize=(10, 8))
        fig.suptitle(f"Wind Interpolation using {interp_label} on U-V Components")
        ax1 = fig.add_subplot(111, projection='polar')
        sc2 = ax1.scatter(bearing1, colors1, c=r1, cmap='winter', alpha=0.75, s=2)
        ax1.title.set_text(f"{forecast_type} 3D Windrose for ({self.LAT}, {self.LON}) on {self.new_timestamp}")
        ax1.set_xticks(ax1.get_xticks())
        ax1.set_xticklabels(['E', '', 'N', '', 'W', '', 'S', ''])
        plt.colorbar(sc2, ax=ax1, label="Wind Velocity (m/s)")


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
