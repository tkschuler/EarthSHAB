import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import netCDF4
from scipy.interpolate import CubicSpline
import fluids
import config_earth
from datetime import datetime, timedelta


'''
This is a rough update of windmap to include ERA5 support.  However hour indecies are hardcoded
right now.  I think I should be able to make object instances of ERA5.py and GFS.py an call the
variables from those objects instead of what I'm doing here.  That's tomorow's project. 

'''

class Windmap:
    def __init__(self):

        if config_earth.forecast_type == "GFS":
            self.file = netCDF4.Dataset(config_earth.netcdf_gfs["nc_file"])
        else:
            self.file = netCDF4.Dataset("forecasts/" + config_earth.netcdf_era5["filename"])
            lla = self.file.variables['u'][0, 0, :, :]
            self.latlon_index_range(lla)
            self.lat = self.file.variables['latitude'][self.lat_top_idx:self.lat_bot_idx]
            self.lon = self.file.variables['longitude'][self.lon_left_idx:self.lon_right_idx]


        print(self.file)
        print("hey")


        #FIX THE DEFINING OF THESE Variables
        #*******************
        self.res = config_earth.netcdf_gfs['res'] #fix this
        self.nc_start = config_earth.netcdf_gfs["nc_start"] #fix this

        self.start_time = config_earth.simulation["start_time"]
        self.coord = config_earth.simulation['start_coord']
        self.LAT = self.coord["lat"]
        self.LON = self.coord["lon"]

        hours = None
        self.new_timestamp = None


        if config_earth.forecast_type == "GFS":
            self.lat  = self.file.variables['lat'][:]
            self.lon  = self.file.variables['lon'][:]
            self.levels = self.file.variables['lev'][:]

            self.vgrdprs = self.file.variables['vgrdprs']
            self.ugrdprs = self.file.variables['ugrdprs']
            self.hgtprs = self.file.variables['hgtprs']
            self.tmpprs = self.file.variables['tmpprs']

            self.hour_index, self.new_timestamp = self.getHourIndex(self.start_time, self.nc_start)

        if config_earth.forecast_type == "ERA5":
            self.lat  = self.file.variables['latitude'][:]
            self.lon  = self.file.variables['longitude'][:]
            self.levels = self.file.variables['level'][:]

            self.vgrdprs = self.file.variables['v']
            self.ugrdprs = self.file.variables['u']
            self.hgtprs = self.file.variables['z']
            self.tmpprs = self.file.variables['t']

            self.hour_index = 0 #WILL THIS MESS THINGS UP FOR GFS?

    def closest(self,arr, k):
        return min(range(len(arr)), key = lambda i: abs(arr[i]-k))


    def getHourIndex(self,start_time, nc_start):
        #ok What does this function do???
        hour_indices = np.linspace(0, 240, num=81)
        hours = (start_time - nc_start).days * 24 + (start_time - nc_start).seconds / 3600
        new_timestamp = start_time + timedelta(hours=hours)
        hour_index = self.closest(hour_indices, hours)
        new_timestamp = nc_start + timedelta(hours=hour_indices[hour_index])

        return hour_index, new_timestamp


    def getNearestLat(self,lat):
        arr = np.arange(start=-90, stop=90.01, step=self.res)
        i = self.closest(arr, lat)
        return i

    def getNearestLon(self,lon):
        lon = lon % 360 #convert from -180-180 to 0-360
        arr = np.arange(start=0, stop=360, step=self.res)

        i = self.closest(arr, lon)
        return i

    #-------------------------
    #New from ERA5

    def latlon_index_range(self, lla):
        """
        Determine the dimensions of actual data. If you have columns or rows with missing data
        as indicated by NaN, then this function will return your actual shape size so you can
        resize your data being used.
        """
        #logging.debug("Scraping given netcdf4 forecast file\n (" + str(self.file) + "\n for subset size")
        #logging.info("Scraping given netcdf4 forecast file (" + str(self.file) + " for subset size")
        #logging.info(colored('...', 'white', attrs=['blink']))

        #logging.debug(f"Shape of u data compent: {lla.shape}")
        results = np.all(~lla.mask)
        if results == False:
            #logging.debug("Found missing data inside the latitude, longitude, grid will determine range and set new latitude, longitude boundary indexes.")
            rows, columns = np.nonzero(~lla.mask)
            #logging.debug('Row values :', (rows.min(), rows.max()))  # print the min and max rows
            #logging.debug('Column values :', (columns.min(), columns.max()))  # print the min and max columns
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










    #--------------------------------------


    def getWindPlotData(self,hour_index,lat_i,lon_i):
        # Extract relevant u/v wind velocity, and altitude

        if config_earth.forecast_type == "ERA5":
            hour_index = 0 #not right, hard_coded for now

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
        if config_earth.forecast_type == "ERA5":
            u = np.flip(u)
            v = np.flip(v)
            h = np.flip(h)

        if config_earth.forecast_type == "GFS":
            cs_u = CubicSpline(h, u)
            cs_v = CubicSpline(h, v)
            h_new = np.arange(0, 50000, 10) # New altitude range
            u = cs_u(h_new)
            v = cs_v(h_new)


        # Forecast data is sparse, so use a cubic spline to add more points
        if config_earth.forecast_type == "ERA5":
            cs_u = CubicSpline(h, u)
            cs_v = CubicSpline(h, v)
            h_new = np.arange(0, 550000, 10) # New altitude range
            u = cs_u(h_new)
            v = cs_v(h_new)



        # Calculate altitude
        bearing = np.arctan2(v,u)
        bearing = np.unwrap(bearing)
        r = np.power((np.power(u,2)+np.power(v,2)),.5)

        # Set up Color Bar
        colors = h_new
        cmap=mpl.colors.ListedColormap(colors)

        return [bearing, r , colors, cmap]

    def plotWindVelocity(self,hour_index,lat,lon):

        if config_earth.forecast_type == "ERA5":
            hour_index = 0 #not right, hard_coded for now
            self.hour_index = 0 #Fix this later,  match self, and not self

        # Find location in data
        if config_earth.forecast_type == "GFS":
            lat_i = self.getNearestLat(lat)
            lon_i = self.getNearestLon(lon)
        elif config_earth.forecast_type == "ERA5":
            lat_i = self.getNearestLatIdx(lat, self.lat_top_idx, self.lat_bot_idx)
            lon_i = self.getNearestLonIdx(lon, self.lon_left_idx, self.lon_right_idx)

        bearing0, r0 , colors0, cmap0 = self.getWindPlotData(hour_index,lat_i,lon_i)
        bearing1, r1 , colors1, cmap1 = self.getWindPlotData(1,lat_i,lon_i)
        bearing2, r2 , colors2, cmap2 = self.getWindPlotData(2,lat_i,lon_i)
        bearing3, r3 , colors3, cmap3 = self.getWindPlotData(3,lat_i,lon_i)

        # Plot figure and legend
        fig = plt.figure(figsize=(10, 8))
        ax1 = fig.add_subplot(111, projection='polar')
        if config_earth.forecast_type == "GFS":
            sc1 = ax1.scatter(bearing0, colors0, c=r0, cmap='rainbow', alpha=0.75, s = 2)
            ax1.title.set_text("GFS 3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))


        elif config_earth.forecast_type == "ERA5":
            #This is hardcoded for now
            sc1 = ax1.scatter(bearing0[0:50000], colors0[0:50000], c=r0[0:50000], cmap='rainbow', alpha=0.75, s = 2)
            ax1.title.set_text("ERA5 3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))

        ax1.set_xticklabels(['E', '', 'N', '', 'W', '', 'S', ''])
        plt.colorbar(sc1, ax=ax1, label=" Wind Velocity (m/s)")
        #ax1.title.set_text("3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))

        '''
        ax2 = fig.add_subplot(222, projection='polar')
        sc2 = ax2.scatter(bearing1, colors1, c=r1, cmap='rainbow', alpha=0.75, s = 2)
        ax3 = fig.add_subplot(223, projection='polar')
        sc3 = ax3.scatter(bearing2, colors2, c=r2, cmap='rainbow', alpha=0.75, s = 2)
        ax4 = fig.add_subplot(224, projection='polar')
        sc4 = ax4.scatter(bearing3, colors3, c=r3, cmap='rainbow', alpha=0.75, s = 2)
        plt.colorbar(sc1, ax=ax1, label=" Wind Velocity (m/s)")
        plt.colorbar(sc2, ax=ax2, label=" Wind Velocity (m/s)")
        plt.colorbar(sc3, ax=ax3, label=" Wind Velocity (m/s)")
        plt.colorbar(sc4, ax=ax4, label=" Wind Velocity (m/s)")

        # Additional Formatting
        ax1.set_xticklabels(['E', '','N', '', 'W', '', 'S', ''])
        ax2.set_xticklabels(['E', '','N', '', 'W', '', 'S', ''])
        ax3.set_xticklabels(['E', '','N', '', 'W', '', 'S', ''])
        ax4.set_xticklabels(['E', '','N', '', 'W', '', 'S', ''])

        ax1.title.set_text('0600 UTC')
        ax2.title.set_text('1200 UTC')
        ax3.title.set_text('1800 UTC')
        ax4.title.set_text('0000(+1) UTC')
        '''
        #plt.title('Wind Velocity (m/s) as function of Altitudes over Tucson, AZ on ' + date +' ' + str(t) + '00 UTC')

    def plotTempAlt(self,hour_index,lat,lon):
        if config_earth.forecast_type == "ERA5":
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

    def makePlots(self):
        self.plotWindVelocity(self.hour_index,self.LAT,self.LON)
        #self.plotTempAlt(self.hour_index,self.LAT,self.LON)
        wind.file.close()
        plt.show()


wind = Windmap()

wind.makePlots()
