"""
This is a rough update of windmap to include ERA5 support.  However hour indecies are hardcoded
right now.  I think I should be able to make object instances of ERA5.py and GFS.py an call the
variables from those objects instead of what I'm doing here.  That's tomorow's project.

testing

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

        if config_earth.forecast_type == "GFS":
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
        if config_earth.forecast_type == "GFS":
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

        if config_earth.forecast_type == "ERA5":

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


        #Verify this function,  maybe include it in the other files era5 and gfs
        #Can probably make this like the ERA5 one to clean up the code a bit
    '''
    def getHourIndex_gfs(self,start_time, nc_start):
        print(self.gfs.time_convert)
        hour_indices = np.linspace(0, 240, num=81)
        hours = (start_time - nc_start).days * 24 + (start_time - nc_start).seconds / 3600
        new_timestamp = start_time + timedelta(hours=hours)
        hour_index = self.closest(hour_indices, hours)
        new_timestamp = nc_start + timedelta(hours=hour_indices[hour_index])

        print(colored("ERA5 simulation start time " + str(self.start_time) + " is not within range of " + str(hours[0]) + " - " + str(hours[-1]) , "red"))


        return hour_index, new_timestamp
    '''

    def time_in_range(self,start, end, x):
        """Return true if x is in the range [start, end]"""
        if start <= end:
            return start <= x <= end
        else:
            return start <= x or x <= end

    #THIS IS ASSUMING THE CORRECT ERA5 TIME PERIOD HAS BEEN DOWNLOADED
    def getHourIndex(self,start_time):
        if config_earth.forecast_type == "GFS":
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
        # Calculate altitude
        bearing = np.arctan2(v,u)
        bearing = np.unwrap(bearing)
        r = np.power((np.power(u,2)+np.power(v,2)),.5)

        # Set up Color Bar
        colors = h
        cmap=mpl.colors.ListedColormap(colors)

        return [bearing, r , colors, cmap]


    def getWind(self,hour_index,lat_i,lon_i):

        g = 9.80665 # gravitation constant used to convert geopotential height to height


        #hard coded hour index to 0 for now?  Update this later.
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
            h = h / g #have to do this for ERA5

        #Fix this interpolation method later, espcially for ERA5
        cs_u = CubicSpline(h, u)
        cs_v = CubicSpline(h, v)
        if config_earth.forecast_type == "GFS":
            h_new = np.arange(0, 50000, 10) # New altitude range
        elif config_earth.forecast_type == "ERA5":
            h_new = np.arange(0, 50000, 10) # New altitude range
        u = cs_u(h_new)
        v = cs_v(h_new)

        return self.windVectorToBearing(u, v, h_new)


    def plotWindVelocity(self,hour_index,lat,lon):

        # Find location in data
        if config_earth.forecast_type == "GFS":
            lat_i = self.gfs.getNearestLat(lat, -90, 90.01 ) #I think instead of min max this is because download netcdf downloads the whole world, but many of the spots are empty.
            lon_i = self.gfs.getNearestLon(lon, 0, 360 )
            bearing1, r1 , colors1, cmap1 = self.getWind(hour_index,lat_i,lon_i)

        elif config_earth.forecast_type == "ERA5":
            lat_i = self.era5.getNearestLatIdx(lat, self.era5.lat_top_idx, self.era5.lat_bot_idx)
            lon_i = self.era5.getNearestLonIdx(lon, self.era5.lon_left_idx, self.era5.lon_right_idx)
            bearing1, r1 , colors1, cmap1 = self.getWind(hour_index,lat_i,lon_i)

        # Plot figure and legend
        fig = plt.figure(figsize=(10, 8))
        ax1 = fig.add_subplot(111, projection='polar')
        if config_earth.forecast_type == "GFS":
            #sc1 = ax1.scatter(bearing0, colors0, c=r0, cmap='rainbow', alpha=0.75, s = 2)
            sc2 = ax1.scatter(bearing1[0:3000], colors1[0:3000], c=r1[0:3000], cmap='winter', alpha=0.75, s = 2)
            ax1.title.set_text("GFS 3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))


        elif config_earth.forecast_type == "ERA5":
            #This is hardcoded for now
            #sc1 = ax1.scatter(bearing0[0:50000], colors0[0:50000], c=r0[0:50000], cmap='rainbow', alpha=0.75, s = 2)
            sc2 = ax1.scatter(bearing1[0:3000], colors1[0:3000], c=r1[0:3000], cmap='winter', alpha=0.75, s = 2)
            ax1.title.set_text("ERA5 3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))

        #ax1.set_xticks([0,9,90,180, 180,270,270,360]) #Fixes the FixedLocator Warning for the line below
        ax1.set_xticks(ax1.get_xticks())
        ax1.set_xticklabels(['E', '', 'N', '', 'W', '', 'S', ''])

        plt.colorbar(sc2, ax=ax1, label=" Wind Velocity (m/s)")


        plt.figure()
        plt.plot(self.ugrdprs[hour_index,:,lat_i,lon_i], self.hgtprs[hour_index,:,lat_i,lon_i])
        plt.plot(self.vgrdprs[hour_index,:,lat_i,lon_i], self.hgtprs[hour_index,:,lat_i,lon_i])
        plt.title("U V Wind Plot")

        #ax1.title.set_text("3D Windrose for (" + str(self.LAT) + ", " + str(self.LON) + ") on " + str(self.new_timestamp))


        #Update this
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
        #self.getWind_ERA5(self.hour_index,self.LAT,self.LON)
        #self.plotTempAlt(self.hour_index,self.LAT,self.LON)
        wind.file.close()
        plt.show()


#wind = Windmap()

#wind.makePlots()
