import numpy as np
import netCDF4
from termcolor import colored
import math
from geographiclib.geodesic import Geodesic
import sys
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()   #Hacky solution for Python 3.6 to use ISO format Strings

import config_earth

""" GFS.py extracts meteorological data from a NOAA netcdf (.nc file) data set.  To speed up predicting trajectories,
saveNETCDF.py should be run before main.py. This way, the large data set doesn't need to be redownloaded each time the trajectory is run.
For now, only wind velocity is used for the simulation prediction.  Atmposheric properties such as temperarature and pressure are based off of
the U.S. Standard Atmosphere tables from 1976, and the fluids library is used for these, which can be seen in radiation3.py
"""

class GFS:
    def __init__(self, centered_coord):
        # import variables from configuration file
        self.centered_coord = centered_coord
        self.min_alt = config_earth.simulation['min_alt']
        self.start_time = config_earth.simulation['start_time']
        self.hours3 = config_earth.netcdf['hours3']

        self.file = netCDF4.Dataset(config_earth.netcdf["nc_file"])  # Only accepting manual uploads for now
        self.gfs_time = config_earth.netcdf['nc_start']
        self.res = config_earth.netcdf['res']

        self.geod = Geodesic.WGS84

        # Initialize min/max lat/lon index values from netcdf4 subset
        self.lat_low = None
        self.lon_low = None
        self.lat_high = None
        self.lon_high = None

        #Determine Index values from netcdf4 subset
        lla = self.file.variables['ugrdprs'][0,0,:,:]
        self.latlonrange(lla)

        # smaller array of downloaded forecast subset
        self.lat  = self.file.variables['lat'][self.lat_low:self.lat_high]
        self.lon  = self.file.variables['lon'][self.lon_low:self.lon_high]

        # min/max lat/lon degree values from netcdf4 subset
        self.LAT_LOW  = self.file.variables['lat'][self.lat_low]
        self.LON_LOW  = self.file.variables['lon'][self.lon_low]
        self.LAT_HIGH = self.file.variables['lat'][self.lat_high]
        self.LON_HIGH = self.file.variables['lon'][self.lon_high]


        print("LAT RANGE - min:" + str(self.LAT_LOW), " max: " + str(self.LAT_HIGH) + " size: " + str(self.lat_high-self.lat_low+1))
        print("LON RANGE - min:" + str(self.LON_LOW), " max: " + str(self.LON_HIGH) + " size: " + str(self.lon_high-self.lon_low+1))

        hour_index = 0  # start simulating at GFS file start time

        # Import the netcdf4 subset to speed up table lookup in this script
        self.levels = self.file.variables['lev'][:]
        self.ugdrps0 = self.file.variables['ugrdprs'][hour_index:hour_index+self.hours3,:,self.lat_low:self.lat_high,self.lon_low:self.lon_high]
        self.vgdrps0 = self.file.variables['vgrdprs'][hour_index:hour_index+self.hours3,:,self.lat_low:self.lat_high,self.lon_low:self.lon_high]
        self.hgtprs  = self.file.variables['hgtprs'][hour_index:hour_index+self.hours3,:,self.lat_low:self.lat_high,self.lon_low:self.lon_high]

        print("Data downloaded.\n\n")

    def latlonrange(self,lla):
        """
        Determine the lat/lon min/max index values from netcdf4 forecast subset.
        """
        print("Scraping given netcdf4 forecast file for subset size")
        print(colored('...', 'white', attrs=['blink']))


        for i in range (0,len(lla)):
            for j in range (0,len(lla[i])):
                if type(lla[i][j]) is np.float32 and self.lon_low is None:
                    self.lon_low = j

                if type(lla[i][j]) is np.float32 and type(lla[i][j+1]) is not np.float32 and self.lon_high is None:
                    self.lon_high = j

            if self.lon_low is not None and self.lat_low is None:
                self.lat_low = i

            if type(lla[i][self.lon_low]) is np.float32 and type(lla[i+1][self.lon_low]) is not np.float32 and self.lat_high is None:
                self.lat_high = i
                break



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

    def wind_alt_Interpolate(self, coord):
        """Performs a 2 step linear interpolation to determine horizontal wind velocity.  First the altitude is interpolated between the two nearest
        .25 degree lat/lon areas.  Once altitude is matched, a second linear interpolation is performed with respect to time.
        :param coord: Coordinate of balloon
        :type coord: dict
        :returns: [x_wind_vel, y_wind_vel]
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

        try:
            x_wind_vel,y_wind_vel = self.wind_alt_Interpolate(coord) #for now hour index is 0
        except:
            print("Simulation Timestamp:    ", coord["timestamp"])
            print("GFS Forecast Start time: ", self.gfs_time)
            print(colored("Mismatch with simulation and forecast timstamps. Check simulation start time and/or download a new forecast.", "red"))
            sys.exit(1)
        bearing = math.degrees(math.atan2(y_wind_vel,x_wind_vel))
        bearing = 90 - bearing # perform 90 degree rotation for bearing from wind data
        d = math.pow((math.pow(y_wind_vel,2)+math.pow(x_wind_vel,2)),.5) * dt #dt multiplier
        g = self.geod.Direct(coord["lat"], coord["lon"], bearing, d)

        if coord["alt"] <= self.min_alt:
            # Balloon should remain stationary if it's reached the minimum altitude
            return [coord['lat'],coord['lon'],x_wind_vel,y_wind_vel,bearing, self.lat[i], self.lon[j], self.hgtprs[0,z,i,j]] # hgtprs doesn't matter here so is set to 0
        else:
            return [g['lat2'],g['lon2'],x_wind_vel,y_wind_vel,bearing, self.lat[i], self.lon[j], self.hgtprs[0,z,i,j]]
