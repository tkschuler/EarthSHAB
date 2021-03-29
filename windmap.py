import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import netCDF4
from scipy.interpolate import CubicSpline
import fluids
import config_earth
from datetime import datetime, timedelta

file = netCDF4.Dataset(config_earth.netcdf["nc_file"])
nc_start = config_earth.netcdf["nc_start"]
start_time = config_earth.simulation["start_time"]
res = config_earth.netcdf['res']
coord = config_earth.simulation['start_coord']
LAT = coord["lat"]
LON = coord["lon"]

hours = None
new_timestamp = None

def closest(arr, k):
    return min(range(len(arr)), key = lambda i: abs(arr[i]-k))

lat  = file.variables['lat'][:]
lon  = file.variables['lon'][:]
levels = file.variables['lev'][:]

vgrdprs = file.variables['vgrdprs']
ugrdprs = file.variables['ugrdprs']
hgtprs = file.variables['hgtprs']
tmpprs = file.variables['tmpprs']

def getHourIndex(start_time, nc_start):
    hour_indices = np.linspace(0, 240, num=81)
    hours = (start_time - nc_start).days * 24 + (start_time - nc_start).seconds / 3600
    new_timestamp = start_time + timedelta(hours=hours)
    hour_index = closest(hour_indices, hours)
    new_timestamp = nc_start + timedelta(hours=hour_indices[hour_index])

    return hour_index, new_timestamp

hour_index, new_timestamp = getHourIndex(start_time, nc_start)



def getNearestLat(lat):
    arr = np.arange(start=-90, stop=90.01, step=res)
    i = closest(arr, lat)
    return i

def getNearestLon(lon):
    lon = lon % 360 #convert from -180-180 to 0-360
    arr = np.arange(start=0, stop=360, step=res)

    i = closest(arr, lon)
    return i


def getWindPlotData(hour_index,lat_i,lon_i):
    # Extract relevant u/v wind velocity, and altitude
    u = ugrdprs[hour_index,:,lat_i,lon_i]
    v = vgrdprs[hour_index,:,lat_i,lon_i]
    h = hgtprs[hour_index,:,lat_i,lon_i]

    # Remove missing data
    u = u.filled(np.nan)
    v = v.filled(np.nan)
    nans = ~np.isnan(u)
    u= u[nans]
    v= v[nans]
    h = h[nans]

    # Forecast data is sparse, so use a cubic spline to add more points
    cs_u = CubicSpline(h, u)
    cs_v = CubicSpline(h, v)
    h_new = np.arange(0, 50000, 10) # New altitude range
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

def plotWindVelocity(hour_index,lat,lon):

    # Find location in data
    lat_i = getNearestLat(lat)
    lon_i = getNearestLon(lon)

    bearing0, r0 , colors0, cmap0 = getWindPlotData(hour_index,lat_i,lon_i)
    bearing1, r1 , colors1, cmap1 = getWindPlotData(1,lat_i,lon_i)
    bearing2, r2 , colors2, cmap2 = getWindPlotData(2,lat_i,lon_i)
    bearing3, r3 , colors3, cmap3 = getWindPlotData(3,lat_i,lon_i)

    # Plot figure and legend
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(111, projection='polar')
    sc1 = ax1.scatter(bearing0, colors0, c=r0, cmap='rainbow', alpha=0.75, s = 2)
    ax1.set_xticklabels(['E', '', 'N', '', 'W', '', 'S', ''])
    plt.colorbar(sc1, ax=ax1, label=" Wind Velocity (m/s)")
    ax1.title.set_text("3D Windrose for (" + str(LAT) + ", " + str(LON) + ") on " + str(new_timestamp))

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

def plotTempAlt(hour_index,lat,lon):
        # Find nearest lat/lon in ncdf4 resolution
        plt.figure(figsize=(10, 8))
        lat_i = getNearestLat(lat)
        lon_i = getNearestLon(lon)

        # Extract relevant u/v wind velocity, and altitude
        T = tmpprs[hour_index,:,lat_i,lon_i]
        h = hgtprs[hour_index,:,lat_i,lon_i]

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
        plt.title('Atmospheric Temperature Profile for (' + str(LAT) + ", " + str(LON) + ") on " + str(new_timestamp))

        plt.plot(T,h, label = "GFS Forecast")
        plt.plot(T_atm,el, label = "ISA Model")
        plt.legend(loc='upper right')

'''
plotWindVelocity(hour_index,LAT,LON)
plotTempAlt(hour_index,LAT,LON)
file.close()
plt.show()
'''



