""" This file shows an example of how to predict solar balloon trajectories and produces several plots
as well as an html trajectory map that uses Google maps and can be opened in an internet browser.

run saveNETCDF.py before running this file to download a forecast from NOAA.

Maybe convert to this new library later https://unidata.github.io/python-training/workshop/Bonus/downloading-gfs-with-siphon/

"""

import math
import solve_states
import GFS
import ERA5
from termcolor import colored
import matplotlib.pyplot as plt
import fluids
import gmplot
import time as tm
import pandas as pd
import os
import numpy as np
import re
import copy
import sys

import seaborn as sns
import xarray as xr
from netCDF4 import Dataset

import radiation
import config_earth
import windmap

if not os.path.exists('trajectories'):
    os.makedirs('trajectories')

scriptstartTime = tm.time()

from simulate import sim_state

# Recreate all the original top-level variables in this module:
globals().update(sim_state)

sns.set_palette("muted")
fig, ax = plt.subplots()
ax.plot(ttt,el, label = "reforecasted simulation")
plt.xlabel('Datetime (MST)')
plt.ylabel('Elevation (m)')
if balloon_trajectory != None:
    ax.plot(df["time"],df["altitude"],label = "trajectory")

    if forecast_type == "GFS":
        gmap1.plot(lat_aprs_gps, lon_aprs_gps,'cyan', edge_width = 2.5) #Trajectory using Altitude balloon data with forecast data
        gmap1.text(coord["lat"]-.3, coord["lon"]-.2, trajectory_name + " Alt + " + forecast_type + " Wind Data" , color='cyan')
    elif forecast_type == "ERA5":
        gmap1.plot(lat_aprs_gps, lon_aprs_gps,'orange', edge_width = 2.5) #Trajectory using Altitude balloon data with forecast data
        gmap1.text(coord["lat"]-.3, coord["lon"]-.2, trajectory_name + " Alt + " + forecast_type + " Wind Data" , color='orange')

fig2, ax2 = plt.subplots()
ax2.plot(ttt,T_s,label="Surface Temperature")
ax2.plot(ttt,T_i,label="Internal Temperature")
ax2.plot(ttt,T_atm,label="Atmospheric Temperature")
#ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M',tz=pytz.timezone(coord['timezone'])))
plt.xlabel('Datetime (MST)')
plt.ylabel('Temperature (K)')
plt.legend(loc='upper right')
plt.title('Solar Balloon Temperature - Earth')


def windVectorToBearing(u, v):
    bearing = np.arctan2(v,u)
    speed = np.power((np.power(u,2)+np.power(v,2)),.5)
    return [bearing, speed]

#'''
plt.figure()
plt.plot(ttt, x_winds_new, label = "X winds New", color = "blue")
plt.plot(ttt, x_winds_old, label = "X winds Old", color = "cyan")

plt.plot(ttt, y_winds_new, label = "Y winds New", color = "red")
plt.plot(ttt, y_winds_old, label = "Y winds Old", color = "orange")
#'''

plt.legend(loc='upper right')
plt.title('Wind Interpolation Comparison')

#Winds Figure
plt.figure()
if any(x_winds_old):
    plt.plot(ttt, np.degrees(windVectorToBearing(x_winds_old, y_winds_old)[0]), label = "Bearing old", color = "blue")
plt.plot(ttt, np.degrees(windVectorToBearing(x_winds_new, y_winds_new)[0]), label = "Bearing New", color = "red")
plt.legend(loc='upper right')
plt.title('Wind Interpolation Comparison')




def draw_bounding_box(min_lat, min_lon, max_lat, max_lon):
    """
    Draws a rectangular bounding box on a Google Map and saves it as an HTML file.

    Args:
        min_lat (float): Minimum latitude (southern boundary).
        min_lon (float): Minimum longitude (western boundary).
        max_lat (float): Maximum latitude (northern boundary).
        max_lon (float): Maximum longitude (eastern boundary).
    """
    
    # Build polygon points (bounding box)
    lats = []
    lons = []

    # Bottom edge
    for lon in range(int(min_lon), int(max_lon) + 1, 10):
        lats.append(min_lat)
        lons.append(lon)

    # Right edge
    for lat in range(int(min_lat), int(max_lat) + 1, 10):
        lats.append(lat)
        lons.append(max_lon)

    # Top edge
    for lon in range(int(max_lon), int(min_lon) - 1, -10):
        lats.append(max_lat)
        lons.append(lon)

    # Left edge
    for lat in range(int(max_lat), int(min_lat) - 1, -10):
        lats.append(lat)
        lons.append(min_lon)

    # Draw polygon
    gmap1.polygon(lats, lons, 'cornflowerblue', edge_width=5, alpha= .2)

    # Save map



# Outline Downloaded NOAA forecast subset:

if forecast_type == "GFS":
    '''
    region= zip(*[
        (gfs.LAT_LOW, 0),
        (gfs.LAT_HIGH, 0),
        (gfs.LAT_HIGH, wrap_lon(gfs.lon).max()),
        (gfs.LAT_LOW, wrap_lon(gfs.lon).max())
    ])
    print(wrap_lon(gfs.lon).min(), wrap_lon(gfs.lon).max())
    print(gfs.LAT_LOW, wrap_lon(gfs.LON_LOW), gfs.LAT_HIGH, wrap_lon(gfs.LON_HIGH))
    print("hello")
    
    gmap1.polygon(*region, color='cornflowerblue', edge_width=5, alpha= .2) #plot region
    '''

    gmap1.plot(lat, lon,'blue', edge_width = 2.5) # Simulated Trajectory
    gmap1.text(coord["lat"]-.2, coord["lon"]-.2, 'Simulated Trajectory with GFS Forecast', color='blue')
    draw_bounding_box(gfs.LAT_LOW, wrap_lon(gfs.lon).min(), gfs.LAT_HIGH, wrap_lon(gfs.lon).max())

elif forecast_type == "ERA5":
    region= zip(*[
        (gfs.LAT_LOW, gfs.LON_LOW),
        (80, gfs.LON_LOW),
        (80, gfs.LON_HIGH),
        (gfs.LAT_LOW, gfs.LON_HIGH)
    ])
    gmap1.plot(lat, lon,'red', edge_width = 2.5) # Simulated Trajectory
    gmap1.text(coord["lat"]-.2, coord["lon"]-.2, 'Simulated Trajectory with ERA5 Reanalysis', color='red')
    gmap1.polygon(*region, color='orange', edge_width=1, alpha= .15) #plot region


year = str(tm.localtime()[0])
month = str(tm.localtime()[1]).zfill(2)
day = str(tm.localtime()[2]).zfill(2)

if balloon_trajectory != None:
    if forecast_type == "GFS":
        gmap1.draw("trajectories/" + trajectory_name +"_GFS_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html" )

    elif forecast_type == "ERA5":
        gmap1.draw("trajectories/" + trajectory_name +"_ERA5_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html" )
else:
    if forecast_type == "GFS":
        gmap1.draw("trajectories/PREDICTION_GFS_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html" )

    elif forecast_type == "ERA5":
        gmap1.draw("trajectories/PREDICTION_ERA5_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html" )


executionTime = (tm.time() - scriptstartTime)
print('\nSimulation executed in ' + str(executionTime) + ' seconds.')

windmap = windmap.Windmap()
#windmap.plotWindVelocity(windmap.hour_index,windmap.LAT,windmap.LON)
windmap.plotWind2(windmap.hour_index,windmap.LAT,windmap.LON)
#windmap.plotWindOLD(windmap.hour_index,windmap.LAT,windmap.LON)

plt.show()
