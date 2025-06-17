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

GMT = 7 #0 # UTC MST
dt = config_earth.simulation['dt']
coord = config_earth.simulation['start_coord']
t = config_earth.simulation['start_time']
start = t
nc_start = config_earth.netcdf_gfs["nc_start"]
min_alt = config_earth.simulation['min_alt']
alt_sp = config_earth.simulation['alt_sp']
v_sp = config_earth.simulation['v_sp']
sim_time = config_earth.simulation['sim_time'] * int(3600*(1/dt))
lat = [coord["lat"]]
lon = [coord["lon"]]
GFSrate = config_earth.forecast['GFSrate']
hourstamp = config_earth.netcdf_gfs['hourstamp']
balloon_trajectory = config_earth.simulation['balloon_trajectory']
forecast_type = config_earth.forecast['forecast_type']
atm = fluids.atmosphere.ATMOSPHERE_1976(min_alt)

def wrap_lon(lon):
    """Convert longitude from 0-360 to -180 to 180"""
    #lon = lon % 360
    return np.where(lon > 180, lon - 360, lon)


#Get trajectory name from config file for Google Maps:
if balloon_trajectory != None:
    trajectory_name = copy.copy(balloon_trajectory)
    replacements=[("balloon_data/", ""), (".csv", "")]
    for pat,repl in replacements:
        trajectory_name = re.sub(pat, repl, trajectory_name)
    print (trajectory_name)


# Variables for Simulation and Plotting
T_s = [atm.T]
T_i = [atm.T]
T_atm = [atm.T]
el = [min_alt]
v= [0.]
coords = [coord]

x_winds_old = [0]
y_winds_old = [0]
x_winds_new = [0]
y_winds_new = [0]

ttt = [t - pd.Timedelta(hours=GMT)] #Just for visualizing plot better]
data_loss = False
burst = False
gmap1 = gmplot.GoogleMapPlotter(coord["lat"],coord["lon"], 9) #9 is how zoomed in the map starts, the lower the number the more zoomed out

e = solve_states.SolveStates()

if forecast_type == "GFS":
    gfs = GFS.GFS(coord)
else:
    gfs = ERA5.ERA5(coord)

lat_aprs_gps = [coord["lat"]]
lon_aprs_gps = [coord["lon"]]
ttt_aprs = [t - pd.Timedelta(hours=GMT)]
coords_aprs = [coord]

for i in range(0,sim_time):
    T_s_new,T_i_new,T_atm_new,el_new,v_new, q_rad, q_surf, q_int = e.solveVerticalTrajectory(t,T_s[i],T_i[i],el[i],v[i],coord,alt_sp,v_sp)

    T_s.append(T_s_new)
    T_i.append(T_i_new)
    el.append(el_new)
    v.append(v_new)
    T_atm.append(T_atm_new)
    t = t + pd.Timedelta(hours=(1/3600*dt))
    ttt.append(t - pd.Timedelta(hours=GMT)) #Just for visualizing plot better


    if i % GFSrate == 0:
        lat_new,lon_new,x_wind_vel,y_wind_vel, x_wind_vel_old, y_wind_vel_old, bearing,nearest_lat, nearest_lon, nearest_alt = gfs.getNewCoord(coords[i],dt*GFSrate)  #(coord["lat"],coord["lon"],0,0,0,0,0,0)
    coord_new  =	{
                      "lat": lat_new,                # (deg) Latitude
                      "lon": lon_new,                # (deg) Longitude
                      "alt": el_new,                 # (m) Elevation
                      "timestamp": t,                # Timestamp
                    }

    coords.append(coord_new)
    lat.append(lat_new)
    lon.append(lon_new)

    x_winds_old.append(x_wind_vel_old)
    y_winds_old.append(y_wind_vel_old)
    x_winds_new.append(x_wind_vel)
    y_winds_new.append(y_wind_vel)

    rad = radiation.Radiation()
    zen = rad.get_zenith(t, coord_new)

    if i % 360*(1/dt) == 0:
        print(str(t - pd.Timedelta(hours=GMT)) #Just for visualizing better
         +  " el " + str("{:.4f}".format(el_new))
         + " v " + str("{:.4f}".format(v_new))
         #+ " accel " + str("{:.4f}".format(dzdotdt))
         + " T_s " + str("{:.4f}".format(T_s_new))
         + " T_i " + str("{:.4f}".format(T_i_new))
         + " zen " + str(math.degrees(zen))
        )

        print(colored(("U wind speed: " + str(x_wind_vel) + " V wind speed: " + str(y_wind_vel) + " Bearing: " + str(bearing)),"yellow"))
        print(colored(("Lat: " + str(lat_new) + " Lon: " + str(lon_new)),"green"))
        print(colored(("Nearest Lat: " + str(nearest_lat) + " Nearest Lon: " + str(nearest_lon) +
                        " Nearest Alt: " + str(nearest_alt)),"cyan"))

df = None
#Plots
'''
Add code to check if the trajectory name exists before running simulations

Can I just switch order???
'''

if balloon_trajectory != None:
    df = pd.read_csv(balloon_trajectory)
    df["time"] = pd.to_datetime(df['time'])
    df["time"] = df['time'] - pd.to_timedelta(7, unit='h') #Convert to MST
    df["dt"] = df["time"].diff().apply(lambda x: x/np.timedelta64(1, 's')).fillna(0).astype('int64')
    gmap1.plot(df['lat'], df['lng'],'white', edge_width = 2.5) # Actual Trajectory
    gmap1.text(coord["lat"]-.1, coord["lon"]-.2, trajectory_name + " True Trajectory", color='white')

#Reforecasting
if balloon_trajectory != None:
    alt_aprs = df["altitude"].to_numpy()
    time_aprs = df["time"].to_numpy()
    dt_aprs = df["dt"].to_numpy()
    t = config_earth.simulation['start_time']
    #t = config_earth.simulation['start_time']

    for i in range(0,len(alt_aprs)-1):

        lat_new,lon_new,x_wind_vel,y_wind_vel, x_wind_vel_old, y_wind_vel_old, bearing,nearest_lat, nearest_lon, nearest_alt = gfs.getNewCoord(coords_aprs[i],dt_aprs[i])

        t = t + pd.Timedelta(seconds=dt_aprs[i+1])
        ttt_aprs.append(t - pd.Timedelta(hours=GMT))


        coord_new  =	{
                          "lat": lat_new,                # (deg) Latitude
                          "lon": lon_new,                # (deg) Longitude
                          "alt": alt_aprs[i],                 # (m) Elevation
                          "timestamp": t,                # Timestamp
                        }

        print(ttt_aprs[i], dt_aprs[i])

        coords_aprs.append(coord_new)
        lat_aprs_gps.append(lat_new)
        lon_aprs_gps.append(lon_new)

        print(colored(("El: " + str(alt_aprs[i]) + " Lat: " + str(lat_new) + " Lon: " + str(lon_new) + " Bearing: " + str(bearing)),"green"))


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

'''
plt.figure()
plt.plot(ttt, x_winds_new, label = "X winds New", color = "blue")
plt.plot(ttt, x_winds_old, label = "X winds Old", color = "cyan")

plt.plot(ttt, y_winds_new, label = "Y winds New", color = "red")
plt.plot(ttt, y_winds_old, label = "Y winds Old", color = "orange")
'''

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
