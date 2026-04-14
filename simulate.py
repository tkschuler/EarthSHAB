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

# --- single line that captures all simulation variables for downstream plotting ---

sim_state = locals()