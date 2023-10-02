import math
import solve_states
import GFS
from termcolor import colored
import matplotlib.pyplot as plt
import fluids
import gmplot
import time as tm
import pandas as pd
import os
import windmap

import radiation
import config_earth

""" This file shows an example of how to predict solar balloon trajectories and produces several plots
as well as an html trajectory map that uses Google maps and can be opened in an internet browser.

run saveNETCDF.py before running this file to download a forecast from NOAA.
"""

if not os.path.exists('trajectories'):
    os.makedirs('trajectories')

scriptstartTime = tm.time()


GMT = 7 # MST
dt = config_earth.dt
coord = config_earth.simulation['start_coord']
t = config_earth.simulation['start_time']
start = t
nc_start = config_earth.netcdf["nc_start"]
min_alt = config_earth.simulation['min_alt']
alt_sp = config_earth.simulation['alt_sp']
v_sp = config_earth.simulation['v_sp']
sim_time = config_earth.simulation['sim_time'] * int(3600*(1/dt))
lat = [coord["lat"]]
lon = [coord["lon"]]
GFSrate = config_earth.GFS['GFSrate']
hourstamp = config_earth.netcdf['hourstamp']


atm = fluids.atmosphere.ATMOSPHERE_1976(min_alt)

# Variables for Simulation and Plotting
T_s = [atm.T]
T_i = [atm.T]
T_atm = [atm.T]
el = [min_alt]
v= [0.]
coords = [coord]

ttt = [t - pd.Timedelta(hours=GMT)] #Just for visualizing plot better]
data_loss = False
burst = False
gmap1 = gmplot.GoogleMapPlotter(coord["lat"],coord["lon"], 10)

e = solve_states.SolveStates()
gfs = GFS.GFS(coord)

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
        lat_new,lon_new,x_wind_vel,y_wind_vel,bearing,nearest_lat, nearest_lon, nearest_alt = gfs.getNewCoord(coords[i],dt*GFSrate)  #(coord["lat"],coord["lon"],0,0,0,0,0,0)

    coord_new  =	{
                      "lat": lat_new,                # (deg) Latitude
                      "lon": lon_new,                # (deg) Longitude
                      "alt": el_new,                 # (m) Elevation
                      "timestamp": t,                # Timestamp
                    }

    coords.append(coord_new)
    lat.append(lat_new)
    lon.append(lon_new)

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

#Plots
plt.style.use('seaborn-pastel')
fig, ax = plt.subplots()
ax.plot(ttt,el)
plt.xlabel('Datetime (MST)')
plt.ylabel('Elevation (m)')

fig2, ax2 = plt.subplots()
ax2.plot(ttt,T_s,label="Surface Temperature")
ax2.plot(ttt,T_i,label="Internal Temperature")
ax2.plot(ttt,T_atm,label="Atmospheric Temperature")
#ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M',tz=pytz.timezone(coord['timezone'])))
plt.xlabel('Datetime (MST)')
plt.ylabel('Temperature (K)')
plt.legend(loc='upper right')
plt.title('Solar Balloon Temperature - Earth')

# Outline Downloaded NOAA forecast subset:

region= zip(*[
    (gfs.LAT_LOW, gfs.LON_LOW),
    (gfs.LAT_HIGH, gfs.LON_LOW),
    (gfs.LAT_HIGH, gfs.LON_HIGH),
    (gfs.LAT_LOW, gfs.LON_HIGH)
])

# Google Plotting of Trajectory
gmap1.plot(lat, lon,'red', edge_width = 2.5)
gmap1.polygon(*region, color='cornflowerblue', edge_width=1, alpha= .2)

year = str(tm.localtime()[0])
month = str(tm.localtime()[1]).zfill(2)
day = str(tm.localtime()[2]).zfill(2)
gmap1.draw("trajectories/SHAB_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html" )

executionTime = (tm.time() - scriptstartTime)
print('\nSimulation executed in ' + str(executionTime) + ' seconds.')


plt.style.use('default')
hour_index, new_timestamp = windmap.getHourIndex(start, nc_start)
windmap.plotWindVelocity(hour_index,coord["lat"],coord["lon"])
windmap.plotTempAlt(hour_index,coord["lat"],coord["lon"])
plt.show()
