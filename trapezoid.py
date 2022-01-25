from termcolor import colored
import math
import gmplot
import pandas as pd
import matplotlib.pyplot as plt
import os

import GFS
import radiation
import config_earth


if not os.path.exists('trajectories'):
    os.makedirs('trajectories')

"""
This file shows an example of using a manual altitude profile to generate a balloon trajectory.

This particular altitude profile is trapezoidal in shape with an ascent/descent velocity 2 and 3 m/s respectively and a float altitude
that is specified in config_earth.
"""

GMT = 7

# Import configuration file variables
coord = config_earth.simulation['start_coord']
start = config_earth.simulation['start_time']
t = start
min_alt = config_earth.simulation['min_alt']
float = config_earth.simulation['float']
dt = config_earth.dt
sim = config_earth.simulation["sim_time"]
GFSrate = config_earth.GFS["GFSrate"]

simulation_time = sim*int(3600*(1/dt)) # Simulation time in seconds

# Initialize trajectroy variables
el = [min_alt]  #9000
el_new = min_alt
coords = [coord]
lat = [coord["lat"]]
lon = [coord["lon"]]


gfs = GFS.GFS(coord)
burst = False
gmap1 = gmplot.GoogleMapPlotter(coord["lat"],coord["lon"],8)

ttt=[t]

#sunset = datetime.fromisoformat("2020-02-02 01:30:00")

for i in range(0,simulation_time):
    t = t + pd.Timedelta(hours=(1/3600*dt))

    if i % GFSrate == 0:
        lat_new,lon_new,x_wind_vel,y_wind_vel,bearing,nearest_lat, nearest_lon, nearest_alt = gfs.getNewCoord(coords[i], dt*GFSrate)

    coord_new  =	{
                      "lat": lat_new,                # (deg) Latitude
                      "lon": lon_new,                # (deg) Longitude
                      "alt": el_new,                 # (m) Elevation
                      "timestamp": t,                # Timestamp
                    }

    rad = radiation.Radiation()
    zen = rad.get_zenith(t, coord_new)

    # Trapezoidal Trajectories for faster simulation

    if el_new < float and zen < math.radians(90):
        el_new += 2

    if el_new >= float and zen < math.radians(90):
        el_new = float

    if zen > math.radians(90):
        el_new -= 3


    if el_new < min_alt:
        el_new = min_alt

    el.append(el_new)
    coords.append(coord_new)
    lat.append(lat_new)
    lon.append(lon_new)
    ttt.append(t)

    # Burst Marker
    #if el_new >= 5000:
    #    gmap1.marker(lat_new, lon_new, color='cornflowerblue')


    if i % 360*(1/dt) == 0:
        print(str(t - pd.Timedelta(hours=GMT)) #Just for visualizing better
         +  " el " + str("{:.4f}".format(el_new))
         + " zen " + str(math.degrees(zen))
        )

        print(colored(("U wind speed: " + str(x_wind_vel) + " V wind speed: " + str(y_wind_vel) + " Bearing: " + str(bearing)),"yellow"))
        print(colored(("Lat: " + str(lat_new) + " Lon: " + str(lon_new)),"green"))
        print(colored(("Nearest Lat:" + str(nearest_lat) + " Nearest Lon:" + str(nearest_lon) + " (" + str(360-nearest_lon) +
                        ") Nearest Alt: " + str(nearest_alt)),"cyan"))


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
gmap1.draw( "trajectories/SHAB_trapezoid_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + "_" + str(int(float)) + ".html" )



plt.style.use('seaborn-pastel')
fig, ax = plt.subplots(figsize=(12,10))
ax.plot(ttt,el)
plt.xlabel('Datetime (MST)')
plt.ylabel('Altitude (m)')
plt.title('Altitude Profile for Solar Balloon')

plt.show()
