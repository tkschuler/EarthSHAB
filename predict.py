""" predict.py creates a family of predictions at different float altitudes.  Float altitudes are adjusted by
changing the payload mass in .25 increments.

"""

import math
from termcolor import colored
import matplotlib.pyplot as plt
import fluids
import gmplot
import config_earth
import pandas as pd
from matplotlib.pyplot import cm
import matplotlib as mpl
import os
import sys
import seaborn as sns

import solve_states
import GFS
import radiation
import windmap

if config_earth.forecast_type == "ERA5":
    print(colored("WARNING: Switch forecast type to GFS. This example is for predicting future flights", "yellow"))
    sys.exit()

if not os.path.exists('trajectories'):
    os.makedirs('trajectories')

coord = config_earth.simulation['start_coord']
nc_start = config_earth.netcdf_gfs["nc_start"]
gfs = GFS.GFS(coord)
gmap1 = gmplot.GoogleMapPlotter(coord["lat"], coord["lon"], 10)
hourstamp = config_earth.netcdf_gfs['hourstamp']

masses = [0, .25, .5, .75, 1, 1.25, 1.5, 1.75, 2]

color = cmap = cm.get_cmap('rainbow_r', len(masses))

sns.set_style("darkgrid")
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(1, 1, figsize=(12,8))


for j in range(0,len(masses)):

    print(colored("---------------------------" + str(masses[j]) + "kg-------------------------------","magenta"))

    #Reset Config Values
    GMT = 7  # MST
    coord = config_earth.simulation['start_coord']
    t = config_earth.simulation['start_time']
    start = t
    lat = math.radians(coord['lat'])
    Ls = t.timetuple().tm_yday
    min_alt = config_earth.simulation['min_alt']
    alt_sp = config_earth.simulation['alt_sp']
    v_sp = config_earth.simulation['v_sp']
    dt = config_earth.dt
    atm = fluids.atmosphere.ATMOSPHERE_1976(min_alt)

    GFSrate = config_earth.GFS['GFSrate']

    # Variables for Simulation and Plotting
    T_s = [atm.T]
    T_i = [atm.T]
    T_atm = [atm.T]
    el = [min_alt]
    v = [0.]
    coords = [coord]
    lat = [coord["lat"]]
    lon = [coord["lon"]]
    ttt = [t - pd.Timedelta(hours=GMT)]  # Just for visualizing plot better]
    data_loss = False
    simulation_time = config_earth.simulation["sim_time"] * int(3600 * (1 / dt))  # seconds
    burst = False

    # Set new payload mass to simulate different float altitude
    config_earth.balloon_properties['mp'] = masses[j]

    e = solve_states.SolveStates()

    descent = False
    for i in range(0,simulation_time):
        T_s_new,T_i_new,T_atm_new,el_new,v_new, q_rad, q_surf, q_int = e.solveVerticalTrajectory(t,T_s[i],T_i[i],el[i],v[i],coord,alt_sp,v_sp)

        # Correct for the infrared affects with low masses
        if v_new < -3.0 and el_new > 15000:
            descent = True

        if descent:
            v_new = -3
            el_new = el[i] + v_new * dt

            if el_new < min_alt:
                el_new = min_alt
                v_new = 0

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
    plt.plot(ttt,el, mpl.colors.rgb2hex(color(j)))

    # Google Plotting of Trajectory
    gmap1.plot(lat, lon,mpl.colors.rgb2hex(color(j)), edge_width = 2.5)


# Plotting

plt.xlabel('Datetime (MST)')
plt.ylabel('Elevation (m)')
ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax.grid(visible=True, which='major', color='w', linewidth=1.0)
ax.grid(visible=True, which='minor', color='w', linewidth=0.5)

region= zip(*[
    (gfs.LAT_LOW, gfs.LON_LOW),
    (gfs.LAT_HIGH, gfs.LON_LOW),
    (gfs.LAT_HIGH, gfs.LON_HIGH),
    (gfs.LAT_LOW, gfs.LON_HIGH)
])

gmap1.polygon(*region, color='cornflowerblue', edge_width=1, alpha= .2)
gmap1.draw("trajectories/" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + "_trajectories.html" )

plt.style.use('default')
windmap = windmap.Windmap()
windmap.plotWindVelocity(windmap.hour_index,windmap.LAT,windmap.LON)
#hour_index, new_timestamp = windmap.getHourIndex(start, nc_start)
#windmap.plotWindVelocity(hour_index,coord["lat"],coord["lon"])
#windmap.plotTempAlt(hour_index,coord["lat"],coord["lon"])
plt.show()
