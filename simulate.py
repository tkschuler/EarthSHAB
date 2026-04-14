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


class BalloonSimulation:
    def __init__(self):
        if not os.path.exists('trajectories'):
            os.makedirs('trajectories')

        self.scriptstartTime = tm.time()

        self.GMT = 7  # 0 # UTC MST
        self.dt = config_earth.simulation['dt']
        self.coord = config_earth.simulation['start_coord']
        self.t = config_earth.simulation['start_time']
        self.start = self.t
        self.nc_start = config_earth.netcdf_gfs["nc_start"]
        self.min_alt = config_earth.simulation['min_alt']
        self.alt_sp = config_earth.simulation['alt_sp']
        self.v_sp = config_earth.simulation['v_sp']
        self.sim_time = config_earth.simulation['sim_time'] * int(3600 * (1 / self.dt))
        self.lat = [self.coord["lat"]]
        self.lon = [self.coord["lon"]]
        self.GFSrate = config_earth.forecast['GFSrate']
        self.hourstamp = config_earth.netcdf_gfs['hourstamp']
        self.balloon_trajectory = config_earth.simulation['balloon_trajectory']
        self.forecast_type = config_earth.forecast['forecast_type']
        self.atm = fluids.atmosphere.ATMOSPHERE_1976(self.min_alt)

        self.trajectory_name = None
        self.df = None

        # Get trajectory name from config file for Google Maps:
        if self.balloon_trajectory is not None:
            self.trajectory_name = copy.copy(self.balloon_trajectory)
            replacements = [("balloon_data/", ""), (".csv", "")]
            for pat, repl in replacements:
                self.trajectory_name = re.sub(pat, repl, self.trajectory_name)
            print(self.trajectory_name)

        # Variables for Simulation and Plotting
        self.T_s = [self.atm.T]
        self.T_i = [self.atm.T]
        self.T_atm = [self.atm.T]
        self.el = [self.min_alt]
        self.v = [0.0]
        self.coords = [self.coord]

        self.x_winds_old = [0]
        self.y_winds_old = [0]
        self.x_winds_new = [0]
        self.y_winds_new = [0]

        self.ttt = [self.t - pd.Timedelta(hours=self.GMT)]  # Just for visualizing plot better
        self.data_loss = False
        self.burst = False
        self.gmap1 = gmplot.GoogleMapPlotter(
            self.coord["lat"], self.coord["lon"], 9
        )  # 9 is how zoomed in the map starts, the lower the number the more zoomed out

        self.e = solve_states.SolveStates()

        if self.forecast_type == "GFS":
            self.gfs = GFS.GFS(self.coord)
        else:
            self.gfs = ERA5.ERA5(self.coord)

        self.lat_aprs_gps = [self.coord["lat"]]
        self.lon_aprs_gps = [self.coord["lon"]]
        self.ttt_aprs = [self.t - pd.Timedelta(hours=self.GMT)]
        self.coords_aprs = [self.coord]

        self.sim_state = None

    def wrap_lon(self, lon):
        """Convert longitude from 0-360 to -180 to 180"""
        # lon = lon % 360
        return np.where(lon > 180, lon - 360, lon)

    def run_simulation(self):
        for i in range(0, self.sim_time):
            T_s_new, T_i_new, T_atm_new, el_new, v_new, q_rad, q_surf, q_int = self.e.solveVerticalTrajectory(
                self.t,
                self.T_s[i],
                self.T_i[i],
                self.el[i],
                self.v[i],
                self.coord,
                self.alt_sp,
                self.v_sp,
            )

            self.T_s.append(T_s_new)
            self.T_i.append(T_i_new)
            self.el.append(el_new)
            self.v.append(v_new)
            self.T_atm.append(T_atm_new)
            self.t = self.t + pd.Timedelta(hours=(1 / 3600 * self.dt))
            self.ttt.append(self.t - pd.Timedelta(hours=self.GMT))  # Just for visualizing plot better

            if i % self.GFSrate == 0:
                (
                    lat_new,
                    lon_new,
                    x_wind_vel,
                    y_wind_vel,
                    x_wind_vel_old,
                    y_wind_vel_old,
                    bearing,
                    nearest_lat,
                    nearest_lon,
                    nearest_alt,
                ) = self.gfs.getNewCoord(self.coords[i], self.dt * self.GFSrate)  # (coord["lat"],coord["lon"],0,0,0,0,0,0)

            coord_new = {
                "lat": lat_new,        # (deg) Latitude
                "lon": lon_new,        # (deg) Longitude
                "alt": el_new,         # (m) Elevation
                "timestamp": self.t,   # Timestamp
            }

            self.coords.append(coord_new)
            self.lat.append(lat_new)
            self.lon.append(lon_new)

            self.x_winds_old.append(x_wind_vel_old)
            self.y_winds_old.append(y_wind_vel_old)
            self.x_winds_new.append(x_wind_vel)
            self.y_winds_new.append(y_wind_vel)

            rad = radiation.Radiation()
            zen = rad.get_zenith(self.t, coord_new)

            if i % 360 * (1 / self.dt) == 0:
                print(
                    str(self.t - pd.Timedelta(hours=self.GMT))  # Just for visualizing better
                    + " el " + str("{:.4f}".format(el_new))
                    + " v " + str("{:.4f}".format(v_new))
                    # + " accel " + str("{:.4f}".format(dzdotdt))
                    + " T_s " + str("{:.4f}".format(T_s_new))
                    + " T_i " + str("{:.4f}".format(T_i_new))
                    + " zen " + str(math.degrees(zen))
                )

                print(colored(("U wind speed: " + str(x_wind_vel) + " V wind speed: " + str(y_wind_vel) + " Bearing: " + str(bearing)), "yellow"))
                print(colored(("Lat: " + str(lat_new) + " Lon: " + str(lon_new)), "green"))
                print(colored(("Nearest Lat: " + str(nearest_lat) + " Nearest Lon: " + str(nearest_lon) +
                               " Nearest Alt: " + str(nearest_alt)), "cyan"))

    def plot_and_reforecast(self):
        # Plots
        '''
        Add code to check if the trajectory name exists before running simulations

        Can I just switch order???
        '''

        if self.balloon_trajectory is not None:
            self.df = pd.read_csv(self.balloon_trajectory)
            self.df["time"] = pd.to_datetime(self.df['time'])
            self.df["time"] = self.df['time'] - pd.to_timedelta(7, unit='h')  # Convert to MST
            self.df["dt"] = self.df["time"].diff().apply(lambda x: x / np.timedelta64(1, 's')).fillna(0).astype('int64')
            self.gmap1.plot(self.df['lat'], self.df['lng'], 'white', edge_width=2.5)  # Actual Trajectory
            self.gmap1.text(
                self.coord["lat"] - .1,
                self.coord["lon"] - .2,
                self.trajectory_name + " True Trajectory",
                color='white'
            )

        # Reforecasting
        if self.balloon_trajectory is not None:
            alt_aprs = self.df["altitude"].to_numpy()
            time_aprs = self.df["time"].to_numpy()
            dt_aprs = self.df["dt"].to_numpy()
            self.t = config_earth.simulation['start_time']
            # self.t = config_earth.simulation['start_time']

            for i in range(0, len(alt_aprs) - 1):
                (
                    lat_new,
                    lon_new,
                    x_wind_vel,
                    y_wind_vel,
                    x_wind_vel_old,
                    y_wind_vel_old,
                    bearing,
                    nearest_lat,
                    nearest_lon,
                    nearest_alt,
                ) = self.gfs.getNewCoord(self.coords_aprs[i], dt_aprs[i])

                self.t = self.t + pd.Timedelta(seconds=dt_aprs[i + 1])
                self.ttt_aprs.append(self.t - pd.Timedelta(hours=self.GMT))

                coord_new = {
                    "lat": lat_new,          # (deg) Latitude
                    "lon": lon_new,          # (deg) Longitude
                    "alt": alt_aprs[i],      # (m) Elevation
                    "timestamp": self.t,     # Timestamp
                }

                print(self.ttt_aprs[i], dt_aprs[i])

                self.coords_aprs.append(coord_new)
                self.lat_aprs_gps.append(lat_new)
                self.lon_aprs_gps.append(lon_new)

                print(colored(("El: " + str(alt_aprs[i]) + " Lat: " + str(lat_new) + " Lon: " + str(lon_new) + " Bearing: " + str(bearing)), "green"))

    def run(self):
        self.run_simulation()
        self.plot_and_reforecast()

        self.sim_state = {
            "scriptstartTime": self.scriptstartTime,
            "GMT": self.GMT,
            "dt": self.dt,
            "coord": self.coord,
            "t": self.t,
            "start": self.start,
            "nc_start": self.nc_start,
            "min_alt": self.min_alt,
            "alt_sp": self.alt_sp,
            "v_sp": self.v_sp,
            "sim_time": self.sim_time,
            "lat": self.lat,
            "lon": self.lon,
            "GFSrate": self.GFSrate,
            "hourstamp": self.hourstamp,
            "balloon_trajectory": self.balloon_trajectory,
            "forecast_type": self.forecast_type,
            "atm": self.atm,
            "trajectory_name": self.trajectory_name,
            "T_s": self.T_s,
            "T_i": self.T_i,
            "T_atm": self.T_atm,
            "el": self.el,
            "v": self.v,
            "coords": self.coords,
            "x_winds_old": self.x_winds_old,
            "y_winds_old": self.y_winds_old,
            "x_winds_new": self.x_winds_new,
            "y_winds_new": self.y_winds_new,
            "ttt": self.ttt,
            "data_loss": self.data_loss,
            "burst": self.burst,
            "gmap1": self.gmap1,
            "e": self.e,
            "gfs": self.gfs,
            "lat_aprs_gps": self.lat_aprs_gps,
            "lon_aprs_gps": self.lon_aprs_gps,
            "ttt_aprs": self.ttt_aprs,
            "coords_aprs": self.coords_aprs,
            "df": self.df,
        }


