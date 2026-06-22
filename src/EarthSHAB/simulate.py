import math
from termcolor import colored
import fluids
import gmplot
import time as tm
import pandas as pd
import numpy as np
import os
import re
import sys
import copy

import EarthSHAB.solve_states as solve_states
from EarthSHAB.Forecast import Forecast
import EarthSHAB.radiation as radiation
import EarthSHAB.config as config_earth


def _load_aprs(path: str) -> tuple:
    """Load and normalize an APRS CSV to standard column names.

    Detects format by column names and returns (df, format_name).

    Standard APRS.fi format columns: time, lat, lng, altitude, comment
    Raw APRS logger format columns:  Date, Time(UTC), Latitude, Longitude, Altitude(m),
                                     Internal Temp(C), Pressure(Pa), ...

    Both are normalized so downstream code sees: time, lat, lng, altitude, comment, dt.
    Time is converted to local (MST = UTC-7) in both cases.
    """
    if not os.path.isfile(path):
        print(colored(f"Unable to locate balloon_trajectory file: {path}", "red"))
        sys.exit(1)
    df = pd.read_csv(path)
    cols = set(df.columns)

    if "time" in cols and "lat" in cols and "lng" in cols and "altitude" in cols:
        df["time"] = pd.to_datetime(df["time"])
        df["time"] = df["time"] - pd.to_timedelta(7, unit="h")
        df["dt"] = (
            df["time"].diff()
            .apply(lambda x: x / np.timedelta64(1, "s"))
            .fillna(0)
            .astype("int64")
        )
        return df, "aprs_fi"

    if "Date" in cols and "Time(UTC)" in cols and "Latitude" in cols \
            and "Longitude" in cols and "Altitude(m)" in cols:
        df["time"] = pd.to_datetime(
            df["Date"].astype(str) + " " + df["Time(UTC)"].astype(str)
        )
        df["time"] = df["time"] - pd.to_timedelta(7, unit="h")
        df = df.rename(columns={
            "Latitude":    "lat",
            "Longitude":   "lng",
            "Altitude(m)": "altitude",
        })

        def _make_comment(row):
            parts = []
            t = row.get("Internal Temp(C)")
            p = row.get("Pressure(Pa)")
            if pd.notna(t):
                parts.append(f"{t}C")
            if pd.notna(p):
                parts.append(f"{p}Pa")
            return " ".join(parts)

        df["comment"] = df.apply(_make_comment, axis=1)
        df["dt"] = (
            df["time"].diff()
            .apply(lambda x: x / np.timedelta64(1, "s"))
            .fillna(0)
            .astype("int64")
        )
        return df, "raw"

    if "UTC_Date_Time" in cols and "Lon" in cols and "Lat" in cols and "Altitude_m" in cols:
        df["time"] = pd.to_datetime(df["UTC_Date_Time"], format="%m/%d/%Y %H:%M:%S")
        df["time"] = df["time"] - pd.to_timedelta(7, unit="h")
        # Some files in this format have Lat/Lon data swapped relative to the header;
        # detect by checking whether the "Lon" column holds a positive value (true lat for NM)
        if df["Lon"].iloc[0] > 0:
            df = df.rename(columns={"Lon": "lat", "Lat": "lng", "Altitude_m": "altitude"})
        else:
            df = df.rename(columns={"Lat": "lat", "Lon": "lng", "Altitude_m": "altitude"})
        # StratoTrack files label the pressure column hPa but store Pa; detect by magnitude
        pressure_in_pa = "Pressure_hPa" in df.columns and df["Pressure_hPa"].median() > 1200

        def _make_lightaprs_comment(row):
            parts = []
            t = row.get("Temperature_C")
            p = row.get("Pressure_hPa")
            if pd.notna(t):
                parts.append(f"{t}C")
            if pd.notna(p):
                p_pa = float(p) if pressure_in_pa else float(p) * 100
                parts.append(f"{p_pa:.1f}Pa")
            return " ".join(parts)
        df["comment"] = df.apply(_make_lightaprs_comment, axis=1)
        df["dt"] = (
            df["time"].diff()
            .apply(lambda x: x / np.timedelta64(1, "s"))
            .fillna(0)
            .astype("int64")
        )
        return df, "lightaprs"

    raise ValueError(f"Unrecognized APRS file format. Columns found: {sorted(cols)}")


def _assert_trajectory_within_forecast(df, model_start_datetime, model_end_datetime, gmt_offset_hours=7):
    """Fatal-exit if the balloon trajectory's recorded time span falls outside the forecast.

    ``df['time']`` holds the trajectory's real timestamps in local time — the
    APRS loader shifts UTC down by ``gmt_offset_hours``. Shift them back to UTC
    and require the full ``[first, last]`` span to lie within the forecast's
    ``[model_start_datetime, model_end_datetime]`` coverage.

    This catches the common misconfiguration of pairing a trajectory with a
    forecast from a different day/year (e.g. a 2022 flight against a 2026
    forecast). The reforecast re-anchors the track to ``simulation['start_time']``
    and steps by each APRS ``dt``, so a date mismatch does not crash — it just
    silently applies the wrong winds to the truth altitude profile — which is
    exactly why it has to be checked explicitly here.
    """
    offset = pd.Timedelta(hours=gmt_offset_hours)
    traj_start = (df["time"].min() + offset).to_pydatetime()
    traj_end = (df["time"].max() + offset).to_pydatetime()

    if traj_start < model_start_datetime or traj_end > model_end_datetime:
        print(colored(
            "FATAL: balloon_trajectory recorded time span "
            f"[{traj_start} -> {traj_end}] UTC falls outside the forecast time "
            f"range [{model_start_datetime} -> {model_end_datetime}].\n"
            "The trajectory was recorded outside the downloaded forecast's "
            "coverage (e.g. a flight from a different day/year). Point "
            "balloon_trajectory at a flight the forecast covers, or download a "
            "forecast spanning the flight's dates.",
            "red"))
        sys.exit(1)


class BalloonSimulation:
    def __init__(self):
        self.scriptstartTime = tm.time()

        self.GMT = 7  # 0 # UTC MST
        self.dt = config_earth.simulation['dt']
        self.coord = config_earth.simulation['start_coord']
        self.t = config_earth.simulation['start_time']
        self.start = self.t
        self.min_alt = config_earth.simulation['min_alt']
        self.alt_sp = config_earth.simulation['alt_sp']
        self.v_sp = config_earth.simulation['v_sp']
        self.sim_time = config_earth.simulation['sim_time'] * int(3600 * (1 / self.dt))
        self.lat = [self.coord["lat"]]
        self.lon = [self.coord["lon"]]
        self.forecast_update_interval = config_earth.forecast['forecast_update_interval']
        self.balloon_trajectory = config_earth.simulation['balloon_trajectory']
        self.atm = fluids.atmosphere.ATMOSPHERE_1976(self.min_alt)

        self.trajectory_name = None
        self.df = None
        self.aprs_format = None

        # Get trajectory name from config file for Google Maps:
        if self.balloon_trajectory is not None:
            self.trajectory_name = os.path.splitext(os.path.basename(self.balloon_trajectory))[0]
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

        self.time_local = [self.t - pd.Timedelta(hours=self.GMT)]  # Just for visualizing plot better
        self.data_loss = False
        self.burst = False
        self.gmap1 = gmplot.GoogleMapPlotter(
            self.coord["lat"], self.coord["lon"], 9
        )  # 9 is how zoomed in the map starts, the lower the number the more zoomed out

        self.e = solve_states.SolveStates()

        self.forecast = Forecast(self.coord)
        self.forecast_type = self.forecast.source

        # Load the balloon trajectory up front (when configured) and verify it
        # was recorded within the forecast's time coverage. Doing this in the
        # constructor makes a misconfigured trajectory/forecast pair fail fast —
        # before any simulation runs — for every entry point (main, predict, ...),
        # not just at the reforecast step.
        if self.balloon_trajectory is not None:
            self.df, self.aprs_format = _load_aprs(self.balloon_trajectory)
            _assert_trajectory_within_forecast(
                self.df,
                self.forecast.model_start_datetime,
                self.forecast.model_end_datetime,
                gmt_offset_hours=self.GMT,
            )

        self.lat_aprs_gps = [self.coord["lat"]]
        self.lon_aprs_gps = [self.coord["lon"]]
        self.time_local_aprs = [self.t - pd.Timedelta(hours=self.GMT)]
        self.coords_aprs = [self.coord]

        self.sim_state = None

    def wrap_lon(self, lon):
        """Convert longitude from 0-360 to -180 to 180"""
        # lon = lon % 360
        return np.where(lon > 180, lon - 360, lon)

    def run_simulation(self, run_reforecast=False, descent_correction=False):
        if not run_reforecast:
            descent = False
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

                if descent_correction:
                    if v_new < -3.0 and el_new > 15000:
                        descent = True
                    if descent:
                        v_new = -3
                        el_new = self.el[i] + v_new * self.dt
                        if el_new < self.min_alt:
                            el_new = self.min_alt
                            v_new = 0

                self.T_s.append(T_s_new)
                self.T_i.append(T_i_new)
                self.el.append(el_new)
                self.v.append(v_new)
                self.T_atm.append(T_atm_new)
                self.t = self.t + pd.Timedelta(hours=(1 / 3600 * self.dt))
                self.time_local.append(self.t - pd.Timedelta(hours=self.GMT))

                if i % self.forecast_update_interval == 0:
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
                    ) = self.forecast.getNewCoord(self.coords[i], self.dt * self.forecast_update_interval)

                coord_new = {
                    "lat": lat_new,
                    "lon": lon_new,
                    "alt": el_new,
                    "timestamp": self.t,
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
                        str(self.t - pd.Timedelta(hours=self.GMT))
                        + " el " + str("{:.4f}".format(el_new))
                        + " v " + str("{:.4f}".format(v_new))
                        + " T_s " + str("{:.4f}".format(T_s_new))
                        + " T_i " + str("{:.4f}".format(T_i_new))
                        + " zen " + str(math.degrees(zen))
                    )

                    print(colored(("U wind speed: " + str(x_wind_vel) + " V wind speed: " + str(y_wind_vel) + " Bearing: " + str(bearing)), "yellow"))
                    print(colored(("Lat: " + str(lat_new) + " Lon: " + str(lon_new)), "green"))
                    print(colored(("Nearest Lat: " + str(nearest_lat) + " Nearest Lon: " + str(nearest_lon) +
                                " Nearest Alt: " + str(nearest_alt)), "cyan"))
        else:
            if self.balloon_trajectory is not None:
                # df was loaded and validated against the forecast in __init__.
                self.gmap1.plot(self.df["lat"], self.df["lng"], "white", edge_width=2.5)
                self.gmap1.text(
                    self.coord["lat"] - 0.1,
                    self.coord["lon"] - 0.2,
                    self.trajectory_name + " True Trajectory",
                    color="white",
                )

                alt_aprs = self.df["altitude"].to_numpy()
                time_aprs = self.df["time"].to_numpy()
                dt_aprs = self.df["dt"].to_numpy()
                self.t = config_earth.simulation["start_time"]

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
                    ) = self.forecast.getNewCoord(self.coords_aprs[i], dt_aprs[i])

                    self.t = self.t + pd.Timedelta(seconds=dt_aprs[i + 1])
                    self.time_local_aprs.append(self.t - pd.Timedelta(hours=self.GMT))

                    # Forward-fill NaN altitudes: a corrupted APRS packet yields NaN
                    # altitude, which would cascade — NaN alt → getNewCoord returns
                    # NaN lat/lon → next iteration also returns NaN, propagating
                    # indefinitely.  Use the previous coord's altitude instead.
                    alt_i = alt_aprs[i]
                    if np.isnan(alt_i):
                        alt_i = self.coords_aprs[-1]["alt"]

                    coord_new = {
                        "lat": lat_new,
                        "lon": lon_new,
                        "alt": alt_i,
                        "timestamp": self.t,
                    }

                    print(self.time_local_aprs[i], dt_aprs[i])

                    self.coords_aprs.append(coord_new)
                    self.lat_aprs_gps.append(lat_new)
                    self.lon_aprs_gps.append(lon_new)

                    print(colored(("El: " + str(alt_aprs[i]) + " Lat: " + str(lat_new) + " Lon: " + str(lon_new) + " Bearing: " + str(bearing)), "green"))
        self.process_sim_state()

    def process_sim_state(self):
        self.sim_state = {
            "scriptstartTime": self.scriptstartTime,
            "GMT": self.GMT,
            "dt": self.dt,
            "coord": self.coord,
            "t": self.t,
            "start": self.start,
            "min_alt": self.min_alt,
            "alt_sp": self.alt_sp,
            "v_sp": self.v_sp,
            "sim_time": self.sim_time,
            "lat": self.lat,
            "lon": self.lon,
            "forecast_update_interval": self.forecast_update_interval,
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
            "time_local": self.time_local,
            "data_loss": self.data_loss,
            "burst": self.burst,
            "gmap1": self.gmap1,
            "e": self.e,
            "forecast": self.forecast,
            "lat_aprs_gps": self.lat_aprs_gps,
            "lon_aprs_gps": self.lon_aprs_gps,
            "time_local_aprs": self.time_local_aprs,
            "coords_aprs": self.coords_aprs,
            "df": self.df,
            "aprs_format": self.aprs_format,
        }