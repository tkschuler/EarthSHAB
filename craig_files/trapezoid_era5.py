"""
Trapezoid program designed to work with ECMWF ERA5 reanalysis data. These data reverse the latitude, longitude, and height coordinates
used by NOAA's GFS data (that data grabbed from NOMAD server). So for example, GFS data starts from the bottom (e.g., 1000mb) and each
new level is higher; in contrast ERA5 data starts from the top with the highest level available in the dataset. So for example, your first
level in ERA5 might be 10mb.

Also, latitude, longitude order is reversed from that of GFS

version 1.2 (02/11/2022)

see also: trapezoid.py
"""
import logging
import sys
import time
from termcolor import colored
import math
import gmplot
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from datetime import datetime
import datetime as dt
from dateutil import parser
import numpy as np
import ERA5 as ERA
import properties.radiation as radiation
import geojson_save as geo
import csv_save as csv
import mathtime_utilities as mt
import telemetry_retriever as tm

logging.basicConfig(filename='trapezoid_era5.log',
                    level=logging.WARNING, # Changed from Info for Demo
                    format='%(asctime)s | %(levelname)s | %(message)s',
                    datefmt="%Y-%m-%d %H:%M")

from properties import simulation_plan as sim_prop
import load_config as lc

def log_configuration(rc):
    """Logs parameters to start the simulation

    Parameters:
    rc (run_config object): Class that stores input paremeters

    Returns:
    None

   """
    readme = "Running trapezoid_era5 version 1.2 (03-01-2022)"
    logging.info(readme)
    logging.info("=======================================================================================================")
    logging.info("Runtime settings are taking from simulation_plan, flight_plans, and to a lesser extent era5_properties and balloon_properties")
    logging.info("We start the simulation as follows: ")
    logging.info(f"flight_id: {rc.flight_id}")
    logging.info(f"NetCDF data will be obtained from {rc.input_file}")
    logging.info(f"GMT offset: {rc.GMT}")
    if rc.compare_to_telemetry:
        logging.info(f"telemetry_dir: {rc.telemetry_file}")
    logging.info(f"gfs_rate_s in seconds: {rc.GFS_RATE_SEC}, sampling rate or dt_sec: {rc.DT_SEC}")

    logging.info(f"balloon_type: {rc.balloon_type}")
    logging.info(f"Balloon start and end datatime:{rc.BALLOON_START_DATETIME}, {rc.BALLOON_END_DATETIME}")
    logging.info(f"launch_latitude: {rc.launch_latitude}, launch_longitude: {rc.launch_longitude}, launch_altitude: {rc.launch_altitude}")
    logging.info(f"last_latitude: {rc.last_latitude}, last_longitude: {rc.last_longitude}, last_altitude: {rc.last_altitude}")
    logging.info(f"Minimum alt: {rc.MIN_ALT_M}")
    logging.info(" ")
    logging.info(f"Simulation time in hours: {rc.SIM_TIME_HOURS} and iterations: {rc.SIM_TIME_ITERATIONS}")
    logging.info(f"ascent_rate_ms: {rc.ASCENT_RATE_MS}, descent_rate_ms: {rc.descent_rate_ms}")


    if rc.strato_boundary_m != None:
        logging.info(f"strato_boundary_m: {rc.strato_boundary_m}")
        logging.info(f"strato_ascent_rate_ms: {rc.strato_ascent_rate_ms}, strato_descent_rate_ms: {rc.strato_descent_rate_ms}")

    logging.info(f"burst_hgt_m: {rc.burst_hgt_m}")


def get_era5(rc):
    """
    Call to get ERA5 data once and store to speed access

    :parameter sim_config: class object created in sim_config which are our control parameters.
                abbreviated as rc for "run" "config

    Returns
    -------
        ERA5 class populated with model data for time and area specified in simulation_plan.py and flight_plans.py
    """

    model = ERA.ERA5(rc.start_coord, rc.end_coord, rc.input_file, rc.DT_SEC,
                     rc.SIM_TIME_HOURS, rc.use_time_dimension)
    return model


def get_height(hgt):
    '''
    Uses balloon ascent, descent, float heights and burst heights to get height estimates

    :param hgt: current height in meters
    :return: new height in meters
    '''
    el_new = hgt["elevation_m"]
    if hgt["burst_occurred"] == False:
        if el_new < hgt["min_alt_m"]:
            el_new = hgt["min_alt_m"]
        elif hgt["ceiling_m"] != None and hgt["burst"] == False:
            if el_new < hgt["ceiling_m"]:
                el_new += hgt["ascent_ms"] * hgt["dt_sec"]
            else:
                el_new = hgt["ceiling_m"]
        elif el_new < hgt["burst_m"]:
            # elif el_new < burst_hgt_m and zen < math.radians(90):
            el_new += hgt["ascent_ms"] * hgt["dt_sec"]
        elif el_new >= hgt["burst_m"] and hgt["burst"]:
            hgt["burst_occurred"] = True
    else:
        el_new += hgt["descent_ms"] * hgt["dt_sec"]

    hgt["elevation_m"] = el_new
    return hgt


def get_new_height(h_m, epoch, timestamp):
    """Here we pass the current timestamp, then we perform linear interpolation to get the index of epoch array, that
    is closest to the the timestamp. Lastly, we use the index to grab the height in meters from h_m.
    :param h_m: An array of Height of balloon in meters
    :param epoch: An array of epoch time values
    :param timestamp: timestamp to use in interpolating values from epoch
    :returns: height associated with where the balloon is in time
    """
    #Given an ordered array and a value, determines the index of the closest item contained in the array.

    closest_idx =  min(range(len(epoch)), key=lambda i: abs(epoch[i] - timestamp))
    height_m = h_m[closest_idx]

    return height_m


def compare_to_telemetry(df, gmap, source):
    '''

    :param df: dataframe from csv or excel
    :param gmap: google map object
    :param source: string that tells us where data came from, helps in parsing
    :return:
    '''

    if source == "NRL":

        gps_hgt_m = df['altitude'].tolist()
        gps_lat = df['lat'].tolist()
        gps_lon = df['lng'].tolist()
        gps_lat_arr = np.array(gps_lat)
        gps_lon_arr = np.array(gps_lon)

        # NRL dataset has a lot of missing numbers which must be removed before plotting (Navy doesn't)
        nan_array = np.isnan(gps_lat_arr)
        not_nan_array = ~ nan_array
        gpslat_deg = gps_lat_arr[not_nan_array]
        nan_array = np.isnan(gps_lon_arr)
        not_nan_array = ~ nan_array
        gpslon_deg = gps_lon_arr[not_nan_array]
    elif source == 'APRSFI' or source == 'APRS' or source == 'APRS2':
        gps_hgt_m = df['altitude'].tolist()
        gpslat_deg = df['latitude'].tolist()
        gpslon_deg = df['longitude'].tolist()

    gmap.plot(gpslat_deg, gpslon_deg, 'blue', edge_width=2.5)
    gmap.marker(gpslat_deg[-1], gpslon_deg[-1], 'blue', title="End Observed Location")


def main(flight_id, rc):
    if rc == None:
        savenetcdf = False
        rc = lc.sim_config(savenetcdf, flight_id)
        logging.info(f"Loading config from {flight_id}") # switch to info level logging when tests are done
    else:
        logging.error('Config provided') # switch to info level logging when tests are done
    readme = "Running trapezoid_era version 1.3 (03-01-2022)"
    logging.debug(readme)
    if rc.compare_to_telemetry:
        df, rc = tm.get_flight_plan_by_telemetry(rc)
        h_m = df['altitude']
        times = df['Time']
        try:
            df['date'] = pd.to_datetime(df['Time'])
            epoch = (df['date'] - dt.datetime(1970, 1, 1)).dt.total_seconds()
            logging.debug(pd.to_datetime(epoch[0], unit='s'))
        except:
            logging.debug("ugh, unable to convert timestamp to unix epoch, will take nap, bye")
            sys.exit(-1)

    log_configuration(rc)
    model = get_era5(rc)
    logging.info(readme)
    logging.info(f"Balloon will be launched at {rc.start_coord['lat']}, {rc.start_coord['lon']}")
    logging.info(f"Starting at {rc.start_coord['timestamp']}")
    GMT = 7

    # Import configuration file variables
    coord = rc.start_coord
    start = rc.BALLOON_START_DATETIME
    end = rc.BALLOON_END_DATETIME
    first_sec = start
    t = start
    min_alt_m = rc.MIN_ALT_M
    # float_hgt_m = sim_prop.simulation['altitude_setpt_m']
    burst_hgt_m = rc.burst_hgt_m
    burst_balloon = rc.burst_balloon
    ceiling_hgt_m = rc.ceiling_hgt_m
    dt_sec = rc.DT_SEC
    SIM_TIME_HOURS = rc.SIM_TIME_HOURS
    GFSrate_sec = rc.GFS_RATE_SEC
    balloon_name = rc.BALLOON_NAME
    balloon_ascent_ms = rc.ASCENT_RATE_MS
    balloon_descent_ms = rc.descent_rate_ms

    # determine how many simulations we will do as a function of our flight time and integration seconds
    if rc.use_telemetry_set:
        print("Telemetry file provided the following:")
        sim_time_sec = rc.sim_time_sec
        SIM_TIME_ITERATIONS = rc.SIM_TIME_ITERATIONS
        print(f"Simulation seconds {sim_time_sec} and {SIM_TIME_ITERATIONS} iterations")
        logging.info("Telemetry file provided the following:")
        logging.info(f"Simulation seconds {sim_time_sec} and {SIM_TIME_ITERATIONS} iterations")
        SIM_TIME_HOURS = sim_time_sec/3600.0
        print(f"This will run simulation for approximately {SIM_TIME_HOURS} hours.")
        logging.info(f"This will run simulation for approximately {SIM_TIME_HOURS} hours.")
    else:
        if SIM_TIME_HOURS != None:
            sim_time_sec = SIM_TIME_HOURS * 3600.0
            SIM_TIME_ITERATIONS = int((SIM_TIME_HOURS * 3600.0)/ dt_sec)  # number of simulation points
        else:
            sim_time_sec = mt.datetime_epoch(end) - mt.datetime_epoch(start)
            SIM_TIME_ITERATIONS = int(sim_time_sec / dt_sec)

    simulation_points = SIM_TIME_ITERATIONS
    logging.debug(f"Simulation time in hours: {SIM_TIME_HOURS} and iterations: {SIM_TIME_ITERATIONS}")
    # Initialize trajectory variables
    el = [min_alt_m]  # 9000
    coord_new = coord

    lat_deg = [coord["lat"]]  # won't be used in trapezoid routines
    lon_deg = [coord["lon"]]  # """
    balloon_datetime = [coord["timestamp"]]
    height_m = el.copy()

    ttt = [t - pd.Timedelta(hours=GMT)]  # Just for visualizing plot better]
    gfscount = 0

    first_burst = False

    # the following segment is used when we float a balloon and then
    # tell it to land. We need to know how long we need to start descending
    if ceiling_hgt_m != None and burst_balloon == False:
        logging.debug("determine how many seconds you will need to descend your balloon")
        distance = ceiling_hgt_m - min_alt_m
        time_in_sec = abs(distance/balloon_descent_ms)
        stop_sec = sim_time_sec - time_in_sec

    if rc.compare_to_telemetry == True:
        logging.debug("Do simulation using telemetry based heights... ")
        for i in range(0, simulation_points):
            t = t + pd.Timedelta(hours=(1 / 3600 * dt_sec))
            epoch_t = mt.datetime_epoch(t)

         # Store coordinates that change every iteration

            if i % GFSrate_sec == 0:
                el_new = get_new_height(h_m, epoch, epoch_t)
                coord_new['alt_m'] = el_new
                coord_new['timestamp'] = t
                if rc.use_time_dimension:
                    lat_new, lon_new, x_wind_vel, y_wind_vel, bearing, nearest_lat, nearest_lon, nearest_alt = model.getNewCoord(coord_new, dt_sec * GFSrate_sec)
                else:
                    lat_new, lon_new, x_wind_vel, y_wind_vel, bearing, nearest_lat, nearest_lon, nearest_alt = model.getNewCoordNoTime(coord_new, dt_sec * GFSrate_sec)
                gfscount += 1
                if lat_new == None:
                    logging.debug("Program likely came to the end of model data at epoch {t}")
                    break

                lat_deg.append(lat_new)
                lon_deg.append(lon_new)
                height_m.append(el_new)
                balloon_datetime.append(t)
                # store coordinates that change at a lower resolution
                coord_new['lat'] = lat_new
                coord_new['lon'] = lon_new
                el.append(el_new)

    else:
        logging.debug("Do simulation using modeled elevation and heights... ")
        height_parameters = dict(
            elevation_m=min_alt_m,
            burst=burst_balloon,
            burst_occurred=False,
            ceiling_m=ceiling_hgt_m,
            burst_m=burst_hgt_m,
            ascent_ms=balloon_ascent_ms,
            descent_ms=balloon_descent_ms,
            min_alt_m=min_alt_m,  # (m) Elevation
            dt_sec=dt_sec

        )

        for i in range(0, simulation_points):
            t = t + pd.Timedelta(hours=(1 / 3600 * dt_sec))
            hr_diff = t - first_sec
            timedelta_seconds = hr_diff.total_seconds()
            if burst_balloon == False:
                if timedelta_seconds >= stop_sec:
                    if first_burst == False:
                        logging.debug(f"Time to start balloon descent at {t}")
                        height_parameters["burst_occurred"] = True
            rad = radiation.Radiation()
            zen = rad.get_zenith(t, coord_new)

            height_parameters = get_height(height_parameters)
            el_new = height_parameters["elevation_m"]

            if height_parameters["burst_occurred"] == True and first_burst == False:
                first_burst = True
                coord_new['alt_m'] = el_new
                coord_new['timestamp'] = t
                if rc.use_time_dimension:
                    lat_new, lon_new, x_wind_vel, y_wind_vel, bearing, nearest_lat, nearest_lon, nearest_alt = model.getNewCoord(coord_new, dt_sec * GFSrate_sec)
                else:
                    lat_new, lon_new, x_wind_vel, y_wind_vel, bearing, nearest_lat, nearest_lon, nearest_alt = model.getNewCoordNoTime(coord_new, dt_sec * GFSrate_sec)
                gfscount += 1
                lat_deg.append(lat_new)
                lon_deg.append(lon_new)
                height_m.append(el_new)
                balloon_datetime.append(t)
                # store coordinates that change at a lower resolution
                coord_new['lat'] = lat_new
                coord_new['lon'] = lon_new
                logging.debug(f"Balloon burst at {t} after {i} iterations")

            if el_new < height_parameters["min_alt_m"]:
                logging.debug(f"Balloon hit grount at {t} after {i} iterations")
                break
            # Store coordinates that change every iteration

            coord_new['alt_m'] = el_new
            coord_new['timestamp'] = t

            if i % GFSrate_sec == 0:
                if rc.use_time_dimension:
                    lat_new, lon_new, x_wind_vel, y_wind_vel, bearing, nearest_lat, nearest_lon, nearest_alt = model.getNewCoord(coord_new, dt_sec * GFSrate_sec)
                else:
                    lat_new, lon_new, x_wind_vel, y_wind_vel, bearing, nearest_lat, nearest_lon, nearest_alt = model.getNewCoordNoTime(coord_new, dt_sec * GFSrate_sec)



                if lat_new == None:
                    break

                gfscount += 1
                lat_deg.append(lat_new)
                lon_deg.append(lon_new)
                height_m.append(el_new)
                balloon_datetime.append(t)
                # store coordinates that change at a lower resolution
                coord_new['lat'] = lat_new
                coord_new['lon'] = lon_new

            # if zen > math.radians(90): uncommented out these two lines in the original
            #    el_new -= 3

            el.append(el_new)
            ttt.append(t - pd.Timedelta(hours=GMT))  # Just for visualizing plot better

            if i % 1000 * (1 / dt_sec) == 0:
                logging.debug(str(t - pd.Timedelta(hours=GMT))  # Just for visualizing better
                  + " el " + str("{:.4f}".format(el_new))
                  + " zen " + str(math.degrees(zen))
                  )

                logging.debug(colored(
                ("U wind speed: " + str(x_wind_vel) + " V wind speed: " + str(y_wind_vel) + " Bearing: " + str(
                    bearing)),
                "yellow"))
                logging.debug(colored(("Lat: " + str(lat_new) + " Lon: " + str(lon_new)), "green"))
                logging.debug(colored(
                ("Nearest Lat:" + str(nearest_lat) + " Nearest Lon:" + str(nearest_lon) + " (" + str(
                    360 - nearest_lon) +
                 ") Nearest Alt: " + str(nearest_alt)), "cyan"))

            if rc.BALLOON_END_DATETIME != None and t >= rc.BALLOON_END_DATETIME:
                logging.debug(f"Ending simulation reached end of balloon flight time")
                logging.debug(t)
                break

    # Store the last point of data after loops finished
    lat_deg.append(coord_new['lat'])
    lon_deg.append(coord_new['lon'])
    height_m.append(coord_new['alt_m'])
    balloon_datetime.append(coord_new['timestamp'])
    max_hgt_m = max(height_m)

    # Outline Downloaded forecast subset:
    region = zip(*[
        (model.lat_bot_deg, model.lon_left_deg),
        (model.lat_top_deg, model.lon_left_deg),
        (model.lat_top_deg, model.lon_right_deg),
        (model.lat_bot_deg, model.lon_right_deg)
    ])
    logging.debug(f"{min(lat_deg)}, {max(lat_deg)}, {min(lon_deg)},  {max(lon_deg)}, {ttt[0]} {ttt[-1]}")
    logging.info(f"Balloon flew to a maximum height of {max_hgt_m}")

    # Google Plotting of Trajectory
    if rc.html_map:
        gmap1 = gmplot.GoogleMapPlotter(coord["lat"], coord["lon"], 8)
        gmap1.plot(lat_deg, lon_deg, 'red', edge_width=2.5)
        gmap1.marker(lat_deg[-1], lon_deg[-1], 'red', title="End Forecasted Position")
        gmap1.polygon(*region, color='cornflowerblue', edge_width=1, alpha=.2)
        gmap1.marker(lat_deg[0], lon_deg[0], 'black', title="Start Location")
        if rc.compare_to_telemetry == True:
            compare_to_telemetry(df, gmap1, rc.telemetry_source)

        my_datetime = datetime.now().strftime("%Y-%m-%d_%H%M")
        model_datetime = str(t.year) + "-" + str(t.month) + "-" + str(start.day)
        full_file_path = sim_prop.output_dir+"TrapezoidERA" + my_datetime+"_" + model_datetime + ".html"

        gmap1.draw(full_file_path)

# Plot balloon height versus time
    if rc.balloon_height:
        plt.style.use('seaborn-pastel')
        fig, ax = plt.subplots(figsize=(8, 8))
        height_ft = np.asarray(height_m) * 3.28084
        ax.plot(balloon_datetime, height_m)
        ax.set_ylabel('Altitude (m)')
        ax2 = ax.twinx()
        ax2.plot(balloon_datetime, height_ft)
        ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
        ax2.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
        ax2.set_ylabel('Altitude (ft)')
        ax.set_xlabel('Datetime (GMT)')
        plt.title('Altitude Profile for Solar Balloon')

        model_datetime = str(t.year) + "-" + str(t.month) + "-" + str(start.day)
        full_file_path = rc.OUTPUT_DIR+'HAB-height_'+model_datetime+".png"
        plt.savefig(full_file_path)
        #plt.show()

# Create geojson output
    if rc.geojson_file:
        number_of_balloon_runs = 1  # do simulation for single balloon in trapezoid
        my_datetime = datetime.now().strftime("%Y-%m-%d_%H%M")
        file_identifier = my_datetime + "_model-" + model.model_start_datetime.strftime("%Y-%m-%d_%H")
        output_file = geo.j_init(rc.OUTPUT_DIR, file_identifier)
        number_of_features = number_of_balloon_runs - 1
        geo.j_write(output_file, lat_deg, lon_deg, height_m, number_of_features)
        geo.j_close(output_file, number_of_balloon_runs)

# Create a csv file output
    if rc.csv_file:
        my_datetime = datetime.now().strftime("%Y-%m-%d_%H%M")
        file_identifier = my_datetime + "_ERA5-model-" + model.model_start_datetime.strftime("%Y-%m-%d_%H")
        column_headers = ["count","latitude","longitude","height","datetime"]
        output_file = csv.init(column_headers,rc.OUTPUT_DIR, file_identifier)
        csv.write(output_file, lat_deg, lon_deg, height_m, balloon_datetime)


    #t1_mod = time.time()
    #return t1_mod - t0_mod
    # These are placeholder values, I need to see whats going on with these varibales
    min_dist_datetime = balloon_datetime[0]
    min_dist_km = 50.0

    
    return lat_deg, lon_deg, height_m, balloon_datetime, min_dist_datetime, min_dist_km

if __name__ == "__main__":
    total_mod = main(flight_id=sim_prop.flight_id, rc=None)
    logging.info(f"Total run time: {total_mod}")