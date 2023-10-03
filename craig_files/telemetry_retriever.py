import numpy as np
import mathtime_utilities as mt
import glob
import pandas as pd
import datetime as dt
from datetime import datetime
import sys

import logging

logging.basicConfig(level=logging.WARNING)

def get_telemetry(df, r_config):
    """
    From telemetry data stored in the dataframe, get position data of the  GPS telemetry
    and return in a data frame

    Note, need to migrate more parsers from discovery_balloon_prop.py

    Parameters
    ----------
    df: dataframe
        Dataframe from parsing the telemetry data

    r_config : class obj
        Object that contains all the parameters of our simulation

    Returns
    -------
    df: dataframe
        dataframe containing position and height information of observed tracks from


    Examples
    --------

    """
    if r_config.telemetry_source == "NRL":
        df_nrl = pd.read_csv("./balloon_data/" + r_config.telemetry_file, index_col=False)
        df_nrl['Time'] = (pd.to_datetime(df['Time'].str.strip(), format='%Y-%m-%d %H:%M:%S'))
        # reduce data sample size to second level precision

        df_nrl['Time'] = df['Time'].values.astype('<M8[s]')
        df_nrl.drop_duplicates(subset="Time", inplace=True, keep='last')
        gps_hgt_m = df_nrl['Adafruit_Alt'].tolist()
        gps_lat = df_nrl['Adafruit_Lat'].tolist()
        gps_lon = df_nrl['Adafruit_Lon'].tolist()
        gps_lat_arr = np.array(gps_lat)
        gps_lon_arr = np.array(gps_lon)

        # NRL dataset has a lot of missing numbers which must be removed before plotting (Navy doesn't)
        nan_array = np.isnan(gps_lat_arr)
        not_nan_array = ~ nan_array
        gpslat_deg = gps_lat_arr[not_nan_array]
        nan_array = np.isnan(gps_lon_arr)
        not_nan_array = ~ nan_array
        gpslon_deg = gps_lon_arr[not_nan_array]
    elif r_config.telemetry_source == "NAVY":

        gps_hgt_m = df['gpsAltitude'].tolist()
        gpslat_deg = df['gpsLatitude'].tolist()
        gpslon_deg = df['gpsLongitude'].tolist()
        flight_phase = df['flightPhase'].tolist()
        logging.info(df)
        df['datetime'] = pd.to_datetime(df['Time'])
        dt = df['datetime'].tolist()

        dictionary = {
            'hgt_m': gps_hgt_m,
            'lat': gpslat_deg,
            'lon': gpslon_deg,
            'time' :dt,
            'phase' : flight_phase
        }
        df_gps = pd.DataFrame(dictionary)

    return df_gps


def get_new_height(h_m, epoch, timestamp):
    """Performs a 2 step linear interpolation to determine horizontal wind velocity.  First the altitude is interpolated between the two nearest
    .25 degree lat/lon areas.  Once altitude is matched, a second linear interpolation is performed with respect to time.
    :param coord: Coordinate of balloon
    :type coord: dict
    :returns: [x_wind_vel, y_wind_vel]
    :rtype: array
    """
    #Given an ordered array and a value, determines the index of the closest item contained in the array.

    closest_idx =  min(range(len(epoch)), key=lambda i: abs(epoch[i] - timestamp))
    height_m = h_m[closest_idx]

    return height_m

def load_dataset(telemetry_file, work_sheet_idx):
    logging.info(f"Read in the telemetry flight data from {telemetry_file}")
    filename = "./balloon_data/" + telemetry_file
    path = "C:/Users/cmotell/PycharmProjects/ehab22/earthhab2/balloon_data/"+telemetry_file
    path = filename

    logging.info(f"Search for pattern match at {path}")
    for filename in glob.glob(path):
        filename_lower = filename.lower()
        logging.info(f"Filename found matching our file pattern: {filename}")
        if "csv" in filename_lower:
            df = pd.read_csv("./balloon_data/" + telemetry_file, index_col=False)
        else:
            df = pd.read_excel(filename, sheet_name=work_sheet_idx)
            break

    try:
        logging.info(f"Data frame size: {df.size}")
    except:
        logging.info(f"./balloon_data/{telemetry_file} was not found, aborting program.")
        sys.exit(-1)

    #for (columnName, columnData) in df.iteritems():
    #    logging.info('Colunm Name : ', columnName)
    #    logging.info('Column Contents : ', columnData.values)

    return df

def parse_navy(df, r_config):
    '''
    :parm df: contents of Excel or CSV file
    :type: dataframe
    :param r_config:
    :return:
    '''
    df['Time'] = df['timestampEpoch'].apply(lambda d: dt.datetime.utcfromtimestamp(int(d)).strftime("%Y-%m-%d %H:%M:%S"))

    flight_phase = df["flightPhase"]
    phase = flight_phase.loc[0]
    if phase == "Ascent":
        id_start = 0
    else:
        for i in range(len(flight_phase)):
            phase = flight_phase.loc[i]
            logging.info(f"flight phase: {flight_phase.loc[i]}")
            if phase == "Ascent":
                id_start = i -1
                break

    id_end = len(flight_phase) - 1
    balloon_start_datetime = datetime.fromisoformat(df['Time'].loc[id_start])
    balloon_end_datetime = datetime.fromisoformat(df['Time'].loc[id_end])
    r_config.balloon_start_datetime = balloon_start_datetime
    r_config.balloon_end_datetime = balloon_end_datetime

    r_config.launch_altitude = df['gpsAltitude'].loc[id_start]
    r_config.launch_latitude = df['gpsLatitude'].loc[id_start]
    r_config.launch_longitude = df['gpsLongitude'].loc[id_start]

    r_config.last_altitude = df['gpsAltitude'].loc[id_end]
    r_config.last_latitude = df['gpsLatitude'].loc[id_end]
    r_config.last_longitude = df['gpsLongitude'].loc[id_end]

    r_config.start_coord = dict({
        "lat": df['gpsLatitude'].loc[id_start],
        "lon": df['gpsLongitude'].loc[id_start],
        "alt_m": df['gpsAltitude'].loc[id_start],
        "timestamp": balloon_start_datetime,  # timestamp
    })
    r_config.end_coord = dict({
        "lat": df['gpsLatitude'].loc[id_end],
        "lon": df['gpsLongitude'].loc[id_end],
        "alt_m": df['gpsAltitude'].loc[id_end],
        "end_datetime": balloon_end_datetime,  # timestamp
    })

    r_config.sim_time_sec = mt.datetime_epoch(df['Time'].loc[id_end]) - mt.datetime_epoch(df['Time'].loc[id_start])
    r_config.sim_time_iterations = int(r_config.sim_time_sec / r_config.dt_sec)
    return r_config


def parse_aprs(df):
    '''
    :parm df: contents of Excel or CSV file
    :type: dataframe
    :param r_config: sim_config object containing run-time parameters
    :return: r_config: populated
    '''
    time = "Time(UTC)"
    date = "Date"
    latitude = "Latitude"
    longitude = "Longitude"
    altitude = "Altitude(m)"
    df['Time'] = (pd.to_datetime(df[date] + " " + df[time].str.strip(), format='%m/%d/%Y %H:%M:%S'))
    df['Time2'] = df['Time'].apply(lambda d: (d - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's'))
    first_epoch = df['Time2'].iloc[0]
    df['epoch'] = df['Time2'].apply(lambda d: int(d - first_epoch))
    df['latitude'] = df[latitude]
    df['longitude'] = df[longitude]
    df['altitude'] = df[altitude]
    return df


def parse_aprsif(df):
    '''
    :parm df: contents of Excel or CSV file
    :type: dataframe
    :param r_config: sim_config object containing run-time parameters
    :return: r_config: populated
    '''

    #df['Time'] = (pd.to_datetime(df['time'], format='%Y-%m-%d %H:%M:%S'))
    df['Time'] = df['time']
    df['Time2'] = df['Time'].apply(lambda d: (d - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's'))
    first_epoch = df['Time2'].iloc[0]
    df['epoch'] = df['Time2'].apply(lambda d: int(d - first_epoch))
    df['latitude'] = df['lat']
    df['longitude'] = df['lng']
    #df['altitude'] = df['altitude']
    return df

def parse_aprs2(df):
    '''
    :parm df: contents of Excel or CSV file
    :type: dataframe
    :param r_config: sim_config object containing run-time parameters
    :return: r_config: populated
    '''

    df['Time'] = (pd.to_datetime(df['time'], format='%Y-%m-%d %H:%M:%S'))
    #df['Time'] = df['time']

    df['latitude'] = df['lat']
    df['longitude'] = df['lng']
    #df['altitude'] = df['altitude']
    return df


def get_flight_plan_by_telemetry(r_config):
    """
    Uses the telemetry data to set all our simulation parameters such as starting latitude,
    longitude, altitude and time. Normally, we set these in flight_plans, but here
    we are usually simulating a flight where we already have this information via a
    telemetry file or intend to control a balloon to match our predicted telemetry.

    Parameters
    ----------
    r_config : class obj
        Object that contains all the parameters of our simulation

    Returns
    -------
    df : dataframe
        dataframe containing telemetry information from either csv or excel file
    r_config: sim_config class
        See simulation class, contains all our simulation parameters

    Examples
    --------

    """
    df = load_dataset(r_config.telemetry_file, r_config.work_sheet_index)
    logging.info(f"Read in the telemetry flight data from {r_config.telemetry_file}")
    telem_src = r_config.telemetry_source.upper()
    if telem_src == 'NAVY':
        r_config = parse_navy(df, r_config)
    elif telem_src == 'APRSFI':
        df = parse_aprsif(df)
    elif telem_src == 'APRS':
        df = parse_aprs(df)
    elif telem_src == 'APRS2':
        df = parse_aprs2(df)
    else:
        logging.info(f"Not yet support telemetry format {r_config.telemetry_file}")
        sys.exit(1)
    # your telemetry may have useless data at the start or end, you can remove these data in the file or  as below
    id_start = 0 # if need you can delete unused rows at the start or end of your file or change these indices manually
    id_end = -1

    if r_config.telemetry_source != 'NAVY':
        balloon_start_datetime = df['Time'].iloc[id_start]
        balloon_end_datetime = df['Time'].iloc[id_end]
        r_config.balloon_start_datetime = balloon_start_datetime
        r_config.balloon_end_datetime = balloon_end_datetime
        r_config.min_alt_m = df['altitude'].iloc[id_start]

        r_config.launch_altitude = df['altitude'].iloc[id_start]
        r_config.launch_latitude = df['latitude'].iloc[id_start]
        r_config.launch_longitude = df['longitude'].iloc[id_start]

        r_config.last_altitude = df['altitude'].iloc[id_end]
        r_config.last_latitude = df['latitude'].iloc[id_end]
        r_config.last_longitude = df['longitude'].iloc[id_end]

        r_config.start_coord = dict({
            "lat": df['latitude'].iloc[id_start],
            "lon": df['longitude'].iloc[id_start],
            "alt_m": df['altitude'].iloc[id_start],
            "timestamp": balloon_start_datetime,  # timestamp
        })
        r_config.end_coord = dict({
            "lat": df['latitude'].iloc[id_end],
            "lon": df['longitude'].iloc[id_end],
            "alt_m": df['altitude'].iloc[id_end],
            "end_datetime": balloon_end_datetime,  # timestamp
        })

        r_config.sim_time_sec = mt.datetime_epoch(df['Time'].iloc[id_end]) - mt.datetime_epoch(df['Time'].iloc[id_start])
        r_config.sim_time_iterations = int(r_config.sim_time_sec / r_config.dt_sec)
        r_config.SIM_TIME_ITERATIONS = r_config.sim_time_iterations #when we get rid of upper case will remove these two lines
        r_config.SIM_TIME_HOURS = r_config.sim_time_sec/3600.0
        r_config.BALLOON_START_DATETIME = balloon_start_datetime
        r_config.BALLOON_END_DATETIME = balloon_end_datetime
        r_config.MIN_ALT_M = df['altitude'].iloc[id_start]

    return df, r_config

def set_estimated(lat, lon, hgt, dt):
    '''
    Simple constructor, that uses array information to set or create a dataframe
    containing the model estimated location. Its just like get_telemetry but in this
    case we are creating dataframe to be compared to the actual dataframe using model produced
    results

    returns:
        df: datafram
            With position data estimated from model data
    '''
    dictionary = {
        'hgt_m': hgt,
        'lat': lat,
        'lon': lon,
        'time': dt
    }
    df_est = pd.DataFrame(dictionary)
    return df_est