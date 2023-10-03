from math import radians, degrees, sin, cos, asin, acos, atan2, pi, sqrt
import datetime as dt
from datetime import datetime
from pytz import timezone

import logging

logging.basicConfig(level=logging.WARNING)

def great_circle(lon1, lat1, lon2, lat2):
    if abs(lat2 - lat1) < .01 and abs(lon2 - lon1) < .01:
        return 0.0
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    try:
        x = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
        gc = 6371.0 * x
        return gc
    except:
        logging.info("Problem doing calculation, with the following values:")
        logging.info(f"lat1={lat1}, lon1={lon1}, lat2={lat2}, lon2={lon2}")
        return None

def move_point(latitude, longitude, distance_metres, bearing):
    bearing_rad = radians(bearing)
    lat_rad = radians(latitude)
    lon_rad = radians(longitude)
    earth_radius_metres = 6371000
    dist_frac = distance_metres / earth_radius_metres

    latitude_result = asin(sin(lat_rad) * cos(dist_frac) + cos(lat_rad) * sin(dist_frac) * cos(bearing_rad));
    a = atan2(sin(bearing_rad) * sin(dist_frac) * cos(lat_rad), cos(dist_frac) - sin(lat_rad) * sin(latitude_result));
    longitude_result = (lon_rad + a + 3 * pi) % (2 * pi) - pi
    return latitude_result, longitude_result

def isValidNumber(num):
    if float('-inf') < float(num) < float('inf'):
        return True
    else:
        return False

'''
By default, Python, when use timestamp for some unexplained reason will give an opoch with an offset
for local time. If you run the code in Hawaii the epoch time will be different from say 3 hours from
PST or 3 X 60 X 60. 

So to force python to give you true GMT you need to take time wasting extra steps. The code is confirmed
using the web page: https://www.unixtimestamp.com/index.php
'''
def datetime_epoch(timedate_obj):
    gmt = timezone('GMT')
    try:
        dt = timedate_obj.strftime("%Y-%m-%d %H:%M:%S")
        naive_ts = datetime.strptime(dt, "%Y-%m-%d %H:%M:%S")
    except:
        naive_ts = datetime.strptime(timedate_obj, "%Y-%m-%d %H:%M:%S")

    local_ts = gmt.localize(naive_ts)
    epoch_ts = local_ts.timestamp()

    return int(epoch_ts)

def str_to_datetime(time_str):
    '''
    Convert string to datetime object
    :param time_str:
    :return:
    '''
    format = '%Y-%m-%d %H:%M:%S'  # The format
    try:
        datetime_obj = datetime.strptime(time_str, format)
    except:
        logging.info(f"Unable to convert string {time_str} to datetime object")

    return datetime_obj

def num_to_month(month):
    return {
            1: 'January',
            2: 'February',
            3: 'March',
            4: 'April',
            5: 'May',
            6:'June',
            7:'July',
            8:'August',
            9:'September',
            10:'October',
            11:'November',
            12:'December'
    }[month]
