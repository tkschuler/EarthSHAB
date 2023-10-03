'''
Flight plans is a file where we store the details about a proposed simulation. Details include launch and optional
landing coordinates and altitudes, balloon ascent and descent rates. This file is written to be at least a temporary
database, where new data is just appended and accessed via the flight_id. This file is meant to be a data
repository for flights.

Programs will access data here, using the file "simulation_plan.py". Simulation_plan.py is modified for each new
simulation and serves as a temporary starting location. In "simulation_plan.py" you will see something like
the following:

flight_id = 'NRL20211217'  # this must be changed for different simulations
GMT = 7  # local time offset
compare_to_telemetry = True   # if set to true will compare to a telemetry file.
etc

This is where you enter the key to grab data from here.


Version 1.4 (15-May-2022)


Optional parameters can be replaced with None:

compare_to_telemetry  When set to False, how height calculations are set using model, when set true, use height from
                    telemetry to get balloon height, and also compare model produced trajectories to telemetry position.
                    Note, you can in fact, have use_telemetry set True and have use_telemetry_set False, but would
                    not advise these settings.
use_telemetry_set   When set to False or None, we will use settings below to (1) determine balloon start and end times,
                    starting latitude, longitude, altitude, and number of simulation points. But when set to True
                    all settings in the flight plan our ignored and all simulation is run off the telemetry data you 
                    provide.

use_time_dimension  When set true, we will assume that our netcdf has a time variable, if it doesn't program will fail.
                    When set false, we will be using monthly or yearly climatological data that has no time variable
                    or dimension.

model_datetime     Take the model datetime start from the netcdf_properties.py file used in saveNETCDF.py to create
                    the model. Note, you might think we just use netcdf_properties here and not bother, but netcdf_properties
                    is used and then the properties are changed but not save. Flight_plans archives what our parameters
                    were

data_source         Confirms your model file is GFS or ERA5. Note, the code will check to see the source from
                    model_file. So this is more of check to make sure you know you are doing.

model_file          GFS or ERA5 netcdf file containing the data.

optimize_by         Our current options are 'static', 'dynamic', 'distance', 'ensemble', 'psuedo_gradient'

telemetry_file      If using telemetry give enough of a path to make your file unique either csv, or excel

telemetry_source    Currently either, NAVY, NRL, or APRSFI, note different telemetry files have different formats.
                    You will almost surely need to add a new format to match yours and modify code in
                    telemetry_retriever.py to properly parse your data

work_sheet_index    If we are using telemetry from an excel file, you need to specify the work sheet number, if
                    only one work sheet then 0 (its 0 based), if not using Excel files, ignore. Excel files are
                    only used when we are comparing our predicted flight to known positions in the telemetry file

last_time           When set the simulation will stop at this time, otherwise (2) will end when model data runs out,
                    or (3) when balloon hits ground or (4) when exceeded by flight_hours.
                    Note, if we set flight_hours to a number that ends up before last_time, flight_hours will
                    be used over last_time or vice versa. If you want to do a full simulation and you know the time a balloon
                    hit ground, usually better to use last_time and set flight_hours = None. Note, it is required
                    that either last_time or flight_hours be set (both can't be None).

last_latitude       Can be used to set a target for say station keeping, or say the last latitude from a telemetry file. In the
                    first case, its the location you want to get to, and the later case it's the observed last location.

last_longitude      See above

last_altitude       When set (not None), flight will terminate when balloon falls below this height. Default minimum
                    altitude is given by the launch altitude. Say you release a balloon on a hill, then you will
                    want to set this altitude to something lower, like sea-level or 0 meters.

burst               When set to False balloon will not burst and its altitude will be controlled (see burst_hgt_m). When
                    set to True, balloon will burst at its burst_hgt_m

burst_hgt_m         When burst is set to True, balloon will not exceed burst height and will fall when it reached this height

max_hgt_m           When burst is set to False, balloon will level out and fly at this height.

strato_boundary_m   When this is not None we assume different ascent and descent rates above this value.
                    If balloon has constant ascent and descent rates set to None.

strato_ascent_m     Only used when strato_boundary_m is not None. When set assumes balloon ascent and descent rates
                    are different in the area above the boundary. Note, it is called strato_boundary_m but it 
                    can be set to any arbitrary boundary where the balloon behavior changes (meters per sec).

strato_descent_m    Note, we assume at most two different ascent/descent rates, we could extend this to n ascent/descent rates (meters per sec)

flight_hours        Normally set to None, but if you want the flight simulation to end after a specific number of hours
                    then use. Note, this option will soon be depreciated. We currently have too many ways of controlling
                    the number of simulation points, such as using balloon start and end times, or the start and end times
                    in the telemetry file or this variable flight_hours. Note, previously, this was referred to as sim_time_hours.
'''
flight = {
    'DEFAULT': {
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'Weather Balloon',
        'title': 'Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        'use_time_dimension': False,
        'compare_to_telemetry': False,
        'data_source': 'era',
        'model_file': 'westpac_jan_2003-2021_mean.nc',
        'model_datetime': '2022-01-01 12:00:00',
        'launch_time': '2022-01-01 11:30:00',  # GMT time we want the simulation to start out
        'last_time': "2022-01-02 14:30:00",  # this is pulled from the defualt for troubleshooing
        'launch_latitude': 15.134,
        'launch_longitude': 145.73,
        'launch_altitude': 74,
        'last_latitude': 17.42,
        'last_longitude': 163.619,
        'last_altitude': 16000,
        'flight_hours': 27,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 20000,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'optimize_by' : 'psuedo_gradient',
        'burst': False,
        'stratosphere_boundary_ave_m': 30000.0,
        'strato_boundary_m': 30000.0,
        'ascent_rate_ms': 2,
        'descent_rate_ms': -14.68,
        'tropo_ascent_rate_ms': 3.628,  # Added from Superpressue for dev
        'tropo_descent_rate_ms': -6.3489,  # Added from Superpressue for dev
        'strato_ascent_rate_ms': 0.2822,  # Added from Superpressue for dev
        'strato_descent_rate_ms': -0.2822  # Added from Superpressue for dev
    },
    'SAMPLE_STATIC': {
        'dt_sec': 3,
        'balloon_type': 'SPB',
        'title': 'Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        'use_time_dimension': False,
        'compare_to_telemetry': False,
        'data_source': 'gfs',
        'model_file': 'gfs_0p25_20220525_12.nc',
        'model_datetime': '2022-05-25 12:00:00',
        'launch_time': '2022-05-25 12:30:00',  # GMT time we want the simulation to start out
        'last_time': None,  # this is pulled from the defualt for troubleshooing
        'launch_latitude': 37.775,   # San Francisco CA
        'launch_longitude': -122.419,
        'launch_altitude': 100,
        'last_latitude': 39.53,  # last_latitude, etc, can also be used as target location, for Reno Nevada
        'last_longitude': -119.8143,
        'last_altitude': 1373,
        'flight_hours': 30,  # Optional, if set will either stop when simulation gets to last_time (if not set to None)
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 28042,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'optimize_by': 'static', # static, if fly whole path along a given altitude and compare results
        'burst': False,
        #'float_hgt_m': [14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000],
        'float_hgt_m': [14000, 15000, 16000],
        'leg_time_hr': 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'path_direction': 'forward',
        'tolerance_level': 0.9,
        'target_range': 250,
        'hover_minute_interval': 30,
        'strato_boundary_m': 14000.0,
        'ascent_rate_ms': 2,
        'descent_rate_ms': -14.68,
        'tropo_ascent_rate_ms': 3.628,  # Added from Superpressure for dev
        'tropo_descent_rate_ms': -6.3489,  # Added from Superpressrue for dev
        'strato_ascent_rate_ms': 0.2822,  # Added from Superpressure for dev
        'strato_descent_rate_ms': -0.2822  # Added from Superpressure for dev
    },
    'SAMPLE_DIRECTION': {
        'dt_sec': 3,
        'balloon_type': 'SPB',
        'title': 'Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        'use_time_dimension': False,
        'compare_to_telemetry': False,
        'data_source': 'gfs',
        'model_file': 'gfs_0p25_20220525_12.nc',
        'model_datetime': '2022-05-25 12:00:00',
        'launch_time': '2022-05-25 12:30:00',  # GMT time we want the simulation to start out
        'last_time': None,  # this is pulled from the defualt for troubleshooing
        'launch_latitude': 37.775,  # San Francisco CA
        'launch_longitude': -122.419,
        'launch_altitude': 100,
        'last_latitude': 39.53,  # last_latitude, etc, can also be used as target location, for Reno Nevada
        'last_longitude': -119.8143,
        'last_altitude': 1373,
        'flight_hours': 30,  # Optional, if set will either stop when simulation gets to last_time (if not set to None)
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 28042,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'optimize_by': 'direction',
        'burst': False,
        'float_hgt_m': [14000, 19000, 20000, 22000, 24000],
        'leg_time_hr': 1,  # (hours) Simulation time in hours (for trapezoid.py)
        'path_direction': 'forward',
        'tolerance_level': 0.9,
        'target_range': 250,
        'hover_minute_interval': 30,
        'strato_boundary_m': 14000.0,
        'ascent_rate_ms': 2,
        'descent_rate_ms': -14.68,
        'tropo_ascent_rate_ms': 3.628,  # Added from Superpressure for dev
        'tropo_descent_rate_ms': -6.3489,  # Added from Superpressrue for dev
        'strato_ascent_rate_ms': 0.2822,  # Added from Superpressure for dev
        'strato_descent_rate_ms': -0.2822  # Added from Superpressure for dev
    },
    'SAMPLE_DISTANCE': {
        'dt_sec' : 3,
        'balloon_type': 'SPB',
        'title': 'Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        'use_time_dimension': False,
        'compare_to_telemetry': False,
        'data_source': 'gfs',
        'model_file': 'gfs_0p25_20220525_12.nc',
        'model_datetime': '2022-05-25 12:00:00',
        'launch_time': '2022-05-25 12:30:00',  # GMT time we want the simulation to start out
        'last_time': None,  # this is pulled from the defualt for troubleshooing
        'launch_latitude': 37.775,   # San Francisco CA
        'launch_longitude': -122.419,
        'launch_altitude': 100,
        'last_latitude': 39.53,  # last_latitude, etc, can also be used as target location, for Reno Nevada
        'last_longitude': -119.8143,
        'last_altitude': 1373,
        'flight_hours': 30,  # Optional, if set will either stop when simulation gets to last_time (if not set to None)
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 28042,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'optimize_by': 'distance',
        'burst': False,
        'float_hgt_m' : [14000, 19000, 20000, 22000, 24000],
        'leg_time_hr' : 1,  # (hours) Simulation time in hours (for trapezoid.py)
        'path_direction' : 'forward',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'strato_boundary_m': 14000.0,
        'ascent_rate_ms': 2,
        'descent_rate_ms': -14.68,
        'tropo_ascent_rate_ms': 3.628,  # Added from Superpressure for dev
        'tropo_descent_rate_ms': -6.3489,  # Added from Superpressrue for dev
        'strato_ascent_rate_ms': 0.2822,  # Added from Superpressure for dev
        'strato_descent_rate_ms': -0.2822  # Added from Superpressure for dev
    },
    'SAMPLE_ENSEMBLE': {
        'balloon_type': 'SPB',
        'title': 'Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        'use_time_dimension': False,
        'compare_to_telemetry': False,
        'data_source': 'gfs',
        'model_file': 'gfs_0p25_20210913_12.nc',
        'model_datetime': '2021-09-13 12:00:00',
        'launch_time': '2021-09-13 11:32:00',  # GMT time we want the simulation to start out
        'last_time': None,  # this is pulled from the defualt for troubleshooing
        'launch_latitude': 37.775,
        'launch_longitude': -122.419,
        'launch_altitude': 100,
        'last_latitude': 40.7608,  # last_latitude, etc, can also be used as target location
        'last_longitude': -111.8910,
        'last_altitude': 1673,
        'flight_hours': 30,  # Optional, if set will either stop when simulation gets to last_time (if not set to None)
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 28042,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'optimize_by': 'ensemble',
        'burst': False,
        'float_hgt_m': [14000, 18000, 22000, 26000],
        'leg_time_hr': 6,  # (hours) Simulation time in hours (for trapezoid.py)
        'path_direction': 'forward',
        'dt_sec' : 6,
        'tolerance_level': 0.9,
        'target_range': 250,
        'hover_minute_interval': 30,
        'strato_boundary_m': 14000.0,
        'ascent_rate_ms': 2,
        'descent_rate_ms': -14.68,
        'tropo_ascent_rate_ms': 3.628,  # Added from Superpressure for dev
        'tropo_descent_rate_ms': -6.3489,  # Added from Superpressrue for dev
        'strato_ascent_rate_ms': 0.2822,  # Added from Superpressure for dev
        'strato_descent_rate_ms': -0.2822  # Added from Superpressure for dev
    },
    'DEFAULT_GFS': { # for some reason this has the era datafile as the data source? why?
        'data_source': 'gfs',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'Weather Balloon',
        'title': 'Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        'use_time_dimension': False,
        'compare_to_telemetry': False,
        'model_file': 'gfs_0p25_20220301_12.nc',
        'model_datetime': '2022-03-01 12:00:00',
        'launch_time': '2022-03-01 11:30:00',  # GMT time we want the simulation to start out
        'last_time': "2022-03-02 14:30:00",  # this is pulled from the defualt for troubleshooing
        'launch_latitude': 15.134,
        'launch_longitude': 145.73,
        'launch_altitude': 74,
        'last_latitude': 17.42,
        'last_longitude': 163.619,
        'last_altitude': 16000,
        'flight_hours': 27,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 20000,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'optimize_by': 'psuedo_gradient',
        'burst': False,
        'stratosphere_boundary_ave_m': 30000.0,
        'strato_boundary_m': 30000.0,
        'ascent_rate_ms': 2,
        'descent_rate_ms': -14.68,
        'tropo_ascent_rate_ms': 3.628,  # Added from Superpressue for dev
        'tropo_descent_rate_ms': -6.3489,  # Added from Superpressue for dev
        'strato_ascent_rate_ms': 0.2822,  # Added from Superpressue for dev
        'strato_descent_rate_ms': -0.2822  # Added from Superpressue for dev
    },
    'SAMPLE': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'Weather Balloon',
        'title': 'Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        'use_time_dimension': False,
        'compare_to_telemetry': False,
        'model_file': 'westpac_jan_2003-2021_mean.nc',
        'model_datetime': '2022-02-17 12:00:00',
        'launch_time': '2020-01-17 14:01:32',  # GMT time we want the simulation to start out
        'last_time': None,  # we don't know how long this
        'launch_latitude': 15.18,
        'launch_longitude': 145.74,
        'launch_altitude': 0,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': 4,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 20000,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': 4.68,
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
        'tropo_ascent_rate_ms': 3.628,  # Added from Superpressue for dev
        'tropo_descent_rate_ms': -6.3489  # Added from Superpressue for dev
    },
    'GFS_SAMPLE': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'Weather Balloon',
        'title': 'GFS Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        'use_time_dimension': False,
        'compare_to_telemetry': False,
        'model_file': 'gfs_0p25_20220301_12.nc',
        'model_datetime': '2022-03-01 12:00:00',
        'launch_time': '2022-03-01 12:30:00',  # GMT time we want the simulation to start out
        'last_time': '2022-03-02 15:30:00',  # we don't know how long this
        'launch_latitude': 15.18,
        'launch_longitude': 145.74,
        'launch_altitude': 0,
        'last_latitude': 17.42,
        'last_longitude': 163.619,
        'last_altitude': None,
        'flight_hours': 27,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 20000,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'optimize_by': 'psuedo_gradient',
        'burst': False,
        'strato_boundary_m': 30000.0,
        'ascent_rate_ms': 4.68,
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
        'tropo_ascent_rate_ms': 3.628,  # Added from Superpressue for dev
        'tropo_descent_rate_ms': -6.3489  # Added from Superpressue for dev
    },
    'SAMPLE_ORIG': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'Weather Balloon',
        'title': 'Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        # 'use_time_dimension': True,
        'use_time_dimension': False,  # Mod for t/s
        'compare_to_telemetry': False,
        # 'model_file': 'gfs_0p25_20220217_12.nc',
        'model_file': 'era.nc',  # This was changed to test the file.
        'model_datetime': '2022-02-17 12:00:00',
        'launch_time': '2022-02-17 14:01:32',  # GMT time we want the simulation to start out
        'last_time': None,  # we don't know how long this
        'launch_latitude': 39.2735,
        'launch_longitude': -80.736,
        'launch_altitude': 579.73,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': 4,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': 20068.0,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': True,
        'strato_boundary_m': None,
        'ascent_rate_ms': 4.68,
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'SAMPLE0': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'Weather Balloon',
        'title': 'Sample Flight Plan',
        'telemetry_file': None,
        'telemetry_source': None,
        'use_telemetry_set': False,
        # 'use_time_dimension': True,
        'use_time_dimension': False,  # Mod for t/s
        'compare_to_telemetry': False,
        # 'model_file': 'gfs_0p25_20220217_12.nc',
        'model_file': 'era.nc',  # This was changed to test the file.
        'model_datetime': '2022-02-17 12:00:00',
        'launch_time': '2022-02-17 12:01:32',  # GMT time we want the simulation to start out
        'last_time': None,  # we don't know how long this
        'launch_latitude': 39.2735,
        'launch_longitude': -80.736,
        'launch_altitude': 579.73,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': 4,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': 20068.0,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 20500,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': 4.68,
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'SAMPLE2': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        # Duruisseau et al: Assessment of the ERA-Interim Winds Using High-Altitude Stratospheric Balloons JAS 2017
        'balloon_type': 'Weather Balloon',
        'title': 'Sample of Using Telemetry to get Balloon Height',
        'telemetry_file': 'aprsfi_export_W0OAH*.xlsx',
        'work_sheet_index': 0,
        'telemetry_source': 'APRSFI',
        'compare_to_telemetry': True,
        'use_telemetry_set': True,
        # 'use_time_dimension': True,
        'use_time_dimension': False,  # Mod for t/s
        # 'model_file': 'gfs_0p25_20210512_12.nc',
        'model_file': 'era.nc',  # This was changed to test the file.
        'model_datetime': '2021-05-12 12:00:00',
        'launch_time': None,  # GMT time we want the simulation to start out
        'last_time': None,  # we don't know how long this
        'launch_latitude': None,
        'launch_longitude': None,
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': 4.68,
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'shab10v': {
        'notes': 'Flight plan used to test against shab10v-aprs.csv using telemetry and era5 model data',
        'data_source': 'era',
        'tolerance_level': 0.9,  # na
        'target_range': 250,  # na
        'hover_minute_interval': 30,  # ignore
        'path_direction': 'forward',  # ignore
        'leg_time_hr': 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m': [None],
        'dt_sec': 3,
        # Duruisseau et al: Assessment of the ERA-Interim Winds Using High-Altitude Stratospheric Balloons JAS 2017
        'balloon_type': 'Weather Balloon',
        'title': 'Sample Flight Plan for Readme',
        'telemetry_file': 'SHAB10V-APRS.CSV',
        'work_sheet_index': 0,
        'telemetry_source': 'APRS2',
        'compare_to_telemetry': True, # compare against a telemetry file the position
        'use_telemetry_set' : True, #when set to true use telemetry to get simulation time, false use  flight_hours or balloon start and end times
        'use_time_dimension': True,
        'model_file': 'shab10_era_2022-04-09to2022-04-10.nc',
        'model_datetime': '2022-04-09 00:00:00',  # ignore for era5 data
        'launch_time': None,  # GMT time we want the simulation to start out
        'last_time': None,  # we don't know how long this
        'launch_latitude': None,
        'launch_longitude': None,
        'launch_altitude': None,  # units are meters
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'optimize_by': None,
        'burst': None,  # balloon will pop
        'strato_boundary_m': None,
        'ascent_rate_ms': None,  # units are meters/second
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
        'tropo_ascent_rate_ms': None,
        'tropo_descent_rate_ms': None
    },
    'era_sample_2022-06-01': {
        'notes': 'Flight plan used to fly a weather balloon that burst at a given altitude, tested 6/1/2022 with trapezoid_era5',
        'data_source': 'era',
        'tolerance_level' : 0.9, #na
        'target_range' : 250, #na
        'hover_minute_interval' : 30, #ignore
        'path_direction' : 'forward', #ignore
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [17000],
        'dt_sec' : 3,
        # Duruisseau et al: Assessment of the ERA-Interim Winds Using High-Altitude Stratospheric Balloons JAS 2017
        'balloon_type': 'Weather Balloon',
        'title': 'Sample Flight Plan for Readme',
        'telemetry_file': None,
        'work_sheet_index': 0,
        'telemetry_source': None,
        'compare_to_telemetry': False,
        'use_telemetry_set': False,
        'use_time_dimension': True,
        'model_file': 'kauai_2021_June_01.nc',
        'model_datetime': '2021-06-01 12:00:00', #ignore for era5 data
        'launch_time': '2021-06-01 18:17:35',  # GMT time we want the simulation to start out
        'last_time': None,  # we don't know how long this
        'launch_latitude': 21.99,
        'launch_longitude': -159.34,
        'launch_altitude': 74.0,  # units are meters
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': 16,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': 28000,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'optimize_by' : None,
        'burst':True,  # balloon will pop
        'strato_boundary_m': None,
        'ascent_rate_ms': 5.00,  # units are meters/second
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
        'tropo_ascent_rate_ms': None,
        'tropo_descent_rate_ms' : None
    },
    'era_sample1': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        # Duruisseau et al: Assessment of the ERA-Interim Winds Using High-Altitude Stratospheric Balloons JAS 2017
        'balloon_type': 'Weather Balloon',
        'title': 'Sample Flight Plan for Readme',
        'telemetry_file': None,
        'work_sheet_index': 0,
        'telemetry_source': None,
        'compare_to_telemetry': False,
        'use_telemetry_set': False,
        'use_time_dimension': True,
        'model_file': 'kauai_2021_jun.nc',
        'launch_time': '2021-06-01 18:17:35',  # GMT time we want the simulation to start out
        'last_time': None,  # we don't know how long this
        'launch_latitude': 21.99,
        'launch_longitude': -159.34,
        'launch_altitude': 74.0,  # units are meters
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': 6,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': 35000,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 20500,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,  # balloon will pop
        'strato_boundary_m': None,
        'ascent_rate_ms': 5.00,  # units are meters/second
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'era_sample2': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        # Duruisseau et al: Assessment of the ERA-Interim Winds Using High-Altitude Stratospheric Balloons JAS 2017
        'balloon_type': 'Weather Balloon',
        'title': 'Sample Flight Plan for Readme using Telemetry',
        'work_sheet_index': 0,
        'telemetry_source': 'APRSFI',
        'telemetry_file': "aprsfi_201015-201016.xlsx",
        # must be in the ./balloon_data with respect to where this program is
        'compare_to_telemetry': True,
        'use_telemetry_set': True,
        'use_time_dimension': True,
        'model_file': 'havasu_2020_oct.nc',
        'launch_time': None,  # GMT time we want the simulation to start out
        'last_time': None,  # we don't know how long this
        'launch_latitude': None,
        'launch_longitude': None,
        'launch_altitude': None,  # units are meters
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': 6,  # Optional, if set will either stop when simulation gets to last_time (if not set to None
        'burst_hgt_m': 35000,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'ceiling_hgt_m': 20500,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,  # balloon will pop
        'strato_boundary_m': None,
        'ascent_rate_ms': 5.00,  # units are meters/second
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    '10_1_20_APRS_RAW': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'Solar',
        'title': 'Solar Balloon from Tucson Tech Park',
        'telemetry_file': "10_1_20_APRS_RAW.csv",
        'telemetry_source': 'APRS',
        'era_file': '20200919-22.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'NRL20210920': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'NRL',
        'title': 'West Virginia Launch',
        'telemetry_file': "gps_analysis_2.csv",
        'telemetry_source': 'NRL',
        'use_time_dimension': True,
        'use_telemetry_set': True,
        'model_file': 'era_20210920_160000.nc',
        'launch_time': '2021-09-20 16:17:35',
        'last_time': '2021-09-20 18:32:15',
        'launch_latitude': 39.0461,
        'launch_longitude': -78.7574,  # -negative degrees are degrees west
        'launch_altitude': 406.,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': 3,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': 22227.4,  # when burst_hgt is not None, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'strato_boundary_m': None,
        'ascent_rate_ms': 4.68,
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': 2.67,
        'strato_descent_rate_ms': -14.68,
    },
    'NRL20210920_Original': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'NRL',
        'title': 'West Virginia Launch',
        'telemetry_file': "gps_analysis_2.csv",
        'telemetry_source': 'NRL',
        'model_file': 'era_20210920_160000.nc',
        'launch_time': '2021-09-20 16:17:35',
        'last_time': '2021-09-20 18:32:15',
        'launch_latitude': 39.0461,
        'launch_longitude': -78.7574,  # -negative degrees are degrees west
        'launch_altitude': 406.,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': 22227.4,  # when burst_hgt is not None, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'strato_boundary_m': 5518.0,
        'ascent_rate_ms': 4.68,
        'descent_rate_ms': -14.68,
        'strato_ascent_rate_ms': 2.67,
        'strato_descent_rate_ms': -14.68,
    },
    'NRL20211217': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'NRL',
        'title': 'West Virginia Launch',
        'telemetry_file': "nrl_20211223_150740.csv",
        'telemetry_source': 'NRL',
        'model_file': '2021-12-21_17--23.nc',
        'launch_time': '2021-12-21 17:06:22',
        'last_time': '2021-12-21 21:40:20',
        'launch_latitude': 39.27350,
        'launch_longitude': -80.73600,  # -negative degrees are degrees west
        'launch_altitude': 579.73,
        'last_latitude': 40.76500,
        'last_longitude': -74.90400,
        'last_altitude': 470.92,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': 33241.49,  # when burst_hgt is not None, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'strato_boundary_m': None,
        'ascent_rate_ms': 2.4798,
        'descent_rate_ms': -16.6603,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0458': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'TW VS-20 CFT-1 Yankee-4',
        'telemetry_file': "HBAL0458*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': '20200919-22.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0448': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0448  -- TW VS-20 CFT-1 Yankee-4',
        'telemetry_file': "HBAL0448*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': '2020_09_18-19.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0445': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0445',
        'telemetry_file': "HBAL0445*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': '2020_09_21-26.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0446': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0446',
        'telemetry_file': "HBAL0446*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': '2020_09_14-15.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0447': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0447',
        'telemetry_file': "HBAL0447*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': '2020_09_18-19b.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0449': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0449',
        'telemetry_file': "HBAL0449*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': '2020_09_18-20.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0450': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0450',
        'telemetry_file': "HBAL0450*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': 'westpac_sep_full.nc',
        'use_telemetry_set': True,
        'use_time_dimension': False,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0450_PLAN': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0450_Plan',
        'telemetry_file': "None",
        'telemetry_source': 'NAVY',
        'era_file': 'westpac_sep_full.nc',
        'use_telemetry_set': True,
        'use_time_dimension': False,
        # optional parameters follow
        'launch_time': '2022-06-21 17:06:22',
        'last_time': '2022-06-24 21:40:20',
        'launch_latitude': 13.44,
        'launch_longitude': 144.0,  # -negative degrees are degrees west
        'launch_altitude': 0,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': 3.71,
        'descent_rate_ms': -8.65,
        'flight_levels': 17500,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0451': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0451',
        'telemetry_file': "HBAL0451*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': '2020_09_16-20.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0452': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0452',
        'telemetry_file': "HBAL0452*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': '2020_09_17-17.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
    'HBAL0448': {
        'data_source': 'era',
        'tolerance_level' : 0.9,
        'target_range' : 250,
        'hover_minute_interval' : 30,
        'path_direction' : 'forward',
        'leg_time_hr' : 3,  # (hours) Simulation time in hours (for trapezoid.py)
        'float_hgt_m' : [20000],
        'dt_sec' : 3,
        'balloon_type': 'TH-208',
        'title': 'HBAL0448',
        'telemetry_file': "HBAL0448*.xlsx",
        'telemetry_source': 'NAVY',
        'era_file': 'era_2020_09_19_00-06.nc',
        'use_telemetry_set': True,
        'use_time_dimension': True,
        # optional parameters follow
        'launch_time': None,
        'last_time': None,
        'launch_latitude': None,
        'launch_longitude': None,  # -negative degrees are degrees west
        'launch_altitude': None,
        'last_latitude': None,
        'last_longitude': None,
        'last_altitude': None,
        'flight_hours': None,  # Set to None if want to either use all model data or last_time
        'burst_hgt_m': None,  # when burst_hgt is greater than 0, implies balloon gets to that height in pops
        'max_hgt_m': None,  # when ceiling hgt is greater than 0, implies balloon gets to that height in floats
        'burst': False,
        'strato_boundary_m': None,
        'ascent_rate_ms': None,
        'descent_rate_ms': None,
        'strato_ascent_rate_ms': None,
        'strato_descent_rate_ms': None,
    },
}
