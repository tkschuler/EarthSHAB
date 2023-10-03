'''
This is to build a common config between all of the methods.
'''
# One thing I want to do is to pass the 
# sim_prop.flight_id 
# to the clase so that we can load based on that.
# Otherwise the config can be adjusted by individual params
# in the respective code 
import sys
from datetime import datetime

import properties.simulation_plan as sim_prop  # getting rid of simulation_properties
import properties.balloon_properties as balloon_prop
import properties.mission_properties as mission_prop
import properties.netcdf_properties as net_prop
import properties.demo_properties as demo_prop
import properties.flight_plans as flt_prop
import properties.era5_properties as era_prop
import logging
import mathtime_utilities as mt

logging.basicConfig(filename='load_config.log',
                    level=logging.WARNING,
                    format='%(asctime)s | %(levelname)s | %(message)s',
                    datefmt="%Y-%m-%d %H:%M")


class sim_config:
    def __init__(self, savenetcdf=False, flight_id=sim_prop.flight_id):
        if savenetcdf:
            # For saveNETCDF only, creates a superset of data to be used in further simulations
            self.snc_lat_range = net_prop.netcdf['lat_range']
            self.snc_lon_range = net_prop.netcdf['lon_range']
            self.snc_hours3 = net_prop.netcdf['hours3']
            self.snc_hourstamp = net_prop.netcdf['hourstamp']
            self.snc_res = net_prop.netcdf['res']
            self.snc_nc_start = net_prop.netcdf['nc_start_datetime']

            self.snc_nc_file = net_prop.netcdf['nc_file']
            self.snc_start_coord = net_prop.netcdf['start_coord']
        else:
            # For the flight plan page
            # This will be added later perhaps
            # for now it seems to have to be dynamic based ont a front end
            # selection.
            self.flight_id = flight_id

            # From Trapezoid
            self.model_datetime = flt_prop.flight[flight_id]["model_datetime"]

            self.forecast_dir = sim_prop.forecast_dir
            self.output_dir = sim_prop.output_dir
            self.GMT = sim_prop.GMT
            self.MIN_ALT_M = flt_prop.flight[flight_id]['launch_altitude']
            self.DT_SEC = flt_prop.flight[flight_id][
                "dt_sec"]  # important time for numerical integration of the physics models
            self.dt_sec  = self.DT_SEC
            self.SIM_TIME_HOURS = flt_prop.flight[flight_id]["flight_hours"]
            self.BALLOON_NAME = flt_prop.flight[flight_id]["balloon_type"]
            self.ASCENT_RATE_MS = flt_prop.flight[flight_id]['ascent_rate_ms']
            self.GFS_RATE_SEC = sim_prop.GFSrate_sec
            self.OUTPUT_DIR = sim_prop.output_dir
            self.FLOAT_HGT_M = flt_prop.flight[flight_id]['float_hgt_m']


            try:
                self.LEG_TIME_HOURS = flt_prop.flight[flight_id]["leg_time_hr"]
            except:
                self.LEG_TIME_HOURS = None
                self.LEG_TIME_ITERATIONS = None
                self.LEG_TIME_SEGMENTS = None
            if self.LEG_TIME_HOURS == None:
                self.LEG_TIME_ITERATIONS = None
            else:
                self.LEG_TIME_ITERATIONS = self.LEG_TIME_HOURS * int(3600 * (1 / self.DT_SEC))  # number of simulation points
            try:
                self.LEG_TIME_SEGMENTS = self.SIM_TIME_HOURS / self.LEG_TIME_HOURS
            except:
                self.LEG_TIME_SEGMENTS = None

            self.calc_methods = sim_prop.calc_methods
            self.optimize_by = flt_prop.flight[flight_id]['optimize_by']
            self.OUTPUT = sim_prop.output
            self.STRATO_BOUNDARY = flt_prop.flight[flight_id]['strato_boundary_m']
            self.TROPO_ASCENT_MS = flt_prop.flight[flight_id]['tropo_ascent_rate_ms']
            self.TROPO_DESCENT_MS = flt_prop.flight[flight_id]['tropo_descent_rate_ms']
            self.STRATO_ASCENT_MS = flt_prop.flight[flight_id]['tropo_ascent_rate_ms']
            self.STRATO_DESCENT_MS = flt_prop.flight[flight_id]['tropo_descent_rate_ms']
            self.CEILING_M = flt_prop.flight[flight_id]['ceiling_hgt_m']
            self.PATH_DIRECTION = flt_prop.flight[flight_id]['path_direction']

            # This is done for date changes
            self.SIM_GFS = self.model_datetime
            self.NET_GFS = self.model_datetime

            # New stuff for the hover stage
            self.TOLERANCE_LEVEL = flt_prop.flight[flight_id][
            'tolerance_level']  # this will be a percentage of target range and adjustable
            self.TARGET_RANGE = flt_prop.flight[flight_id]['target_range']  # This will be pulled from the mission params
            self.HOVER_BOUNDRY = self.TOLERANCE_LEVEL * self.TARGET_RANGE  # This is calculated based on the above
            self.HOVER_TIME_MINUTES = flt_prop.flight[flight_id]['hover_minute_interval']
            self.HOVER_ITERATIONS = self.HOVER_TIME_MINUTES * int(60 * (1 / self.DT_SEC))
            # For the Graph
            self.ALTITUDE_DIR = sim_prop.alt_graph_dir

            # For the datasource
            self.data_source = flt_prop.flight[flight_id]["data_source"]
            if self.data_source == None:
                logging.error("You must define data_source as either 'era', 'gfs', or 'navgem'")
                sys.exit(2)

            # new era options
            self.telemetry_dir = sim_prop.telemetry_dir
            # These are outputs (from simulation_plan)
            self.html_map = sim_prop.output['html']  # this is a htmlmap output
            self.leaflet_map = sim_prop.output['leaflet_map']  # this is a leaflet map output
            self.geojson_file = sim_prop.output['geojson']  # this is a geojson output
            self.csv_file = sim_prop.output['csv']  # this is a csv output
            self.cartopy_map = sim_prop.output['cartopy_map']  # this is a cartopy output
            self.balloon_height = sim_prop.output['balloon_height']  # this is a blloon hieght output

            self.compare_to_telemetry = flt_prop.flight[flight_id]["compare_to_telemetry"]
            self.telemetry_file = None
            self.telemetry_source = None
            self.use_telemetry_set = flt_prop.flight[flight_id]["use_telemetry_set"]
            self.use_time_dimension = flt_prop.flight[flight_id]["use_time_dimension"]
            if self.compare_to_telemetry:
                self.telemetry_file = flt_prop.flight[flight_id]["telemetry_file"]
                self.telemetry_source = flt_prop.flight[flight_id]["telemetry_source"]
                self.work_sheet_index = flt_prop.flight[flight_id]["work_sheet_index"]
            else:
                self.telemetry_file = None
                self.telemetry_source = None
                self.use_telemetry_set = False

            self.balloon_type = flt_prop.flight[flight_id]["balloon_type"]
            self.start = flt_prop.flight[flight_id]['launch_time']
            self.end = flt_prop.flight[flight_id]['last_time']
            if self.start == None:
                self.BALLOON_START_DATETIME = None
            else:
                self.BALLOON_START_DATETIME = datetime.fromisoformat(self.start)

            if self.end != None:
                self.BALLOON_END_DATETIME = datetime.fromisoformat(self.end)
            else:
                self.BALLOON_END_DATETIME = None

            self.balloon_end_datetime = self.BALLOON_END_DATETIME
            self.input_file = flt_prop.flight[flight_id]["model_file"]

            self.launch_latitude = flt_prop.flight[flight_id]["launch_latitude"]
            self.launch_longitude = flt_prop.flight[flight_id]['launch_longitude']
            self.launch_altitude = flt_prop.flight[flight_id]['launch_altitude']

            self.last_latitude = flt_prop.flight[flight_id]["last_latitude"]
            self.last_longitude = flt_prop.flight[flight_id]['last_longitude']
            self.last_altitude = flt_prop.flight[flight_id]['last_altitude']
            # min_alt_m = launch_altitude # matched

            # sim_time_hrs = flt_prop.flight[flight_id]["flight_hours"] #Matched
            self.sim_time_hours = flt_prop.flight[flight_id]["flight_hours"]
            self.SIM_TIME_HOURS = self.sim_time_hours
            if self.sim_time_hours == None and self.compare_to_telemetry == False:
                if self.balloon_end_datetime == None:
                    logging.error(
                    "In your flight time you must either set the number of hours to simulate via 'flight_hours' or ")
                    logging.error("you must set the 'last_time' or being using telemetry else unknown length of simulation")
                    sys.exit(1)
            elif self.compare_to_telemetry == True:
                logging.info("Simulation time is taken from the first and last points of your telemetry file")

            # ascent_rate_ms = flt_prop.flight[flight_id]['ascent_rate_ms'] # matched
            self.descent_rate_ms = flt_prop.flight[flight_id]['descent_rate_ms']
            self.strato_boundary_m = flt_prop.flight[flight_id]['strato_boundary_m']
            if self.strato_boundary_m == None:
                self.strato_ascent_rate_ms = None
                self.strato_descent_rate_ms = None
            else:
                self.strato_ascent_rate_ms = flt_prop.flight[flight_id]['strato_ascent_rate_ms']
                self.strato_descent_rate_ms = flt_prop.flight[flight_id]['strato_descent_rate_ms']

            self.burst_hgt_m = flt_prop.flight[flight_id]['burst_hgt_m']
            if self.burst_hgt_m == None:
                self.burst_balloon = False
            else:
                self.burst_balloon = flt_prop.flight[flight_id]['burst']
            if self.burst_balloon == True and self.burst_hgt_m == None:
                print("Your flight plan you set 'burst=True', then you must have a burst height but you set to None")
                logging.error("Your flight plan you set 'burst=True', then you must have a burst height but you set to None")
                sys.exit(2)
            elif self.burst_balloon == False and self.burst_hgt_m != None:
                print(f"Your flight plan you set 'burst=False', but you gave a burst height {self.burst_hgt_m}, will override and set Burst=True")
                logging.warning(f"Your flight plan you set 'burst=False', but you gave a burst height {self.burst_hgt_m}, will override and set Burst=True")
                self.burst_balloon = True

            self.ceiling_hgt_m = flt_prop.flight[flight_id]['ceiling_hgt_m']
            if self.burst_balloon == True and self.ceiling_hgt_m != None:
                logging.debug(
                "Warning you say you want the balloon to burst, but you also gave us a ceiling or float height, results may not be as expected")

            # hrs2sec = 3600

            # start_timestamp = mt.datetime_epoch(balloon_start_datetime)
            self.start_timestamp = self.BALLOON_START_DATETIME
            #if self.BALLOON_END_DATETIME != None:
            self.end_timestamp = self.BALLOON_END_DATETIME

            self.start_coord = dict({
                "lat": flt_prop.flight[flight_id]['launch_latitude'],
                "lon": flt_prop.flight[flight_id]['launch_longitude'],
                "alt_m": flt_prop.flight[flight_id]['launch_altitude'],
                "timestamp": self.BALLOON_START_DATETIME,  # timestamp
            })
            self.end_coord = dict({
                "lat": flt_prop.flight[flight_id]['last_latitude'],
                "lon": flt_prop.flight[flight_id]['last_longitude'],
                "alt_m": flt_prop.flight[flight_id]['last_altitude'],
                "end_datetime": self.BALLOON_END_DATETIME,  # timestamp
            })

            if self.SIM_TIME_HOURS != None:
                self.SIM_TIME_ITERATIONS = self.SIM_TIME_HOURS * int(
                3600 * (1 / self.DT_SEC))  # number of simulation points
            elif self.end_timestamp != None and self.start_timestamp != None:
                self.sim_time_sec = mt.datetime_epoch(self.end_timestamp) - mt.datetime_epoch(self.start_timestamp)
                self.SIM_TIME_ITERATIONS = int(self.sim_time_sec / self.DT_SEC)

            self.TRACK_DATETIME = self.BALLOON_START_DATETIME
            self.TRACK_ALTITUDE = self.MIN_ALT_M
            print("Finished loading configuration properties...")

        def load_new(self, flight_id):
            GMT = sim_prop.GMT
            MIN_ALT_M = flt_prop.flight[flight_id]['launch_altitude']
            DT_SEC = sim_prop.dt_sec  # important time for numerical integration of the physics models
            SIM_TIME_HOURS = flt_prop.flight[flight_id]["flight_hours"]
            BALLOON_NAME = sim_prop.simulation["balloon_name"]
            ASCENT_RATE_MS = flt_prop.flight[flight_id]['ascent_rate_ms']
            GFS_RATE_SEC = sim_prop.GFS["GFSrate_sec"]
            OUTPUT_DIR = sim_prop.output_dir
            FLOAT_HGT_M = sim_prop.simulation['float_hgt_m']
            START_COORD = dict(sim_prop.simulation['start_coord'])
            END_COORD = dict(sim_prop.simulation['end_coord'])

            # New Variable to add to properites
            LEG_TIME_HOURS = sim_prop.simulation["leg_time_hr"]
            LEG_TIME_ITERATIONS = LEG_TIME_HOURS * int(3600 * (1 / DT_SEC))  # number of simulation points
            LEG_TIME_SEGMENTS = SIM_TIME_HOURS / LEG_TIME_HOURS
            calc_methods = sim_prop.simulation["calc_methods"]
            OUTPUT = sim_prop.output
            STRATO_BOUNDARY = flt_prop.flight[flight_id]['strato_boundary_m']
            TROPO_ASCENT_MS = flt_prop.flight[flight_id]['tropo_ascent_rate_ms']
            TROPO_DESCENT_MS = flt_prop.flight[flight_id]['tropo_descent_rate_ms']
            STRATO_ASCENT_MS = flt_prop.flight[flight_id]['tropo_ascent_rate_ms']
            STRATO_DESCENT_MS = flt_prop.flight[flight_id]['tropo_descent_rate_ms']
            CEILING_M = sim_prop.simulation['ceiling_alt_m']
            PATH_DIRECTION = sim_prop.simulation['path_direction']

            # This is done for date changes
            SIM_GFS = sim_prop.gfs
            NET_GFS = net_prop.gfs

            # New stuff for the hover stage
            TOLERANCE_LEVEL = sim_prop.simulation[
                'tolerance_level']  # this will be a percentage of target range and adjustable
            TARGET_RANGE = sim_prop.simulation['target_range']  # This will be pulled from the mission params
            HOVER_BOUNDRY = TOLERANCE_LEVEL * TARGET_RANGE  # This is calculated based on the above

            # For the Graph
            ALTITUDE_DIR = sim_prop.simulation['alt_graph_dir']

            # For the daatsource
            DATA_SOURCE = sim_prop.datasource

            # new era options
            telemetry_dir = sim_prop.telemetry_dir

            # These are outputs
            html_map = sim_prop.output['html']  # this is a htmlmap output
            leaflet_map = sim_prop.output['leaflet_map']  # this is a leaflet map output
            geojson_file = sim_prop.output['geojson']  # this is a geojson output
            csv_file = sim_prop.output['csv']  # this is a csv output
            cartopy_map = sim_prop.output['cartopy_map']  # this is a cartopy output
            balloon_height = sim_prop.output['balloon_height']  # this is a blloon hieght output

            compare_to_telemetry = flt_prop.flight[flight_id]["compare_to_telemetry"]
            telemetry_file = None
            telemetry_source = None
            use_telemetry_set = flt_prop.flight[flight_id]["use_telemetry_set"]
            use_time_dimension = flt_prop.flight[flight_id]["use_time_dimension"]
            if compare_to_telemetry or use_telemetry_set:
                telemetry_file = flt_prop.flight[flight_id]["telemetry_file"]
                telemetry_source = flt_prop.flight[flight_id]["telemetry_source"]
                work_sheet_index = flt_prop.flight[flight_id]["work_sheet_index"]

            balloon_type = flt_prop.flight[flight_id]["balloon_type"]
            start = flt_prop.flight[flight_id]['launch_time']
            end = flt_prop.flight[flight_id]['last_time']
            if start == None:
                BALLOON_START_DATETIME = None
            else:
                BALLOON_START_DATETIME = datetime.fromisoformat(start)

            if end != None:
                BALLOON_END_DATETIME = datetime.fromisoformat(end)
            else:
                BALLOON_END_DATETIME = None

            input_file = flt_prop.flight[flight_id]["model_file"]

            launch_latitude = flt_prop.flight[flight_id]["launch_latitude"]
            launch_longitude = flt_prop.flight[flight_id]['launch_longitude']
            launch_altitude = flt_prop.flight[flight_id]['launch_altitude']

            last_latitude = flt_prop.flight[flight_id]["last_latitude"]
            last_longitude = flt_prop.flight[flight_id]['last_longitude']
            last_altitude = flt_prop.flight[flight_id]['last_altitude']

            if SIM_TIME_HOURS == None:
                if BALLOON_END_DATETIME == None:
                    logging.debug("In your flight time you must either set the number of hours to simulate via 'flight_hours' or ")
                    logging.debug("you must set the 'last_time' or being using telemetry else unknown length of simulation")

            descent_rate_ms = flt_prop.flight[flight_id]['descent_rate_ms']
            strato_boundary_m = flt_prop.flight[flight_id]['strato_boundary_m']
            if strato_boundary_m == None:
                strato_ascent_rate_ms = None
                strato_descent_rate_ms = None
            else:
                strato_ascent_rate_ms = flt_prop.flight[flight_id]['strato_ascent_rate_ms']
                strato_descent_rate_ms = flt_prop.flight[flight_id]['strato_descent_rate_ms']

            burst_hgt_m = flt_prop.flight[flight_id]['burst_hgt_m']
            if burst_hgt_m == None:
                burst_balloon = False
            else:
                burst_balloon = flt_prop.flight[flight_id]['burst']

            ceiling_hgt_m = flt_prop.flight[flight_id]['ceiling_hgt_m']
            if burst_balloon == True and ceiling_hgt_m != None:
                logging.debug("Warning you say you want the balloon to burst, but you also gave us a ceiling or float height, results may not be as expected")

            start_timestamp = BALLOON_START_DATETIME
            if BALLOON_END_DATETIME != None:
                end_timestamp = BALLOON_END_DATETIME

            start_coord = dict({
                "lat": flt_prop.flight[flight_id]['launch_latitude'],
                "lon": flt_prop.flight[flight_id]['launch_longitude'],
                "alt_m": flt_prop.flight[flight_id]['launch_altitude'],
                "timestamp": BALLOON_START_DATETIME,  # timestamp
            })
            end_coord = dict({
                "lat": flt_prop.flight[flight_id]['last_latitude'],
                "lon": flt_prop.flight[flight_id]['last_longitude'],
                "alt_m": flt_prop.flight[flight_id]['last_altitude'],
                "end_datetime": BALLOON_END_DATETIME,  # timestamp
            })

            if SIM_TIME_HOURS != None:
                SIM_TIME_ITERATIONS = SIM_TIME_HOURS * int(3600 * (1 / DT_SEC))  # number of simulation points
            else:
                sim_time_sec = mt.datetime_epoch(end_timestamp) - mt.datetime_epoch(start_timestamp)
                SIM_TIME_ITERATIONS = int(sim_time_sec / DT_SEC)

            TRACK_DATETIME = BALLOON_START_DATETIME
            TRACK_ALTITUDE = MIN_ALT_M

            return 1

            '''
            This will be added later to remove it and standardize it, 
            We will need to have the version be passed someway, maybe thats the thing passed?
            def log_configuration(self):
                """Logs parameters to start the simulation

                Parameters:
                rc (run_config object): Class that stores input paremeters

                Returns:
                None

                readme = "Running trapezoid version 1.4 (03-01-2022)"
                logging.info(readme)
                logging.info("=======================================================================================================")
                logging.info("Runtime settings are taking from simulation_plan, flight_plans, and to a lesser extent era5_properties and balloon_properties")
                logging.info("We start the simulation as follows: ")
                logging.info(f"flight_id: {self.flight_id}")
                logging.info(f"NetCDF data will be obtained from {self.input_file}")
                logging.info(f"GMT offset: {rc.GMT}")
                if rc.compare_to_telemetry:
                    logging.info(f"telemetry_dir: {rc.telemetry_file}")
                logging.info(f"gfs_rate_s in seconds: {rc.gfs_rate_s}, sampling rate or dt_sec: {rc.dt_sec}")

                logging.info(f"balloon_type: {rc.balloon_type}")
                logging.info(f"Balloon start and end datatime:{rc.balloon_start_datetime}, {rc.balloon_end_datetime}")
                logging.info(f"launch_latitude: {rc.launch_latitude}, launch_longitude: {rc.launch_longitude}, launch_altitude: {rc.launch_altitude}")
                logging.info(f"last_latitude: {rc.last_latitude}, last_longitude: {rc.last_longitude}, last_altitude: {rc.last_altitude}")
                logging.info(f"Minimum alt: {rc.min_alt_m}")
                logging.info(" ")
                logging.info(f"Simulation time in hours: {rc.sim_time_hrs}")
                logging.info(f"ascent_rate_ms: {rc.ascent_rate_ms}, descent_rate_ms: {rc.descent_rate_ms}")


                if rc.strato_boundary_m != None:
                    logging.info(f"strato_boundary_m: {rc.strato_boundary_m}")
                    logging.info(f"strato_ascent_rate_ms: {rc.strato_ascent_rate_ms}, strato_descent_rate_ms: {rc.strato_descent_rate_ms}")

                logging.info(f"burst_hgt_m: {rc.burst_hgt_m}")
                '''


def main(flight_id):
    default_config = sim_config(flight_id)

    return default_config


if __name__ == "__main__":
    main(flight_id=sim_prop.flight_id)

