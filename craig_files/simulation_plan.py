'''
This is a very simple file that just tells us what flight plan we want to do Its based on.
what we give it as flight_id.

We separate flight plans from simulation so that we can create a database of previous flights, if we
want to go back. So here we just point to a flight plan (based on flight_id) and tell the simulation
what output options we want.

Data dependency note:
This program can compare computed flight paths with actual flight paths. It should be noted that currently,
we assume a specific format for the data as given by our first sample sets of NRL data. In the future,
we will build in support for different telemetry types.

See also: flight_simulator, flight_plans.py
'''
#flight_id = 'NRL20210920'
#flight_id = 'HBAL0445'
#flight_id = 'HBAL0446'
#flight_id = 'HBAL0447'
#flight_id = 'HBAL0449'
#flight_id = 'HBAL0448'

flight_id = 'era_sample_2022-06-01'
#flight_id = 'SAMPLE_DISTANCE'
#flight_id = 'SAMPLE_DIRECTION'
#flight_id = 'era_sample2'
#flight_id = 'HBAL0458'
flight_id = 'shab10v'
#flight_id = 'NRL20211217'  # this must be changed for different simulations

GMT = 7  # local time offset

output = dict(
    html = True,           # create a google.map representation of the optimized track
    leaflet_map = False,   # create a leaflet map representation of the optimized track
    geojson = False,       # output optimized track to geojson
    cartopy_map = True,    # create a mapplotlib
    balloon_height = True, # will map balloon height on a graph as a function of time
    csv = False,
    altitude = 'all',
    track = 'all' 
)
forecast_dir = "./forecasts"
output_dir = "./trajectories/"
telemetry_dir = "./balloon_data"
alt_graph_dir = './alt_graph'
data_sample_rate_s = 60  # (s) How often New Wind speeds are looked up
dt_sec = 3.0  # (s) Time Step for integrating (If error's occur, use a lower step size)
GFSrate_sec = 60 # (s) How often New Wind speeds are looked up
tmp_calc_methods = ['trapezoid_optimizer', 'ensemble_cc', 'psuedo_gradient', 'best_track'],
calc_methods = list(tmp_calc_methods[0]) # this is a quick workaround, need to solve the underlying issue
# for some reason the above is cast as a tuple. so this is forcing it to not be, theres a better way we an address later. 
# optional leaflet map configuration don't normally use this
leaflet_config = {
    "center": [39.05,-78.75],       # initial map center point
    "zoom": 6,                      # initial zoom level
    "show_zoom_control": False,     # hide / show zoom control
    "tile_provider": 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}',
    "attribution": 'Earth HAB',     # map attribution
    "show_attribution": False       # show / hide map attribution
}
