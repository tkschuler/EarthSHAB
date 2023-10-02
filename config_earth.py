from datetime import datetime
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()     # Hacky solution for Python 3.6 to use ISO format Strings

balloon_properties = dict(
    shape = 'sphere',
    d = 5.8,                          # (m) Diameter of Sphere Balloon
    mp = 1.15,                        # (kg) Mass of Payload
    areaDensityEnv = 939.*7.87E-6,    # (Kg/m^2) rhoEnv*envThickness
    mEnv = 1.3,                       # (kg) Mass of Envelope - SHAB1
    cp = 2000.,                       # (J/(kg K)) Specific heat of envelope material
    absEnv = .93,                     # Absorbiviy of envelope material
    emissEnv = .92,                   # Emisivity of enevelope material
    Upsilon = 2.5,                    # Ascent Resistance coefficient
)

gfs = "2021-03-29 12:00:00" # Forecast start time, should match a downloaded forecast
start_time = datetime.fromisoformat("2021-03-29 11:32:00") # Simulation start time. The end time needs to be within the downloaded forecast

#These parameters are for both downloading new forecasts, and running simulations with downloaded forecasts.
netcdf = dict(
    nc_file = ("forecasts/gfs_0p25_" + gfs[0:4] + gfs[5:7] + gfs[8:10] + "_" + gfs[11:13] + ".nc"),  # file structure for downloading .25 resolution NOAA forecast data.
    nc_start = datetime.fromisoformat(gfs),    # Start time of the downloaded netCDF file
    hourstamp = gfs[11:13],  # parsed from gfs timestamp

    res = 0.25,       # (deg) Do not change
    lat_range = 40,  # (.25 deg)
    lon_range= 60,   # (.25 deg)
    hours3 = 8,      # (1-80) In intervals of 3 hours.  hour_index of 8 is 8*3=24 hours
)


simulation = dict(
    start_time = start_time,    # (UTC) Simulation Start Time, updated above
    sim_time = 16,              # (hours) Simulation time in hours (for trapezoid.py)

    vent = 0.0,                 # (kg/s) Vent Mass Flow Rate  (Do not have an accurate model of the vent yet, this is innacurate)
    alt_sp = 15000.0,           # (m) Altitude Setpoint
    v_sp = 0.,                  # (m/s) Altitude Setpoint, Not Implemented right now
    start_coord =	{
                      "lat": 39.828,           # (deg) Latitude
                      "lon": -98.5795,         # (deg) Longitude
                      "alt": 408.,             # (m) Elevation
                      "timestamp": start_time, # timestamp
                    },
    min_alt = 408.,             # starting altitude. Generally the same as initial coordinate
    float = 11500.              # Maximum float altitude for simple trapezoidal trajectories
)

GFS = dict(
    GFSrate = 60,               # (#) After how many iterated dt steps are new wind speeds are looked up
)

earth_properties = dict(
    Cp_air0 = 1003.8,           # (J/Kg*K)  Specifc Heat Capacity, Constant Pressure
    Cv_air0 = 716.,             # (J/Kg*K)  Specifc Heat Capacity, Constant Volume
    Rsp_air = 287.058,          # (J/Kg*K) Gas Constant
    P0 = 101325.0,              # (Pa) Pressure @ Surface Level
    emissGround = .95,          # assumption
    albedo = 0.17,              # assumption
)

dt = 3.0 # (s) Time Step for integrating (If error's occur, use a lower step size)


