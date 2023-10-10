from datetime import datetime
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()     # Hacky solution for Python 3.6 to use ISO format Strings

balloon_properties = dict(
    shape = 'sphere',
    d = 5.8,                          # (m) Diameter of Sphere Balloon
    mp = 1.7,                        # (kg) Mass of Payload
    areaDensityEnv = 939.*7.87E-6,    # (Kg/m^2) rhoEnv*envThickness
    mEnv = 1.3,                       # (kg) Mass of Envelope - SHAB1
    cp = 2000.,                       # (J/(kg K)) Specific heat of envelope material
    absEnv = .98,                     # Absorbiviy of envelope material
    emissEnv = .95,                   # Emisivity of enevelope material
    Upsilon = 4.5,                    # Ascent Resistance coefficient
)

# Original GitHub EarthShab Test Data
#gfs = "2021-03-29 12:00:00" # Forecast start time, should match a downloaded forecast
#start_time = datetime.fromisoformat("2021-03-29 11:32:00") # Simulation start time. The end time needs to be within the downloaded forecast
#balloon_trajectory = None

#SHAB10
#gfs = "2022-04-09 12:00:00" # Forecast start time, should match a downloaded forecast
#start_time = datetime.fromisoformat("2022-04-09 18:14:00") # Simulation start time. The end time needs to be within the downloaded forecast
#balloon_trajectory = "balloon_data/SHAB10V-APRS.csv"  # Only Accepting Files in the Standard APRS.fi format for now

#SHAB3
#gfs = "2020-11-20 06:00:00" # Forecast start time, should match a downloaded forecast in the forecasts directory
#start_time = datetime.fromisoformat("2020-11-20 15:47:00") # Simulation start time. The end time needs to be within the downloaded forecast
#balloon_trajectory = "balloon_data/SHAB3V-APRS.csv"  # Only Accepting Files in the Standard APRS.fi format for now

#SHAB5
#gfs = "2021-05-12 12:00:00" # Forecast start time, should match a downloaded forecast in the forecasts directory
#start_time = datetime.fromisoformat("2021-05-12 14:01:00") # Simulation start time. The end time needs to be within the downloaded forecast
#balloon_trajectory = "balloon_data/SHAB5V_APRS_Processed.csv"  # Only Accepting Files in the Standard APRS.fi format for now

#SHAB12/15
gfs =  "2022-08-22 12:00:00" # Forecast start time, should match a downloaded forecast in the forecasts directory
start_time = datetime.fromisoformat("2022-08-22 14:21:00") # Simulation start time. The end time needs to be within the downloaded forecast
balloon_trajectory = None #"balloon_data/SHAB15V-APRS.csv"  # Only Accepting Files in the Standard APRS.fi format for now


#GFS AND ERA5 CAN BE DIFFERENT START TIMES???


forecast_type = "GFS" # GFS or ERA5

#These parameters are for both downloading new forecasts, and running simulations with downloaded forecasts.
netcdf_gfs = dict(
    nc_file = ("forecasts/gfs_0p25_" + gfs[0:4] + gfs[5:7] + gfs[8:10] + "_" + gfs[11:13] + ".nc"),  # file structure for downloading .25 resolution NOAA forecast data.
    nc_start = datetime.fromisoformat(gfs),    # Start time of the downloaded netCDF file
    hourstamp = gfs[11:13],  # parsed from gfs timestamp

    res = 0.25,       # (deg) Do not change
    lat_range = 40,  # (.25 deg)
    lon_range= 60,   # (.25 deg)
    hours3 = 8,      # (1-80) In intervals of 3 hours.  hour_index of 8 is 8*3=24 hours

)

netcdf_era5 = dict(
    #Need to add in start time stuff for ERA5 forecasts


    #filename = "SHAB3V_era_20201120_20201121.nc", #SHAB3
    #filename = "SHAB5V-ERA5_20210512_20210513.nc", #SHAB5V
    #filename = "shab10_era_2022-04-09to2022-04-10.nc", #SHAB10V
    #filename = "SHAB12V_ERA5_20220822_20220823.nc", #SHAB12V
    filename = "SHAB15V_ERA5_20220822_20220823.nc", #SHAB15V
    resolution_hr = 1
    )


simulation = dict(
    start_time = start_time,    # (UTC) Simulation Start Time, updated above
    sim_time = 14, #8,              # (hours) Simulation time in hours (for trapezoid.py)

    vent = 0.0,                 # (kg/s) Vent Mass Flow Rate  (Do not have an accurate model of the vent yet, this is innacurate)
    alt_sp = 15000.0,           # (m) Altitude Setpoint
    v_sp = 0.,                  # (m/s) Altitude Setpoint, Not Implemented right now
    start_coord =	{
                      "lat": 34.60,    #34.60,   32.44,    33.66,    # (deg) Latitude
                      "lon": -106.80,    #-106.80, -111.61, -114.22    # (deg) Longitude
                      "alt": 480., #720.,             # (m) Elevation
                      "timestamp": start_time, # timestamp
                    },
    min_alt = 480, #720.,             # starting altitude. Generally the same as initial coordinate
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

dt = 1.0 # (s) Time Step for integrating (If error's occur, use a lower step size)
