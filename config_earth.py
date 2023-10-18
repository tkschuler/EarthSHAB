from datetime import datetime
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()     # Hacky solution for Python 3.6 to use ISO format Strings

balloon_properties = dict(
    shape = 'sphere',
    d = 6,                            # (m) Diameter of Sphere Balloon
    mp = 0.7,                         # (kg) Mass of Payload
    areaDensityEnv = 939.*7.87E-6,    # (Kg/m^2) rhoEnv*envThickness
    mEnv = 2,                         # (kg) Mass of Envelope - SHAB1
    cp = 2000.,                       # (J/(kg K)) Specific heat of envelope material
    absEnv = .98,                     # Absorbiviy of envelope material
    emissEnv = .95,                   # Emisivity of enevelope material
    Upsilon = 4.5,                    # Ascent Resistance coefficient
)

#SHAB14-V Example for EarthSHAB software
forecast_start_time =  "2022-08-22 12:00:00" # Forecast start time, should match a downloaded forecast in the forecasts directory
start_time = datetime.fromisoformat("2022-08-22 14:36:00") # Simulation start time. The end time needs to be within the downloaded forecast
balloon_trajectory = "balloon_data/SHAB14V-APRS.csv"  # Only Accepting Files in the Standard APRS.fi format for now

forecast = dict(
    forecast_type = "ERA5",      # GFS or ERA5
    forecast_start_time = forecast_start_time, # Forecast start time, should match a downloaded forecast in the forecasts directory
    GFSrate = 60,               # (s) After how many iterated dt steps are new wind speeds are looked up
)

#These parameters are for both downloading new forecasts, and running simulations with downloaded forecasts.
netcdf_gfs = dict(
    #DO NOT CHANGE
    nc_file = ("forecasts/gfs_0p25_" + forecast['forecast_start_time'][0:4] + forecast['forecast_start_time'][5:7] + forecast['forecast_start_time'][8:10] + "_" + forecast['forecast_start_time'][11:13] + ".nc"),  # DO NOT CHANGE -  file structure for downloading .25 resolution NOAA forecast data.
    nc_start = datetime.fromisoformat(forecast['forecast_start_time']),    # DO NOT CHANGE - Start time of the downloaded netCDF file
    hourstamp = forecast['forecast_start_time'][11:13],  # parsed from gfs timestamp

    res = 0.25,        # (deg) DO NOT CHANGE

    #The following values are for savenetcdf.py for forecast downloading and saving
    lat_range = 40,    # (.25 deg)
    lon_range= 60,     # (.25 deg)
    download_days = 1, # (1-10) Number of days to download for forecast This value is only used in saveNETCDF.py
)

netcdf_era5 = dict(
    filename = "SHAB14V_ERA5_20220822_20220823.nc", #SHAB14V ERA5 download from Copernicus
    resolution_hr = 1
    )

simulation = dict(
    start_time = start_time,    # (UTC) Simulation Start Time, updated above
    sim_time = 15,              # (hours) Number of hours to simulate

    vent = 0.0,                 # (kg/s) Vent Mass Flow Rate  (Do not have an accurate model of the vent yet, this is innacurate)
    alt_sp = 15000.0,           # (m) Altitude Setpoint
    v_sp = 0.,                  # (m/s) Altitude Setpoint, Not Implemented right now
    start_coord =	{
                      "lat": 34.60,             # (deg) Latitude
                      "lon": -106.80,           # (deg) Longitude
                      "alt": 1600.,             # (m) Elevation
                      "timestamp": start_time,  # current timestamp
                    },
    min_alt = 1600.,            # starting altitude. Generally the same as initial coordinate
    float = 23000,              # for simulating in trapezoid.py
    dt = 1.0,                   # (s) Integration timestep for simulation (If error's occur, use a lower step size)

    balloon_trajectory = balloon_trajectory # Default is None. Only accepting trajectories in aprs.fi csv format.
)

earth_properties = dict(
    Cp_air0 = 1003.8,           # (J/Kg*K)  Specifc Heat Capacity, Constant Pressure
    Cv_air0 = 716.,             # (J/Kg*K)  Specifc Heat Capacity, Constant Volume
    Rsp_air = 287.058,          # (J/Kg*K) Gas Constant
    P0 = 101325.0,              # (Pa) Pressure @ Surface Level
    emissGround = .95,          # assumption
    albedo = 0.17,              # assumption
)
