from datetime import datetime
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()   # Hacky solution for Python 3.6 to use ISO format Strings

balloon_properties = dict(
    shape = 'sphere',
    d = 6.0,                          # (m) Diameter of Sphere Balloon
    mp = 1.15,                        # (kg) Mass of Payload
    areaDensityEnv = 939.*7.87E-6,    # (Kg/m^2) rhoEnv*envThickness
    mEnv = 1.25,                      # (kg) Mass of Envelope - SHAB1
    cp = 2000.,                       #(J/(kg K)) Specific heat of envelope material
    absEnv = 0.93,                    # Absorbiviy of envelope material
    absEnvIR = 0.9,                   # Absorbiviy of envelope material
    emissEnv = 0.9,                   # Emisivity of enevelope material
    Beta = 4,                         # Additional Drag coefficient for upward velocity.
)

# Make sure nc_date and gfs_time match
nc_date = "gfs_0p25_20210223_12.nc"                         # file structure for downloading .25 resolution NOAA forecast data.
gfs_time = datetime.fromisoformat("2021-02-23 12:00:00")    # Start time of the downloaded netCDF file
start_time = datetime.fromisoformat("2021-02-23 15:00:00")  # UTC
sim = 12                                                    # Simulation time in hours (for trajectory.py)

GNC = dict(
    vent = 0.0,                 # (kg/s) Vent Mass Flow Rate  (Do not have an accurate model of the vent yet, this is innacurate)
    alt_sp = 15000.0,           # (m) Altitude Setpoint
    v_sp = 0.,                  # (m/s) Altitude Setpoint, Not Implemented right now
    start_time = start_time,    # Start Time
    start_coord =	{
                      "lat": 32.2226,          # (deg) Latitude
                      "lon": -110.9747,        # (deg) Longitude
                      "alt": 408.,             # (m) Elevation
                      "timestamp": start_time, # timestamp
                    },
    min_alt = 408.,             # starting altitude. Generally the same as initial coordinate
    float = 11500.              # Maximum float altitude for simple trapezoidal trajectories
)

GFS = dict(
    res = .25,                  # (deg) Do not change
    lat_range = 25,             # (.25 deg)
    lon_range = 50,             # (.25 deg)
    hour_index = 8,             # (1-80) In intervals of 3 hours.  hour_index of 8 is 8*3=24 hours
    file = nc_date,             # forecast file specified above
)

earth_properties = dict(
    Cp_air0 = 1003.8,           # (J/Kg*K)  Specifc Heat Capacity, Constant Pressure
    Cv_air0 = 716.,             # (J/Kg*K)  Specifc Heat Capacity, Constant Volume
    Rsp_air = 287.058,          # (J/Kg*K) Gas Constant
    P0 = 101325.0,              # (Pa) Pressure @ Surface Level
    emissGround = .95,          # assumption
    albedo = 0.17,              # assumption
)

dt = 1.0 # (s) Time Step for integrating (If error's occur, use a lower step size, However may not work with GFS trajectories)
