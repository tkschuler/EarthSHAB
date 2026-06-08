from datetime import datetime, timedelta
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()     # Hacky solution for Python 3.6 to use ISO format Strings

balloon_properties = dict(
    shape = 'sphere',
    d = 5.8,                          # (m) Diameter of Sphere Balloon
    mp = 0.9,                         # (kg) Mass of Payload
    areaDensityEnv = 939.*7.87E-6,    # (Kg/m^2) rhoEnv*envThickness
    mEnv = 2.1,                       # (kg) Mass of Envelope - SHAB1
    cp = 2000.,                       # (J/(kg K)) Specific heat of envelope material
    absEnv = .93,                     # Absorbiviy of envelope material
    emissEnv = .92,                   # Emisivity of enevelope material
    Upsilon = 4.5,                    # Ascent Resistance coefficient
)

parent_dir = "src/EarthSHAB/"

# SHAB14-V Example for EarthSHAB software. Runs main.py against the bundled
# SHAB14-V flight (GFS + APRS truth track) so the prediction can be compared to
# the real trajectory; the same flight is registered in evaluation/launches.json
# for the evaluation suite. To run a current-day prediction instead, set a recent
# forecast_start_time (cycle hour 00/06/12/18 UTC), set balloon_trajectory = None,
# and download the forecast first with `python -m EarthSHAB.saveNETCDF`.
forecast_start_time =  "2022-08-22 12:00:00" # Forecast start time, should match a downloaded forecast in the forecasts directory
start_time = datetime.fromisoformat("2022-08-22 14:36:00") # Simulation start time. The end time needs to be within the downloaded forecast
balloon_trajectory = parent_dir + "balloon_data/SHAB14V-APRS.csv"  # Only Accepting Files in the Standard APRS.fi format for now

# Single forecast file path. The reader (EarthSHAB.Forecast.Forecast) opens
# this file regardless of whether it came from GFS or ERA5 — source is read
# from the file's `institution` global attribute. To run an ERA5-based
# simulation, point `file` at an ERA5 .nc instead.
_gfs_res = 0.25   # (deg) GFS grid resolution: 0.25, 0.5, or 1.0
_gfs_res_token = ("%.2f" % _gfs_res).replace(".", "p")   # 0.25->0p25, 0.5->0p50, 1.0->1p00
_gfs_step_hours = 3   # (h) temporal step between forecast hours (see netcdf_gfs below)
_gfs_step_token = f"{int(_gfs_step_hours)}h"             # 1->1h, 3->3h
# Cycles older than NOAA's ~9-day NOMADS live-retention window can't be fetched
# live; saveNETCDF auto-switches to the AWS archive, which writes an `_archive`-
# marked file (see saveNETCDF_archive._archive_output_path). Mirror that marker
# here so `file` below points at whatever the downloader actually produces.
_GFS_RETENTION_DAYS = 9
_oldest_live_cycle = (datetime.utcnow() - timedelta(days=_GFS_RETENTION_DAYS)).replace(
    hour=0, minute=0, second=0, microsecond=0)
_is_archive_cycle = (
    datetime.fromisoformat(forecast_start_time).replace(minute=0, second=0, microsecond=0)
    < _oldest_live_cycle
)
_archive_suffix = "_archive" if _is_archive_cycle else ""

# Filename encodes both the spatial resolution and the temporal step, e.g.
# gfs_0p25_3h_20260605_06.nc, so files of different grids/cadences don't collide.
# Past-retention cycles get an `_archive` suffix: gfs_0p25_3h_20220822_12_archive.nc.
_default_gfs_file = (
    parent_dir + "forecasts/gfs_" + _gfs_res_token + "_" + _gfs_step_token + "_"
    + forecast_start_time[0:4] + forecast_start_time[5:7] + forecast_start_time[8:10]
    + "_" + forecast_start_time[11:13] + _archive_suffix + ".nc"
)

forecast = dict(
    file = _default_gfs_file,
    forecast_start_time = forecast_start_time, # used to build the default file path above
    forecast_update_interval = 60,               # (s) After how many iterated dt steps are new wind speeds are looked up

    # Wind interpolation method used inside Forecast.wind_alt_Interpolate2:
    #   'linear_neighbors' - (default, historical) bearing+speed linearly
    #                        interpolated between the 2 nearest pressure
    #                        levels, with angle-wrap correction.
    #   'linear_full'      - np.interp on u and v independently across the
    #                        full altitude profile.
    #   'spline_full'      - scipy CubicSpline on u and v across the full
    #                        altitude profile; extrapolate=False, with
    #                        np.interp fallback when alt is outside the
    #                        profile bounds (avoids spline overshoot above
    #                        the highest pressure level).
    wind_interpolation = 'linear_full',
)

# GFS downloader configuration. Consumed only by saveNETCDF.py — the Forecast
# reader does not look at this dict.
netcdf_gfs = dict(
    nc_file = _default_gfs_file,     # output path for the downloaded forecast
    nc_start = datetime.fromisoformat(forecast_start_time),  # forecast cycle to fetch
    hourstamp = forecast_start_time[11:13],

    res = _gfs_res,    # (deg) GFS grid resolution: 0.25, 0.5, or 1.0. Must match
                       #     the value used to build nc_file above (_gfs_res).

    lat_range = 40,    # (deg) bounding-box height to download
    lon_range= 60,     # (deg) bounding-box width to download
    download_days = 1, # (1-10) forecast horizon in days

    step_hours = _gfs_step_hours,  # (h) temporal resolution of forecast steps to
                       #     download; also encoded in nc_file above (_gfs_step_hours).
                       #     1 = hourly. NOTE: only the 0.25° grid provides hourly
                       #     steps (out to f120 / 5 days; 3-hourly beyond). The
                       #     0.5° and 1.0° grids are 3-hourly only, so set
                       #     step_hours = 3 (or a multiple) for those — finer
                       #     steps that don't exist on the server are skipped.
)

simulation = dict(
    start_time = start_time,    # (UTC) Simulation Start Time, updated above
    sim_time = 15,              # (int) (hours) Number of hours to simulate

    vent = 0.0,                 # (kg/s) Vent Mass Flow Rate  (Do not have an accurate model of the vent yet, this is innacurate)
    alt_sp = 15000.0,           # (m) Altitude Setpoint
    v_sp = 0.,                  # (m/s) Altitude Setpoint, Not Implemented right now
    start_coord =	{
                      "lat": 34.60,             # (deg) Latitude
                      "lon": -106.80,           # (deg) Longitude
                      "alt": 1000.,             # (m) Elevation
                      "timestamp": start_time,  # current timestamp
                    },
    min_alt = 1000.,            # starting altitude. Generally the same as initial coordinate
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
