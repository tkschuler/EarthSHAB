'''
Note, currently this program must be manually meet to match how the user made a data download from copernicus dataset.
Unlike use of GFS data where the file downloaded is automatically named appropriately, you must manually rename it
to fit the format shown below. So for example, if our time for runtime = "2021-09-20 16:00:00" rename your downloaded
file as "era_20210920_160000.nc"
'''
from datetime import datetime
runtime = "2022-03-01 12:00:00"  # Forecast start time, should match a downloaded forecast
#filename = "era_" + runtime[0:4] + runtime[5:7] + runtime[8:10] + "_" + runtime[11:13] + runtime[14:16]+runtime[17:19]+".nc"
filename = "era.nc"

forecast_dir = "./forecasts/"
netcdf = dict(
    nc_file=(filename),
    # file structure for downloading .25 resolution NOAA forecast data.
    nc_start_datetime=datetime.fromisoformat(runtime),  # Start time of the downloaded netCDF file
    hourstamp=runtime[11:13],  # parsed from gfs timestamp
    res =-0.25,  # (deg) Do not change
    lat_range =40,  # (index into number .25 deg coordinates plus or minus from the start coordinates)
    lon_range =60,  # (.25 deg)
    hours1 =5,  # originally was 8 in (1-80) In intervals of 3 hours.  hour_index of 8 is 8*3=24 hours
    lat_bot =47.0, # How is this figured out? How does this relate to the index, for soem reason I am leaning to have these be index value?
    lat_top =27.0, # How is this figured out?
    lon_left =137.0, # How is this figured out?
    lon_right =165.0, # How is this figured out?
    resolution_hr = 1, # Added this
    res_deg =0.25 # added this as well
)
