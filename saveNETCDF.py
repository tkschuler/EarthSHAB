"""
saveNETCDF.py downloads a portion of the netCDF file for the specified day and
area specified in earth_config.py.

Archive forecast dataset retirval is not currently supported.
"""

import netCDF4 as nc4
import config_earth
from termcolor import colored
import numpy as np
import sys

import os

coord = config_earth.simulation['start_coord']
lat_range = config_earth.netcdf_gfs['lat_range']
lon_range = config_earth.netcdf_gfs['lon_range']
download_days = config_earth.netcdf_gfs['download_days']
hourstamp = config_earth.netcdf_gfs['hourstamp']
res = config_earth.netcdf_gfs['res']

nc_start = config_earth.netcdf_gfs['nc_start']

if not os.path.exists('forecasts'):
    os.makedirs('forecasts')

def closest(arr, k):
    """ Given an ordered array and a value, determines the index of the closest item
    contained in the array.
    """
    return min(range(len(arr)), key = lambda i: abs(arr[i]-k))

def getNearestLat(lat,min,max):
    """ Determines the nearest lattitude (to .25 degrees)
    """
    arr = np.arange(start=min, stop=max, step=res)
    i = closest(arr, lat)
    return i

def getNearestLon(lon,min,max):
    """ Determines the nearest longitude (to .25 degrees)
    """
    lon = lon % 360 #convert from -180-180 to 0-360
    arr = np.arange(start=min, stop=max, step=res)

    i = closest(arr, lon)
    return i

lat_i = getNearestLat(coord["lat"],-90,90.01)
lon_i = getNearestLon(coord["lon"],0,360)

coords = ['time', 'lat', 'lon', 'lev']
vars_out = ['ugrdprs', 'vgrdprs', 'hgtprs', 'tmpprs']

# parse GFS forecast start time from config file
year = str(nc_start.year)
month = str(nc_start.month).zfill(2)
day = str(nc_start.day).zfill(2)

print("Downloading data from")
url = "https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs" + year + month + day + "/gfs_0p25_" + str(hourstamp) + "z"
print(colored(url,"cyan"))

# Open input file in read (r), and output file in write (w) mode:
try:
    nc_in = nc4.Dataset(url)

except:
    print(colored("NOAA DODS Server error with timestamp " + str(nc_start) + ". Data not downloaded.", "red"))
    sys.exit()


#Some print statistics:
# Get dimensions directly from dataset
num_lats = len(nc_in.dimensions['lat'])
num_lons = len(nc_in.dimensions['lon'])

#print netcdf structre
print("GFS NetCDF structure:")

print(f"Number of latitude points: {num_lats}")
print(f"Number of longitude points: {num_lons}")

lat_vals = nc_in.variables['lat'][:]
lon_vals = nc_in.variables['lon'][:]

print(f"Lat min: {lat_vals.min()}, Lat max: {lat_vals.max()}")
print(f"Lon min: {lon_vals.min()}, Lon max: {lon_vals.max()}\n\n")

print("Download/Fill Data dimesnions")
lat_start_i = max(0, lat_i - lat_range)
lat_end_i = min(num_lats, lat_i + lat_range)

lon_start_i = max(0, lon_i - lon_range)
lon_end_i = min(num_lons, lon_i + lon_range)

print(f"Clipped lat indices: {lat_start_i}-{lat_end_i}, size: {lat_end_i - lat_start_i}")
print(f"Clipped lon indices: {lon_start_i}-{lon_end_i}, size: {lon_end_i - lon_start_i}")

# Print corresponding lat/lon degree bounds
print(f"Lat degrees: {lat_vals[lat_start_i]} to {lat_vals[lat_end_i-1]}")
print(f"Lon degrees: {lon_vals[lon_start_i]} to {lon_vals[lon_end_i-1]}\n\n")



# Create output netCDF file
nc_out = nc4.Dataset(config_earth.netcdf_gfs['nc_file'], 'w')

# Create dimensions in output file
for name, dimension in nc_in.dimensions.items():
    nc_out.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

# Create coordinate variables in output file
for name, variable in nc_in.variables.items():
    if name in coords:
        x = nc_out.createVariable(name, variable.datatype, variable.dimensions)
        v = nc_in.variables[name][:]
        nc_out.variables[name][:] = v
        print ("Downloaded " + name)

    if name in vars_out:
        print ("Downloading " + name)
        x = nc_out.createVariable(name, variable.datatype, variable.dimensions, zlib=True) #Without zlib the file will be MASSIVE

        #Download only a chunk of the data
        for i in range(0,download_days*8+1):  #In intervals of 3 hours. hour_index of 8 is 8*3=24 hours. Add one more index to get full day range
            data = nc_in.variables[name][i,0:34,lat_start_i:lat_end_i,lon_start_i:lon_end_i] #This array can only have a maximum of  536,870,912 elements, Need to dynamically add.
            nc_out.variables[name][i,0:34,lat_start_i:lat_end_i,lon_start_i:lon_end_i] = data
            #print("Downloaded and added to output file ", name, ' hour index - ', i, ' time - ', i*3)
            print(f"{name} hour index {i}: downloaded data shape {data.shape}")
