from pyproj import CRS, Transformer
import math
import numpy as np


def wrap_lon(lon):
    """Wrap longitude into [-180, 180)."""
    return ((lon + 180.0) % 360.0) - 180.0


def clamp_lat(lat):
    """Clamp latitude to the valid open range (avoids the poles exactly)."""
    return float(np.clip(lat, -89.999999, 89.999999))


def latlon_to_meters_merc(origin_lat, origin_lon, target_lat, target_lon):
    # Define CRS for WGS84 (lat/lon) and UTM (appropriate zone)
    crs_latlon = CRS("epsg:4326")  # WGS84
    # These 2 should be equivelent
    #crs_merc = CRS(proj="merc", lon_0=0)  # merc?
    crs_merc = CRS("epsg:3857")  # Web Mercator (Projected Coordinate System in meters)

    #crs_utm = CRS("epsg:32612")  # UTM Zone 12N (for Northern Hemisphere)
    #crs_stateplane = CRS("epsg:32140")  # Example for New Mexico State Plane

    # Create a Transformer for converting lat/lon to UTM in meters
    transformer = Transformer.from_crs(crs_latlon, crs_merc, always_xy=True)

    # Transform origin and target coordinates to UTM
    origin_x, origin_y = transformer.transform(origin_lon, origin_lat)
    target_x, target_y = transformer.transform(target_lon, target_lat)

    # Calculate relative position in meters
    relative_x = target_x - origin_x
    relative_y = target_y - origin_y

    # If lat or lon is nan from the transform step being too close to the central coordinate, revert back to central coordinate distance
    if np.isnan(relative_x):
        relative_x = 0
    if np.isnan(relative_y):
        relative_y = 0

    return relative_x, relative_y

def latlon_to_meters_spherical(origin_lat, origin_lon, target_lat, target_lon):
    # Convert degrees to radians
    origin_lat_rad = math.radians(origin_lat)
    #target_lat_rad = math.radians(target_lat) # unused now

    # Calculate distance for latitude change
    lat_distance = (target_lat - origin_lat) * 111320  # meters per degree of latitude

    # Calculate distance for longitude change
    #longitude_distance = (target_lon - origin_lon) * 111000 * math.cos((origin_lat_rad + target_lat_rad) / 2) #OLD OUTDATED TRANSFROM.  THIS FAILED OUTSIDE of 150 km due to Transform errors
    longitude_distance = (target_lon - origin_lon) * 111320 * math.cos(origin_lat_rad)

    # If lat or lon is nan from the transform step being too close to the central coordinate, revert back to central coordinate distance
    if np.isnan(longitude_distance):
        longitude_distance = 0
    if np.isnan(lat_distance):
        lat_distance = 0

    # Latitude distance will be positive or negative depending on direction
    return longitude_distance, lat_distance


def meters_to_latlon_spherical(origin_lat, origin_lon, x_m, y_m):
    # Convert degrees to radians
    origin_lat_rad = math.radians(origin_lat)

    # Convert meter distances back to latitude and longitude changes
    # Latitude change
    lat_change = y_m / 111_320.0  # meters to degrees

    # Longitude change
    lon_change = x_m / (111_320.0 * math.cos(origin_lat_rad))  # meters to degrees

    # Calculate target latitude and longitude
    target_lat = origin_lat + lat_change
    target_lon = origin_lon + lon_change


    # If lat or lon is nan from the transform step being too close to the central coordinate, revert back to central coordinate. 
    if np.isnan(target_lat):
        target_lat = origin_lat
    if np.isnan(target_lon):
        target_lon = origin_lon


    return target_lat, target_lon


def mercator_to_latlon(origin_lat, origin_lon, x_m, y_m):
    # Convert degrees to radians for origin latitude
    origin_lat_rad = math.radians(origin_lat)

    # Mercator projection constants
    R = 6378137  # Radius of Earth in meters (WGS84)

    # Calculate longitude change in degrees
    lon_change = x_m / (R * math.cos(origin_lat_rad)) * (180 / math.pi)

    # Calculate latitude change in degrees using inverse Mercator formula
    # y_m is measured from the origin (equator) in Mercator projection
    lat_change = (2 * math.atan(math.exp(y_m / R)) - math.pi / 2) * (180 / math.pi)

    # Calculate target latitude and longitude
    target_lat = origin_lat + lat_change
    target_lon = origin_lon + lon_change

    # If lat or lon is nan from the transform step being too close to the central coordinate, revert back to central coordinate. 
    if np.isnan(target_lat):
        target_lat = origin_lat
    if np.isnan(target_lon):
        target_lon = origin_lon

    return target_lat, target_lon


if __name__ == '__main__':

    origin_lat = 12.75
    origin_lon = 121.5
    x_m = 150000
    y_m = 150000

    print(meters_to_latlon_spherical(origin_lat, origin_lon, x_m, y_m))



    '''
    # Example usage
    origin_lat, origin_lon = 36, -105
    target_lat = 36  # Adding 1 degree to latitude
    target_lon = -105  # Longitude unchanged


    #Show Example Transformations using mercator versus spherical

    print("Spherical")
    x_m, y_m = latlon_to_meters_spherical(origin_lat, origin_lon, target_lat+1, target_lon)
    print(f"Relative position in meters: x = {x_m:.2f}, y = {y_m:.2f}")

    target_lon_inv, target_lat_inv = meters_to_latlon_spherical(origin_lat, origin_lon, x_m, y_m)
    print(f"Target Coordinates: lat = {target_lat_inv:.2f}, lon = {target_lon_inv:.2f}")

    print()


    x_m, y_m = latlon_to_meters_spherical(origin_lat, origin_lon, target_lat, target_lon+1)
    print(f"Relative position in meters: x = {x_m:.2f}, y = {y_m:.2f}")

    target_lon_inv, target_lat_inv = meters_to_latlon_spherical(origin_lat, origin_lon, x_m, y_m)
    print(f"Target Coordinates: lat = {target_lat_inv:.2f}, lon = {target_lon_inv:.2f}")

    ##########################
    print()
    print()

    print("Mercator")
    x_m, y_m = latlon_to_meters_merc(origin_lat, origin_lon, target_lat+1, target_lon)
    print(f"Relative position in meters: x = {x_m:.2f}, y = {y_m:.2f}")

    target_lon_inv, target_lat_inv = mercator_to_latlon(origin_lat, origin_lon, x_m, y_m)
    print(f"Target Coordinates: lat = {target_lat_inv:.2f}, lon = {target_lon_inv:.2f}")

    print()

    x_m, y_m = latlon_to_meters_merc(origin_lat, origin_lon, target_lat, target_lon+1)
    print(f"Relative position in meters: x = {x_m:.2f}, y = {y_m:.2f}")

    target_lon_inv, target_lat_inv = mercator_to_latlon(origin_lat, origin_lon, x_m, y_m)
    print(f"Target Coordinates: lat = {target_lat_inv:.2f}, lon = {target_lon_inv:.2f}")

    '''