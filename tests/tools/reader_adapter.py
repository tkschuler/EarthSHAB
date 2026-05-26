"""Thin normalized API over GFS and ERA5 reader classes.

The two readers have differently-named methods and different signatures for
the same physical operations. Tests call this adapter so a single test body
can run against both readers without if/else branches everywhere.

When the consolidation refactor happens, only this adapter changes — the
tests stay the same. That property is the consolidation safety net.

Mapping:

    nearest_lat(reader, lat)              -> int index
    nearest_lon(reader, lon)              -> int index
    nearest_alt(reader, hour_idx, lat, lon, alt_m)  -> int index
    wind_interp(reader, coord)            -> (u, v, u_diag, v_diag)
    closest_idx(reader, arr, k)           -> int index (pure utility)

``coord`` is the same dict shape that GFS uses natively:
    {"lat": float, "lon": float, "alt": float, "timestamp": datetime}

For ERA5, the adapter converts internally to (alt_m, diff_time, lat_idx,
lon_idx) by reusing the reader's own model_start_datetime and lookup helpers.
"""

from __future__ import annotations

from EarthSHAB.ERA5 import ERA5
from EarthSHAB.GFS import GFS


def is_gfs(reader) -> bool:
    return isinstance(reader, GFS)


def is_era5(reader) -> bool:
    return isinstance(reader, ERA5)


def closest_idx(reader, arr, k):
    """Pure nearest-index search. GFS exposes `closest`, ERA5 exposes `closestIdx`."""
    if is_gfs(reader):
        return reader.closest(arr, k)
    return reader.closestIdx(arr, k)


def nearest_lat(reader, lat):
    """Index into the reader's lat array nearest to `lat` (degrees)."""
    if is_gfs(reader):
        return reader.getNearestLat(lat, reader.LAT_LOW, reader.LAT_HIGH)
    return reader.getNearestLatIdx(lat, reader.lat_max_idx, reader.lat_min_idx)


def nearest_lon(reader, lon):
    """Index into the reader's lon array nearest to `lon` (degrees).

    GFS normalizes lon to [0, 360) internally; ERA5 expects lon in the same
    convention as the file (synthetic dummies use [0, 2] which is valid in
    both 0-360 and -180-180 styles).
    """
    if is_gfs(reader):
        return reader.getNearestLon(lon, reader.LON_LOW, reader.LON_HIGH)
    return reader.getNearestLonIdx(lon, reader.lon_min_idx, reader.lon_max_idx)


def nearest_alt(reader, hour_idx, lat, lon, alt_m):
    """Index into the reader's altitude column nearest to alt_m at (hour, lat, lon).

    GFS takes (hour_idx, lat, lon, alt). ERA5 takes
    (int_hr_idx, lat_idx, lon_idx, alt_m). The adapter does the lookup of
    lat_idx and lon_idx for ERA5 internally.
    """
    if is_gfs(reader):
        return reader.getNearestAlt(hour_idx, lat, lon, alt_m)
    lat_i = nearest_lat(reader, lat)
    lon_i = nearest_lon(reader, lon)
    return reader.getNearestAltbyIndex(int(hour_idx), lat_i, lon_i, alt_m)


def wind_interp(reader, coord):
    """Run the composite altitude+time wind interpolation.

    Both readers expose ``wind_alt_Interpolate2`` but with different
    signatures:
        GFS:  wind_alt_Interpolate2(coord)         -> [u, v, u_lf, v_lf]
        ERA5: wind_alt_Interpolate2(alt_m, diff_time, lat_idx, lon_idx)
                                                    -> [u, v, u_lf, v_lf]

    The adapter takes the GFS-style dict and translates as needed.
    """
    if is_gfs(reader):
        return reader.wind_alt_Interpolate2(coord)

    diff_time = coord["timestamp"] - reader.model_start_datetime
    lat_idx = nearest_lat(reader, coord["lat"])
    lon_idx = nearest_lon(reader, coord["lon"])
    return reader.wind_alt_Interpolate2(
        coord["alt"], diff_time, lat_idx, lon_idx
    )


def get_new_coord(reader, coord, dt):
    """Top-level propagation. Same name and signature on both readers."""
    return reader.getNewCoord(coord, dt)
