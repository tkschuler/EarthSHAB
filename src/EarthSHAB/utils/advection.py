"""Trajectory advection: advance (lat, lon) under a horizontal wind over dt.

Two models, both vectorized over members:

  geodesic      -- WGS84 ellipsoid step (geographiclib), re-anchored at the
                   balloon each step. Loops member-wise because geographiclib
                   is not array-aware.
  tangent_plane -- spherical tangent-plane step about a FIXED anchor (the launch
                   coordinate); closed-form array ops since cos(central_lat) is
                   constant.

Single source of truth shared by ``Forecast._advect`` (the scalar/batch reader)
and the vectorized Monte Carlo runner, so an advection change lands in exactly
one place. Kept byte-identical between EarthSHAB and HAB-COM.

Each function returns ``(lat_new, lon_new, x_disp, y_disp)``. The displacement
fields mirror the MC runner's diagnostic convention (they differ by model):

  geodesic      -> (u*dt, v*dt)    per-step incremental displacement (m)
  tangent_plane -> (x_new, y_new)  cumulative displacement from the anchor (m)

``M_PER_DEG`` is the single meters-per-degree constant for the tangent plane.
Forward and inverse conversions MUST use the same value or the per-step
lat/lon<->meters round-trip about the fixed anchor stops being an identity and
drifts the trajectory toward the anchor (the 111000/111320 bug).
"""
import math

import numpy as np

# meters per degree of latitude (and of longitude at the equator).
M_PER_DEG = 111320.0


def step_geodesic(lats, lons, u, v, dt, geod):
    """WGS84 ellipsoid step, re-anchored at each member's position."""
    M = len(lats)
    lat_new = np.empty(M, dtype=np.float64)
    lon_new = np.empty(M, dtype=np.float64)
    for m in range(M):
        bearing = 90.0 - math.degrees(math.atan2(v[m], u[m]))
        d = math.hypot(u[m], v[m]) * dt
        g = geod.Direct(float(lats[m]), float(lons[m]), bearing, d)
        lat_new[m] = g['lat2']
        lon_new[m] = g['lon2']
    return lat_new, lon_new, u * dt, v * dt


def step_tangent_plane(lats, lons, u, v, dt, central_lat, central_lon, cos_central):
    """Spherical tangent-plane step about the fixed launch anchor."""
    cx = (lons - central_lon) * M_PER_DEG * cos_central
    cy = (lats - central_lat) * M_PER_DEG
    x_new = cx + u * dt
    y_new = cy + v * dt
    lat_new = np.clip(central_lat + y_new / M_PER_DEG, -89.999999, 89.999999)
    lon_new = ((central_lon + x_new / (M_PER_DEG * cos_central) + 180.0) % 360.0) - 180.0
    return lat_new, lon_new, x_new, y_new


def advect(method, lats, lons, u, v, dt, *, geod=None,
           central_lat=None, central_lon=None, cos_central=None):
    """Dispatch the advection step for ``method``.

    :param method: 'geodesic' or 'tangent_plane'
    :param lats/lons: 1-D arrays (length M) of current positions (deg)
    :param u/v: 1-D arrays of east/north wind (m/s)
    :param dt: integration interval (s)
    :param geod: geographiclib Geodesic instance (required for 'geodesic')
    :param central_lat/central_lon/cos_central: fixed anchor (required for
        'tangent_plane')
    :returns: (lat_new, lon_new, x_disp, y_disp); see module docstring.
    """
    if method == 'tangent_plane':
        return step_tangent_plane(lats, lons, u, v, dt,
                                  central_lat, central_lon, cos_central)
    if method == 'geodesic':
        return step_geodesic(lats, lons, u, v, dt, geod)
    raise ValueError(
        f"unknown advection: {method!r}. Expected 'geodesic' or 'tangent_plane'.")
