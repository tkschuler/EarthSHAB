"""Shared wind-interpolation primitives used by both the scalar Forecast.py
reader and the vectorized Monte Carlo runner.

The point of this module is *exactness*: the altitude/time interpolation math
here is a faithful copy of Forecast.py's original inline implementations, so the
fast vectorized runner and the canonical Forecast.py engine compute the same
wind for the same method. The only intended differences between the two engines
are (a) float32 vs float64 storage of the forecast cube and (b) the horizontal
sampling stencil (nearest cell vs bilinear), both of which are method choices.

Four wind-interpolation methods are supported, each = (horizontal stencil) +
(altitude/time scheme):

    linear_neighbors : nearest lat/lon cell + bearing/speed altitude & time interp
    linear_full      : nearest lat/lon cell + linear u/v in altitude, linear in time
    spline_full      : nearest lat/lon cell + cubic-spline u/v in altitude, linear in time
    bilinear         : bilinear lat/lon     + linear u/v in altitude, linear in time

All ``apply_*`` functions take pre-sampled vertical columns (one at the lower
forecast hour t0, one at the upper t1) plus a fractional ``hour_index`` and the
target altitude, and return scalar ``(u, v)`` wind components (east, north).
"""
import numpy as np
from scipy.interpolate import CubicSpline


# ---------------------------------------------------------------------------
# 1-D helpers
# ---------------------------------------------------------------------------

def fill_missing_data(data):
    """Linearly interpolate over missing entries in a 1-D profile column.

    Verbatim copy of Forecast.fill_missing_data so both engines treat NaN /
    masked levels identically.
    """
    data = data.filled(np.nan) if np.ma.isMaskedArray(data) else np.asarray(data, dtype=float)
    nans, x = np.isnan(data), lambda z: z.nonzero()[0]
    if nans.any() and (~nans).any():
        data[nans] = np.interp(x(nans), x(~nans), data[~nans])
    return data


def find_bracket_extrap(arr, x):
    """Bracketing indices + fractional weight for x in a monotone-increasing arr.

    Matches utils.vectorization._find_bracketing_extrap: inside the bounds the
    weight is in [0, 1]; outside, the first/last segment is reused so the weight
    goes <0 or >1 (linear extrapolation). Returns (i0, i1, w) with
    value = v[i0] + w*(v[i1]-v[i0]).
    """
    n = len(arr)
    if n < 2:
        return 0, 0, 0.0
    if x <= arr[0]:
        i0, i1 = 0, 1
    elif x >= arr[-1]:
        i0, i1 = n - 2, n - 1
    else:
        i1 = int(np.searchsorted(arr, x, side="right"))
        i0 = i1 - 1
    w = (x - arr[i0]) / (arr[i1] - arr[i0])
    return i0, i1, float(w)


def spline_uv(alt_m, h_ascending, u_ascending, v_ascending):
    """CubicSpline on u and v across the altitude profile, with np.interp fallback.

    Verbatim copy of Forecast._spline_uv (kept here to avoid a circular import).
    """
    h_arr = np.asarray(h_ascending)
    u_arr = np.asarray(u_ascending)
    v_arr = np.asarray(v_ascending)
    keep = np.concatenate(([True], np.diff(h_arr) > 0))
    h_s, u_s, v_s = h_arr[keep], u_arr[keep], v_arr[keep]

    h_lo, h_hi = h_s[0], h_s[-1]
    if alt_m < h_lo or alt_m > h_hi or len(h_s) < 4:
        return (
            float(np.interp(alt_m, h_s, u_s)),
            float(np.interp(alt_m, h_s, v_s)),
        )
    cs_u = CubicSpline(h_s, u_s, extrapolate=False)
    cs_v = CubicSpline(h_s, v_s, extrapolate=False)
    return float(cs_u(alt_m)), float(cs_v(alt_m))


# ---------------------------------------------------------------------------
# Bearing/speed helpers (linear_neighbors) — verbatim from Forecast.py
# ---------------------------------------------------------------------------

def wind_vector_to_bearing(u, v):
    bearing = np.arctan2(v, u)
    speed = np.power((np.power(u, 2) + np.power(v, 2)), .5)
    return bearing, speed


def _two_nearest_alt_idxs(h, alt_m):
    """Two enclosing altitude indices (Forecast.get2NearestAltIdxs)."""
    h_nearest = int(np.argmin(np.abs(np.asarray(h) - alt_m)))
    if alt_m > h[h_nearest]:
        return h_nearest, h_nearest + 1
    return h_nearest - 1, h_nearest


def interpolate_bearing(h, u, v, alt_m):
    """Linear interp between adjacent pressure levels on bearing+speed, with
    0/360 axis-crossover correction. Verbatim from Forecast.interpolateBearing.
    """
    h_idx0, h_idx1 = _two_nearest_alt_idxs(h, alt_m)
    if h_idx0 == -1:
        h_idx0 = 0
    if h_idx1 == 0:
        h_idx1 = -1

    bearing0, speed0 = wind_vector_to_bearing(u[h_idx0], v[h_idx0])
    bearing1, speed1 = wind_vector_to_bearing(u[h_idx1], v[h_idx1])

    angle1 = np.degrees(bearing0) % 360
    angle2 = np.degrees(bearing1) % 360
    if abs(angle2 - angle1) > 180:
        if angle2 > angle1:
            angle1 += 360
        else:
            angle2 += 360

    interp_speed = np.interp(alt_m, [h[h_idx0], h[h_idx1]], [speed0, speed1])
    interp_dir_deg = np.interp(alt_m, [h[h_idx0], h[h_idx1]], [angle1, angle2]) % 360
    return interp_dir_deg, interp_speed


def interpolate_bearing_time(bearing0, speed0, bearing1, speed1, hour_index):
    """Interpolate bearing+speed in time. Verbatim from Forecast.interpolateBearingTime."""
    angle1 = bearing0 % 360
    angle2 = bearing1 % 360
    if abs(angle2 - angle1) > 180:
        if angle2 > angle1:
            angle1 += 360
        else:
            angle2 += 360
    fp = [int(hour_index), int(hour_index) + 1]
    interp_speed = np.interp(hour_index, fp, [speed0, speed1])
    interp_dir_deg = np.interp(hour_index, fp, [angle1, angle2]) % 360
    return interp_dir_deg, interp_speed


# ---------------------------------------------------------------------------
# Altitude + time appliers (operate on already-sampled columns)
# ---------------------------------------------------------------------------

def apply_linear_full(h0, u0, v0, h1, u1, v1, alt_m, hour_index):
    """Linear u/v in altitude on each forecast hour, then linear in time."""
    u0_lf = np.interp(alt_m, h0, u0)
    v0_lf = np.interp(alt_m, h0, v0)
    u1_lf = np.interp(alt_m, h1, u1)
    v1_lf = np.interp(alt_m, h1, v1)
    fp = [int(hour_index), int(hour_index) + 1]
    u = np.interp(hour_index, fp, [u0_lf, u1_lf])
    v = np.interp(hour_index, fp, [v0_lf, v1_lf])
    return float(u), float(v)


def apply_linear_neighbors(h0, u0, v0, h1, u1, v1, alt_m, hour_index):
    """Bearing/speed altitude interp at each hour, then bearing/speed in time."""
    bearing_t0, speed_t0 = interpolate_bearing(h0, u0, v0, alt_m)
    bearing_t1, speed_t1 = interpolate_bearing(h1, u1, v1, alt_m)
    bearing_i, speed_i = interpolate_bearing_time(
        bearing_t0, speed_t0, bearing_t1, speed_t1, hour_index)
    u = speed_i * np.cos(np.radians(bearing_i))
    v = speed_i * np.sin(np.radians(bearing_i))
    return float(u), float(v)


def apply_spline_full(h0, u0, v0, h1, u1, v1, alt_m, hour_index):
    """Cubic-spline u/v in altitude on each forecast hour, then linear in time."""
    u0_sp, v0_sp = spline_uv(alt_m, h0, u0, v0)
    u1_sp, v1_sp = spline_uv(alt_m, h1, u1, v1)
    fp = [int(hour_index), int(hour_index) + 1]
    u = np.interp(hour_index, fp, [u0_sp, u1_sp])
    v = np.interp(hour_index, fp, [v0_sp, v1_sp])
    return float(u), float(v)


_APPLIERS = {
    "linear_full": apply_linear_full,
    "linear_neighbors": apply_linear_neighbors,
    "spline_full": apply_spline_full,
    # 'bilinear' uses the linear_full altitude/time scheme; the difference is the
    # horizontal stencil, applied when sampling the columns (see sample_columns).
    "bilinear": apply_linear_full,
}


def apply_method(method, h0, u0, v0, h1, u1, v1, alt_m, hour_index):
    """Dispatch the vertical-altitude + time applier for ``method``.

    Used by backends that fetch the (h, u, v) columns themselves (e.g. the
    xarray backend) and only need the shared alt/time math.
    """
    return _APPLIERS[method](h0, u0, v0, h1, u1, v1, alt_m, hour_index)


# ---------------------------------------------------------------------------
# Column samplers (horizontal stencil) — operate on ascending lat/lon cubes
# ---------------------------------------------------------------------------

def nearest_columns(z_t, u_t, v_t, lat_vals, lon_vals, lat, lon, u_bias=None, v_bias=None):
    """Nearest lat/lon cell vertical columns from a single time slice.

    z_t/u_t/v_t are [level, lat, lon]; lat_vals/lon_vals are 1-D coordinate
    arrays. ``argmin(abs(...))`` reproduces Forecast.closestIdx, so the same
    physical cell is selected regardless of ascending/descending coordinate
    order.
    """
    lat_i = int(np.argmin(np.abs(np.asarray(lat_vals) - lat)))
    lon_i = int(np.argmin(np.abs(np.asarray(lon_vals) - lon)))
    h = fill_missing_data(z_t[:, lat_i, lon_i])
    u = fill_missing_data(u_t[:, lat_i, lon_i])
    v = fill_missing_data(v_t[:, lat_i, lon_i])
    if u_bias is not None:
        u = u + u_bias
    if v_bias is not None:
        v = v + v_bias
    return h, u, v


def bilinear_columns(z_t, u_t, v_t, lat_vals, lon_vals, lat, lon, u_bias=None, v_bias=None):
    """Bilinear (lat, lon) vertical columns from a single time slice.

    lat_vals/lon_vals MUST be ascending (binary search + linear extrapolation
    outside bounds). Mirrors utils.vectorization.interp4d_alt_batch_match_xr's
    per-level bilinear blend.
    """
    la0, la1, a_la = find_bracket_extrap(lat_vals, lat)
    lo0, lo1, a_lo = find_bracket_extrap(lon_vals, lon)

    def _bilerp(arr):
        a = arr[:, la0, lo0] + a_lo * (arr[:, la0, lo1] - arr[:, la0, lo0])
        b = arr[:, la1, lo0] + a_lo * (arr[:, la1, lo1] - arr[:, la1, lo0])
        return a + a_la * (b - a)

    h = fill_missing_data(_bilerp(z_t))
    u = fill_missing_data(_bilerp(u_t))
    v = fill_missing_data(_bilerp(v_t))
    if u_bias is not None:
        u = u + u_bias
    if v_bias is not None:
        v = v + v_bias
    return h, u, v


def sample_uv(method, lat_vals, lon_vals, z_np, u_np, v_np,
              i_hr, hour_index, lat, lon, alt_m,
              u_bias=None, v_bias=None):
    """Sample (u, v) for one member at (lat, lon, alt) and fractional hour_index.

    z_np/u_np/v_np are [time, level, lat, lon] cubes (z already in meters). For
    'bilinear' lat_vals/lon_vals must be ascending. The two bracketing forecast
    hours are i_hr and i_hr+1 (matching Forecast.py's int(hour_index) scheme).
    """
    sampler = bilinear_columns if method == "bilinear" else nearest_columns
    h0, u0, v0 = sampler(z_np[i_hr], u_np[i_hr], v_np[i_hr],
                         lat_vals, lon_vals, lat, lon, u_bias, v_bias)
    h1, u1, v1 = sampler(z_np[i_hr + 1], u_np[i_hr + 1], v_np[i_hr + 1],
                         lat_vals, lon_vals, lat, lon, u_bias, v_bias)
    return _APPLIERS[method](h0, u0, v0, h1, u1, v1, alt_m, hour_index)
