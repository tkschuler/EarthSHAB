"""Analytical-truth functions for synthetic forecast scenarios.

For each scenario in make_dummy_forecasts.SCENARIOS, this module provides
a ``truth(scenario, t_dt, alt_m, lat, lon)`` function returning the analytical
(u, v) that a correct reader+interpolator should produce at that query.

The truth functions assume the synthetic grid defined in make_dummy_forecasts
(GRID_SHAPE, LAT_VALUES_ASC, LON_VALUES, LEVELS_HPA_GFS, ALT_M_GFS_ORDER,
TIMES_DT). Off-grid queries return the linearly-interpolated truth where the
scenario field is linear in the relevant index; for non-linear scenarios
(e.g. bearing_wrap) the off-grid truth is computed by the same bearing/speed
linear interpolation that ``linear_neighbors`` performs.

For scenarios that vary with a particular index, fractional indices are
computed via np.interp against the stored coordinate arrays. This gives
the analytical truth for any query inside the grid.
"""

from __future__ import annotations

import math
from datetime import datetime

import numpy as np

from tests.tools.make_dummy_forecasts import (
    ALL_STATIC_U,
    ALL_STATIC_V,
    ALT_M_GFS_ORDER,
    LAT_VALUES_ASC,
    LON_VALUES,
    TIMES_DT,
    TIME_STEP_HOURS,
)


def _frac_time_idx(t_dt: datetime) -> float:
    """Return fractional index into TIMES_DT for the given datetime."""
    hrs = (t_dt - TIMES_DT[0]).total_seconds() / 3600.0
    return hrs / TIME_STEP_HOURS


def _frac_lat_idx_gfs_order(lat: float) -> float:
    """Fractional lat index in GFS ascending-lat ordering."""
    return float(np.interp(lat, LAT_VALUES_ASC, np.arange(len(LAT_VALUES_ASC))))


def _frac_lon_idx(lon: float) -> float:
    """Fractional lon index. Synthetic lon frame is [0, 2] in 0-360 style;
    queries are expected within that range (GFS reader applies % 360 first).
    """
    lon_norm = float(lon) % 360.0
    # If query lon falls in [0, 360) but outside [0, 2] synthetic range,
    # clip to nearest endpoint — matches np.interp clamping.
    return float(np.interp(lon_norm, LON_VALUES, np.arange(len(LON_VALUES))))


def _frac_level_idx_by_alt(alt_m: float) -> float:
    """Fractional level index given altitude in meters.

    ALT_M_GFS_ORDER is ascending (index 0 = lowest altitude), matching
    the canonical in-memory layout used by the generator. np.interp
    needs ascending xp, which is what we have.
    """
    idx = np.arange(len(ALT_M_GFS_ORDER))
    return float(np.interp(alt_m, ALT_M_GFS_ORDER, idx))


# ---------------------------------------------------------------------------
# Per-scenario truth
# ---------------------------------------------------------------------------

def _truth_all_static(t_dt, alt_m, lat, lon):
    return (ALL_STATIC_U, ALL_STATIC_V)


def _truth_static_by_level(t_dt, alt_m, lat, lon):
    return (_frac_level_idx_by_alt(alt_m), 0.0)


def _truth_linear_ramp_all_dims(t_dt, alt_m, lat, lon):
    u = (
        _frac_time_idx(t_dt)
        + _frac_level_idx_by_alt(alt_m)
        + _frac_lat_idx_gfs_order(lat)
        + _frac_lon_idx(lon)
    )
    return (u, 0.0)


def _truth_time_ramp_only(t_dt, alt_m, lat, lon):
    return (_frac_time_idx(t_dt), 0.0)


def _truth_lat_ramp_only(t_dt, alt_m, lat, lon):
    return (_frac_lat_idx_gfs_order(lat), 0.0)


def _truth_lon_ramp_only(t_dt, alt_m, lat, lon):
    return (_frac_lon_idx(lon), 0.0)


def _truth_altitude_ramp_only(t_dt, alt_m, lat, lon):
    return (_frac_level_idx_by_alt(alt_m), 0.0)


def _truth_bearing_wrap(t_dt, alt_m, lat, lon):
    """Analytical bearing-wrap-aware linear interpolation between adjacent
    levels, matching what wind_alt_Interpolate2's 'linear_neighbors' branch
    computes internally. Speed is constant 10 m/s; bearing rotates with k.
    """
    speed = 10.0
    frac_k = _frac_level_idx_by_alt(alt_m)  # ascending: 0 = lowest altitude
    k0 = int(np.floor(frac_k))
    k1 = k0 + 1
    n_lev = len(ALT_M_GFS_ORDER)
    if k0 < 0:
        k0, k1 = 0, 1
    if k1 >= n_lev:
        k0, k1 = n_lev - 2, n_lev - 1
    alpha = (frac_k - k0)  # 0 at k0, 1 at k1
    bearing0 = (-90.0 + 20.0 * k0) % 360.0
    bearing1 = (-90.0 + 20.0 * k1) % 360.0
    # Wrap correction matching interpolateBearing.
    a1, a2 = bearing0, bearing1
    if abs(a2 - a1) > 180.0:
        if a2 > a1:
            a1 += 360.0
        else:
            a2 += 360.0
    bearing_interp = ((1.0 - alpha) * a1 + alpha * a2) % 360.0
    u = speed * math.cos(math.radians(bearing_interp))
    v = speed * math.sin(math.radians(bearing_interp))
    return (u, v)


def _truth_nan_at_top(t_dt, alt_m, lat, lon):
    """After fill_missing_data, the masked top levels are linearly extrapolated
    from the valid range, then np.interp clamps anything past the top valid
    level to the last valid value. For most query altitudes within the valid
    column, truth is identical to static_by_level (u = frac_level_idx).

    Queries above the original top valid altitude resolve to whatever
    fill_missing_data extrapolated to — by inspection of GFS.fill_missing_data,
    that's a linear extrapolation in INDEX SPACE over the level dimension,
    i.e. the masked indices get u-values continuing the level_idx ramp.
    So truth here is still u = frac_level_idx_by_alt(alt_m), provided the
    query altitude maps to a fractional level index in [0, n_lev-1].

    Returns (u, 0.0). Used only for sanity-checking inside the valid column;
    tests of fill_missing_data itself use the raw masked array directly.
    """
    return (_frac_level_idx_by_alt(alt_m), 0.0)


def _truth_geopotential_vs_height(t_dt, alt_m, lat, lon):
    return (7.5, -3.25)


TRUTH_FNS = {
    "all_static": _truth_all_static,
    "static_by_level": _truth_static_by_level,
    "linear_ramp_all_dims": _truth_linear_ramp_all_dims,
    "time_ramp_only": _truth_time_ramp_only,
    "lat_ramp_only": _truth_lat_ramp_only,
    "lon_ramp_only": _truth_lon_ramp_only,
    "altitude_ramp_only": _truth_altitude_ramp_only,
    "bearing_wrap": _truth_bearing_wrap,
    "nan_at_top": _truth_nan_at_top,
    "geopotential_vs_height": _truth_geopotential_vs_height,
}


def truth(scenario: str, t_dt: datetime, alt_m: float, lat: float, lon: float):
    """Analytical (u, v) for the named scenario at the given query.

    Assumes the synthetic grid from make_dummy_forecasts. Returns a tuple of
    floats. Off-grid queries are interpolated linearly per the scenario's
    construction; bearing_wrap uses bearing-aware linear interp.
    """
    if scenario not in TRUTH_FNS:
        raise ValueError(f"Unknown scenario {scenario!r}")
    return TRUTH_FNS[scenario](t_dt, alt_m, lat, lon)
