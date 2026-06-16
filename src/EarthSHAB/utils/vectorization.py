import math

import numpy as np
import numba as nb

@nb.njit(cache=True)
def _find_bracketing_pad_backfill(time_unix_s: np.ndarray, t_s: np.int64):
    """
    REPLACES: xarray time selection in xr_interp_4d_alt:
        t0 = ds.time.sel(time=timestamp, method="pad")
        t1 = ds.time.sel(time=timestamp, method="backfill")
        alpha = (tv - t0) / (t1 - t0)  (0 if equal)

    Behavior:
      - If t is outside time range: clamp to edge with t0=t1 and alpha=0
      - Otherwise: choose nearest lower (pad) and nearest higher (backfill)
    """
    n = time_unix_s.size
    if t_s <= time_unix_s[0]:
        return 0, 0, 0.0
    if t_s >= time_unix_s[n - 1]:
        return n - 1, n - 1, 0.0

    # Find first index where time >= t_s (standard lower_bound)
    left = 0
    right = n - 1
    while left <= right:
        mid = (left + right) // 2
        if time_unix_s[mid] < t_s:
            left = mid + 1
        else:
            right = mid - 1

    t1 = left
    t0 = t1 - 1
    denom = time_unix_s[t1] - time_unix_s[t0]
    a = 0.0 if denom == 0 else (t_s - time_unix_s[t0]) / denom
    return t0, t1, a


@nb.njit(cache=True)
def _find_bracketing_extrap(arr: np.ndarray, x: float):
    """
    REPLACES: xarray .interp(..., kwargs={"fill_value":"extrapolate"}) behavior for 1D coords.

    - For x within bounds: return bracketing indices and fractional weight a in [0,1]
    - For x outside bounds: return first/last segment and allow a < 0 or a > 1 (linear extrapolation)

    Assumes arr is monotone increasing.
    """
    n = arr.size
    if n < 2:
        return 0, 0, 0.0

    if x <= arr[0]:
        a0 = arr[0]
        a1 = arr[1]
        a = (x - a0) / (a1 - a0)
        return 0, 1, a

    if x >= arr[n - 1]:
        a0 = arr[n - 2]
        a1 = arr[n - 1]
        a = (x - a0) / (a1 - a0)
        return n - 2, n - 1, a

    left = 0
    right = n - 1
    while right - left > 1:
        mid = (left + right) // 2
        if arr[mid] <= x:
            left = mid
        else:
            right = mid
    a0 = arr[left]
    a1 = arr[right]
    a = (x - a0) / (a1 - a0)
    return left, right, a


@nb.njit(cache=True)
def _interp_vertical_like_user(zcol: np.ndarray, ucol: np.ndarray, vcol: np.ndarray, alt_m: float):
    """
    REPLACES: your Python interpolate_wind() called from xr_interp_4d_alt()

    Your original interpolate_wind():
      - If z is descending, reverse z/u/v
      - Clamp alt to [zmin, zmax]
      - np.interp in z

    This numba version performs the same logic (including clamp),
    without allocating reversed arrays if possible.
    """
    L = zcol.size

    # If z descending with height, treat a virtual reversed view
    if zcol[0] > zcol[L - 1]:
        zmin = zcol[L - 1]
        zmax = zcol[0]
        alt = alt_m
        if alt < zmin:
            alt = zmin
        elif alt > zmax:
            alt = zmax

        # binary search in virtual increasing array: z_inc[i] = zcol[L-1-i]
        if alt <= zcol[L - 1]:
            k0 = k1 = L - 1
            a = 0.0
        elif alt >= zcol[0]:
            k0 = k1 = 0
            a = 0.0
        else:
            left = 0
            right = L - 1
            while right - left > 1:
                mid = (left + right) // 2
                if zcol[L - 1 - mid] <= alt:
                    left = mid
                else:
                    right = mid
            k0v = left
            k1v = right
            z0 = zcol[L - 1 - k0v]
            z1 = zcol[L - 1 - k1v]
            a = 0.0 if z1 == z0 else (alt - z0) / (z1 - z0)
            k0 = L - 1 - k0v
            k1 = L - 1 - k1v

        u = ucol[k0] + a * (ucol[k1] - ucol[k0])
        v = vcol[k0] + a * (vcol[k1] - vcol[k0])
        return u, v

    # z ascending
    zmin = zcol[0]
    zmax = zcol[L - 1]
    alt = alt_m
    if alt < zmin:
        alt = zmin
    elif alt > zmax:
        alt = zmax

    if alt <= zcol[0]:
        k0 = k1 = 0
        a = 0.0
    elif alt >= zcol[L - 1]:
        k0 = k1 = L - 1
        a = 0.0
    else:
        left = 0
        right = L - 1
        while right - left > 1:
            mid = (left + right) // 2
            if zcol[mid] <= alt:
                left = mid
            else:
                right = mid
        k0 = left
        k1 = right
        z0 = zcol[k0]
        z1 = zcol[k1]
        a = 0.0 if z1 == z0 else (alt - z0) / (z1 - z0)

    u = ucol[k0] + a * (ucol[k1] - ucol[k0])
    v = vcol[k0] + a * (vcol[k1] - vcol[k0])
    return u, v


@nb.njit(cache=True)
def interp4d_alt_batch_match_xr(
    time_unix_s: np.ndarray,
    lat_vals: np.ndarray,
    lon_vals: np.ndarray,
    z_np: np.ndarray,
    u_np: np.ndarray,
    v_np: np.ndarray,
    t_s: np.int64,
    lats: np.ndarray,
    lons: np.ndarray,
    alts: np.ndarray,
    u_bias: np.ndarray,
    v_bias: np.ndarray,
):
    """
    REPLACES: xr_interp_4d_alt() (the slow xarray-based interpolator)

    Matches your original "altitude-first then time" conceptual pipeline:
      1) Time bracketing: t0 pad, t1 backfill
      2) At each time:
           bilinear interp in (lat, lon) to get vertical profiles z(level), u(level), v(level)
           with fill_value="extrapolate" -> linear extrap outside bounds
      3) Vertical interp to altitude (like interpolate_wind):
           reverse if z descending; clamp altitude to z-range
      4) Linear blend in time between t0 and t1 winds

    Returns:
      u_out[m], v_out[m] for each member m
    """
    M = lats.size
    T, L, La, Lo = u_np.shape
    u_out = np.empty(M, np.float32)
    v_out = np.empty(M, np.float32)

    t0, t1, at = _find_bracketing_pad_backfill(time_unix_s, t_s)

    for m in range(M):
        lat = lats[m]
        lon = lons[m]
        alt = alts[m]

        la0, la1, a_la = _find_bracketing_extrap(lat_vals, lat)
        lo0, lo1, a_lo = _find_bracketing_extrap(lon_vals, lon)

        # Build vertical columns at t0 and t1 for this member.
        # (L is usually small-ish, this is fine; Numba makes it fast.)
        zcol0 = np.empty(L, np.float32)
        ucol0 = np.empty(L, np.float32)
        vcol0 = np.empty(L, np.float32)
        zcol1 = np.empty(L, np.float32)
        ucol1 = np.empty(L, np.float32)
        vcol1 = np.empty(L, np.float32)

        for lev in range(L):
            # ---- t0 bilinear ----
            z00 = z_np[t0, lev, la0, lo0]; z01 = z_np[t0, lev, la0, lo1]
            z10 = z_np[t0, lev, la1, lo0]; z11 = z_np[t0, lev, la1, lo1]
            u00 = u_np[t0, lev, la0, lo0]; u01 = u_np[t0, lev, la0, lo1]
            u10 = u_np[t0, lev, la1, lo0]; u11 = u_np[t0, lev, la1, lo1]
            v00 = v_np[t0, lev, la0, lo0]; v01 = v_np[t0, lev, la0, lo1]
            v10 = v_np[t0, lev, la1, lo0]; v11 = v_np[t0, lev, la1, lo1]

            z0a = z00 + a_lo * (z01 - z00)
            z0b = z10 + a_lo * (z11 - z10)
            zcol0[lev] = z0a + a_la * (z0b - z0a)

            u0a = u00 + a_lo * (u01 - u00)
            u0b = u10 + a_lo * (u11 - u10)
            ucol0[lev] = (u0a + a_la * (u0b - u0a)) + u_bias[m, lev]

            v0a = v00 + a_lo * (v01 - v00)
            v0b = v10 + a_lo * (v11 - v10)
            vcol0[lev] = (v0a + a_la * (v0b - v0a)) + v_bias[m, lev]

            # ---- t1 bilinear ----
            z00 = z_np[t1, lev, la0, lo0]; z01 = z_np[t1, lev, la0, lo1]
            z10 = z_np[t1, lev, la1, lo0]; z11 = z_np[t1, lev, la1, lo1]
            u00 = u_np[t1, lev, la0, lo0]; u01 = u_np[t1, lev, la0, lo1]
            u10 = u_np[t1, lev, la1, lo0]; u11 = u_np[t1, lev, la1, lo1]
            v00 = v_np[t1, lev, la0, lo0]; v01 = v_np[t1, lev, la0, lo1]
            v10 = v_np[t1, lev, la1, lo0]; v11 = v_np[t1, lev, la1, lo1]

            z1a = z00 + a_lo * (z01 - z00)
            z1b = z10 + a_lo * (z11 - z10)
            zcol1[lev] = z1a + a_la * (z1b - z1a)

            u1a = u00 + a_lo * (u01 - u00)
            u1b = u10 + a_lo * (u11 - u10)
            ucol1[lev] = (u1a + a_la * (u1b - u1a)) + u_bias[m, lev]

            v1a = v00 + a_lo * (v01 - v00)
            v1b = v10 + a_lo * (v11 - v10)
            vcol1[lev] = (v1a + a_la * (v1b - v1a)) + v_bias[m, lev]

        # vertical interpolation to altitude
        u0, v0 = _interp_vertical_like_user(zcol0, ucol0, vcol0, alt)
        u1, v1 = _interp_vertical_like_user(zcol1, ucol1, vcol1, alt)

        # time blend
        u_out[m] = u0 + at * (u1 - u0)
        v_out[m] = v0 + at * (v1 - v0)

    return u_out, v_out


# ===========================================================================
# Numba batch implementations of the nearest-cell wind methods.
#
# These reproduce utils.wind_interp's pure-numpy math (nearest lat/lon cell,
# fill_missing_data, then the per-method altitude/time scheme) so the fast
# vectorized runner produces the same trajectories it did on the numpy path —
# just ~40-70x faster on the wind lookup. 'spline_full' is intentionally NOT
# here: scipy.CubicSpline can't run under numba, so it stays on the numpy path.
#
# Time bracketing matches the runner: i_hr = clamp(floor(hour_index), 0, T-2),
# linear blend with weight (hour_index - i_hr). Bias (wind noise) is added
# AFTER fill, per level, exactly like wind_interp.nearest_columns.
# ===========================================================================

@nb.njit(cache=True)
def _nearest_idx(arr, x):
    """argmin(abs(arr - x)), first occurrence (matches Forecast.closestIdx)."""
    best = 0
    bestd = abs(arr[0] - x)
    for i in range(1, arr.size):
        d = abs(arr[i] - x)
        if d < bestd:
            bestd = d
            best = i
    return best


@nb.njit(cache=True)
def _fill_missing_1d(col):
    """In-place-on-copy linear fill over NaNs by index (matches np.interp fill)."""
    n = col.size
    out = col.copy()
    first_valid = -1
    last_valid = -1
    for i in range(n):
        if not np.isnan(out[i]):
            if first_valid == -1:
                first_valid = i
            last_valid = i
    if first_valid == -1:
        return out  # all NaN -> leave as is (matches: no valid points to interp from)
    for i in range(0, first_valid):       # leading NaNs -> first valid (np.interp clamp)
        out[i] = out[first_valid]
    for i in range(last_valid + 1, n):    # trailing NaNs -> last valid
        out[i] = out[last_valid]
    i = first_valid
    while i <= last_valid:                 # interior NaNs -> linear in index
        if np.isnan(out[i]):
            j = i + 1
            while np.isnan(out[j]):
                j += 1
            lo = i - 1
            for k in range(i, j):
                w = (k - lo) / (j - lo)
                out[k] = out[lo] + w * (out[j] - out[lo])
            i = j + 1
        else:
            i += 1
    return out


@nb.njit(cache=True)
def _interp_clamped(x, xp, fp):
    """np.interp for monotone-increasing xp (clamps to endpoints)."""
    n = xp.size
    if x <= xp[0]:
        return fp[0]
    if x >= xp[n - 1]:
        return fp[n - 1]
    lo = 0
    hi = n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xp[mid] <= x:
            lo = mid
        else:
            hi = mid
    return fp[lo] + (x - xp[lo]) / (xp[hi] - xp[lo]) * (fp[hi] - fp[lo])


@nb.njit(cache=True)
def _interp2(x, x0, x1, y0, y1):
    """np.interp on a 2-point increasing table (clamps)."""
    if x <= x0:
        return y0
    if x >= x1:
        return y1
    return y0 + (x - x0) / (x1 - x0) * (y1 - y0)


@nb.njit(cache=True)
def _mod360(x):
    """x mod 360 in [0, 360) — matches numpy/Python % for positive divisor."""
    return x - 360.0 * math.floor(x / 360.0)


@nb.njit(cache=True)
def _time_bracket(time_unix_s, t_s):
    """(i_hr, hour_index) matching the runner's floor scheme, clamped to [0, T-2]."""
    T = time_unix_s.size
    res = time_unix_s[1] - time_unix_s[0] if T >= 2 else 1.0
    hour_index = (t_s - time_unix_s[0]) / res
    i_hr = int(math.floor(hour_index))
    if i_hr > T - 2:
        i_hr = T - 2
    if i_hr < 0:
        i_hr = 0
    return i_hr, hour_index


@nb.njit(cache=True)
def _nearest_cols(z_np, u_np, v_np, i_hr, lat_i, lon_i, ubias_m, vbias_m):
    """Filled vertical columns at t0=i_hr and t1=i_hr+1 for one cell, +bias."""
    L = u_np.shape[1]
    h0 = np.empty(L, np.float64); u0 = np.empty(L, np.float64); v0 = np.empty(L, np.float64)
    h1 = np.empty(L, np.float64); u1 = np.empty(L, np.float64); v1 = np.empty(L, np.float64)
    for lev in range(L):
        h0[lev] = z_np[i_hr, lev, lat_i, lon_i]
        u0[lev] = u_np[i_hr, lev, lat_i, lon_i]
        v0[lev] = v_np[i_hr, lev, lat_i, lon_i]
        h1[lev] = z_np[i_hr + 1, lev, lat_i, lon_i]
        u1[lev] = u_np[i_hr + 1, lev, lat_i, lon_i]
        v1[lev] = v_np[i_hr + 1, lev, lat_i, lon_i]
    h0 = _fill_missing_1d(h0); u0 = _fill_missing_1d(u0); v0 = _fill_missing_1d(v0)
    h1 = _fill_missing_1d(h1); u1 = _fill_missing_1d(u1); v1 = _fill_missing_1d(v1)
    for lev in range(L):          # bias added AFTER fill (matches wind_interp)
        u0[lev] += ubias_m[lev]; v0[lev] += vbias_m[lev]
        u1[lev] += ubias_m[lev]; v1[lev] += vbias_m[lev]
    return h0, u0, v0, h1, u1, v1


@nb.njit(cache=True)
def interp_nearest_linear_full_batch(time_unix_s, lat_vals, lon_vals,
                                     z_np, u_np, v_np, t_s,
                                     lats, lons, alts, u_bias, v_bias):
    """Nearest cell + linear u/v in altitude + linear in time (matches wind_interp)."""
    M = lats.size
    u_out = np.empty(M, np.float64)
    v_out = np.empty(M, np.float64)
    i_hr, hour_index = _time_bracket(time_unix_s, t_s)
    tf = hour_index - i_hr
    if tf < 0.0:
        tf = 0.0
    elif tf > 1.0:
        tf = 1.0
    for m in range(M):
        lat_i = _nearest_idx(lat_vals, lats[m])
        lon_i = _nearest_idx(lon_vals, lons[m])
        h0, u0, v0, h1, u1, v1 = _nearest_cols(z_np, u_np, v_np, i_hr, lat_i, lon_i,
                                               u_bias[m], v_bias[m])
        alt = alts[m]
        u0a = _interp_clamped(alt, h0, u0); v0a = _interp_clamped(alt, h0, v0)
        u1a = _interp_clamped(alt, h1, u1); v1a = _interp_clamped(alt, h1, v1)
        u_out[m] = u0a + tf * (u1a - u0a)
        v_out[m] = v0a + tf * (v1a - v0a)
    return u_out, v_out


@nb.njit(cache=True)
def _interp_bearing(h, u, v, alt):
    """Bearing(deg)+speed at altitude (matches wind_interp.interpolate_bearing)."""
    L = h.size
    nearest = _nearest_idx(h, alt)
    if alt > h[nearest]:
        i0 = nearest; i1 = nearest + 1
    else:
        i0 = nearest - 1; i1 = nearest
    if i0 == -1:
        i0 = 0
    if i1 == 0:
        i1 = -1
    # Resolve the original code's -1 (= last element) and clamp to valid range.
    if i1 < 0:
        i1 = L - 1
    if i1 > L - 1:
        i1 = L - 1
    if i0 < 0:
        i0 = 0
    if i0 > L - 1:
        i0 = L - 1

    b0 = math.atan2(v[i0], u[i0]); s0 = math.sqrt(u[i0] * u[i0] + v[i0] * v[i0])
    b1 = math.atan2(v[i1], u[i1]); s1 = math.sqrt(u[i1] * u[i1] + v[i1] * v[i1])
    a1 = _mod360(math.degrees(b0))
    a2 = _mod360(math.degrees(b1))
    if abs(a2 - a1) > 180.0:
        if a2 > a1:
            a1 += 360.0
        else:
            a2 += 360.0
    sp = _interp2(alt, h[i0], h[i1], s0, s1)
    dr = _mod360(_interp2(alt, h[i0], h[i1], a1, a2))
    return dr, sp


@nb.njit(cache=True)
def interp_nearest_linear_neighbors_batch(time_unix_s, lat_vals, lon_vals,
                                          z_np, u_np, v_np, t_s,
                                          lats, lons, alts, u_bias, v_bias):
    """Nearest cell + bearing/speed in altitude & time (matches wind_interp)."""
    M = lats.size
    u_out = np.empty(M, np.float64)
    v_out = np.empty(M, np.float64)
    i_hr, hour_index = _time_bracket(time_unix_s, t_s)
    tf = hour_index - i_hr
    if tf < 0.0:
        tf = 0.0
    elif tf > 1.0:
        tf = 1.0
    for m in range(M):
        lat_i = _nearest_idx(lat_vals, lats[m])
        lon_i = _nearest_idx(lon_vals, lons[m])
        h0, u0, v0, h1, u1, v1 = _nearest_cols(z_np, u_np, v_np, i_hr, lat_i, lon_i,
                                               u_bias[m], v_bias[m])
        alt = alts[m]
        dr0, sp0 = _interp_bearing(h0, u0, v0, alt)
        dr1, sp1 = _interp_bearing(h1, u1, v1, alt)
        a1 = _mod360(dr0)
        a2 = _mod360(dr1)
        if abs(a2 - a1) > 180.0:
            if a2 > a1:
                a1 += 360.0
            else:
                a2 += 360.0
        sp = sp0 + tf * (sp1 - sp0)
        dr = _mod360(a1 + tf * (a2 - a1))
        rad = math.radians(dr)
        u_out[m] = sp * math.cos(rad)
        v_out[m] = sp * math.sin(rad)
    return u_out, v_out
