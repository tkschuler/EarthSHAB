"""Hot-swappable wind-lookup backends for Forecast.py.

Three tiers with an identical interface, selected by ``config.forecast['backend']``
(or the ``backend=`` constructor kwarg). The tier swaps **only** the wind lookup
(data access + horizontal stencil); the vertical altitude/time math and advection
live elsewhere (``utils.wind_interp`` appliers; Forecast's ``_advect``):

    XarrayBackend : slow / original. Queries the live ``xr.Dataset`` every call
                    (no preload). Bilinear via ``.interp``, nearest methods via
                    ``.sel(method="nearest")``, then the shared appliers.
    NumpyBackend  : preloaded contiguous cache, pure-numpy per-member sampling
                    (``utils.wind_interp.sample_uv``).
    NumbaBackend  : preloaded cache, JIT batch kernels (``utils.vectorization``).
                    ``spline_full`` falls back to the numpy path — scipy's
                    CubicSpline cannot run under ``@njit``.

Every backend exposes one method used by the Forecast skeleton::

    sample_winds(method, lats, lons, alts, hour_index) -> (u, v)

``lats``/``lons``/``alts`` are 1-D float arrays (length M), ``hour_index`` is the
fractional forecast-hour index (scalar), and the returned ``u``/``v`` are float64
arrays (east, north) of length M. ``method`` is one of ``bilinear`` /
``linear_neighbors`` / ``linear_full`` / ``spline_full``.

All three should agree within float32-cache tolerance for a given
(method, advection); ``verify_backends.py`` checks that, using xarray as the
per-method reference.
"""
import numpy as np

import EarthSHAB.utils.wind_interp as wind_interp
import EarthSHAB.utils.vectorization as vectorization

GRAVITY = 9.80665  # geopotential (m^2 s^-2) -> geopotential height (m)


def _ascending_cache(lat, lon, hgt_m, u, v):
    """Flip lat/lon to ascending order (binary-search sampling requires it)."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    z = np.asarray(hgt_m, dtype=float)
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        z = z[:, :, ::-1, :]; u = u[:, :, ::-1, :]; v = v[:, :, ::-1, :]
    if lon[0] > lon[-1]:
        lon = lon[::-1]
        z = z[:, :, :, ::-1]; u = u[:, :, :, ::-1]; v = v[:, :, :, ::-1]
    return lat, lon, z, u, v


class NumpyBackend:
    """Preloaded contiguous cache; pure-numpy per-member interpolation."""

    name = "numpy"

    def __init__(self, *, lat, lon, hgt_m, u, v, **_ignored):
        self.lat, self.lon, self.z, self.u, self.v = _ascending_cache(lat, lon, hgt_m, u, v)
        self.T = self.z.shape[0]

    def sample_winds(self, method, lats, lons, alts, hour_index):
        M = len(lats)
        i_hr = min(int(np.floor(hour_index)), self.T - 2)
        u_out = np.empty(M, dtype=np.float64)
        v_out = np.empty(M, dtype=np.float64)
        for m in range(M):
            u_out[m], v_out[m] = wind_interp.sample_uv(
                method, self.lat, self.lon, self.z, self.u, self.v,
                i_hr, hour_index, float(lats[m]), float(lons[m]), float(alts[m]))
        return u_out, v_out

    def close(self):
        pass


class NumbaBackend:
    """Preloaded cache; JIT batch kernels (spline_full -> numpy fallback)."""

    name = "numba"

    def __init__(self, *, lat, lon, hgt_m, u, v, time_unix_s, **_ignored):
        la, lo, z, uu, vv = _ascending_cache(lat, lon, hgt_m, u, v)
        self.lat = la.astype(np.float32)
        self.lon = lo.astype(np.float32)
        self.z = z.astype(np.float32)
        self.u = uu.astype(np.float32)
        self.v = vv.astype(np.float32)
        self.time_unix_s = np.asarray(time_unix_s, dtype=np.int64)
        self.T, self.L = self.z.shape[0], self.z.shape[1]
        self.t0 = float(self.time_unix_s[0])
        self.res_s = float(self.time_unix_s[1] - self.time_unix_s[0]) if self.T >= 2 else 1.0

    def _t_s(self, hour_index):
        # Real sim timestamps are integer seconds, so this round-trips exactly:
        # the kernel recomputes hour_index = (t_s - t0) / res_s.
        return np.int64(round(self.t0 + hour_index * self.res_s))

    def sample_winds(self, method, lats, lons, alts, hour_index):
        M = len(lats)
        lats64 = np.ascontiguousarray(lats, dtype=np.float64)
        lons64 = np.ascontiguousarray(lons, dtype=np.float64)
        alts64 = np.ascontiguousarray(alts, dtype=np.float64)
        zb = np.zeros((M, self.L), dtype=np.float32)

        if method == "spline_full":
            # scipy CubicSpline can't run under numba -> numpy path (per design).
            i_hr = min(int(np.floor(hour_index)), self.T - 2)
            u_out = np.empty(M, dtype=np.float64)
            v_out = np.empty(M, dtype=np.float64)
            for m in range(M):
                u_out[m], v_out[m] = wind_interp.sample_uv(
                    "spline_full", self.lat, self.lon, self.z, self.u, self.v,
                    i_hr, hour_index, float(lats64[m]), float(lons64[m]), float(alts64[m]))
            return u_out, v_out

        t_s = self._t_s(hour_index)
        if method == "bilinear":
            u, v = vectorization.interp4d_alt_batch_match_xr(
                self.time_unix_s, self.lat, self.lon, self.z, self.u, self.v, t_s,
                lats64.astype(np.float32), lons64.astype(np.float32),
                alts64.astype(np.float32), zb, zb)
            return u.astype(np.float64), v.astype(np.float64)
        if method == "linear_full":
            return vectorization.interp_nearest_linear_full_batch(
                self.time_unix_s, self.lat, self.lon, self.z, self.u, self.v, t_s,
                lats64, lons64, alts64, zb, zb)
        if method == "linear_neighbors":
            return vectorization.interp_nearest_linear_neighbors_batch(
                self.time_unix_s, self.lat, self.lon, self.z, self.u, self.v, t_s,
                lats64, lons64, alts64, zb, zb)
        raise ValueError(f"unknown wind_interpolation: {method!r}")

    def close(self):
        pass


class XarrayBackend:
    """Slow / original: query the live ``xr.Dataset`` each call (no preload)."""

    name = "xarray"

    def __init__(self, *, path, **_ignored):
        import xarray as xr
        self.ds = xr.open_dataset(path)
        self.T = int(self.ds.sizes["valid_time"])

    def _columns(self, method, lat, lon, i_hr):
        """(h, u, v) profile columns at forecast hours i_hr and i_hr+1."""
        n0 = self.ds.isel(valid_time=i_hr)
        n1 = self.ds.isel(valid_time=i_hr + 1)
        if method == "bilinear":
            n0 = n0.interp(latitude=lat, longitude=lon)
            n1 = n1.interp(latitude=lat, longitude=lon)
        else:
            n0 = n0.sel(latitude=lat, longitude=lon, method="nearest")
            n1 = n1.sel(latitude=lat, longitude=lon, method="nearest")
        f = wind_interp.fill_missing_data
        return (f(n0.z.values / GRAVITY), f(n0.u.values), f(n0.v.values),
                f(n1.z.values / GRAVITY), f(n1.u.values), f(n1.v.values))

    def sample_winds(self, method, lats, lons, alts, hour_index):
        M = len(lats)
        i_hr = min(int(np.floor(hour_index)), self.T - 2)
        u_out = np.empty(M, dtype=np.float64)
        v_out = np.empty(M, dtype=np.float64)
        for m in range(M):
            h0, u0, v0, h1, u1, v1 = self._columns(method, float(lats[m]), float(lons[m]), i_hr)
            u_out[m], v_out[m] = wind_interp.apply_method(
                method, h0, u0, v0, h1, u1, v1, float(alts[m]), hour_index)
        return u_out, v_out

    def close(self):
        try:
            self.ds.close()
        except Exception:
            pass


_BACKENDS = {"xarray": XarrayBackend, "numpy": NumpyBackend, "numba": NumbaBackend}


def make_backend(name, **kwargs):
    """Instantiate a backend by name, passing it whatever state it needs.

    All backends accept the same kwargs (``path``, ``lat``, ``lon``, ``hgt_m``,
    ``u``, ``v``, ``time_unix_s``) and ignore the ones they don't use, so the
    Forecast skeleton constructs any of them the same way.
    """
    try:
        cls = _BACKENDS[name]
    except KeyError:
        raise ValueError(
            f"unknown backend: {name!r}. Expected one of {sorted(_BACKENDS)}.")
    return cls(**kwargs)
