"""evaluate.py - Compares an EarthSHAB simulation against balloon flight ground truth.

Structure mirrors main.py. BalloonEvaluator accepts optional config overrides so
the same class can be reused for batch testing across flights and parameter sweeps.

Usage (single flight):
    python -m EarthSHAB.Evaluation.evaluate

Batch example (future):
    for mp in [0.5, 0.7, 1.0]:
        ev = BalloonEvaluator(config_overrides={'balloon_properties': {'mp': mp}})
        ev.run()
        results[mp] = ev.compute_metrics()
"""

import csv
import os
import re
import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
import fluids
from geographiclib.geodesic import Geodesic
from termcolor import colored
import matplotlib.pyplot as plt

import EarthSHAB.config_earth as config_earth
from EarthSHAB.simulate import BalloonSimulation
from EarthSHAB.Evaluation.plot_evaluation import plot_comparison
from EarthSHAB.Plotting.plot_trajectory_map import plot_map
from EarthSHAB.Plotting.plot_windmap import plot_windmap


# ── Data classes ─────────────────────────────────────────────────────────────

@dataclass
class FlightMetrics:
    """Summary statistics for one flight (simulation or ground truth)."""
    float_alt_mean: float = math.nan
    float_alt_std:  float = math.nan
    ascent_rate_mean: float = math.nan
    ascent_rate_std:  float = math.nan
    descent_rate_mean: float = math.nan
    descent_rate_std:  float = math.nan
    elapsed_time_min: float = math.nan
    landing_lat:  float = math.nan
    landing_lon:  float = math.nan
    landing_time: Optional[object] = None  # pd.Timestamp or datetime


@dataclass
class EvaluationResult:
    """Cross-comparison metrics between simulation and ground truth."""
    sim:   FlightMetrics = field(default_factory=FlightMetrics)
    truth: FlightMetrics = field(default_factory=FlightMetrics)
    landing_distance_km:  float = math.nan
    landing_time_diff_min: float = math.nan
    temp_mae_k:     float = math.nan
    pressure_mae_pa: float = math.nan
    # GFS Forecast applied to truth altitude profile (reforecast trajectory)
    gfs_truth_landing_lat:  float = math.nan
    gfs_truth_landing_lon:  float = math.nan
    gfs_truth_landing_dist_m: float = math.nan


# ── Evaluator ────────────────────────────────────────────────────────────────

class BalloonEvaluator:
    GMT = 7  # MST offset
    OUTPUT_DIR = "src/EarthSHAB/Evaluation/"

    def __init__(self, config_overrides: dict = None):
        """
        Parameters
        ----------
        config_overrides : dict, optional
            Nested dict of config_earth attributes to override before running.
            Keys are attribute names on the config_earth module; values are
            either replacement dicts (for dict-typed config sections) or
            scalar replacements.  Example::

                {
                    'balloon_properties': {'mp': 0.5},
                    'simulation': {'sim_time': 12},
                }
        """
        if config_overrides:
            self._apply_overrides(config_overrides)

        self.geod = Geodesic.WGS84
        self.sim: Optional[BalloonSimulation] = None
        self.sim_state: Optional[dict] = None
        self.df_truth: Optional[pd.DataFrame] = None
        self.result: Optional[EvaluationResult] = None

    def _apply_overrides(self, overrides: dict):
        for attr, value in overrides.items():
            cfg = getattr(config_earth, attr)
            if isinstance(cfg, dict):
                cfg.update(value)
            else:
                setattr(config_earth, attr, value)

    # ── Run ──────────────────────────────────────────────────────────────────

    def run(self):
        """Run the simulation (mirrors main.py) and prepare ground truth data."""
        traj_path = config_earth.simulation['balloon_trajectory']
        if traj_path is None:
            raise ValueError(
                "balloon_trajectory must be set in config_earth to run an evaluation"
            )

        self.sim = BalloonSimulation()
        self.sim.run_simulation(run_reforecast=False)
        self.sim.run_simulation(run_reforecast=True)
        self.sim_state = self.sim.sim_state

        # Enrich the APRS dataframe (already loaded by run_reforecast=True)
        df = self.sim_state["df"].copy()
        df["temp_k"]      = df["comment"].apply(self._parse_temp_k)
        df["pressure_pa"] = df["comment"].apply(self._parse_pressure_pa)
        # Vertical velocity from altitude finite-differences (m/s)
        df["v_truth"] = (
            df["altitude"].diff() / df["dt"].replace(0, np.nan)
        ).fillna(0)
        self.df_truth = df

    # ── APRS parsers ─────────────────────────────────────────────────────────

    @staticmethod
    def _parse_temp_k(comment: str) -> float:
        m = re.search(r'(-?\d+(?:\.\d+)?)C', str(comment))
        return float(m.group(1)) + 273.15 if m else math.nan

    @staticmethod
    def _parse_pressure_pa(comment: str) -> float:
        m = re.search(r'(\d+(?:\.\d+)?)Pa', str(comment))
        return float(m.group(1)) if m else math.nan

    # ── Phase detection ───────────────────────────────────────────────────────

    @staticmethod
    def _detect_phases(el, v, min_alt, alt_fraction=0.90, v_float=1.0, v_linear=1.0):
        """Return (ascent_mask, float_mask, descent_mask) boolean arrays.

        Strategy:
        - **Float**: altitude >= alt_fraction * el_max AND |v| < v_float.
          The combined condition excludes the curved transition regions where
          the balloon is decelerating into or accelerating out of float.
        - **Ascent**: v > v_linear, strictly before the balloon enters the
          high-altitude region (i.e. before the deceleration curve begins).
        - **Descent**: v < -v_linear, strictly after the balloon exits the
          high-altitude region (i.e. after the initial slow post-float drop).
        """
        el = np.asarray(el, dtype=float)
        v  = np.asarray(v,  dtype=float)
        n  = len(el)

        el_max        = el.max()
        alt_threshold = alt_fraction * el_max

        # Find where the balloon enters and exits the high-altitude region
        high_indices = np.where(el >= alt_threshold)[0]
        if len(high_indices) > 0:
            i_enter = int(high_indices[0])
            i_exit  = int(high_indices[-1])
        else:
            i_peak  = int(np.argmax(el))
            i_enter = i_exit = i_peak

        # Float: slow motion while at high altitude (excludes rounded transitions)
        float_mask = (el >= alt_threshold) & (np.abs(v) < v_float)

        # Ascent: linear climb before the deceleration into float
        ascent_mask = np.zeros(n, dtype=bool)
        ascent_mask[:i_enter] = v[:i_enter] > v_linear

        # Descent: linear fall after the balloon leaves the float altitude band
        descent_mask = np.zeros(n, dtype=bool)
        descent_mask[i_exit:] = v[i_exit:] < -v_linear

        return ascent_mask, float_mask, descent_mask, i_enter, i_exit

    @staticmethod
    def _find_landing(el, time_arr, min_alt):
        """Return (time, index) of first post-peak descent to near min_alt."""
        el = np.asarray(el, dtype=float)
        i_peak = int(np.argmax(el))
        threshold = min_alt + 200
        for i in range(i_peak, len(el)):
            if el[i] <= threshold:
                return time_arr[i], i
        return time_arr[-1], len(el) - 1

    # ── Metrics ───────────────────────────────────────────────────────────────

    def compute_metrics(self) -> EvaluationResult:
        ss  = self.sim_state
        df  = self.df_truth
        min_alt = config_earth.simulation['min_alt']

        # ── Simulation ───────────────────────────────────────────────────────
        el_sim = np.array(ss["el"])
        v_sim  = np.array(ss["v"])
        t_sim  = ss["time_local"]

        asc_m, flt_m, des_m, i_enter_sim, i_exit_sim = self._detect_phases(el_sim, v_sim, min_alt)

        def _mean(arr, mask): return float(np.mean(arr[mask])) if mask.any() else math.nan
        def _std(arr, mask):  return float(np.std(arr[mask]))  if mask.any() else math.nan

        land_t_sim, land_i_sim = self._find_landing(el_sim, t_sim, min_alt)

        elapsed_sim = (
            pd.Timestamp(land_t_sim) - pd.Timestamp(t_sim[0])
        ).total_seconds() / 60.0

        sim_metrics = FlightMetrics(
            float_alt_mean=_mean(el_sim, flt_m),
            float_alt_std=_std(el_sim,  flt_m),
            ascent_rate_mean=_mean(v_sim, asc_m),
            ascent_rate_std=_std(v_sim,  asc_m),
            descent_rate_mean=_mean(v_sim, des_m),
            descent_rate_std=_std(v_sim,  des_m),
            elapsed_time_min=elapsed_sim,
            landing_lat=ss["lat"][land_i_sim],
            landing_lon=ss["lon"][land_i_sim],
            landing_time=land_t_sim,
        )

        # ── Ground truth ─────────────────────────────────────────────────────
        el_truth = df["altitude"].to_numpy(dtype=float)
        v_truth  = df["v_truth"].to_numpy(dtype=float)
        t_truth  = df["time"].tolist()

        asc_t, flt_t, des_t, i_enter_truth, i_exit_truth = self._detect_phases(el_truth, v_truth, min_alt)

        elapsed_truth = (
            pd.Timestamp(t_truth[-1]) - pd.Timestamp(t_truth[0])
        ).total_seconds() / 60.0

        truth_metrics = FlightMetrics(
            float_alt_mean=_mean(el_truth, flt_t),
            float_alt_std=_std(el_truth,  flt_t),
            ascent_rate_mean=_mean(v_truth, asc_t),
            ascent_rate_std=_std(v_truth,  asc_t),
            descent_rate_mean=_mean(v_truth, des_t),
            descent_rate_std=_std(v_truth,  des_t),
            elapsed_time_min=elapsed_truth,
            landing_lat=float(df["lat"].iloc[-1]),
            landing_lon=float(df["lng"].iloc[-1]),
            landing_time=t_truth[-1],
        )

        # ── Cross metrics ─────────────────────────────────────────────────────
        geo = self.geod.Inverse(
            sim_metrics.landing_lat, sim_metrics.landing_lon,
            truth_metrics.landing_lat, truth_metrics.landing_lon,
        )
        landing_dist_km = geo['s12'] / 1000.0

        if sim_metrics.landing_time is not None and truth_metrics.landing_time is not None:
            landing_dt_min = (
                pd.Timestamp(sim_metrics.landing_time)
                - pd.Timestamp(truth_metrics.landing_time)
            ).total_seconds() / 60.0
        else:
            landing_dt_min = math.nan

        # ── GFS Forecast + Truth Altitude (reforecast) ────────────────────────
        gfs_lat_arr = ss.get("lat_aprs_gps", [])
        gfs_lon_arr = ss.get("lon_aprs_gps", [])
        if gfs_lat_arr and gfs_lon_arr:
            gfs_land_lat = float(gfs_lat_arr[-1])
            gfs_land_lon = float(gfs_lon_arr[-1])
            geo_gfs = self.geod.Inverse(
                gfs_land_lat, gfs_land_lon,
                truth_metrics.landing_lat, truth_metrics.landing_lon,
            )
            gfs_land_dist_m = geo_gfs['s12']
        else:
            gfs_land_lat = gfs_land_lon = gfs_land_dist_m = math.nan

        self.result = EvaluationResult(
            sim=sim_metrics,
            truth=truth_metrics,
            landing_distance_km=landing_dist_km,
            landing_time_diff_min=landing_dt_min,
            temp_mae_k=self._compute_temp_mae(ss, df),
            pressure_mae_pa=self._compute_pressure_mae(df),
            gfs_truth_landing_lat=gfs_land_lat,
            gfs_truth_landing_lon=gfs_land_lon,
            gfs_truth_landing_dist_m=gfs_land_dist_m,
        )

        # Store phase info for use in plot()
        self._sim_phases = {
            'ascent_mask': asc_m,
            'float_mask':  flt_m,
            'descent_mask': des_m,
            'i_enter': i_enter_sim,
            'i_exit':  i_exit_sim,
            'float_alt_mean': sim_metrics.float_alt_mean,
            'float_alt_std':  sim_metrics.float_alt_std,
        }
        self._truth_phases = {
            'ascent_mask': asc_t,
            'float_mask':  flt_t,
            'descent_mask': des_t,
            'i_enter': i_enter_truth,
            'i_exit':  i_exit_truth,
            'float_alt_mean': truth_metrics.float_alt_mean,
            'float_alt_std':  truth_metrics.float_alt_std,
        }

        return self.result

    @staticmethod
    def _suggest_start_time(df: pd.DataFrame, min_alt: float, gmt: int):
        """Estimate actual launch time by linearly extrapolating the initial APRS
        ascent rate back to min_alt.

        Returns (suggested_utc, first_aprs_utc, v_mean_ms) or (None, None, nan)
        if there is insufficient ascending data.
        """
        # Skip the first row (dt=0 → v_truth=0) and find initial ascending points
        ascending = df.iloc[1:][df.iloc[1:]['v_truth'] > 0.5].head(5)
        if ascending.empty:
            return None, None, math.nan

        v_mean    = float(ascending['v_truth'].mean())
        first_alt = float(df['altitude'].iloc[0])
        first_mst = pd.Timestamp(df['time'].iloc[0])

        # Time to ascend from min_alt to first_alt at the measured initial rate
        dt_s = (first_alt - min_alt) / v_mean
        suggested_mst = first_mst - pd.Timedelta(seconds=dt_s)

        suggested_utc    = suggested_mst + pd.Timedelta(hours=gmt)
        first_aprs_utc   = first_mst    + pd.Timedelta(hours=gmt)
        return suggested_utc, first_aprs_utc, v_mean

    def _compute_temp_mae(self, ss: dict, df: pd.DataFrame) -> float:
        """Interpolate sim T_atm to APRS timestamps, return MAE vs onboard sensor."""
        valid = df["temp_k"].notna()
        if not valid.any():
            return math.nan
        t0 = pd.Timestamp(ss["time_local"][0])
        t_sim_s  = np.array([(pd.Timestamp(t) - t0).total_seconds() for t in ss["time_local"]])
        t_aprs_s = np.array([(pd.Timestamp(t) - t0).total_seconds() for t in df["time"][valid]])
        T_interp = np.interp(t_aprs_s, t_sim_s, np.array(ss["T_atm"]))
        return float(np.mean(np.abs(T_interp - df["temp_k"].to_numpy()[valid])))

    def _compute_pressure_mae(self, df: pd.DataFrame) -> float:
        """Compare ISA model pressure (from sim altitude profile) vs APRS measured pressure."""
        valid = df["pressure_pa"].notna()
        if not valid.any():
            return math.nan
        p_isa = np.array([
            fluids.atmosphere.ATMOSPHERE_1976(alt).P
            for alt in df["altitude"].to_numpy()[valid]
        ])
        return float(np.mean(np.abs(p_isa - df["pressure_pa"].to_numpy()[valid])))

    def _stem(self) -> str:
        """Build a consistent filename stem: {trajectory}_{forecast}_{Y}_{M}_{D}."""
        ss = self.sim_state
        name     = ss.get("trajectory_name") or "Unknown"
        forecast = ss.get("forecast_type", "")
        start    = ss.get("start")
        date_str = f"{start.year}_{start.month}_{start.day}" if start else "unknown"
        return f"{name}_{forecast}_{date_str}"

    # ── Report ────────────────────────────────────────────────────────────────

    def report(self, result: Optional[EvaluationResult] = None):
        if result is None:
            result = self.result

        self._report_body(result)

        os.makedirs(self.OUTPUT_DIR, exist_ok=True)
        csv_path = os.path.join(self.OUTPUT_DIR, f"{self._stem()}.csv")
        self._write_csv(result, csv_path)
        print(f"  Report saved → {csv_path}")

    def _report_body(self, result: EvaluationResult):
        traj_path = config_earth.simulation.get('balloon_trajectory', '') or ''
        name = traj_path.split("/")[-1].replace(".csv", "") or "Unknown"

        W = 68

        def _fmt(v, dec=1):
            return "N/A" if (isinstance(v, float) and math.isnan(v)) else f"{v:.{dec}f}"

        def _row(label, sv, tv, dv, unit):
            return f"  {label:<30} {sv:>9}  {tv:>9}  {dv:>9}  {unit}"

        s, t = result.sim, result.truth

        def _diff(a, b):
            return a - b if not (math.isnan(a) or math.isnan(b)) else math.nan

        rows = [
            ("Float Alt Mean (m)",       s.float_alt_mean,     t.float_alt_mean,     _diff(s.float_alt_mean,     t.float_alt_mean),     "m",   0),
            ("Float Alt Std (m)",        s.float_alt_std,      t.float_alt_std,       _diff(s.float_alt_std,      t.float_alt_std),      "m",   0),
            ("Ascent Rate Mean (m/s)",   s.ascent_rate_mean,   t.ascent_rate_mean,   _diff(s.ascent_rate_mean,   t.ascent_rate_mean),   "m/s", 2),
            ("Ascent Rate Std (m/s)",    s.ascent_rate_std,    t.ascent_rate_std,    _diff(s.ascent_rate_std,    t.ascent_rate_std),    "m/s", 2),
            ("Descent Rate Mean (m/s)",  s.descent_rate_mean,  t.descent_rate_mean,  _diff(s.descent_rate_mean,  t.descent_rate_mean),  "m/s", 2),
            ("Descent Rate Std (m/s)",   s.descent_rate_std,   t.descent_rate_std,   _diff(s.descent_rate_std,   t.descent_rate_std),   "m/s", 2),
            ("Elapsed Time (min)",       s.elapsed_time_min,   t.elapsed_time_min,   _diff(s.elapsed_time_min,   t.elapsed_time_min),   "min", 1),
            ("Landing Lat (°)",          s.landing_lat,        t.landing_lat,        _diff(s.landing_lat,        t.landing_lat),        "°",   4),
            ("Landing Lon (°)",          s.landing_lon,        t.landing_lon,        _diff(s.landing_lon,        t.landing_lon),        "°",   4),
        ]

        land_dist_m = (
            result.landing_distance_km * 1000
            if not math.isnan(result.landing_distance_km) else math.nan
        )

        print()
        print("=" * W)
        print(f"  EarthSHAB Evaluation: {name}")
        print("=" * W)
        print(f"  {'Metric':<30} {'Sim':>9}  {'Truth':>9}  {'Diff':>9}  Unit")
        print("-" * W)
        for label, sv, tv, dv, unit, dec in rows:
            print(_row(label, _fmt(sv, dec), _fmt(tv, dec), _fmt(dv, dec), unit))
        print("-" * W)
        print(_row("Distance Off (m)",   "",  "", _fmt(land_dist_m, 0), "m"))

        t_sim_str   = str(s.landing_time)[:16] if s.landing_time else "N/A"
        t_truth_str = str(t.landing_time)[:16] if t.landing_time else "N/A"
        print(_row("Landing Time (MST)", t_sim_str, t_truth_str,
                   _fmt(result.landing_time_diff_min, 1), "min"))
        print("-" * W)
        print(_row("Temperature MAE",    "", "", _fmt(result.temp_mae_k, 2), "K"))
        print(_row("Pressure MAE",       "", "", _fmt(result.pressure_mae_pa, 0), "Pa"))
        print("-" * W)
        print(f"  GFS Forecast + Truth Altitude (reforecast landing vs truth)")
        print(f"  {'Metric':<30} {'GFS+TA':>9}  {'Truth':>9}  {'Diff':>9}  Unit")
        gfs = result
        print(_row("Landing Lat (°)",
                   _fmt(gfs.gfs_truth_landing_lat, 4),
                   _fmt(t.landing_lat, 4),
                   _fmt(_diff(gfs.gfs_truth_landing_lat, t.landing_lat), 4), "°"))
        print(_row("Landing Lon (°)",
                   _fmt(gfs.gfs_truth_landing_lon, 4),
                   _fmt(t.landing_lon, 4),
                   _fmt(_diff(gfs.gfs_truth_landing_lon, t.landing_lon), 4), "°"))
        print(_row("Distance Off (m)", "", "", _fmt(gfs.gfs_truth_landing_dist_m, 0), "m"))
        print("=" * W)

        # ── Start-time suggestion ──────────────────────────────────────────
        if self.df_truth is not None:
            min_alt = config_earth.simulation['min_alt']
            current_utc = config_earth.simulation['start_time']
            sugg_utc, first_aprs_utc, v_mean = self._suggest_start_time(
                self.df_truth, min_alt, self.GMT
            )
            print()
            print("  Start-time analysis")
            print(f"  Current config start_time : {current_utc} UTC")
            if first_aprs_utc is not None:
                print(f"  First APRS transmission   : {first_aprs_utc} UTC  "
                      f"({self.df_truth['altitude'].iloc[0]:.0f} m)")
            if sugg_utc is not None:
                print(f"  Estimated launch time     : {sugg_utc.strftime('%Y-%m-%d %H:%M:%S')} UTC"
                      f"  (extrapolated at {v_mean:.2f} m/s from first APRS point)")
                print(f"  Suggested start_time      : \"{sugg_utc.strftime('%Y-%m-%d %H:%M:%S')}\"")

    def _write_csv(self, result: EvaluationResult, path: str):
        """Write evaluation metrics to a structured CSV file."""
        s, t = result.sim, result.truth

        def _v(x, dec=4):
            return "" if (isinstance(x, float) and math.isnan(x)) else f"{x:.{dec}f}"

        def _diff(a, b):
            return a - b if not (math.isnan(a) or math.isnan(b)) else math.nan

        land_dist_m = result.landing_distance_km * 1000 if not math.isnan(result.landing_distance_km) else math.nan

        rows = [
            ("Float Alt Mean",        "m",   _v(s.float_alt_mean, 0),    _v(t.float_alt_mean, 0),    _v(_diff(s.float_alt_mean,    t.float_alt_mean),    0)),
            ("Float Alt Std",         "m",   _v(s.float_alt_std,  0),    _v(t.float_alt_std,  0),    _v(_diff(s.float_alt_std,     t.float_alt_std),     0)),
            ("Ascent Rate Mean",      "m/s", _v(s.ascent_rate_mean, 2),  _v(t.ascent_rate_mean, 2),  _v(_diff(s.ascent_rate_mean,  t.ascent_rate_mean),  2)),
            ("Ascent Rate Std",       "m/s", _v(s.ascent_rate_std,  2),  _v(t.ascent_rate_std,  2),  _v(_diff(s.ascent_rate_std,   t.ascent_rate_std),   2)),
            ("Descent Rate Mean",     "m/s", _v(s.descent_rate_mean, 2), _v(t.descent_rate_mean, 2), _v(_diff(s.descent_rate_mean, t.descent_rate_mean), 2)),
            ("Descent Rate Std",      "m/s", _v(s.descent_rate_std,  2), _v(t.descent_rate_std,  2), _v(_diff(s.descent_rate_std,  t.descent_rate_std),  2)),
            ("Elapsed Time",          "min", _v(s.elapsed_time_min, 1),  _v(t.elapsed_time_min, 1),  _v(_diff(s.elapsed_time_min,  t.elapsed_time_min),  1)),
            ("Landing Lat",           "deg", _v(s.landing_lat, 4),       _v(t.landing_lat, 4),       _v(_diff(s.landing_lat,       t.landing_lat),       4)),
            ("Landing Lon",           "deg", _v(s.landing_lon, 4),       _v(t.landing_lon, 4),       _v(_diff(s.landing_lon,       t.landing_lon),       4)),
            ("Landing Time (MST)",    "",    str(s.landing_time)[:16] if s.landing_time else "",
                                             str(t.landing_time)[:16] if t.landing_time else "",
                                             _v(result.landing_time_diff_min, 1)),
            ("Distance Off",          "m",   "", "", _v(land_dist_m, 0)),
            ("Temperature MAE",       "K",   "", "", _v(result.temp_mae_k, 2)),
            ("Pressure MAE",          "Pa",  "", "", _v(result.pressure_mae_pa, 0)),
            ("GFS+TA Landing Lat",    "deg", _v(result.gfs_truth_landing_lat, 4),
                                             _v(t.landing_lat, 4),
                                             _v(_diff(result.gfs_truth_landing_lat, t.landing_lat), 4)),
            ("GFS+TA Landing Lon",    "deg", _v(result.gfs_truth_landing_lon, 4),
                                             _v(t.landing_lon, 4),
                                             _v(_diff(result.gfs_truth_landing_lon, t.landing_lon), 4)),
            ("GFS+TA Distance Off",   "m",   "", "", _v(result.gfs_truth_landing_dist_m, 0)),
        ]

        with open(path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(["Metric", "Unit", "Sim", "Truth", "Diff"])
            w.writerows(rows)

    # ── Plot ─────────────────────────────────────────────────────────────────

    def plot(self):
        ss = self.sim_state
        os.makedirs(self.OUTPUT_DIR, exist_ok=True)
        stem = self._stem()

        fig = plot_comparison(
            time_sim=ss["time_local"],
            el_sim=ss["el"],
            v_sim=ss["v"],
            T_atm_sim=ss["T_atm"],
            df_truth=self.df_truth,
            sim_phases=getattr(self, '_sim_phases', None),
            truth_phases=getattr(self, '_truth_phases', None),
        )
        png_path = os.path.join(self.OUTPUT_DIR, f"{stem}.png")
        fig.savefig(png_path, dpi=150, bbox_inches='tight')
        print(f"  Plot saved     → {png_path}")

        plot_map(
            gmap1=ss["gmap1"],
            coord=ss["coord"],
            forecast_type=ss["forecast_type"],
            balloon_trajectory=ss["balloon_trajectory"],
            trajectory_name=ss["trajectory_name"],
            lat=ss["lat"],
            lon=ss["lon"],
            lat_aprs_gps=ss["lat_aprs_gps"],
            lon_aprs_gps=ss["lon_aprs_gps"],
            df=ss["df"],
            gfs=ss["gfs"],
            sim=self.sim,
            t=ss["t"],
            start=ss["start"],
            html_prefix="EVALUATION_",
            output_dir=self.OUTPUT_DIR,
        )
        print(f"  Map saved      → {os.path.join(self.OUTPUT_DIR, 'EVALUATION_' + stem + '.html')}")
        plot_windmap()


# ── Entry point ───────────────────────────────────────────────────────────────

def evaluate():
    ev = BalloonEvaluator()
    ev.run()
    result = ev.compute_metrics()
    ev.report(result)
    ev.plot()
    plt.show()


if __name__ == "__main__":
    evaluate()
