"""evaluate.py - Compares an EarthSHAB simulation against balloon flight ground truth.

Structure mirrors main.py. BalloonEvaluator accepts optional config overrides so
the same class can be reused for batch testing across flights and parameter sweeps.

Usage (single flight):
    python -m evaluation.evaluate

Batch example (future):
    for mp in [0.5, 0.7, 1.0]:
        ev = BalloonEvaluator(config_overrides={'balloon_properties': {'mp': mp}})
        ev.run()
        results[mp] = ev.compute_metrics()
"""

import os
import re
import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
import fluids
from geographiclib.geodesic import Geodesic
import matplotlib.pyplot as plt

import EarthSHAB.config as config_earth
from EarthSHAB.simulate import BalloonSimulation
from EarthSHAB.radiation import solar_zenith_adjusted
from evaluation.plot_evaluation import plot_comparison
from evaluation import reporting
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
    elapsed_time_min:    float = math.nan
    time_to_float_min:   float = math.nan
    landing_lat:  float = math.nan
    landing_lon:  float = math.nan
    landing_time: Optional[object] = None  # pd.Timestamp or datetime


@dataclass
class EvaluationResult:
    """Cross-comparison metrics between simulation and ground truth."""
    sim:   FlightMetrics = field(default_factory=FlightMetrics)
    truth: FlightMetrics = field(default_factory=FlightMetrics)
    landing_distance_km:  float = math.nan
    landing_time_diff_min:   float = math.nan
    time_to_float_diff_min:  float = math.nan
    time_to_ground_diff_min: float = math.nan
    temp_mae_k:     float = math.nan
    pressure_mae_pa: float = math.nan
    # GFS Forecast applied to truth altitude profile (reforecast trajectory)
    gfs_truth_landing_lat:  float = math.nan
    gfs_truth_landing_lon:  float = math.nan
    gfs_truth_landing_dist_m: float = math.nan


# ── Evaluator ────────────────────────────────────────────────────────────────

class BalloonEvaluator:
    GMT = 7  # MST offset
    OUTPUT_DIR = "evaluation/"

    # Launches whose ascent is non-physical for EarthSHAB's solar-balloon model:
    # carried up by helium ("helium_augmented") or by a weather balloon
    # ("grand_slam"). For these the ascent metrics are blanked out — only the
    # float and descent phases get scored.
    NON_STANDARD_LAUNCH_TYPES = {"helium_augmented", "grand_slam"}

    def __init__(self, config_overrides: dict = None, output_dir: str = None,
                 launch_type: str = "standard"):
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
        output_dir : str, optional
            Directory where reports, plots, and maps are saved.
            Defaults to OUTPUT_DIR ("evaluation/").
        """
        if config_overrides:
            self._apply_overrides(config_overrides)

        self._output_dir = output_dir if output_dir is not None else self.OUTPUT_DIR
        self.launch_type = launch_type
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
    def _detect_phases(el, v, min_alt, alt_fraction=0.90, v_float=1.0, v_linear=1.0,
                       launch_type="standard"):
        """Return (ascent_mask, float_mask, descent_mask, i_enter, i_exit).

        Strategy:
        - **Float**: within the high-altitude region [i_enter, i_exit], find where
          the rolling-mean velocity is approximately zero.  Using a rolling mean
          rather than instantaneous velocity prevents two failure modes:
            * sim: deceleration curve into float has |v| < v_float but mean > v_float
            * truth: APRS noise spikes have |v| > v_float but local mean ≈ 0
        - **Ascent**: v > v_linear strictly before i_enter.
        - **Descent**: v < -v_linear strictly after i_exit.

        Grand Slam launches: the balloon is carried by a weather balloon to
        well above its natural float altitude, then released and *descends*
        into the SHAB float plateau.  ``alt ≥ alt_fraction · max_alt`` would
        clip the search to the brief release peak and miss the actual float —
        so for ``launch_type == "grand_slam"`` we widen the search bracket to
        the entire post-apex region of the trajectory.  Descent is then
        defined as everything after the detected float ends.
        """
        el = np.asarray(el, dtype=float)
        v  = np.asarray(v,  dtype=float)
        n  = len(el)

        # Use nanmax so a single missing altitude reading doesn't collapse detection.
        el_max        = np.nanmax(el)
        alt_threshold = alt_fraction * el_max
        i_peak        = int(np.nanargmax(el))

        if launch_type == "grand_slam":
            # Search for SHAB float anywhere from the release peak onward.
            i_enter = i_peak
            i_exit  = n - 1
        else:
            # Standard / helium-augmented: balloon ascends, peaks at float.
            high_indices = np.where(el >= alt_threshold)[0]
            if len(high_indices) > 0:
                i_enter = int(high_indices[0])
                i_exit  = int(high_indices[-1])
            else:
                i_enter = i_exit = i_peak

        # Rolling-mean velocity within the high-altitude slice only.
        # Window ≈ 1/6 of the high-altitude span, capped so it stays meaningful.
        span = max(i_exit - i_enter, 1)
        win  = max(10, min(span // 6, 600))
        v_hi_roll = (pd.Series(v[i_enter:i_exit + 1])
                     .rolling(window=win, center=True, min_periods=1)
                     .mean().to_numpy())

        candidate = np.abs(v_hi_roll) < v_float

        # Keep only the LARGEST contiguous True block inside the candidate mask.
        # Discard it if it's shorter than 1/4 of the high-altitude span — this
        # rejects spurious "float" detections caused by the velocity zero-crossing
        # at the apex of a non-float trajectory.
        float_mask = np.zeros(n, dtype=bool)
        if candidate.any():
            chg        = np.diff(np.concatenate([[False], candidate, [False]]).astype(int))
            blk_starts = np.where(chg ==  1)[0]
            blk_ends   = np.where(chg == -1)[0]
            lengths    = blk_ends - blk_starts
            best       = int(np.argmax(lengths))
            min_len    = max(10, span // 4)
            if lengths[best] >= min_len:
                bs = blk_starts[best]
                be = blk_ends[best]
                # Trim the block front: skip the ascending approach into float
                # (rolling mean still > 0.3 m/s means balloon is still climbing).
                # Trim the block back: skip rapid descent out of float
                # (rolling mean < -0.5 m/s means balloon has clearly left float).
                v_trim_hi =  v_float * 0.3
                v_trim_lo = -v_float * 0.5
                front_trim = be  # default: skip if threshold never met
                for idx in range(bs, be):
                    if v_hi_roll[idx] < v_trim_hi:
                        front_trim = idx
                        break
                back_trim = front_trim  # default
                for idx in range(be - 1, front_trim - 1, -1):
                    if v_hi_roll[idx] > v_trim_lo:
                        back_trim = idx + 1
                        break
                if back_trim - front_trim >= min_len:
                    float_mask[i_enter + front_trim : i_enter + back_trim] = True

        # Ascent: linear climb before the deceleration into float
        ascent_mask = np.zeros(n, dtype=bool)
        ascent_mask[:i_enter] = v[:i_enter] > v_linear

        # Descent: linear fall after the balloon leaves the float altitude band.
        # For Grand Slam, "after the float" is what matters — `i_exit` was
        # widened to the end of the trajectory, so use the actual float-end
        # index instead.
        descent_mask = np.zeros(n, dtype=bool)
        if launch_type == "grand_slam":
            flt_idx = np.where(float_mask)[0]
            descent_start = int(flt_idx[-1] + 1) if flt_idx.size else i_peak
        else:
            descent_start = i_exit
        descent_mask[descent_start:] = v[descent_start:] < -v_linear

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

        asc_m, flt_m, des_m, i_enter_sim, i_exit_sim = self._detect_phases(
            el_sim, v_sim, min_alt, launch_type=self.launch_type,
        )

        # For non-standard launches (helium-augmented, grand-slam) the balloon is
        # carried up by helium / a weather balloon — EarthSHAB's solar-balloon
        # ascent model doesn't apply. Blank out the ascent mask so all
        # ascent-derived metrics report N/A.
        skip_ascent = self.launch_type in self.NON_STANDARD_LAUNCH_TYPES
        if skip_ascent:
            asc_m = np.zeros_like(asc_m)

        def _mean(arr, mask): return float(np.mean(arr[mask])) if mask.any() else math.nan
        def _std(arr, mask):  return float(np.std(arr[mask]))  if mask.any() else math.nan

        land_t_sim, land_i_sim = self._find_landing(el_sim, t_sim, min_alt)

        elapsed_sim = (
            pd.Timestamp(land_t_sim) - pd.Timestamp(t_sim[0])
        ).total_seconds() / 60.0

        flt_idx_sim = np.where(flt_m)[0]
        time_to_float_sim = (
            (pd.Timestamp(t_sim[flt_idx_sim[0]]) - pd.Timestamp(t_sim[0])).total_seconds() / 60.0
            if flt_idx_sim.size and not skip_ascent else math.nan
        )

        sim_metrics = FlightMetrics(
            float_alt_mean=_mean(el_sim, flt_m),
            float_alt_std=_std(el_sim,  flt_m),
            ascent_rate_mean=_mean(v_sim, asc_m),
            ascent_rate_std=_std(v_sim,  asc_m),
            descent_rate_mean=_mean(v_sim, des_m),
            descent_rate_std=_std(v_sim,  des_m),
            elapsed_time_min=elapsed_sim,
            time_to_float_min=time_to_float_sim,
            landing_lat=ss["lat"][land_i_sim],
            landing_lon=ss["lon"][land_i_sim],
            landing_time=land_t_sim,
        )

        # ── Ground truth ─────────────────────────────────────────────────────
        el_truth = df["altitude"].to_numpy(dtype=float)
        v_truth  = df["v_truth"].to_numpy(dtype=float)
        t_truth  = df["time"].tolist()

        # Smooth before phase detection — raw APRS finite-differences are noisy
        # (altitude quantization causes many zero-velocity readings during ascent/descent)
        v_truth_smooth = (pd.Series(v_truth)
                          .rolling(window=5, center=True, min_periods=1)
                          .median().to_numpy())
        asc_t, flt_t, des_t, i_enter_truth, i_exit_truth = self._detect_phases(
            el_truth, v_truth_smooth, min_alt, launch_type=self.launch_type,
        )

        if skip_ascent:
            asc_t = np.zeros_like(asc_t)

        elapsed_truth = (
            pd.Timestamp(t_truth[-1]) - pd.Timestamp(t_truth[0])
        ).total_seconds() / 60.0

        flt_idx_truth = np.where(flt_t)[0]
        time_to_float_truth = (
            (pd.Timestamp(t_truth[flt_idx_truth[0]]) - pd.Timestamp(t_truth[0])).total_seconds() / 60.0
            if flt_idx_truth.size and not skip_ascent else math.nan
        )

        truth_metrics = FlightMetrics(
            float_alt_mean=_mean(el_truth, flt_t),
            float_alt_std=_std(el_truth,  flt_t),
            ascent_rate_mean=_mean(v_truth, asc_t),
            ascent_rate_std=_std(v_truth,  asc_t),
            descent_rate_mean=_mean(v_truth, des_t),
            descent_rate_std=_std(v_truth,  des_t),
            elapsed_time_min=elapsed_truth,
            time_to_float_min=time_to_float_truth,
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

        def _diff(a, b):
            return a - b if not (math.isnan(a) or math.isnan(b)) else math.nan

        time_to_float_diff  = _diff(sim_metrics.time_to_float_min,  truth_metrics.time_to_float_min)
        time_to_ground_diff = _diff(sim_metrics.elapsed_time_min,   truth_metrics.elapsed_time_min)

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
            time_to_float_diff_min=time_to_float_diff,
            time_to_ground_diff_min=time_to_ground_diff,
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

        traj_path = config_earth.simulation.get('balloon_trajectory', '') or ''
        reporting.print_report(
            result,
            df_truth=self.df_truth,
            traj_path=traj_path,
            current_start_utc=config_earth.simulation['start_time'],
            min_alt=config_earth.simulation['min_alt'],
            gmt=self.GMT,
        )

        os.makedirs(self._output_dir, exist_ok=True)
        csv_path = os.path.join(self._output_dir, f"{self._stem()}.csv")
        reporting.write_metrics_csv(result, csv_path)
        print(f"  Report saved → {csv_path}")

    # ── Plot ─────────────────────────────────────────────────────────────────

    def _compute_sunset_sim(self):
        """Sunset local Timestamp along the simulated trajectory, or None."""
        ss   = self.sim_state
        gmt  = ss["GMT"]
        lats = ss["lat"]
        lons = ss["lon"]
        alts = ss["el"]
        step = 60
        for i_s, t_local in enumerate(ss["time_local"][::step]):
            i = i_s * step
            t_utc = pd.Timestamp(t_local) + pd.Timedelta(hours=gmt)
            if solar_zenith_adjusted(t_utc, lats[i], lons[i], alts[i]) >= math.pi / 2:
                return pd.Timestamp(t_local)
        return None

    def _compute_sunset_truth(self):
        """Sunset local Timestamp along the APRS truth trajectory, or None.

        Scans the truth trajectory in chronological order, using each APRS
        point's actual lat/lng/altitude.  The sort is required because aprs.fi
        exports are newest-first.  The UTC offset (7 h) is the same constant
        used by _load_aprs when converting timestamps — no model data is used.
        """
        _gmt = self.sim_state["GMT"]
        for _, row in self.df_truth.sort_values("time").iterrows():
            t_utc = row["time"] + pd.Timedelta(hours=_gmt)
            if solar_zenith_adjusted(t_utc, row["lat"], row["lng"], float(row["altitude"])) >= math.pi / 2:
                return row["time"]
        return None

    def plot(self):
        ss = self.sim_state
        os.makedirs(self._output_dir, exist_ok=True)
        stem = self._stem()

        fig = plot_comparison(
            time_sim=ss["time_local"],
            el_sim=ss["el"],
            v_sim=ss["v"],
            T_atm_sim=ss["T_atm"],
            df_truth=self.df_truth,
            sim_phases=getattr(self, '_sim_phases', None),
            truth_phases=getattr(self, '_truth_phases', None),
            t_sunset_sim=self._compute_sunset_sim(),
            t_sunset_truth=self._compute_sunset_truth(),
        )
        png_path = os.path.join(self._output_dir, f"{stem}.png")
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
            forecast=ss["forecast"],
            sim=self.sim,
            t=ss["t"],
            start=ss["start"],
            html_prefix="EVALUATION_",
            output_dir=self._output_dir,
        )
        print(f"  Map saved      → {os.path.join(self._output_dir, 'EVALUATION_' + stem + '.html')}")
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
