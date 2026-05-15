"""reporting.py - Reporting and CSV/HTML formatting for the evaluation suite.

All console pretty-printing, per-launch CSV writing, batch-summary CSV row
construction, and the interactive batch-summary HTML are in here so the runner
modules (`evaluate.py` and `run_batch.py`) can stay focused on *running*
simulations.

Public API:
    print_report(result, df_truth, traj_path, current_start_utc, min_alt, gmt)
    write_metrics_csv(result, path)
    suggest_start_time(df, min_alt, gmt)

    SUMMARY_FIELDNAMES
    result_to_summary_row(...)
    write_summary_html(summary_rows, batch_dir, batch_id, note, git)
"""

import csv
import html as _html
import math
import os
from typing import Optional

import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# Per-launch console report + CSV
# ─────────────────────────────────────────────────────────────────────────────

def suggest_start_time(df: pd.DataFrame, min_alt: float, gmt: int):
    """Estimate actual launch time by linearly extrapolating the initial APRS
    ascent rate back to min_alt.

    Returns (suggested_utc, first_aprs_utc, v_mean_ms) or (None, None, NaN)
    if there is insufficient ascending data.
    """
    ascending = df.iloc[1:][df.iloc[1:]['v_truth'] > 0.5].head(5)
    if ascending.empty:
        return None, None, math.nan

    v_mean    = float(ascending['v_truth'].mean())
    first_alt = float(df['altitude'].iloc[0])
    first_mst = pd.Timestamp(df['time'].iloc[0])

    dt_s = (first_alt - min_alt) / v_mean
    suggested_mst = first_mst - pd.Timedelta(seconds=dt_s)

    suggested_utc  = suggested_mst + pd.Timedelta(hours=gmt)
    first_aprs_utc = first_mst    + pd.Timedelta(hours=gmt)
    return suggested_utc, first_aprs_utc, v_mean


def print_report(result, df_truth: Optional[pd.DataFrame],
                 traj_path: str, current_start_utc, min_alt: float, gmt: int):
    """Pretty-print the evaluation metrics block and start-time analysis."""
    name = (traj_path or "").split("/")[-1].replace(".csv", "") or "Unknown"

    W = 68

    def _fmt(v, dec=1):
        return "N/A" if (isinstance(v, float) and math.isnan(v)) else f"{v:.{dec}f}"

    def _row(label, sv, tv, dv, unit):
        return f"  {label:<30} {sv:>9}  {tv:>9}  {dv:>9}  {unit}"

    s, t = result.sim, result.truth

    def _diff(a, b):
        return a - b if not (math.isnan(a) or math.isnan(b)) else math.nan

    rows = [
        ("Float Alt Mean (m)",      s.float_alt_mean,    t.float_alt_mean,    _diff(s.float_alt_mean,    t.float_alt_mean),    "m",   0),
        ("Float Alt Std (m)",       s.float_alt_std,     t.float_alt_std,     _diff(s.float_alt_std,     t.float_alt_std),     "m",   0),
        ("Ascent Rate Mean (m/s)",  s.ascent_rate_mean,  t.ascent_rate_mean,  _diff(s.ascent_rate_mean,  t.ascent_rate_mean),  "m/s", 2),
        ("Ascent Rate Std (m/s)",   s.ascent_rate_std,   t.ascent_rate_std,   _diff(s.ascent_rate_std,   t.ascent_rate_std),   "m/s", 2),
        ("Descent Rate Mean (m/s)", s.descent_rate_mean, t.descent_rate_mean, _diff(s.descent_rate_mean, t.descent_rate_mean), "m/s", 2),
        ("Descent Rate Std (m/s)",  s.descent_rate_std,  t.descent_rate_std,  _diff(s.descent_rate_std,  t.descent_rate_std),  "m/s", 2),
        ("Elapsed Time (min)",      s.elapsed_time_min,  t.elapsed_time_min,  _diff(s.elapsed_time_min,  t.elapsed_time_min),  "min", 1),
        ("Landing Lat (°)",         s.landing_lat,       t.landing_lat,       _diff(s.landing_lat,       t.landing_lat),       "°",   4),
        ("Landing Lon (°)",         s.landing_lon,       t.landing_lon,       _diff(s.landing_lon,       t.landing_lon),       "°",   4),
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
    print(_row("Distance Off (m)", "", "", _fmt(land_dist_m, 0), "m"))

    t_sim_str   = str(s.landing_time)[:16] if s.landing_time else "N/A"
    t_truth_str = str(t.landing_time)[:16] if t.landing_time else "N/A"
    print(_row("Landing Time (MST)", t_sim_str, t_truth_str,
               _fmt(result.landing_time_diff_min, 1), "min"))
    print("-" * W)
    print(_row("Temperature MAE", "", "", _fmt(result.temp_mae_k, 2), "K"))
    print(_row("Pressure MAE",    "", "", _fmt(result.pressure_mae_pa, 0), "Pa"))
    print("-" * W)
    print(f"  GFS Forecast + Truth Altitude (reforecast landing vs truth)")
    print(f"  {'Metric':<30} {'GFS+TA':>9}  {'Truth':>9}  {'Diff':>9}  Unit")
    print(_row("Landing Lat (°)",
               _fmt(result.gfs_truth_landing_lat, 4),
               _fmt(t.landing_lat, 4),
               _fmt(_diff(result.gfs_truth_landing_lat, t.landing_lat), 4), "°"))
    print(_row("Landing Lon (°)",
               _fmt(result.gfs_truth_landing_lon, 4),
               _fmt(t.landing_lon, 4),
               _fmt(_diff(result.gfs_truth_landing_lon, t.landing_lon), 4), "°"))
    print(_row("Distance Off (m)", "", "", _fmt(result.gfs_truth_landing_dist_m, 0), "m"))
    print("=" * W)

    if df_truth is None:
        return

    sugg_utc, first_aprs_utc, v_mean = suggest_start_time(df_truth, min_alt, gmt)
    print()
    print("  Start-time analysis")
    print(f"  Current config start_time : {current_start_utc} UTC")
    if first_aprs_utc is not None:
        print(f"  First APRS transmission   : {first_aprs_utc} UTC  "
              f"({df_truth['altitude'].iloc[0]:.0f} m)")
    if sugg_utc is not None:
        print(f"  Estimated launch time     : {sugg_utc.strftime('%Y-%m-%d %H:%M:%S')} UTC"
              f"  (extrapolated at {v_mean:.2f} m/s from first APRS point)")
        print(f"  Suggested start_time      : \"{sugg_utc.strftime('%Y-%m-%d %H:%M:%S')}\"")


def write_metrics_csv(result, path: str):
    """Write per-launch evaluation metrics to a structured CSV file."""
    s, t = result.sim, result.truth

    def _v(x, dec=4):
        return "" if (isinstance(x, float) and math.isnan(x)) else f"{x:.{dec}f}"

    def _diff(a, b):
        return a - b if not (math.isnan(a) or math.isnan(b)) else math.nan

    land_dist_m = (
        result.landing_distance_km * 1000
        if not math.isnan(result.landing_distance_km) else math.nan
    )

    rows = [
        ("Float Alt Mean",      "m",   _v(s.float_alt_mean, 0),    _v(t.float_alt_mean, 0),    _v(_diff(s.float_alt_mean,    t.float_alt_mean),    0)),
        ("Float Alt Std",       "m",   _v(s.float_alt_std,  0),    _v(t.float_alt_std,  0),    _v(_diff(s.float_alt_std,     t.float_alt_std),     0)),
        ("Ascent Rate Mean",    "m/s", _v(s.ascent_rate_mean, 2),  _v(t.ascent_rate_mean, 2),  _v(_diff(s.ascent_rate_mean,  t.ascent_rate_mean),  2)),
        ("Ascent Rate Std",     "m/s", _v(s.ascent_rate_std,  2),  _v(t.ascent_rate_std,  2),  _v(_diff(s.ascent_rate_std,   t.ascent_rate_std),   2)),
        ("Descent Rate Mean",   "m/s", _v(s.descent_rate_mean, 2), _v(t.descent_rate_mean, 2), _v(_diff(s.descent_rate_mean, t.descent_rate_mean), 2)),
        ("Descent Rate Std",    "m/s", _v(s.descent_rate_std,  2), _v(t.descent_rate_std,  2), _v(_diff(s.descent_rate_std,  t.descent_rate_std),  2)),
        ("Elapsed Time",        "min", _v(s.elapsed_time_min, 1),  _v(t.elapsed_time_min, 1),  _v(_diff(s.elapsed_time_min,  t.elapsed_time_min),  1)),
        ("Landing Lat",         "deg", _v(s.landing_lat, 4),       _v(t.landing_lat, 4),       _v(_diff(s.landing_lat,       t.landing_lat),       4)),
        ("Landing Lon",         "deg", _v(s.landing_lon, 4),       _v(t.landing_lon, 4),       _v(_diff(s.landing_lon,       t.landing_lon),       4)),
        ("Landing Time (MST)",  "",    str(s.landing_time)[:16] if s.landing_time else "",
                                        str(t.landing_time)[:16] if t.landing_time else "",
                                        _v(result.landing_time_diff_min, 1)),
        ("Distance Off",        "m",   "", "", _v(land_dist_m, 0)),
        ("Temperature MAE",     "K",   "", "", _v(result.temp_mae_k, 2)),
        ("Pressure MAE",        "Pa",  "", "", _v(result.pressure_mae_pa, 0)),
        ("GFS+TA Landing Lat",  "deg", _v(result.gfs_truth_landing_lat, 4),
                                        _v(t.landing_lat, 4),
                                        _v(_diff(result.gfs_truth_landing_lat, t.landing_lat), 4)),
        ("GFS+TA Landing Lon",  "deg", _v(result.gfs_truth_landing_lon, 4),
                                        _v(t.landing_lon, 4),
                                        _v(_diff(result.gfs_truth_landing_lon, t.landing_lon), 4)),
        ("GFS+TA Distance Off", "m",   "", "", _v(result.gfs_truth_landing_dist_m, 0)),
    ]

    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["Metric", "Unit", "Sim", "Truth", "Diff"])
        w.writerows(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Batch summary CSV row construction
# ─────────────────────────────────────────────────────────────────────────────

SUMMARY_FIELDNAMES = [
    "batch_id", "launch_id", "forecast_type", "status", "aprs_format",
    "campaign", "organization", "launch_date",
    "payload_weight_kg", "balloon_size_m",
    # sim FlightMetrics
    "sim_float_alt_mean_m", "sim_float_alt_std_m",
    "sim_ascent_rate_mean_ms", "sim_ascent_rate_std_ms",
    "sim_descent_rate_mean_ms", "sim_descent_rate_std_ms",
    "sim_elapsed_time_min",
    "sim_time_to_float_min",
    "sim_landing_lat_deg", "sim_landing_lon_deg", "sim_landing_time",
    # truth FlightMetrics
    "truth_float_alt_mean_m", "truth_float_alt_std_m",
    "truth_ascent_rate_mean_ms", "truth_ascent_rate_std_ms",
    "truth_descent_rate_mean_ms", "truth_descent_rate_std_ms",
    "truth_elapsed_time_min",
    "truth_time_to_float_min",
    "truth_landing_lat_deg", "truth_landing_lon_deg", "truth_landing_time",
    # cross metrics
    "landing_distance_km",
    "landing_time_diff_min",
    "time_to_float_diff_min",
    "time_to_ground_diff_min",
    "temp_mae_k",
    "pressure_mae_pa",
    "reforecast_landing_lat_deg",
    "reforecast_landing_lon_deg",
    "reforecast_landing_dist_m",
]


def _nan(v):
    return "" if (isinstance(v, float) and math.isnan(v)) else v


def result_to_summary_row(batch_id: str, launch_id: str, forecast_type: str,
                          result, status: str, aprs_format: str = "",
                          campaign: str = "", organization: str = "",
                          launch_date: str = "",
                          payload_weight_kg: float = math.nan,
                          balloon_size_m: float = math.nan) -> dict:
    """Build one row of summary.csv for the given evaluation result."""
    if result is None:
        return {f: "" for f in SUMMARY_FIELDNAMES} | {
            "batch_id":          batch_id,
            "launch_id":         launch_id,
            "forecast_type":     forecast_type,
            "status":            status,
            "aprs_format":       aprs_format,
            "campaign":          campaign,
            "organization":      organization,
            "launch_date":       launch_date,
            "payload_weight_kg": _nan(payload_weight_kg),
            "balloon_size_m":    _nan(balloon_size_m),
        }

    s, t = result.sim, result.truth
    return {
        "batch_id":                    batch_id,
        "launch_id":                   launch_id,
        "forecast_type":               forecast_type,
        "status":                      status,
        "aprs_format":                 aprs_format,
        "campaign":                    campaign,
        "organization":                organization,
        "launch_date":                 launch_date,
        "payload_weight_kg":           _nan(payload_weight_kg),
        "balloon_size_m":              _nan(balloon_size_m),
        "sim_float_alt_mean_m":        _nan(s.float_alt_mean),
        "sim_float_alt_std_m":         _nan(s.float_alt_std),
        "sim_ascent_rate_mean_ms":     _nan(s.ascent_rate_mean),
        "sim_ascent_rate_std_ms":      _nan(s.ascent_rate_std),
        "sim_descent_rate_mean_ms":    _nan(s.descent_rate_mean),
        "sim_descent_rate_std_ms":     _nan(s.descent_rate_std),
        "sim_elapsed_time_min":        _nan(s.elapsed_time_min),
        "sim_time_to_float_min":       _nan(s.time_to_float_min),
        "sim_landing_lat_deg":         _nan(s.landing_lat),
        "sim_landing_lon_deg":         _nan(s.landing_lon),
        "sim_landing_time":            str(s.landing_time)[:19] if s.landing_time else "",
        "truth_float_alt_mean_m":      _nan(t.float_alt_mean),
        "truth_float_alt_std_m":       _nan(t.float_alt_std),
        "truth_ascent_rate_mean_ms":   _nan(t.ascent_rate_mean),
        "truth_ascent_rate_std_ms":    _nan(t.ascent_rate_std),
        "truth_descent_rate_mean_ms":  _nan(t.descent_rate_mean),
        "truth_descent_rate_std_ms":   _nan(t.descent_rate_std),
        "truth_elapsed_time_min":      _nan(t.elapsed_time_min),
        "truth_time_to_float_min":     _nan(t.time_to_float_min),
        "truth_landing_lat_deg":       _nan(t.landing_lat),
        "truth_landing_lon_deg":       _nan(t.landing_lon),
        "truth_landing_time":          str(t.landing_time)[:19] if t.landing_time else "",
        "landing_distance_km":         _nan(result.landing_distance_km),
        "landing_time_diff_min":       _nan(result.landing_time_diff_min),
        "time_to_float_diff_min":      _nan(result.time_to_float_diff_min),
        "time_to_ground_diff_min":     _nan(result.time_to_ground_diff_min),
        "temp_mae_k":                  _nan(result.temp_mae_k),
        "pressure_mae_pa":             _nan(result.pressure_mae_pa),
        "reforecast_landing_lat_deg":  _nan(result.gfs_truth_landing_lat),
        "reforecast_landing_lon_deg":  _nan(result.gfs_truth_landing_lon),
        "reforecast_landing_dist_m":   _nan(result.gfs_truth_landing_dist_m),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Batch summary HTML report
# ─────────────────────────────────────────────────────────────────────────────

_GREEN  = "#b7e4b7"
_YELLOW = "#ffeaa0"
_RED    = "#f4a5a5"
_FAIL   = "#e8e8e8"


def _pct(s, t):
    try:
        sv, tv = float(s), float(t)
        if math.isnan(sv) or math.isnan(tv) or tv == 0:
            return math.nan
        return (sv - tv) / abs(tv) * 100.0
    except (ValueError, TypeError):
        return math.nan


def _fv(v, fmt=".1f"):
    try:
        fv = float(v)
        return "—" if math.isnan(fv) else f"{fv:{fmt}}"
    except (ValueError, TypeError):
        return str(v) if v else "—"


def _rv(v):
    """Raw numeric string for data-val; '' if not available."""
    if v in ("", None):
        return ""
    try:
        fv = float(v)
        return "" if math.isnan(fv) else str(fv)
    except (ValueError, TypeError):
        return str(v) if v else ""


def _pct_cell(pct, lo=10, hi=25):
    if math.isnan(pct):
        return "—", "", ""
    color = _GREEN if abs(pct) <= lo else _YELLOW if abs(pct) <= hi else _RED
    sign  = "+" if pct > 0 else ""
    return f"{sign}{pct:.1f}%", f' style="background:{color}"', str(pct)


def _abs_cell(val, lo, hi):
    if math.isnan(val):
        return "—", "", ""
    color = _GREEN if abs(val) <= lo else _YELLOW if abs(val) <= hi else _RED
    return f"{abs(val):.1f}", f' style="background:{color}"', str(val)


def _avg(col, rows_subset):
    """Mean of a numeric column across ok rows, NaN if none."""
    vals = []
    for r in rows_subset:
        if r.get("status") != "ok":
            continue
        v = r.get(col)
        if v in ("", None):
            continue
        try:
            fv = float(v)
            if not math.isnan(fv):
                vals.append(fv)
        except (ValueError, TypeError):
            pass
    return sum(vals) / len(vals) if vals else math.nan


def _td(raw, txt, style=""):
    return f'<td data-val="{_html.escape(str(raw))}"{style}>{txt}</td>'


def _signed(v):
    """Format a signed minute difference as '+X' / '-X'."""
    if math.isnan(v):
        return "—"
    return f"{v:+.0f}"


def _build_data_row(row: dict, orig_idx: int) -> str:
    status  = row.get("status", "")
    lid_raw = str(row.get("launch_id", ""))
    ft_raw  = str(row.get("forecast_type", ""))
    lid     = _html.escape(lid_raw)
    ft      = _html.escape(ft_raw)

    if status != "ok":
        msg = _html.escape(str(status)[:80])
        return (
            f'<tr class="failed-row" data-orig="{orig_idx}" style="background:{_FAIL}">'
            f'<td data-val="{_html.escape(lid_raw)}" style="text-align:left">{lid}</td>'
            f'<td data-val="{_html.escape(ft_raw)}">{ft}</td>'
            f'<td colspan="21" style="text-align:left">FAILED: {msg}</td></tr>'
        )

    def _f(key):
        try:    return float(row.get(key, math.nan))
        except: return math.nan

    pct_alt = _pct(row.get("sim_float_alt_mean_m"),     row.get("truth_float_alt_mean_m"))
    pct_asc = _pct(row.get("sim_ascent_rate_mean_ms"),  row.get("truth_ascent_rate_mean_ms"))
    pct_des = _pct(row.get("sim_descent_rate_mean_ms"), row.get("truth_descent_rate_mean_ms"))

    ttf_diff = _f("time_to_float_diff_min")
    ttg_diff = _f("time_to_ground_diff_min")
    land_km  = _f("landing_distance_km")
    tdiff    = _f("landing_time_diff_min")
    temp_mae = _f("temp_mae_k")
    pres_mae = _f("pressure_mae_pa")

    pa_txt, pa_s, pa_r = _pct_cell(pct_alt)
    pc_txt, pc_s, pc_r = _pct_cell(pct_asc)
    pd_txt, pd_s, pd_r = _pct_cell(pct_des)
    tf_txt, tf_s, tf_r = _abs_cell(abs(ttf_diff) if not math.isnan(ttf_diff) else math.nan, 15, 45)
    tg_txt, tg_s, tg_r = _abs_cell(abs(ttg_diff) if not math.isnan(ttg_diff) else math.nan, 30, 90)
    ld_txt, ld_s, ld_r = _abs_cell(land_km, 20, 50)
    td_txt, td_s, td_r = _abs_cell(abs(tdiff) if not math.isnan(tdiff) else math.nan, 30, 90)
    tm_txt, tm_s, tm_r = _abs_cell(temp_mae, 5, 15)
    pm_txt, pm_s, pm_r = _abs_cell(pres_mae, 500, 2000)

    return (
        f'<tr data-orig="{orig_idx}">'
        f'<td data-val="{_html.escape(lid_raw)}" style="text-align:left">{lid}</td>'
        f'<td data-val="{_html.escape(ft_raw)}">{ft}</td>'
        + _td(_rv(row.get("payload_weight_kg")), _fv(row.get("payload_weight_kg"), ".2f"))
        + _td(_rv(row.get("balloon_size_m")),    _fv(row.get("balloon_size_m"),    ".1f"))
        + _td(_rv(row.get("sim_float_alt_mean_m")),    _fv(row.get("sim_float_alt_mean_m"),    ".0f"))
        + _td(_rv(row.get("truth_float_alt_mean_m")),  _fv(row.get("truth_float_alt_mean_m"),  ".0f"))
        + _td(pa_r, pa_txt, pa_s)
        + _td(_rv(row.get("sim_ascent_rate_mean_ms")),  _fv(row.get("sim_ascent_rate_mean_ms"),  ".2f"))
        + _td(_rv(row.get("truth_ascent_rate_mean_ms")), _fv(row.get("truth_ascent_rate_mean_ms"), ".2f"))
        + _td(pc_r, pc_txt, pc_s)
        + _td(_rv(row.get("sim_descent_rate_mean_ms")),  _fv(row.get("sim_descent_rate_mean_ms"),  ".2f"))
        + _td(_rv(row.get("truth_descent_rate_mean_ms")), _fv(row.get("truth_descent_rate_mean_ms"), ".2f"))
        + _td(pd_r, pd_txt, pd_s)
        + _td(_rv(row.get("sim_time_to_float_min")),   _fv(row.get("sim_time_to_float_min"),   ".0f"))
        + _td(_rv(row.get("truth_time_to_float_min")), _fv(row.get("truth_time_to_float_min"), ".0f"))
        + _td(tf_r, _signed(ttf_diff), tf_s)
        + _td(_rv(row.get("sim_elapsed_time_min")),    _fv(row.get("sim_elapsed_time_min"),    ".0f"))
        + _td(_rv(row.get("truth_elapsed_time_min")),  _fv(row.get("truth_elapsed_time_min"),  ".0f"))
        + _td(tg_r, _signed(ttg_diff), tg_s)
        + _td(ld_r, ld_txt, ld_s)
        + _td(td_r, td_txt, td_s)
        + _td(tm_r, tm_txt, tm_s)
        + _td(pm_r, pm_txt, pm_s)
        + '</tr>'
    )


def _build_avg_row(label: str, rows_subset: list) -> str:
    """Averages footer row for ok rows in rows_subset."""
    a = lambda col: _avg(col, rows_subset)

    avg_sa  = a("sim_float_alt_mean_m");    avg_ta  = a("truth_float_alt_mean_m")
    avg_sas = a("sim_ascent_rate_mean_ms"); avg_tas = a("truth_ascent_rate_mean_ms")
    avg_sd  = a("sim_descent_rate_mean_ms"); avg_td = a("truth_descent_rate_mean_ms")
    avg_sf  = a("sim_time_to_float_min");   avg_tf  = a("truth_time_to_float_min")
    avg_se  = a("sim_elapsed_time_min");    avg_te  = a("truth_elapsed_time_min")

    pct_alt = _pct(avg_sa,  avg_ta)
    pct_asc = _pct(avg_sas, avg_tas)
    pct_des = _pct(avg_sd,  avg_td)

    avg_ttf_diff = a("time_to_float_diff_min")
    avg_ttg_diff = a("time_to_ground_diff_min")
    avg_land     = a("landing_distance_km")
    avg_tdiff    = a("landing_time_diff_min")
    avg_tmae     = a("temp_mae_k")
    avg_pmae     = a("pressure_mae_pa")

    pa_txt, pa_s, _ = _pct_cell(pct_alt)
    pc_txt, pc_s, _ = _pct_cell(pct_asc)
    pd_txt, pd_s, _ = _pct_cell(pct_des)
    tf_txt, tf_s, _ = _abs_cell(abs(avg_ttf_diff) if not math.isnan(avg_ttf_diff) else math.nan, 15, 45)
    tg_txt, tg_s, _ = _abs_cell(abs(avg_ttg_diff) if not math.isnan(avg_ttg_diff) else math.nan, 30, 90)
    ld_txt, ld_s, _ = _abs_cell(avg_land, 20, 50)
    td_txt, td_s, _ = _abs_cell(abs(avg_tdiff) if not math.isnan(avg_tdiff) else math.nan, 30, 90)
    tm_txt, tm_s, _ = _abs_cell(avg_tmae, 5, 15)
    pm_txt, pm_s, _ = _abs_cell(avg_pmae, 500, 2000)

    avg_payload = a("payload_weight_kg")
    avg_balsize = a("balloon_size_m")

    return (
        f'<tr>'
        f'<td colspan="2" style="text-align:left">{_html.escape(label)}</td>'
        f'<td>{_fv(avg_payload, ".2f")}</td>'
        f'<td>{_fv(avg_balsize, ".1f")}</td>'
        f'<td>{_fv(avg_sa,  ".0f")}</td><td>{_fv(avg_ta,  ".0f")}</td><td{pa_s}>{pa_txt}</td>'
        f'<td>{_fv(avg_sas, ".2f")}</td><td>{_fv(avg_tas, ".2f")}</td><td{pc_s}>{pc_txt}</td>'
        f'<td>{_fv(avg_sd,  ".2f")}</td><td>{_fv(avg_td,  ".2f")}</td><td{pd_s}>{pd_txt}</td>'
        f'<td>{_fv(avg_sf,  ".0f")}</td><td>{_fv(avg_tf,  ".0f")}</td><td{tf_s}>{_signed(avg_ttf_diff)}</td>'
        f'<td>{_fv(avg_se,  ".0f")}</td><td>{_fv(avg_te,  ".0f")}</td><td{tg_s}>{_signed(avg_ttg_diff)}</td>'
        f'<td{ld_s}>{ld_txt}</td>'
        f'<td{td_s}>{td_txt}</td>'
        f'<td{tm_s}>{tm_txt}</td>'
        f'<td{pm_s}>{pm_txt}</td>'
        f'</tr>'
    )


_TABLE_HEADER = (
    '<thead>'
    '<tr>'
    '<th data-col="0"  rowspan="2" style="text-align:left">Launch<span class="si"></span></th>'
    '<th data-col="1"  rowspan="2">Fcst<span class="si"></span></th>'
    '<th data-col="2"  rowspan="2">Payload<br>(kg)<span class="si"></span></th>'
    '<th data-col="3"  rowspan="2">Bal&nbsp;&Oslash;<br>(m)<span class="si"></span></th>'
    '<th colspan="3">Float Alt (m)</th>'
    '<th colspan="3">Ascent (m/s)</th>'
    '<th colspan="3">Descent (m/s)</th>'
    '<th colspan="3">Time to Float (min)</th>'
    '<th colspan="3">Time to Ground (min)</th>'
    '<th data-col="22" rowspan="2">Land Dist<br>(km)<span class="si"></span></th>'
    '<th data-col="23" rowspan="2">|Time &Delta;|<br>(min)<span class="si"></span></th>'
    '<th data-col="24" rowspan="2">Temp MAE<br>(K)<span class="si"></span></th>'
    '<th data-col="25" rowspan="2">Press MAE<br>(Pa)<span class="si"></span></th>'
    '</tr>'
    '<tr>'
    '<th data-col="4"  class="sub">Sim<span class="si"></span></th>'
    '<th data-col="5"  class="sub">Truth<span class="si"></span></th>'
    '<th data-col="6"  class="sub">%&Delta;<span class="si"></span></th>'
    '<th data-col="7"  class="sub">Sim<span class="si"></span></th>'
    '<th data-col="8"  class="sub">Truth<span class="si"></span></th>'
    '<th data-col="9"  class="sub">%&Delta;<span class="si"></span></th>'
    '<th data-col="10" class="sub">Sim<span class="si"></span></th>'
    '<th data-col="11" class="sub">Truth<span class="si"></span></th>'
    '<th data-col="12" class="sub">%&Delta;<span class="si"></span></th>'
    '<th data-col="13" class="sub">Sim<span class="si"></span></th>'
    '<th data-col="14" class="sub">Truth<span class="si"></span></th>'
    '<th data-col="15" class="sub">&Delta;(min)<span class="si"></span></th>'
    '<th data-col="16" class="sub">Sim<span class="si"></span></th>'
    '<th data-col="17" class="sub">Truth<span class="si"></span></th>'
    '<th data-col="18" class="sub">&Delta;(min)<span class="si"></span></th>'
    '</tr>'
    '</thead>'
)


_SORT_JS = """\
<script>
(function () {
  var states = {};
  document.querySelectorAll('th[data-col]').forEach(function (th) {
    th.style.cursor = 'pointer';
    th.title = 'Click to sort';
    th.addEventListener('click', function () {
      var table = th.closest('table');
      var tid   = table.id;
      var col   = parseInt(th.dataset.col, 10);
      var st    = states[tid] || {col: -1, dir: 0};

      table.querySelectorAll('th[data-col] .si').forEach(function (s) { s.textContent = ''; });

      if (st.col === col) {
        st.dir = st.dir === 1 ? -1 : (st.dir === -1 ? 0 : 1);
      } else {
        st.col = col; st.dir = 1;
      }
      states[tid] = st;

      var icon = th.querySelector('.si');
      if (icon) icon.textContent = st.dir === 1 ? ' ↑' : (st.dir === -1 ? ' ↓' : '');

      var tbody = table.querySelector('tbody');
      var rows  = Array.from(tbody.querySelectorAll('tr'));

      if (st.dir === 0) {
        rows.sort(function (a, b) {
          return parseInt(a.dataset.orig || 0, 10) - parseInt(b.dataset.orig || 0, 10);
        });
      } else {
        rows.sort(function (a, b) {
          var fa = a.classList.contains('failed-row');
          var fb = b.classList.contains('failed-row');
          if (fa && fb) return 0;
          if (fa) return 1;
          if (fb) return -1;
          var ca = a.cells[col], cb = b.cells[col];
          var da = ca ? (ca.dataset.val || '') : '';
          var db = cb ? (cb.dataset.val || '') : '';
          var na = parseFloat(da), nb = parseFloat(db);
          var an = isNaN(na), bn = isNaN(nb);
          if (an && bn) return da.localeCompare(db) * st.dir;
          if (an) return 1;
          if (bn) return -1;
          return (na - nb) * st.dir;
        });
      }
      rows.forEach(function (r) { tbody.appendChild(r); });
    });
  });
})();
</script>"""


def write_summary_html(summary_rows: list, batch_dir: str,
                       batch_id: str, note: str, git: dict) -> str:
    """Write summary.html with colored % diff table, averages row, and sortable columns."""
    overall_data = "".join(_build_data_row(r, i) for i, r in enumerate(summary_rows))
    overall_avg  = _build_avg_row("Average", summary_rows)

    def _group_key(row):
        camp = row.get("campaign") or ""
        if camp:
            return camp
        org  = row.get("organization") or ""
        date = row.get("launch_date") or ""
        return f"{org} {date}" if (org or date) else ""

    groups: dict[str, list[dict]] = {}
    for row in summary_rows:
        if row.get("status") != "ok":
            continue
        key = _group_key(row)
        if key:
            groups.setdefault(key, []).append(row)

    campaign_html = ""
    tbl_n = 1
    for camp_key, rows in sorted(groups.items()):
        if len({r["launch_id"] for r in rows}) < 2:
            continue
        tbl_id      = f"tbl-camp-{tbl_n}"
        tbl_n      += 1
        camp_data   = "".join(_build_data_row(r, i) for i, r in enumerate(rows))
        camp_avg    = _build_avg_row("Campaign Average", rows)
        n_flights   = len({r["launch_id"] for r in rows})
        n_forecasts = len(rows)
        campaign_html += (
            f'<div class="campaign-card">'
            f'<h2>Campaign: {_html.escape(camp_key)}</h2>'
            f'<p style="color:#666;font-size:0.85em">'
            f'{n_flights} balloon(s), {n_forecasts} evaluation(s)</p>'
            f'<table id="{tbl_id}">{_TABLE_HEADER}'
            f'<tbody>{camp_data}</tbody>'
            f'<tfoot>{camp_avg}</tfoot>'
            f'</table></div>'
        )

    dirty_note = " <span style='color:orange'>(dirty)</span>" if git.get("git_dirty") else ""
    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>EarthSHAB Batch: {_html.escape(batch_id)}</title>
  <style>
    body      {{ font-family: monospace; margin: 24px; background: #f5f5f5; color: #222; }}
    h1        {{ margin-bottom: 4px; }}
    h2        {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 4px; margin-top: 0; }}
    .meta     {{ margin-bottom: 18px; font-size: 0.9em; color: #555; line-height: 1.7; }}
    .legend   {{ display:flex; gap:10px; flex-wrap:wrap; margin-bottom: 18px; font-size:0.85em; align-items:center; }}
    .legend span {{ padding: 3px 10px; border: 1px solid #bbb; border-radius: 3px; }}
    table     {{ border-collapse: collapse; width: 100%; margin-bottom: 6px; font-size: 0.82em; }}
    th, td    {{ border: 1px solid #bbb; padding: 5px 8px; text-align: center; white-space: nowrap; }}
    th        {{ background: #2c3e50; color: #fff; font-weight: bold; user-select: none; }}
    th.sub    {{ background: #34495e; }}
    th[data-col]:hover {{ filter: brightness(1.2); }}
    .si       {{ font-size: 0.8em; }}
    tr:hover  {{ filter: brightness(0.96); }}
    tfoot tr  {{ font-weight: bold; border-top: 2px solid #555; background: #eef; }}
    tfoot tr:hover {{ filter: none; }}
    .campaign-card {{
      background: #fff; border: 1px solid #ccc; border-radius: 6px;
      padding: 20px 24px; margin-bottom: 32px; box-shadow: 0 2px 6px rgba(0,0,0,0.08);
    }}
    .overall-section {{
      background: #fff; border: 1px solid #ccc; border-radius: 6px;
      padding: 20px 24px; margin-bottom: 32px; box-shadow: 0 2px 6px rgba(0,0,0,0.08);
    }}
  </style>
</head>
<body>
  <h1>EarthSHAB Batch Evaluation</h1>
  <div class="meta">
    <b>Batch ID:</b> {_html.escape(batch_id)}<br>
    <b>Note:</b> {_html.escape(note)}<br>
    <b>Git:</b> {_html.escape(git.get("git_hash",""))} &mdash; {_html.escape(git.get("git_commit_message",""))}{dirty_note}<br>
    <b>Branch:</b> {_html.escape(git.get("git_branch",""))}
  </div>

  <div class="legend">
    <b>Color key (%&Delta; from truth):</b>
    <span style="background:{_GREEN}">&le;&thinsp;10%</span>
    <span style="background:{_YELLOW}">10&ndash;25%</span>
    <span style="background:{_RED}">&gt;&thinsp;25%</span>
    &nbsp;&nbsp;<b>Error metrics:</b>
    <span style="background:{_GREEN}">good</span>
    <span style="background:{_YELLOW}">marginal</span>
    <span style="background:{_RED}">poor</span>
    &nbsp;&nbsp;<span style="color:#666">Click any column header to sort &uarr;&darr;</span>
  </div>

  <div class="overall-section">
    <h2>Overall Summary</h2>
    <table id="tbl-overall">
      {_TABLE_HEADER}
      <tbody>{overall_data}</tbody>
      <tfoot>{overall_avg}</tfoot>
    </table>
  </div>

  {campaign_html}
  {_SORT_JS}
</body>
</html>"""

    path = os.path.join(batch_dir, "summary.html")
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    return path
