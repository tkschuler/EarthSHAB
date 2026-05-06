"""run_batch.py - Run all launches in launches.json against the current codebase.

Each batch is tagged with the current git hash and a timestamp, and outputs are
written to a dedicated folder so results from different code states can be compared.

Usage:
    python -m evaluation.run_batch --note "describe what changed"

Output structure:
    Evaluation/batches/<timestamp>_<git_hash>/
        batch_info.json
        summary.csv
        <org>_<shab>_<date>/
            <stem>_GFS.csv / <stem>_ERA5.csv  (per-forecast metrics)
            <stem>_GFS.png / ...               (comparison plots)
            EVALUATION_<stem>_GFS.html / ...   (trajectory maps)
"""

import argparse
import copy
import csv
import json
import math
import os
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone

import pandas as pd

import EarthSHAB.config_earth as config_earth
from evaluation.evaluate import BalloonEvaluator

# ── Paths ─────────────────────────────────────────────────────────────────────

BALLOON_DATA_DIR = "src/EarthSHAB/balloon_data/"
FORECASTS_DIR    = "src/EarthSHAB/forecasts/"
EVAL_DIR         = "evaluation/"
BATCHES_DIR      = EVAL_DIR + "batches/"
LAUNCHES_JSON    = EVAL_DIR + "launches.json"

REQUIRED_FIELDS = [
    "shab_name", "organization", "launch_time", "sim_time_hr",
    "aprs_file", "launch_lat", "launch_lon", "launch_alt_m",
    "payload_weight_kg", "envelope_weight_kg", "balloon_shape", "balloon_size",
]

# ── Git helpers ───────────────────────────────────────────────────────────────

def _git(cmd: list[str]) -> str:
    try:
        return subprocess.check_output(["git"] + cmd, stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return ""


def _git_info() -> dict:
    hash_short  = _git(["rev-parse", "--short", "HEAD"])
    hash_full   = _git(["rev-parse", "HEAD"])
    branch      = _git(["rev-parse", "--abbrev-ref", "HEAD"])
    commit_msg  = _git(["log", "-1", "--pretty=%s"])
    dirty       = bool(_git(["status", "--porcelain"]))
    return {
        "git_hash":           hash_short,
        "git_hash_full":      hash_full,
        "git_branch":         branch,
        "git_commit_message": commit_msg,
        "git_dirty":          dirty,
    }


def _earthshab_version() -> str:
    try:
        from importlib.metadata import version
        return version("EarthSHAB")
    except Exception:
        return "unknown"

# ── Config helpers ────────────────────────────────────────────────────────────

def _snapshot_config() -> dict:
    """Deep-copy the current config_earth state so we can restore it per launch."""
    return {
        "balloon_properties": copy.deepcopy(config_earth.balloon_properties),
        "earth_properties":   copy.deepcopy(config_earth.earth_properties),
        "simulation":         copy.deepcopy(config_earth.simulation),
        "forecast":           copy.deepcopy(config_earth.forecast),
        "netcdf_gfs":         copy.deepcopy(config_earth.netcdf_gfs),
        "netcdf_era5":        copy.deepcopy(config_earth.netcdf_era5),
    }


def _parse_gfs_filename(filename: str) -> tuple[str, str]:
    """Return (forecast_start_str, hour_str) from a GFS filename.

    Handles both gfs_0p25_YYYYMMDD_HH.nc and gfs_0p25_YYYYMMDD_HH_<suffix>.nc.
    """
    stem  = filename.replace(".nc", "")
    parts = stem.split("_")
    date_part = parts[2]   # e.g. "20220822"
    hour_part = parts[3]   # e.g. "12"
    forecast_start = (
        f"{date_part[:4]}-{date_part[4:6]}-{date_part[6:8]} {hour_part}:00:00"
    )
    return forecast_start, hour_part


def _build_overrides(launch: dict, forecast_type: str, orig: dict) -> dict:
    """Build a complete config_overrides dict for one launch + forecast type.

    Always includes every field so previous launches cannot bleed through.
    """
    launch_time = datetime.fromisoformat(launch["launch_time"])
    aprs_path   = BALLOON_DATA_DIR + launch["aprs_file"]

    # Balloon properties: start from original defaults, apply metadata on top
    balloon_props = dict(orig["balloon_properties"])
    balloon_props.update({
        "shape": launch["balloon_shape"],
        "d":     launch["balloon_size"],
        "mp":    launch["payload_weight_kg"],
        "mEnv":  launch["envelope_weight_kg"],
    })
    for key in ["areaDensityEnv", "cp", "absEnv", "emissEnv", "Upsilon"]:
        if key in launch:
            balloon_props[key] = launch[key]

    # Earth properties: start from original defaults, apply metadata on top
    earth_props = dict(orig["earth_properties"])
    for key in ["Cp_air0", "Cv_air0", "Rsp_air", "P0", "emissGround", "albedo"]:
        if key in launch:
            earth_props[key] = launch[key]

    start_coord = {
        "lat":       launch["launch_lat"],
        "lon":       launch["launch_lon"],
        "alt":       launch["launch_alt_m"],
        "timestamp": launch_time,
    }

    overrides = {
        "balloon_properties": balloon_props,
        "earth_properties":   earth_props,
        "simulation": {
            "start_time":         launch_time,
            "sim_time":           launch["sim_time_hr"],
            "vent":               orig["simulation"]["vent"],
            "alt_sp":             orig["simulation"]["alt_sp"],
            "v_sp":               orig["simulation"]["v_sp"],
            "start_coord":        start_coord,
            "min_alt":            launch["launch_alt_m"],
            "float":              orig["simulation"]["float"],
            "dt":                 orig["simulation"]["dt"],
            "balloon_trajectory": aprs_path,
        },
    }

    if forecast_type == "GFS":
        gfs_file           = launch["gfs_file"]
        gfs_path           = FORECASTS_DIR + gfs_file
        forecast_start, hour_str = _parse_gfs_filename(gfs_file)
        forecast_start_dt  = datetime.fromisoformat(forecast_start)

        overrides["forecast"] = {
            "forecast_type":      "GFS",
            "forecast_start_time": forecast_start,
            "GFSrate":            orig["forecast"]["GFSrate"],
        }
        overrides["netcdf_gfs"] = {
            "nc_file":       gfs_path,
            "nc_start":      forecast_start_dt,
            "hourstamp":     hour_str,
            "res":           orig["netcdf_gfs"]["res"],
            "lat_range":     orig["netcdf_gfs"]["lat_range"],
            "lon_range":     orig["netcdf_gfs"]["lon_range"],
            "download_days": orig["netcdf_gfs"]["download_days"],
        }

    elif forecast_type == "ERA5":
        overrides["forecast"] = {
            "forecast_type":      "ERA5",
            "forecast_start_time": orig["forecast"]["forecast_start_time"],
            "GFSrate":            orig["forecast"]["GFSrate"],
        }
        # ERA5.py prepends "src/EarthSHAB/forecasts/" itself — filename only
        overrides["netcdf_era5"] = {
            "filename":       launch["era5_file"],
            "resolution_hr":  orig["netcdf_era5"]["resolution_hr"],
        }

    return overrides

# ── Summary CSV ───────────────────────────────────────────────────────────────

_SUMMARY_FIELDNAMES = [
    "batch_id", "launch_id", "forecast_type", "status", "aprs_format",
    "campaign", "organization", "launch_date",
    # sim FlightMetrics
    "sim_float_alt_mean_m", "sim_float_alt_std_m",
    "sim_ascent_rate_mean_ms", "sim_ascent_rate_std_ms",
    "sim_descent_rate_mean_ms", "sim_descent_rate_std_ms",
    "sim_elapsed_time_min",
    "sim_landing_lat_deg", "sim_landing_lon_deg", "sim_landing_time",
    # truth FlightMetrics
    "truth_float_alt_mean_m", "truth_float_alt_std_m",
    "truth_ascent_rate_mean_ms", "truth_ascent_rate_std_ms",
    "truth_descent_rate_mean_ms", "truth_descent_rate_std_ms",
    "truth_elapsed_time_min",
    "truth_landing_lat_deg", "truth_landing_lon_deg", "truth_landing_time",
    # cross metrics
    "landing_distance_km",
    "landing_time_diff_min",
    "temp_mae_k",
    "pressure_mae_pa",
    "reforecast_landing_lat_deg",
    "reforecast_landing_lon_deg",
    "reforecast_landing_dist_m",
]


def _nan(v):
    return "" if (isinstance(v, float) and math.isnan(v)) else v


def _result_to_row(batch_id: str, launch_id: str, forecast_type: str,
                   result, status: str, aprs_format: str = "",
                   campaign: str = "", organization: str = "",
                   launch_date: str = "") -> dict:
    if result is None:
        return {f: "" for f in _SUMMARY_FIELDNAMES} | {
            "batch_id": batch_id,
            "launch_id": launch_id,
            "forecast_type": forecast_type,
            "status": status,
            "aprs_format": aprs_format,
            "campaign": campaign,
            "organization": organization,
            "launch_date": launch_date,
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
        "sim_float_alt_mean_m":        _nan(s.float_alt_mean),
        "sim_float_alt_std_m":         _nan(s.float_alt_std),
        "sim_ascent_rate_mean_ms":     _nan(s.ascent_rate_mean),
        "sim_ascent_rate_std_ms":      _nan(s.ascent_rate_std),
        "sim_descent_rate_mean_ms":    _nan(s.descent_rate_mean),
        "sim_descent_rate_std_ms":     _nan(s.descent_rate_std),
        "sim_elapsed_time_min":        _nan(s.elapsed_time_min),
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
        "truth_landing_lat_deg":       _nan(t.landing_lat),
        "truth_landing_lon_deg":       _nan(t.landing_lon),
        "truth_landing_time":          str(t.landing_time)[:19] if t.landing_time else "",
        "landing_distance_km":         _nan(result.landing_distance_km),
        "landing_time_diff_min":       _nan(result.landing_time_diff_min),
        "temp_mae_k":                  _nan(result.temp_mae_k),
        "pressure_mae_pa":             _nan(result.pressure_mae_pa),
        "reforecast_landing_lat_deg":  _nan(result.gfs_truth_landing_lat),
        "reforecast_landing_lon_deg":  _nan(result.gfs_truth_landing_lon),
        "reforecast_landing_dist_m":   _nan(result.gfs_truth_landing_dist_m),
    }

# ── HTML Summary ─────────────────────────────────────────────────────────────

def _write_summary_html(summary_rows: list[dict], batch_dir: str,
                        batch_id: str, note: str, git: dict) -> str:
    """Write summary.html: a colored percent-diff table + campaign sections."""
    import html as _html

    GREEN  = "#b7e4b7"
    YELLOW = "#ffeaa0"
    RED    = "#f4a5a5"
    FAIL   = "#e8e8e8"

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

    def _pct_cell(pct, lo=10, hi=25):
        if math.isnan(pct):
            return "—", ""
        color = GREEN if abs(pct) <= lo else YELLOW if abs(pct) <= hi else RED
        sign  = "+" if pct > 0 else ""
        return f"{sign}{pct:.1f}%", f' style="background:{color}"'

    def _abs_cell(val, lo, hi):
        if math.isnan(val):
            return "—", ""
        color = GREEN if abs(val) <= lo else YELLOW if abs(val) <= hi else RED
        return f"{abs(val):.1f}", f' style="background:{color}"'

    def _build_data_row(row: dict) -> str:
        status = row.get("status", "")
        lid    = _html.escape(str(row.get("launch_id", "")))
        ft     = _html.escape(str(row.get("forecast_type", "")))

        if status != "ok":
            msg = _html.escape(str(status)[:80])
            return (f'<tr style="background:{FAIL}">'
                    f'<td>{lid}</td><td>{ft}</td>'
                    f'<td colspan="12" style="text-align:left">FAILED: {msg}</td></tr>')

        def _float(key):
            try:    return float(row.get(key, math.nan))
            except: return math.nan

        pct_alt = _pct(row.get("sim_float_alt_mean_m"),     row.get("truth_float_alt_mean_m"))
        pct_asc = _pct(row.get("sim_ascent_rate_mean_ms"),  row.get("truth_ascent_rate_mean_ms"))
        pct_des = _pct(row.get("sim_descent_rate_mean_ms"), row.get("truth_descent_rate_mean_ms"))
        pct_ela = _pct(row.get("sim_elapsed_time_min"),     row.get("truth_elapsed_time_min"))

        land_km  = _float("landing_distance_km")
        tdiff    = _float("landing_time_diff_min")
        temp_mae = _float("temp_mae_k")
        pres_mae = _float("pressure_mae_pa")

        pa_txt, pa_s = _pct_cell(pct_alt)
        pc_txt, pc_s = _pct_cell(pct_asc)
        pd_txt, pd_s = _pct_cell(pct_des)
        pe_txt, pe_s = _pct_cell(pct_ela)
        ld_txt, ld_s = _abs_cell(land_km, 20, 50)
        td_txt, td_s = _abs_cell(abs(tdiff) if not math.isnan(tdiff) else math.nan, 30, 90)
        tm_txt, tm_s = _abs_cell(temp_mae, 5, 15)
        pm_txt, pm_s = _abs_cell(pres_mae, 500, 2000)

        sim_alt  = _fv(row.get("sim_float_alt_mean_m"),    ".0f")
        tru_alt  = _fv(row.get("truth_float_alt_mean_m"),  ".0f")
        sim_asc  = _fv(row.get("sim_ascent_rate_mean_ms"), ".2f")
        tru_asc  = _fv(row.get("truth_ascent_rate_mean_ms"), ".2f")
        sim_des  = _fv(row.get("sim_descent_rate_mean_ms"), ".2f")
        tru_des  = _fv(row.get("truth_descent_rate_mean_ms"), ".2f")
        sim_ela  = _fv(row.get("sim_elapsed_time_min"),    ".0f")
        tru_ela  = _fv(row.get("truth_elapsed_time_min"),  ".0f")

        return (
            f'<tr>'
            f'<td style="text-align:left">{lid}</td><td>{ft}</td>'
            f'<td>{sim_alt} / {tru_alt}</td><td{pa_s}>{pa_txt}</td>'
            f'<td>{sim_asc} / {tru_asc}</td><td{pc_s}>{pc_txt}</td>'
            f'<td>{sim_des} / {tru_des}</td><td{pd_s}>{pd_txt}</td>'
            f'<td>{sim_ela} / {tru_ela}</td><td{pe_s}>{pe_txt}</td>'
            f'<td{ld_s}>{ld_txt}</td>'
            f'<td{td_s}>{td_txt}</td>'
            f'<td{tm_s}>{tm_txt}</td>'
            f'<td{pm_s}>{pm_txt}</td>'
            f'</tr>'
        )

    TABLE_HEADER = """
    <table>
      <thead>
        <tr>
          <th rowspan="2" style="text-align:left">Launch</th>
          <th rowspan="2">Fcst</th>
          <th colspan="2">Float Alt (m)</th>
          <th colspan="2">Ascent (m/s)</th>
          <th colspan="2">Descent (m/s)</th>
          <th colspan="2">Elapsed (min)</th>
          <th rowspan="2">Land Dist<br>(km)</th>
          <th rowspan="2">|Time Δ|<br>(min)</th>
          <th rowspan="2">Temp MAE<br>(K)</th>
          <th rowspan="2">Press MAE<br>(Pa)</th>
        </tr>
        <tr>
          <th class="sub">Sim / Truth</th><th class="sub">%Δ</th>
          <th class="sub">Sim / Truth</th><th class="sub">%Δ</th>
          <th class="sub">Sim / Truth</th><th class="sub">%Δ</th>
          <th class="sub">Sim / Truth</th><th class="sub">%Δ</th>
        </tr>
      </thead>
      <tbody>"""

    # ── Overall table ─────────────────────────────────────────────────────────
    overall_rows = "".join(_build_data_row(r) for r in summary_rows)

    # ── Campaign sections ─────────────────────────────────────────────────────
    # Group by explicit campaign field; fall back to (organization, launch_date)
    # when multiple flights share the same org+date.
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

    # Only show campaign sections for groups with ≥2 distinct balloons
    campaign_html = ""
    for camp_key, rows in sorted(groups.items()):
        if len({r["launch_id"] for r in rows}) < 2:
            continue
        camp_rows = "".join(_build_data_row(r) for r in rows)

        # Aggregate averages
        def _avg(col, _rows=rows):
            vals = []
            for r in _rows:
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

        avg_land  = _avg("landing_distance_km")
        avg_tdiff = _avg("landing_time_diff_min")
        avg_tmae  = _avg("temp_mae_k")
        avg_pmae  = _avg("pressure_mae_pa")
        avg_palt  = _pct(_avg("sim_float_alt_mean_m"),     _avg("truth_float_alt_mean_m"))
        avg_pasc  = _pct(_avg("sim_ascent_rate_mean_ms"),  _avg("truth_ascent_rate_mean_ms"))
        avg_pdes  = _pct(_avg("sim_descent_rate_mean_ms"), _avg("truth_descent_rate_mean_ms"))
        avg_pela  = _pct(_avg("sim_elapsed_time_min"),     _avg("truth_elapsed_time_min"))

        def _avg_row(label, p_alt, p_asc, p_des, p_ela, land, tdiff, tmae, pmae):
            pa_t, pa_s = _pct_cell(p_alt)
            pc_t, pc_s = _pct_cell(p_asc)
            pd_t, pd_s = _pct_cell(p_des)
            pe_t, pe_s = _pct_cell(p_ela)
            ld_t, ld_s = _abs_cell(land, 20, 50)
            td_t, td_s = _abs_cell(abs(tdiff) if not math.isnan(tdiff) else math.nan, 30, 90)
            tm_t, tm_s = _abs_cell(tmae, 5, 15)
            pm_t, pm_s = _abs_cell(pmae, 500, 2000)
            return (
                f'<tr style="font-weight:bold; border-top:2px solid #555">'
                f'<td colspan="2" style="text-align:left">{_html.escape(label)}</td>'
                f'<td>—</td><td{pa_s}>{pa_t}</td>'
                f'<td>—</td><td{pc_s}>{pc_t}</td>'
                f'<td>—</td><td{pd_s}>{pd_t}</td>'
                f'<td>—</td><td{pe_s}>{pe_t}</td>'
                f'<td{ld_s}>{ld_t}</td>'
                f'<td{td_s}>{td_t}</td>'
                f'<td{tm_s}>{tm_t}</td>'
                f'<td{pm_s}>{pm_t}</td>'
                f'</tr>'
            )

        avg_row_html = _avg_row("Campaign Average",
                                avg_palt, avg_pasc, avg_pdes, avg_pela,
                                avg_land, avg_tdiff, avg_tmae, avg_pmae)

        n_flights = len({r["launch_id"] for r in rows})
        n_forecasts = len(rows)
        campaign_html += f"""
    <div class="campaign-card">
      <h2>Campaign: {_html.escape(camp_key)}</h2>
      <p style="color:#666;font-size:0.85em">{n_flights} balloon(s), {n_forecasts} evaluation(s)</p>
      {TABLE_HEADER}
        {camp_rows}
        {avg_row_html}
      </tbody></table>
    </div>"""

    # ── Assemble page ─────────────────────────────────────────────────────────
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
    .legend   {{ display:flex; gap:10px; margin-bottom: 18px; font-size:0.85em; }}
    .legend span {{ padding: 3px 10px; border: 1px solid #bbb; border-radius: 3px; }}
    table     {{ border-collapse: collapse; width: 100%; margin-bottom: 6px; font-size: 0.82em; }}
    th, td    {{ border: 1px solid #bbb; padding: 5px 8px; text-align: center; white-space: nowrap; }}
    th        {{ background: #2c3e50; color: #fff; font-weight: bold; }}
    th.sub    {{ background: #34495e; }}
    tr:hover  {{ filter: brightness(0.96); }}
    .campaign-card {{
      background: #fff;
      border: 1px solid #ccc;
      border-radius: 6px;
      padding: 20px 24px;
      margin-bottom: 32px;
      box-shadow: 0 2px 6px rgba(0,0,0,0.08);
    }}
    .overall-section {{
      background: #fff;
      border: 1px solid #ccc;
      border-radius: 6px;
      padding: 20px 24px;
      margin-bottom: 32px;
      box-shadow: 0 2px 6px rgba(0,0,0,0.08);
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
    <span style="background:{GREEN}">&le;&thinsp;10%</span>
    <span style="background:{YELLOW}">10&ndash;25%</span>
    <span style="background:{RED}">&gt;&thinsp;25%</span>
    &nbsp;&nbsp;<b>Error metrics:</b>
    <span style="background:{GREEN}">good</span>
    <span style="background:{YELLOW}">marginal</span>
    <span style="background:{RED}">poor</span>
  </div>

  <div class="overall-section">
    <h2>Overall Summary</h2>
    {TABLE_HEADER}
      {overall_rows}
    </tbody></table>
  </div>

  {campaign_html}
</body>
</html>"""

    path = os.path.join(batch_dir, "summary.html")
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    return path


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_launch(launch: dict) -> list[str]:
    """Return list of validation error strings, empty if the launch is valid."""
    errors = []
    missing = [f for f in REQUIRED_FIELDS if f not in launch]
    if missing:
        errors.append(f"Missing required fields: {missing}")
    if not launch.get("gfs_file") and not launch.get("era5_file"):
        errors.append("Must have at least one of gfs_file or era5_file")
    if launch.get("aprs_file"):
        p = BALLOON_DATA_DIR + launch["aprs_file"]
        if not os.path.exists(p):
            errors.append(f"APRS file not found: {p}")
    if launch.get("gfs_file"):
        p = FORECASTS_DIR + launch["gfs_file"]
        if not os.path.exists(p):
            errors.append(f"GFS forecast not found: {p}")
    if launch.get("era5_file"):
        p = FORECASTS_DIR + launch["era5_file"]
        if not os.path.exists(p):
            errors.append(f"ERA5 forecast not found: {p}")
    return errors

# ── Main batch runner ─────────────────────────────────────────────────────────

def run_batch(note: str):
    if not os.path.exists(LAUNCHES_JSON):
        print(f"ERROR: {LAUNCHES_JSON} not found.")
        print(f"  Copy the example to get started:")
        print(f"    cp evaluation/launches.example.json evaluation/launches.json")
        print(f"  Then add your own launch entries and re-run.")
        sys.exit(1)

    with open(LAUNCHES_JSON) as f:
        launches = json.load(f)["launches"]

    git       = _git_info()
    ts        = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M")
    batch_id  = f"{ts}_{git['git_hash']}"
    batch_dir = BATCHES_DIR + batch_id + "/"
    os.makedirs(batch_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  EarthSHAB Batch Evaluation")
    print(f"  Batch ID : {batch_id}")
    print(f"  Note     : {note}")
    print(f"  Launches : {len(launches)}")
    print(f"{'='*60}\n")

    orig = _snapshot_config()

    summary_rows: list[dict] = []
    attempted:    list[str]  = []
    succeeded:    list[str]  = []
    failed:       dict       = {}
    per_launch_times: list[float] = []

    batch_start = time.monotonic()

    for launch in launches:
        launch_date   = launch.get("launch_time", "")[:10]
        launch_org    = launch.get("organization", "?")
        launch_camp   = launch.get("campaign") or ""
        launch_id     = f"{launch_org}_{launch.get('shab_name','?')}_{launch_date}"
        attempted.append(launch_id)

        print(f"── {launch_id} ──")

        # Validate before running
        errors = _validate_launch(launch)
        if errors:
            msg = "; ".join(errors)
            print(f"  [SKIP] {msg}\n")
            failed[launch_id] = msg
            for ft in ["GFS", "ERA5"]:
                if ft.lower() + "_file" in launch or ft == "GFS" and "gfs_file" in launch:
                    summary_rows.append(_result_to_row(
                        batch_id, launch_id, ft, None, msg,
                        campaign=launch_camp, organization=launch_org,
                        launch_date=launch_date,
                    ))
            continue

        forecast_types = []
        if launch.get("gfs_file"):
            forecast_types.append("GFS")
        if launch.get("era5_file"):
            forecast_types.append("ERA5")

        launch_dir = batch_dir + launch_id + "/"
        os.makedirs(launch_dir, exist_ok=True)

        launch_errors = []
        launch_start  = time.monotonic()

        for ft in forecast_types:
            print(f"  Running {ft}...")
            aprs_fmt = ""
            try:
                overrides = _build_overrides(launch, ft, orig)
                ev = BalloonEvaluator(config_overrides=overrides, output_dir=launch_dir)
                ev.run()
                aprs_fmt = (ev.sim_state or {}).get("aprs_format", "") or ""
                result = ev.compute_metrics()
                ev.report(result)
                ev.plot()
                summary_rows.append(_result_to_row(
                    batch_id, launch_id, ft, result, "ok", aprs_fmt,
                    campaign=launch_camp, organization=launch_org,
                    launch_date=launch_date,
                ))
                print(f"  [{ft}] done\n")
            except Exception:
                tb = traceback.format_exc()
                short_msg = tb.strip().splitlines()[-1]
                print(f"  [{ft}] FAILED: {short_msg}")
                print(f"  Full traceback:\n{tb}")
                summary_rows.append(_result_to_row(
                    batch_id, launch_id, ft, None, short_msg, aprs_fmt,
                    campaign=launch_camp, organization=launch_org,
                    launch_date=launch_date,
                ))
                launch_errors.append(f"{ft}: {short_msg}")

        elapsed = time.monotonic() - launch_start
        per_launch_times.append(elapsed)

        if launch_errors:
            failed[launch_id] = "; ".join(launch_errors)
        else:
            succeeded.append(launch_id)

    total_runtime = time.monotonic() - batch_start
    avg_runtime   = (sum(per_launch_times) / len(per_launch_times)) if per_launch_times else 0.0

    # ── Write summary.csv ────────────────────────────────────────────────────
    summary_path = batch_dir + "summary.csv"
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=_SUMMARY_FIELDNAMES)
        w.writeheader()
        w.writerows(summary_rows)
    print(f"Summary → {summary_path}")

    html_path = _write_summary_html(summary_rows, batch_dir, batch_id, note, git)
    print(f"Report  → {html_path}")

    # ── Write batch_info.json ────────────────────────────────────────────────
    with open(LAUNCHES_JSON) as f:
        launches_snapshot = json.load(f)

    batch_info = {
        "batch_id":                batch_id,
        "note":                    note,
        "timestamp_utc":           datetime.now(timezone.utc).isoformat(),
        "git_hash":                git["git_hash"],
        "git_hash_full":           git["git_hash_full"],
        "git_branch":              git["git_branch"],
        "git_commit_message":      git["git_commit_message"],
        "git_dirty":               git["git_dirty"],
        "earthshab_version":       _earthshab_version(),
        "total_runtime_s":         round(total_runtime, 2),
        "per_launch_avg_runtime_s": round(avg_runtime, 2),
        "launches_attempted":      attempted,
        "launches_succeeded":      succeeded,
        "launches_failed":         failed,
        "launches_json_snapshot":  launches_snapshot,
    }

    info_path = batch_dir + "batch_info.json"
    with open(info_path, "w") as f:
        json.dump(batch_info, f, indent=2, default=str)
    print(f"Batch info → {info_path}")

    # ── Final summary ────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  Batch complete: {len(succeeded)}/{len(attempted)} launches succeeded")
    if failed:
        print(f"  Failed launches:")
        for lid, msg in failed.items():
            print(f"    {lid}: {msg}")
    print(f"  Total runtime: {total_runtime:.1f}s")
    print(f"{'='*60}\n")

    return len(failed) == 0


def main():
    parser = argparse.ArgumentParser(
        description="Run all launches in launches.json against the current codebase."
    )
    parser.add_argument(
        "--note", required=True,
        help="Description of what changed in this batch (e.g. 'tuned vent rate')"
    )
    args = parser.parse_args()
    success = run_batch(args.note)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
