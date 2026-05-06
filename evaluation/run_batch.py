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
    "payload_weight_kg", "balloon_size_m",
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
                   launch_date: str = "",
                   payload_weight_kg: float = math.nan,
                   balloon_size_m: float = math.nan) -> dict:
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
    """Write summary.html with colored % diff table, averages row, and sortable columns."""
    import html as _html

    GREEN  = "#b7e4b7"
    YELLOW = "#ffeaa0"
    RED    = "#f4a5a5"
    FAIL   = "#e8e8e8"

    # Columns: Launch(0) Fcst(1) Payload(2) BalSize(3)
    #          | FA-S(4) FA-T(5) FA-%(6) | Asc-S(7) Asc-T(8) Asc-%(9)
    #          | Des-S(10) Des-T(11) Des-%(12) | Ela-S(13) Ela-T(14) Ela-%(15)
    #          | LandDist(16) |TimeΔ|(17) TempMAE(18) PressMAE(19)

    # ── Value helpers ─────────────────────────────────────────────────────────

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
        color = GREEN if abs(pct) <= lo else YELLOW if abs(pct) <= hi else RED
        sign  = "+" if pct > 0 else ""
        return f"{sign}{pct:.1f}%", f' style="background:{color}"', str(pct)

    def _abs_cell(val, lo, hi):
        if math.isnan(val):
            return "—", "", ""
        color = GREEN if abs(val) <= lo else YELLOW if abs(val) <= hi else RED
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

    # ── Row builders ──────────────────────────────────────────────────────────

    def _td(raw, txt, style=""):
        return f'<td data-val="{_html.escape(str(raw))}"{style}>{txt}</td>'

    def _build_data_row(row: dict, orig_idx: int) -> str:
        status  = row.get("status", "")
        lid_raw = str(row.get("launch_id", ""))
        ft_raw  = str(row.get("forecast_type", ""))
        lid     = _html.escape(lid_raw)
        ft      = _html.escape(ft_raw)

        if status != "ok":
            msg = _html.escape(str(status)[:80])
            return (
                f'<tr class="failed-row" data-orig="{orig_idx}" style="background:{FAIL}">'
                f'<td data-val="{_html.escape(lid_raw)}" style="text-align:left">{lid}</td>'
                f'<td data-val="{_html.escape(ft_raw)}">{ft}</td>'
                f'<td colspan="18" style="text-align:left">FAILED: {msg}</td></tr>'
            )

        def _f(key):
            try:    return float(row.get(key, math.nan))
            except: return math.nan

        pct_alt = _pct(row.get("sim_float_alt_mean_m"),     row.get("truth_float_alt_mean_m"))
        pct_asc = _pct(row.get("sim_ascent_rate_mean_ms"),  row.get("truth_ascent_rate_mean_ms"))
        pct_des = _pct(row.get("sim_descent_rate_mean_ms"), row.get("truth_descent_rate_mean_ms"))
        pct_ela = _pct(row.get("sim_elapsed_time_min"),     row.get("truth_elapsed_time_min"))

        land_km  = _f("landing_distance_km")
        tdiff    = _f("landing_time_diff_min")
        temp_mae = _f("temp_mae_k")
        pres_mae = _f("pressure_mae_pa")

        pa_txt, pa_s, pa_r = _pct_cell(pct_alt)
        pc_txt, pc_s, pc_r = _pct_cell(pct_asc)
        pd_txt, pd_s, pd_r = _pct_cell(pct_des)
        pe_txt, pe_s, pe_r = _pct_cell(pct_ela)
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
            + _td(_rv(row.get("sim_elapsed_time_min")),   _fv(row.get("sim_elapsed_time_min"),   ".0f"))
            + _td(_rv(row.get("truth_elapsed_time_min")), _fv(row.get("truth_elapsed_time_min"), ".0f"))
            + _td(pe_r, pe_txt, pe_s)
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
        avg_se  = a("sim_elapsed_time_min");    avg_te  = a("truth_elapsed_time_min")

        pct_alt = _pct(avg_sa,  avg_ta)
        pct_asc = _pct(avg_sas, avg_tas)
        pct_des = _pct(avg_sd,  avg_td)
        pct_ela = _pct(avg_se,  avg_te)

        avg_land  = a("landing_distance_km")
        avg_tdiff = a("landing_time_diff_min")
        avg_tmae  = a("temp_mae_k")
        avg_pmae  = a("pressure_mae_pa")

        pa_txt, pa_s, _ = _pct_cell(pct_alt)
        pc_txt, pc_s, _ = _pct_cell(pct_asc)
        pd_txt, pd_s, _ = _pct_cell(pct_des)
        pe_txt, pe_s, _ = _pct_cell(pct_ela)
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
            f'<td>{_fv(avg_se,  ".0f")}</td><td>{_fv(avg_te,  ".0f")}</td><td{pe_s}>{pe_txt}</td>'
            f'<td{ld_s}>{ld_txt}</td>'
            f'<td{td_s}>{td_txt}</td>'
            f'<td{tm_s}>{tm_txt}</td>'
            f'<td{pm_s}>{pm_txt}</td>'
            f'</tr>'
        )

    # ── Shared table header (18 columns) ──────────────────────────────────────
    # Sortable <th> elements carry data-col=N and a <span class="si"> for the
    # sort indicator.  Group headers (colspan) are labels only, not sortable.

    TABLE_HEADER = (
        '<thead>'
        '<tr>'
        '<th data-col="0"  rowspan="2" style="text-align:left">Launch<span class="si"></span></th>'
        '<th data-col="1"  rowspan="2">Fcst<span class="si"></span></th>'
        '<th data-col="2"  rowspan="2">Payload<br>(kg)<span class="si"></span></th>'
        '<th data-col="3"  rowspan="2">Bal&nbsp;&Oslash;<br>(m)<span class="si"></span></th>'
        '<th colspan="3">Float Alt (m)</th>'
        '<th colspan="3">Ascent (m/s)</th>'
        '<th colspan="3">Descent (m/s)</th>'
        '<th colspan="3">Elapsed (min)</th>'
        '<th data-col="16" rowspan="2">Land Dist<br>(km)<span class="si"></span></th>'
        '<th data-col="17" rowspan="2">|Time &Delta;|<br>(min)<span class="si"></span></th>'
        '<th data-col="18" rowspan="2">Temp MAE<br>(K)<span class="si"></span></th>'
        '<th data-col="19" rowspan="2">Press MAE<br>(Pa)<span class="si"></span></th>'
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
        '<th data-col="15" class="sub">%&Delta;<span class="si"></span></th>'
        '</tr>'
        '</thead>'
    )

    # ── Overall table ─────────────────────────────────────────────────────────
    overall_data = "".join(_build_data_row(r, i) for i, r in enumerate(summary_rows))
    overall_avg  = _build_avg_row("Average", summary_rows)

    # ── Campaign sections ─────────────────────────────────────────────────────
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
    tbl_n = [1]
    for camp_key, rows in sorted(groups.items()):
        if len({r["launch_id"] for r in rows}) < 2:
            continue
        tbl_id      = f"tbl-camp-{tbl_n[0]}"
        tbl_n[0]   += 1
        camp_data   = "".join(_build_data_row(r, i) for i, r in enumerate(rows))
        camp_avg    = _build_avg_row("Campaign Average", rows)
        n_flights   = len({r["launch_id"] for r in rows})
        n_forecasts = len(rows)
        campaign_html += (
            f'<div class="campaign-card">'
            f'<h2>Campaign: {_html.escape(camp_key)}</h2>'
            f'<p style="color:#666;font-size:0.85em">'
            f'{n_flights} balloon(s), {n_forecasts} evaluation(s)</p>'
            f'<table id="{tbl_id}">{TABLE_HEADER}'
            f'<tbody>{camp_data}</tbody>'
            f'<tfoot>{camp_avg}</tfoot>'
            f'</table></div>'
        )

    # ── Sort JavaScript ───────────────────────────────────────────────────────
    sort_js = """\
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
    <span style="background:{GREEN}">&le;&thinsp;10%</span>
    <span style="background:{YELLOW}">10&ndash;25%</span>
    <span style="background:{RED}">&gt;&thinsp;25%</span>
    &nbsp;&nbsp;<b>Error metrics:</b>
    <span style="background:{GREEN}">good</span>
    <span style="background:{YELLOW}">marginal</span>
    <span style="background:{RED}">poor</span>
    &nbsp;&nbsp;<span style="color:#666">Click any column header to sort &uarr;&darr;</span>
  </div>

  <div class="overall-section">
    <h2>Overall Summary</h2>
    <table id="tbl-overall">
      {TABLE_HEADER}
      <tbody>{overall_data}</tbody>
      <tfoot>{overall_avg}</tfoot>
    </table>
  </div>

  {campaign_html}
  {sort_js}
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
        launch_date    = launch.get("launch_time", "")[:10]
        launch_org     = launch.get("organization", "?")
        launch_camp    = launch.get("campaign") or ""
        launch_id      = f"{launch_org}_{launch.get('shab_name','?')}_{launch_date}"
        launch_payload = float(launch.get("payload_weight_kg") or math.nan)
        launch_balsize = float(launch.get("balloon_size") or math.nan)
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
                        payload_weight_kg=launch_payload,
                        balloon_size_m=launch_balsize,
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
                    payload_weight_kg=launch_payload,
                    balloon_size_m=launch_balsize,
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
                    payload_weight_kg=launch_payload,
                    balloon_size_m=launch_balsize,
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
