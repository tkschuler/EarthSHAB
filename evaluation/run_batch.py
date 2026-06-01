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
import gc
import json
import math
import os
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone

import matplotlib.pyplot as plt
import pandas as pd

import EarthSHAB.config_earth as config_earth
from evaluation.evaluate import BalloonEvaluator
from evaluation.reporting import SUMMARY_FIELDNAMES, result_to_summary_row, write_summary_html

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
    }


def _parse_gfs_filename(filename: str) -> tuple[str, str]:
    """Return (forecast_start_str, hour_str) from a GFS filename.

    Handles both gfs_0p25_YYYYMMDD_HH.nc and gfs_0p25_YYYYMMDD_HH_<suffix>.nc.
    Uses basename so path components with underscores (e.g. "balloon_data/")
    don't corrupt the split.
    """
    stem  = os.path.basename(filename).replace(".nc", "")
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
        gfs_file = launch["gfs_file"]
        gfs_path = FORECASTS_DIR + gfs_file
        forecast_start, hour_str = _parse_gfs_filename(gfs_file)
        forecast_start_dt = datetime.fromisoformat(forecast_start)

        overrides["forecast"] = {
            "file":                gfs_path,
            "forecast_start_time": forecast_start,
            "GFSrate":             orig["forecast"]["GFSrate"],
            "wind_interpolation":  orig["forecast"].get("wind_interpolation", "linear_full"),
        }
        # netcdf_gfs is still used by saveNETCDF.py for downloads; some
        # downstream sim_state fields (nc_start, hourstamp) still read it.
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
        era5_path = FORECASTS_DIR + launch["era5_file"]
        overrides["forecast"] = {
            "file":                era5_path,
            "forecast_start_time": orig["forecast"]["forecast_start_time"],
            "GFSrate":             orig["forecast"]["GFSrate"],
            "wind_interpolation":  orig["forecast"].get("wind_interpolation", "linear_full"),
        }
        # Keep netcdf_gfs in overrides so simulate.py's read of
        # config_earth.netcdf_gfs['nc_start']/['hourstamp'] doesn't trip.
        overrides["netcdf_gfs"] = copy.deepcopy(orig["netcdf_gfs"])

    return overrides

# ── Summary CSV / HTML ────────────────────────────────────────────────────────
# Per-row, per-table, and per-page formatting all live in evaluation.reporting.


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
        launch_type    = launch.get("launch_type", "standard")
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
                    summary_rows.append(result_to_summary_row(
                        batch_id, launch_id, ft, None, msg,
                        campaign=launch_camp, organization=launch_org,
                        launch_date=launch_date,
                        payload_weight_kg=launch_payload,
                        balloon_size_m=launch_balsize,
                        launch_type=launch_type,
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
                # launch_type is optional — defaults to "standard" if absent
                ev = BalloonEvaluator(
                    config_overrides=overrides,
                    output_dir=launch_dir,
                    launch_type=launch_type,
                )
                ev.run()
                aprs_fmt = (ev.sim_state or {}).get("aprs_format", "") or ""
                result = ev.compute_metrics()
                ev.report(result)
                ev.plot()
                summary_rows.append(result_to_summary_row(
                    batch_id, launch_id, ft, result, "ok", aprs_fmt,
                    campaign=launch_camp, organization=launch_org,
                    launch_date=launch_date,
                    payload_weight_kg=launch_payload,
                    balloon_size_m=launch_balsize,
                    launch_type=launch_type,
                ))
                print(f"  [{ft}] done\n")
            except Exception:
                tb = traceback.format_exc()
                short_msg = tb.strip().splitlines()[-1]
                print(f"  [{ft}] FAILED: {short_msg}")
                print(f"  Full traceback:\n{tb}")
                summary_rows.append(result_to_summary_row(
                    batch_id, launch_id, ft, None, short_msg, aprs_fmt,
                    campaign=launch_camp, organization=launch_org,
                    launch_date=launch_date,
                    payload_weight_kg=launch_payload,
                    balloon_size_m=launch_balsize,
                    launch_type=launch_type,
                ))
                launch_errors.append(f"{ft}: {short_msg}")

        elapsed = time.monotonic() - launch_start
        per_launch_times.append(elapsed)

        # Release accumulated matplotlib figures (evaluate.py and windmap.py
        # don't close their figs explicitly). Without this, ~6-7 launches in
        # the process segfaults from accumulated figure + netCDF handle state.
        plt.close('all')
        # Force GC so lingering GFS/ERA5 reader objects (each holds an open
        # netCDF4.Dataset) are reclaimed; netCDF4 closes the file in __del__.
        # Without this, HDF5 chunk-cache memory from raw CDS files accumulates
        # across launches and the process segfaults on the 6th-7th launch.
        gc.collect()

        if launch_errors:
            failed[launch_id] = "; ".join(launch_errors)
        else:
            succeeded.append(launch_id)

    total_runtime = time.monotonic() - batch_start
    avg_runtime   = (sum(per_launch_times) / len(per_launch_times)) if per_launch_times else 0.0

    # ── Write summary.csv ────────────────────────────────────────────────────
    summary_path = batch_dir + "summary.csv"
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_FIELDNAMES)
        w.writeheader()
        w.writerows(summary_rows)
    print(f"Summary → {summary_path}")

    html_path = write_summary_html(summary_rows, batch_dir, batch_id, note, git)
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
        "wind_interpolation":      config_earth.forecast.get(
                                       "wind_interpolation", "linear_neighbors"),
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
