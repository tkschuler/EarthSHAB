"""compare_batches.py - Compare summary metrics between two batch runs.

Usage:
    python -m evaluation.compare_batches <batch_id_a> <batch_id_b>

    batch_id_a is the baseline; batch_id_b is the candidate.
    Positive delta means batch_b improved (lower error) for error metrics.

Example:
    python -m evaluation.compare_batches \\
        2026-04-28T1423_a3f9c12 \\
        2026-04-28T1601_b7d4e21
"""

import argparse
import json
import math
import os
import sys

import pandas as pd

EVAL_DIR    = "evaluation/"
BATCHES_DIR = EVAL_DIR + "batches/"

# Metrics where a smaller absolute value is better (errors).
# Delta = baseline - candidate: positive means candidate improved.
ERROR_METRICS = [
    "landing_distance_km",
    "landing_time_diff_min",
    "temp_mae_k",
    "pressure_mae_pa",
    "reforecast_landing_dist_m",
]

# Metrics shown as raw values for both batches (no sign convention).
INFO_METRICS = [
    "sim_float_alt_mean_m",
    "sim_ascent_rate_mean_ms",
    "sim_descent_rate_mean_ms",
    "sim_elapsed_time_min",
]


def _load_batch(batch_id: str) -> tuple[pd.DataFrame, dict]:
    batch_dir = BATCHES_DIR + batch_id + "/"
    summary   = batch_dir + "summary.csv"
    info_file = batch_dir + "batch_info.json"

    if not os.path.exists(summary):
        print(f"ERROR: summary.csv not found for batch '{batch_id}'")
        print(f"       Expected: {summary}")
        sys.exit(1)

    df = pd.read_csv(summary)
    info = {}
    if os.path.exists(info_file):
        with open(info_file) as f:
            info = json.load(f)
    return df, info


def _fmt(v, width=12) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "N/A".rjust(width)
    if isinstance(v, float):
        return f"{v:+.3f}".rjust(width) if abs(v) < 1e6 else f"{v:+.0f}".rjust(width)
    return str(v).rjust(width)


def compare(batch_id_a: str, batch_id_b: str):
    df_a, info_a = _load_batch(batch_id_a)
    df_b, info_b = _load_batch(batch_id_b)

    print(f"\n{'='*72}")
    print(f"  Batch comparison")
    print(f"  Baseline  (A): {batch_id_a}")
    print(f"             Note: {info_a.get('note', 'N/A')}")
    print(f"             Commit: {info_a.get('git_commit_message', 'N/A')}")
    print(f"  Candidate (B): {batch_id_b}")
    print(f"             Note: {info_b.get('note', 'N/A')}")
    print(f"             Commit: {info_b.get('git_commit_message', 'N/A')}")
    print(f"{'='*72}")

    merge_keys = ["launch_id", "forecast_type"]
    merged = df_a.merge(df_b, on=merge_keys, suffixes=("_a", "_b"), how="outer")

    if merged.empty:
        print("  No matching launch_id + forecast_type pairs found between batches.")
        return

    for _, row in merged.iterrows():
        launch_id    = row["launch_id"]
        forecast_type = row["forecast_type"]
        status_a     = row.get("status_a", "N/A")
        status_b     = row.get("status_b", "N/A")

        print(f"\n  ── {launch_id} [{forecast_type}] ──")
        print(f"     Status A: {status_a}")
        print(f"     Status B: {status_b}")

        if status_a != "ok" or status_b != "ok":
            print("     (skipping metric diff — at least one batch failed this launch)")
            continue

        print(f"\n  {'Metric':<35} {'A':>12} {'B':>12} {'Delta (A-B)':>14}")
        print(f"  {'-'*35} {'-'*12} {'-'*12} {'-'*14}")

        for col in ERROR_METRICS + INFO_METRICS:
            col_a = col + "_a"
            col_b = col + "_b"
            if col_a not in row or col_b not in row:
                continue
            va = row[col_a]
            vb = row[col_b]
            try:
                va_f = float(va)
                vb_f = float(vb)
                delta = va_f - vb_f
            except (ValueError, TypeError):
                va_f = vb_f = delta = float("nan")

            label = col.replace("_", " ")
            print(f"  {label:<35} {_fmt(va_f)} {_fmt(vb_f)} {_fmt(delta)}")

    # Launches only in one batch
    only_a = set(df_a["launch_id"]) - set(df_b["launch_id"])
    only_b = set(df_b["launch_id"]) - set(df_a["launch_id"])
    if only_a:
        print(f"\n  Launches only in A: {sorted(only_a)}")
    if only_b:
        print(f"  Launches only in B: {sorted(only_b)}")

    print(f"\n{'='*72}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Compare summary metrics between two batch runs."
    )
    parser.add_argument("batch_a", help="Baseline batch ID")
    parser.add_argument("batch_b", help="Candidate batch ID")
    args = parser.parse_args()
    compare(args.batch_a, args.batch_b)


if __name__ == "__main__":
    main()
