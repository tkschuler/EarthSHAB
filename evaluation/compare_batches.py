"""compare_batches.py - Compare two batches: A (baseline) vs B (experiment).

Answers "did my code change make trajectory sims better, worse, or different."
Produces an HTML report with per-launch deltas, aggregates broken out by
forecast type / campaign, win/loss counts, and overview bar+scatter plots.

Usage:
    PYTHONPATH=src python -m evaluation.compare_batches <batch_A> <batch_B>

Both arguments accept either a full batch_id (e.g. "2026-05-06T1937_57c5ff9")
or just a short git hash (e.g. "57c5ff9"); ambiguous hashes error out.

Output structure:
    evaluation/comparisons/<timestamp>_<hashA>_vs_<hashB>/
        compare.html
        plot_<metric>_<forecast>_bar.png
        plot_<metric>_<forecast>_scatter.png
        ...

Sign convention: delta = B - A. For all metrics (which are magnitudes of
error), negative delta = improvement = green.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Optional

import pandas as pd


# ── Paths ─────────────────────────────────────────────────────────────────────

EVAL_DIR        = "evaluation/"
BATCHES_DIR     = EVAL_DIR + "batches/"
COMPARISONS_DIR = EVAL_DIR + "comparisons/"


# ── Thresholds ────────────────────────────────────────────────────────────────

# Cells are colored only if BOTH |delta| > abs_floor (per-metric) AND |%| > PCT_FLOOR.
PCT_FLOOR = 5.0


# ── Metric definitions ────────────────────────────────────────────────────────
#
# All metrics here represent "magnitude of error." For signed diff columns
# (landing_time_diff_min, time_to_*_diff_min) we take abs() before comparing so
# that "got closer to zero" reads as improvement uniformly.
#
# Derived metrics (`derived` key set) are computed in-script as |sim - truth|
# from the named pair of columns.

METRICS = [
    {"key": "landing_distance_km",     "label": "Landing dist (km)",       "abs_floor": 1.0,    "fmt": ".1f"},
    {"key": "landing_time_diff_min",   "label": "Landing time err (min)",  "abs_floor": 1.0,    "fmt": ".1f", "use_abs": True},
    {"key": "time_to_float_diff_min",  "label": "Time-to-float err (min)", "abs_floor": 1.0,    "fmt": ".1f", "use_abs": True},
    {"key": "time_to_ground_diff_min", "label": "Time-to-ground err (min)","abs_floor": 1.0,    "fmt": ".1f", "use_abs": True},
    {"key": "temp_mae_k",              "label": "Temp MAE (K)",            "abs_floor": 0.5,    "fmt": ".2f"},
    {"key": "pressure_mae_pa",         "label": "Pressure MAE (Pa)",       "abs_floor": 100.0,  "fmt": ".0f"},
    {"key": "float_alt_abs_err_m",     "label": "Float alt err (m)",       "abs_floor": 50.0,   "fmt": ".0f",
     "derived": ("sim_float_alt_mean_m", "truth_float_alt_mean_m")},
    {"key": "ascent_rate_abs_err_ms",  "label": "Ascent rate err (m/s)",   "abs_floor": 0.1,    "fmt": ".2f",
     "derived": ("sim_ascent_rate_mean_ms", "truth_ascent_rate_mean_ms")},
    {"key": "descent_rate_abs_err_ms", "label": "Descent rate err (m/s)",  "abs_floor": 0.1,    "fmt": ".2f",
     "derived": ("sim_descent_rate_mean_ms", "truth_descent_rate_mean_ms")},
]

REFORECAST_METRIC = {
    "key": "reforecast_landing_dist_m", "label": "Reforecast landing dist (m)",
    "abs_floor": 1000.0, "fmt": ".0f",
}

# Headlines shown in the stdout/HTML Summary block (in order).
PRIMARY_METRIC_KEYS = [
    "landing_distance_km",
    "time_to_ground_diff_min",
    "float_alt_abs_err_m",
]


# ── Pure logic (unit-tested in tests/test_compare_batches.py) ─────────────────


def resolve_batch_id(arg: str, batches_dir: str = BATCHES_DIR) -> str:
    """Resolve a CLI arg (full batch_id OR short git hash) to a full batch_id.

    Raises FileNotFoundError if nothing matches, ValueError if a short hash
    matches more than one folder.
    """
    if not os.path.isdir(batches_dir):
        raise FileNotFoundError(f"Batches directory not found: {batches_dir}")
    candidates = sorted(
        d for d in os.listdir(batches_dir)
        if os.path.isdir(os.path.join(batches_dir, d))
    )
    # Exact match wins.
    if arg in candidates:
        return arg
    # Try as a short hash (matches folders ending with "_" + arg).
    suffix = "_" + arg
    matches = [d for d in candidates if d.endswith(suffix)]
    if not matches:
        raise FileNotFoundError(
            f"No batch matches '{arg}'. Available:\n  " + "\n  ".join(candidates)
        )
    if len(matches) > 1:
        raise ValueError(
            f"Hash '{arg}' is ambiguous. Matches:\n  " + "\n  ".join(matches)
        )
    return matches[0]


def add_derived_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Add |sim - truth| columns for derived metrics. Returns df modified in place."""
    for m in METRICS:
        derived = m.get("derived")
        if not derived:
            continue
        sim_col, truth_col = derived
        if sim_col in df.columns and truth_col in df.columns:
            df[m["key"]] = (
                pd.to_numeric(df[sim_col], errors="coerce")
                - pd.to_numeric(df[truth_col], errors="coerce")
            ).abs()
        else:
            df[m["key"]] = math.nan
    # Apply abs() to use_abs metrics so all values represent error magnitude.
    for m in METRICS:
        if m.get("use_abs") and m["key"] in df.columns:
            df[m["key"]] = pd.to_numeric(df[m["key"]], errors="coerce").abs()
    return df


def delta_and_pct(a: float, b: float, abs_floor: float) -> tuple[float, float]:
    """Return (delta, pct). pct is NaN if |a| < abs_floor (too small to be
    meaningful as a denominator)."""
    if a is None or b is None:
        return math.nan, math.nan
    try:
        av, bv = float(a), float(b)
    except (TypeError, ValueError):
        return math.nan, math.nan
    if math.isnan(av) or math.isnan(bv):
        return math.nan, math.nan
    delta = bv - av
    if abs(av) < abs_floor:
        return delta, math.nan
    pct = delta / abs(av) * 100.0
    return delta, pct


def classify_change(delta: float, pct: float, abs_floor: float) -> str:
    """Classify a delta as 'improved' | 'worsened' | 'unchanged'.

    A change counts only if both:
      - |delta| > abs_floor (passes absolute noise floor), AND
      - |pct|   > PCT_FLOOR (passes relative noise floor) — OR pct is NaN
                 (baseline too small to compute %), in which case abs alone decides.
    """
    if delta is None or (isinstance(delta, float) and math.isnan(delta)):
        return "unchanged"
    if abs(delta) <= abs_floor:
        return "unchanged"
    if not (isinstance(pct, float) and math.isnan(pct)):
        if abs(pct) <= PCT_FLOOR:
            return "unchanged"
    return "improved" if delta < 0 else "worsened"


def win_loss_counts(deltas_pcts: list[tuple[float, float]], abs_floor: float) -> dict:
    """Count improved / unchanged / worsened across a list of (delta, pct) pairs."""
    counts = {"improved": 0, "unchanged": 0, "worsened": 0}
    for delta, pct in deltas_pcts:
        counts[classify_change(delta, pct, abs_floor)] += 1
    return counts


def nan_mean(values) -> float:
    """Mean ignoring NaN; NaN if no valid values."""
    clean = []
    for v in values:
        if v is None:
            continue
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        if not math.isnan(fv):
            clean.append(fv)
    return sum(clean) / len(clean) if clean else math.nan


def find_asymmetries(df_a: pd.DataFrame, df_b: pd.DataFrame) -> dict:
    """Classify (launch_id, forecast_type) keys across the two batches.

    Returns a dict with sets:
        intersection_ok    — present + ok in both (the diff set)
        only_in_a          — present in A only (any status)
        only_in_b          — present in B only
        failed_in_a        — present in both, ok in B, not-ok in A
        failed_in_b        — present in both, ok in A, not-ok in B
        failed_in_both     — present in both, not-ok in both
    """
    def keymap(df):
        if df is None or df.empty:
            return {}
        return {
            (str(r["launch_id"]), str(r["forecast_type"])): str(r.get("status", ""))
            for _, r in df.iterrows()
        }
    a_map, b_map = keymap(df_a), keymap(df_b)
    a_keys, b_keys = set(a_map), set(b_map)

    intersection = a_keys & b_keys
    out = {
        "intersection_ok":  set(),
        "only_in_a":        a_keys - b_keys,
        "only_in_b":        b_keys - a_keys,
        "failed_in_a":      set(),
        "failed_in_b":      set(),
        "failed_in_both":   set(),
    }
    for k in intersection:
        a_ok = a_map[k] == "ok"
        b_ok = b_map[k] == "ok"
        if a_ok and b_ok:
            out["intersection_ok"].add(k)
        elif a_ok and not b_ok:
            out["failed_in_b"].add(k)
        elif b_ok and not a_ok:
            out["failed_in_a"].add(k)
        else:
            out["failed_in_both"].add(k)
    return out


def find_missing_metric_columns(
    df_a: pd.DataFrame, df_b: pd.DataFrame, metric_keys: list[str]
) -> dict:
    """Return {metric_key: 'A' | 'B' | 'both'} for metrics absent from a CSV.

    Derived metrics are excluded — they're computed by add_derived_metrics()
    from sim_*/truth_* columns, so absence is handled at that layer.
    """
    out = {}
    a_cols, b_cols = set(df_a.columns), set(df_b.columns)
    for k in metric_keys:
        miss_a = k not in a_cols
        miss_b = k not in b_cols
        if miss_a and miss_b:
            out[k] = "both"
        elif miss_a:
            out[k] = "A"
        elif miss_b:
            out[k] = "B"
    return out


# ── IO ────────────────────────────────────────────────────────────────────────


@dataclass
class BatchData:
    batch_id: str
    info: dict          # batch_info.json contents
    df: pd.DataFrame    # summary.csv contents


def load_batch(batch_id: str, batches_dir: str = BATCHES_DIR) -> BatchData:
    folder = os.path.join(batches_dir, batch_id)
    summary_csv = os.path.join(folder, "summary.csv")
    info_json   = os.path.join(folder, "batch_info.json")
    if not os.path.isfile(summary_csv):
        raise FileNotFoundError(f"Missing summary.csv in {folder}")
    df = pd.read_csv(summary_csv)
    info = {}
    if os.path.isfile(info_json):
        with open(info_json) as f:
            info = json.load(f)
    return BatchData(batch_id=batch_id, info=info, df=df)


# ── Comparison computation ────────────────────────────────────────────────────


def compute_per_launch_rows(
    df_a: pd.DataFrame, df_b: pd.DataFrame,
    intersection: set, available_metrics: list[dict]
) -> list[dict]:
    """Build one row per (launch_id, forecast_type) in the intersection.

    Each row contains the launch_id, forecast_type, campaign, organization,
    and per-metric dicts {a, b, delta, pct, color_class}.
    """
    a_idx = {(str(r["launch_id"]), str(r["forecast_type"])): r for _, r in df_a.iterrows()}
    b_idx = {(str(r["launch_id"]), str(r["forecast_type"])): r for _, r in df_b.iterrows()}

    rows = []
    for key in sorted(intersection):
        lid, ft = key
        ra, rb = a_idx[key], b_idx[key]
        row = {
            "launch_id":      lid,
            "forecast_type":  ft,
            "campaign":       str(ra.get("campaign", "") or rb.get("campaign", "")),
            "organization":   str(ra.get("organization", "") or rb.get("organization", "")),
            "launch_date":    str(ra.get("launch_date", "") or rb.get("launch_date", "")),
            "metrics":        {},
        }
        for m in available_metrics:
            k = m["key"]
            av = ra.get(k, math.nan) if k in ra.index else math.nan
            bv = rb.get(k, math.nan) if k in rb.index else math.nan
            delta, pct = delta_and_pct(av, bv, m["abs_floor"])
            row["metrics"][k] = {
                "a":     _to_float(av),
                "b":     _to_float(bv),
                "delta": delta,
                "pct":   pct,
                "class": classify_change(delta, pct, m["abs_floor"]),
            }
        rows.append(row)
    return rows


def _to_float(v) -> float:
    try:
        fv = float(v)
        return fv
    except (TypeError, ValueError):
        return math.nan


def compute_aggregate(rows: list[dict], available_metrics: list[dict],
                      group_filter=None, label: str = "Overall") -> dict:
    """One aggregate row: mean(A), mean(B), delta, pct, plus win/loss counts.

    group_filter is a callable(row) -> bool, or None for "all rows".
    """
    subset = [r for r in rows if (group_filter is None or group_filter(r))]
    agg = {"label": label, "n": len(subset), "metrics": {}}
    for m in available_metrics:
        k = m["key"]
        a_mean = nan_mean(r["metrics"][k]["a"] for r in subset)
        b_mean = nan_mean(r["metrics"][k]["b"] for r in subset)
        delta, pct = delta_and_pct(a_mean, b_mean, m["abs_floor"])
        cls = classify_change(delta, pct, m["abs_floor"])
        deltas_pcts = [
            (r["metrics"][k]["delta"], r["metrics"][k]["pct"])
            for r in subset
            if not math.isnan(r["metrics"][k].get("delta", math.nan))
        ]
        counts = win_loss_counts(deltas_pcts, m["abs_floor"])
        agg["metrics"][k] = {
            "a_mean": a_mean, "b_mean": b_mean,
            "delta":  delta,  "pct":    pct,
            "class":  cls,
            "counts": counts,
            "n":      len(deltas_pcts),
        }
    return agg


def compute_aggregates(rows: list[dict], available_metrics: list[dict]) -> list[dict]:
    """Build the standard set of aggregates: overall, per forecast_type, per campaign."""
    out = [compute_aggregate(rows, available_metrics, None, "Overall")]
    for ft in sorted({r["forecast_type"] for r in rows}):
        out.append(compute_aggregate(
            rows, available_metrics,
            (lambda ft_=ft: lambda r: r["forecast_type"] == ft_)(),
            f"Forecast = {ft}",
        ))
    for camp in sorted({r["campaign"] for r in rows if r["campaign"]}):
        out.append(compute_aggregate(
            rows, available_metrics,
            (lambda c_=camp: lambda r: r["campaign"] == c_)(),
            f"Campaign = {camp}",
        ))
    return out


# ── Summary ─────────────────────────────────────────────────────────────────────


def build_summary_lines(rows: list[dict], available_metrics: list[dict]) -> list[str]:
    """Build short headline strings: one per primary metric per forecast type."""
    keys_by = {m["key"]: m for m in available_metrics}
    forecasts = sorted({r["forecast_type"] for r in rows})
    lines = []
    for mkey in PRIMARY_METRIC_KEYS:
        if mkey not in keys_by:
            continue
        m = keys_by[mkey]
        for ft in forecasts:
            subset = [r for r in rows if r["forecast_type"] == ft]
            if not subset:
                continue
            a_mean = nan_mean(r["metrics"][mkey]["a"] for r in subset)
            b_mean = nan_mean(r["metrics"][mkey]["b"] for r in subset)
            delta, pct = delta_and_pct(a_mean, b_mean, m["abs_floor"])
            deltas_pcts = [
                (r["metrics"][mkey]["delta"], r["metrics"][mkey]["pct"])
                for r in subset
                if not math.isnan(r["metrics"][mkey].get("delta", math.nan))
            ]
            counts = win_loss_counts(deltas_pcts, m["abs_floor"])
            n = counts["improved"] + counts["unchanged"] + counts["worsened"]
            pct_str = f"{pct:+.1f}%" if not math.isnan(pct) else "n/a"
            a_str = f"{a_mean:{m['fmt']}}" if not math.isnan(a_mean) else "—"
            b_str = f"{b_mean:{m['fmt']}}" if not math.isnan(b_mean) else "—"
            lines.append(
                f"  {m['label']:32s} [{ft:4s}]  "
                f"{a_str} → {b_str}  ({pct_str}, {counts['improved']}/{n} improved)"
            )
    return lines


# ── CLI ───────────────────────────────────────────────────────────────────────


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Compare two batches (A = baseline, B = experiment).")
    ap.add_argument("batch_a", help="Baseline batch_id or short git hash")
    ap.add_argument("batch_b", help="Experiment batch_id or short git hash")
    args = ap.parse_args(argv)

    # Resolve and load.
    try:
        a_id = resolve_batch_id(args.batch_a)
        b_id = resolve_batch_id(args.batch_b)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    print(f"  A (baseline)   : {a_id}")
    print(f"  B (experiment) : {b_id}")

    batch_a = load_batch(a_id)
    batch_b = load_batch(b_id)

    # Add derived |sim - truth| columns and abs() the use_abs metrics.
    add_derived_metrics(batch_a.df)
    add_derived_metrics(batch_b.df)

    # Detect schema mismatch on the non-derived metric columns.
    non_derived_keys = [m["key"] for m in METRICS if "derived" not in m] + [REFORECAST_METRIC["key"]]
    missing = find_missing_metric_columns(batch_a.df, batch_b.df, non_derived_keys)

    # Build the available-metrics list (drop ones missing in either side).
    available_metrics = [m for m in METRICS if m["key"] not in missing]
    available_reforecast = REFORECAST_METRIC if REFORECAST_METRIC["key"] not in missing else None

    # Asymmetry classification.
    asym = find_asymmetries(batch_a.df, batch_b.df)
    print(f"  Intersection (ok in both) : {len(asym['intersection_ok'])} (launch, forecast) pairs")
    if asym["only_in_a"] or asym["only_in_b"] or asym["failed_in_a"] or asym["failed_in_b"]:
        print(f"  Asymmetries               : "
              f"only_in_A={len(asym['only_in_a'])} "
              f"only_in_B={len(asym['only_in_b'])} "
              f"failed_in_A={len(asym['failed_in_a'])} "
              f"failed_in_B={len(asym['failed_in_b'])}")

    # Per-launch rows + aggregates.
    per_launch = compute_per_launch_rows(batch_a.df, batch_b.df,
                                         asym["intersection_ok"], available_metrics)
    aggregates = compute_aggregates(per_launch, available_metrics)

    # Reforecast — computed for every (launch, forecast_type) pair that has data.
    reforecast_section = None
    if available_reforecast:
        ref_rows = compute_per_launch_rows(batch_a.df, batch_b.df,
                                           asym["intersection_ok"], [available_reforecast])
        ref_aggs = compute_aggregates(ref_rows, [available_reforecast])
        reforecast_section = {"rows": ref_rows, "aggregates": ref_aggs, "metric": available_reforecast}

    # Output directory.
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M")
    short_a = _git_hash_from_id(a_id)
    short_b = _git_hash_from_id(b_id)
    out_dir = os.path.join(COMPARISONS_DIR, f"{ts}_{short_a}_vs_{short_b}")
    os.makedirs(out_dir, exist_ok=True)

    # Plots.
    from evaluation.compare_plots import generate_all_plots
    plot_files = generate_all_plots(
        per_launch, available_metrics,
        reforecast_section, out_dir,
    )

    # Summary.
    summary_lines = build_summary_lines(per_launch, available_metrics)
    print()
    print("Summary:")
    for line in summary_lines:
        print(line)
    print()

    # HTML.
    from evaluation.compare_reporting import write_compare_html
    html_path = write_compare_html(
        out_dir=out_dir,
        batch_a=batch_a, batch_b=batch_b,
        per_launch=per_launch, aggregates=aggregates,
        reforecast_section=reforecast_section,
        asymmetries=asym, missing=missing,
        available_metrics=available_metrics,
        summary_lines=summary_lines, plot_files=plot_files,
    )
    print(f"  Wrote: {html_path}")
    return 0


def _git_hash_from_id(batch_id: str) -> str:
    """Extract the trailing short git hash from a batch_id like '2026-05-08T2002_4afc699'."""
    if "_" in batch_id:
        return batch_id.rsplit("_", 1)[-1]
    return batch_id


if __name__ == "__main__":
    sys.exit(main())
