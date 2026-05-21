"""compare_plots.py - Overview plots (bar + scatter) for batch comparison.

For each metric x {GFS, ERA5} we generate:
  - bar chart of per-launch signed deltas, sorted by signed delta
  - scatter of A vs B with y=x reference line

The GFS and ERA5 bar charts for the same metric share Y-axis bounds so the
side-by-side HTML layout makes magnitudes visually comparable. Scatter plots
share both axes for the same reason.
"""

from __future__ import annotations

import math
import os
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


_GREEN = "#4caf50"
_RED   = "#e57373"
_GRAY  = "#bdbdbd"


def _bar_color(delta: float, abs_floor: float, pct: float, pct_floor: float = 5.0) -> str:
    if delta is None or (isinstance(delta, float) and math.isnan(delta)):
        return _GRAY
    if abs(delta) <= abs_floor:
        return _GRAY
    if not (isinstance(pct, float) and math.isnan(pct)) and abs(pct) <= pct_floor:
        return _GRAY
    return _GREEN if delta < 0 else _RED


def _safe_max(values, default=1.0):
    clean = [abs(v) for v in values if v is not None and not (isinstance(v, float) and math.isnan(v))]
    return max(clean) if clean else default


def _bar_chart(rows, metric: dict, out_path: str,
               y_abs_limit: float, title_extra: str = ""):
    """One bar chart of per-launch signed deltas, sorted by signed delta."""
    items = []
    for r in rows:
        m = r["metrics"][metric["key"]]
        if m["delta"] is None or (isinstance(m["delta"], float) and math.isnan(m["delta"])):
            continue
        items.append((r["launch_id"], m["delta"], m["pct"]))
    items.sort(key=lambda x: x[1])  # signed

    fig, ax = plt.subplots(figsize=(max(8, len(items) * 0.22), 4.5))
    if not items:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
    else:
        labels   = [it[0] for it in items]
        deltas   = [it[1] for it in items]
        colors   = [_bar_color(it[1], metric["abs_floor"], it[2]) for it in items]
        x        = range(len(items))
        ax.bar(x, deltas, color=colors, edgecolor="#444", linewidth=0.4)
        ax.axhline(0, color="#000", linewidth=0.6)
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, rotation=75, ha="right", fontsize=7)
        ax.set_ylabel(f"Δ {metric['label']}  (B − A)")
        ax.set_title(f"{metric['label']} — per-launch delta {title_extra}".strip())
        # Shared y-limits across GFS/ERA5 for this metric (passed in).
        if y_abs_limit > 0:
            ax.set_ylim(-y_abs_limit * 1.05, y_abs_limit * 1.05)
        ax.grid(True, axis="y", linestyle=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def _scatter_chart(rows, metric: dict, out_path: str,
                   axis_limit: float, title_extra: str = ""):
    """Scatter of A vs B with y=x reference. Points below line = improved."""
    xs, ys = [], []
    labels = []
    for r in rows:
        m = r["metrics"][metric["key"]]
        a, b = m["a"], m["b"]
        if any(isinstance(v, float) and math.isnan(v) for v in (a, b)):
            continue
        xs.append(a); ys.append(b); labels.append(r["launch_id"])

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    if not xs:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
    else:
        colors = [_GREEN if b < a else _RED if b > a else _GRAY for a, b in zip(xs, ys)]
        ax.scatter(xs, ys, c=colors, edgecolor="#444", s=35, linewidth=0.4)
        lim = axis_limit if axis_limit > 0 else max(_safe_max(xs), _safe_max(ys), 1.0) * 1.05
        ax.plot([0, lim], [0, lim], color="#888", linestyle="--", linewidth=0.8, label="y = x")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel(f"A:  {metric['label']}")
        ax.set_ylabel(f"B:  {metric['label']}")
        ax.set_title(f"{metric['label']} — A vs B {title_extra}".strip())
        ax.legend(loc="upper left", fontsize=8)
        ax.grid(True, linestyle=":", alpha=0.4)
        ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def generate_all_plots(per_launch, available_metrics, reforecast_section,
                       out_dir: str) -> dict:
    """Generate bar + scatter PNGs for every metric x forecast_type pair.

    Returns dict[(metric_key, forecast_type, kind)] -> filename (relative to out_dir).
    For reforecast: only GFS plots are generated.
    """
    files: dict = {}
    forecasts = sorted({r["forecast_type"] for r in per_launch}) or ["GFS"]

    for m in available_metrics:
        # Compute shared bounds across forecast types so visual comparison is honest.
        all_deltas = [
            r["metrics"][m["key"]]["delta"]
            for r in per_launch
            if not (isinstance(r["metrics"][m["key"]]["delta"], float)
                    and math.isnan(r["metrics"][m["key"]]["delta"]))
        ]
        y_abs_limit = _safe_max(all_deltas)
        all_vals = []
        for r in per_launch:
            for k in ("a", "b"):
                v = r["metrics"][m["key"]][k]
                if isinstance(v, float) and not math.isnan(v):
                    all_vals.append(v)
        axis_limit = (max(all_vals) * 1.05) if all_vals else 1.0

        for ft in forecasts:
            subset = [r for r in per_launch if r["forecast_type"] == ft]
            bar_fname     = f"plot_{m['key']}_{ft}_bar.png"
            scatter_fname = f"plot_{m['key']}_{ft}_scatter.png"
            _bar_chart(subset, m, os.path.join(out_dir, bar_fname),
                       y_abs_limit, title_extra=f"({ft})")
            _scatter_chart(subset, m, os.path.join(out_dir, scatter_fname),
                           axis_limit, title_extra=f"({ft})")
            files[(m["key"], ft, "bar")]     = bar_fname
            files[(m["key"], ft, "scatter")] = scatter_fname

    # Reforecast: GFS-only.
    if reforecast_section:
        m = reforecast_section["metric"]
        rows = reforecast_section["rows"]
        deltas = [r["metrics"][m["key"]]["delta"]
                  for r in rows
                  if not (isinstance(r["metrics"][m["key"]]["delta"], float)
                          and math.isnan(r["metrics"][m["key"]]["delta"]))]
        y_abs_limit = _safe_max(deltas)
        vals = []
        for r in rows:
            for k in ("a", "b"):
                v = r["metrics"][m["key"]][k]
                if isinstance(v, float) and not math.isnan(v):
                    vals.append(v)
        axis_limit = (max(vals) * 1.05) if vals else 1.0
        bar_fname     = f"plot_{m['key']}_GFS_bar.png"
        scatter_fname = f"plot_{m['key']}_GFS_scatter.png"
        _bar_chart(rows, m, os.path.join(out_dir, bar_fname),
                   y_abs_limit, title_extra="(GFS reforecast)")
        _scatter_chart(rows, m, os.path.join(out_dir, scatter_fname),
                       axis_limit, title_extra="(GFS reforecast)")
        files[(m["key"], "GFS", "bar")]     = bar_fname
        files[(m["key"], "GFS", "scatter")] = scatter_fname

    return files
