"""compare_reporting.py - HTML report for batch comparison.

Produces a single compare.html with:
  - Minimal header (batch IDs, timestamps, git hashes, notes for A and B)
  - Schema-mismatch warning banner (if any metrics were dropped)
  - Asymmetries section (only-in-A, only-in-B, failed-in-one)
  - TL;DR block (one line per primary metric x forecast type)
  - Aggregate table (Overall, per forecast_type, per campaign)
  - Win/Loss count table (per metric, per group)
  - Per-launch diff table (3 sub-cells per metric: A | B | Δ)
  - Per-metric overview plots (GFS and ERA5 side-by-side)
  - Reforecast diagnostics section (GFS-only, separate header)

Colored cells indicate improvement (green) or regression (red); gray = within
the noise floor. Columns are sortable (click headers).
"""

from __future__ import annotations

import html as _html
import math
import os
from typing import Optional

from evaluation.compare_batches import BatchData, PCT_FLOOR


_GREEN = "#b7e4b7"
_RED   = "#f4a5a5"
_GRAY  = "#ececec"
_FAIL  = "#e8e8e8"

_CLASS_BG = {"improved": _GREEN, "worsened": _RED, "unchanged": ""}


def _fv(v, fmt=".1f") -> str:
    try:
        fv = float(v)
        return "—" if math.isnan(fv) else f"{fv:{fmt}}"
    except (TypeError, ValueError):
        return "—"


def _signed(v, fmt=".1f") -> str:
    try:
        fv = float(v)
        if math.isnan(fv):
            return "—"
        return f"{fv:+{fmt}}"
    except (TypeError, ValueError):
        return "—"


def _pct(v) -> str:
    try:
        fv = float(v)
        if math.isnan(fv):
            return "—"
        return f"{fv:+.1f}%"
    except (TypeError, ValueError):
        return "—"


def _td(raw, txt, style="", align="center"):
    """One <td> with data-val for sorting."""
    style_attr = f' style="{style}"' if style else ""
    return f'<td data-val="{_html.escape(str(raw))}" style="text-align:{align}"{style_attr}>{txt}</td>'


def _delta_cell(m_cell: dict, fmt: str) -> str:
    delta = m_cell["delta"]
    pct   = m_cell["pct"]
    cls   = m_cell["class"]
    bg    = _CLASS_BG.get(cls, "")
    style = f"background:{bg}" if bg else ""
    if isinstance(delta, float) and math.isnan(delta):
        return _td("", "—", style)
    txt = f"{_signed(delta, fmt)} ({_pct(pct)})"
    return _td(delta, txt, style)


def _three_cells(m_cell: dict, fmt: str) -> str:
    a_str = _fv(m_cell["a"], fmt)
    b_str = _fv(m_cell["b"], fmt)
    return _td(m_cell["a"], a_str) + _td(m_cell["b"], b_str) + _delta_cell(m_cell, fmt)


# ── Header rows for tables ────────────────────────────────────────────────────


def _build_metric_header(metrics: list[dict], leading_cols: list[str]) -> str:
    """Two-row header: row1 = metric labels (colspan=3), row2 = 'A | B | Δ' sub-headers."""
    col_index = 0
    row1 = "<tr>"
    for c in leading_cols:
        row1 += f'<th rowspan="2" data-col="{col_index}">{_html.escape(c)}<span class="si"></span></th>'
        col_index += 1
    for m in metrics:
        row1 += f'<th colspan="3" style="background:#34495e">{_html.escape(m["label"])}</th>'
    row1 += "</tr>"

    row2 = "<tr>"
    for _ in metrics:
        row2 += (
            f'<th class="sub" data-col="{col_index}">A<span class="si"></span></th>'
        )
        col_index += 1
        row2 += (
            f'<th class="sub" data-col="{col_index}">B<span class="si"></span></th>'
        )
        col_index += 1
        row2 += (
            f'<th class="sub" data-col="{col_index}">Δ (Δ%)<span class="si"></span></th>'
        )
        col_index += 1
    row2 += "</tr>"
    return f"<thead>{row1}{row2}</thead>"


# ── Body builders ─────────────────────────────────────────────────────────────


def _build_per_launch_body(rows: list[dict], metrics: list[dict]) -> str:
    body = ""
    for i, r in enumerate(rows):
        cells = (
            _td(r["launch_id"], _html.escape(r["launch_id"]), align="left")
            + _td(r["forecast_type"], _html.escape(r["forecast_type"]))
            + _td(r["campaign"], _html.escape(r["campaign"] or "—"))
        )
        for m in metrics:
            cells += _three_cells(r["metrics"][m["key"]], m["fmt"])
        body += f'<tr data-orig="{i}">{cells}</tr>'
    return f"<tbody>{body}</tbody>"


def _build_aggregate_body(aggregates: list[dict], metrics: list[dict]) -> str:
    body = ""
    for i, agg in enumerate(aggregates):
        cells = (
            _td(agg["label"], _html.escape(agg["label"]), align="left")
            + _td("", "")
            + _td(agg["n"], str(agg["n"]))
        )
        for m in metrics:
            mc = agg["metrics"][m["key"]]
            pseudo = {"a": mc["a_mean"], "b": mc["b_mean"], "delta": mc["delta"],
                      "pct": mc["pct"], "class": mc["class"]}
            cells += _three_cells(pseudo, m["fmt"])
        body += f'<tr data-orig="{i}">{cells}</tr>'
    return f"<tbody>{body}</tbody>"


def _build_winloss_table(aggregates: list[dict], metrics: list[dict]) -> str:
    """Win/loss counts: rows = aggregate groups, columns = metrics."""
    header = "<tr><th>Group</th>"
    for m in metrics:
        header += f'<th>{_html.escape(m["label"])}</th>'
    header += "</tr>"

    body = ""
    for agg in aggregates:
        cells = f"<td style=\"text-align:left\">{_html.escape(agg['label'])}</td>"
        for m in metrics:
            c = agg["metrics"][m["key"]]["counts"]
            n = c["improved"] + c["unchanged"] + c["worsened"]
            if n == 0:
                cells += "<td>—</td>"
                continue
            cells += (
                f"<td><span style='color:#1a7a1a'>{c['improved']}</span> "
                f"/ {c['unchanged']} / "
                f"<span style='color:#a32424'>{c['worsened']}</span>"
                f"<br><span style='font-size:0.78em;color:#666'>(of {n})</span></td>"
            )
        body += f"<tr>{cells}</tr>"
    return f"<table><thead>{header}</thead><tbody>{body}</tbody></table>"


# ── Section builders ──────────────────────────────────────────────────────────


def _build_header_section(batch_a: BatchData, batch_b: BatchData) -> str:
    def info_block(label: str, b: BatchData) -> str:
        info = b.info or {}
        dirty = " <span style='color:orange'>(dirty)</span>" if info.get("git_dirty") else ""
        return (
            f"<div class='batch-card'>"
            f"<h3>{label}</h3>"
            f"<b>Batch ID:</b> {_html.escape(b.batch_id)}<br>"
            f"<b>Timestamp:</b> {_html.escape(str(info.get('timestamp_utc', '')))}<br>"
            f"<b>Git:</b> {_html.escape(info.get('git_hash', ''))}{dirty}<br>"
            f"<b>Note:</b> {_html.escape(str(info.get('note', '')))}"
            f"</div>"
        )
    return (
        "<div class='batch-row'>"
        + info_block("A — Baseline", batch_a)
        + info_block("B — Experiment", batch_b)
        + "</div>"
    )


def _build_tldr_section(lines: list[str]) -> str:
    body = _html.escape("\n".join(lines)) if lines else "(no primary metrics available)"
    return (
        "<div class='section'><h2>TL;DR</h2>"
        f"<pre style='background:#fff;padding:12px;border:1px solid #ccc;border-radius:4px;"
        f"font-size:0.95em;line-height:1.5'>{body}</pre></div>"
    )


def _build_asymmetries_section(asym: dict) -> str:
    def lst(keys):
        if not keys:
            return "<i>none</i>"
        items = sorted(f"{k[0]} ({k[1]})" for k in keys)
        return "<ul style='margin:4px 0 0 0;padding-left:18px'>" \
               + "".join(f"<li>{_html.escape(x)}</li>" for x in items) + "</ul>"
    any_asym = any(asym[k] for k in ("only_in_a", "only_in_b", "failed_in_a",
                                     "failed_in_b", "failed_in_both"))
    body = (
        f"<p><b>Intersection (ok in both):</b> {len(asym['intersection_ok'])} (launch, forecast) pairs</p>"
        f"<p><b>Only in A:</b> {lst(asym['only_in_a'])}</p>"
        f"<p><b>Only in B:</b> {lst(asym['only_in_b'])}</p>"
        f"<p><b>Failed in A only (succeeded in B):</b> {lst(asym['failed_in_a'])}</p>"
        f"<p><b>Failed in B only (succeeded in A):</b> {lst(asym['failed_in_b'])}</p>"
        f"<p><b>Failed in both:</b> {lst(asym['failed_in_both'])}</p>"
    )
    if not any_asym:
        body = (
            f"<p><b>Intersection (ok in both):</b> {len(asym['intersection_ok'])} (launch, forecast) pairs</p>"
            "<p>No asymmetries — both batches have the same launches with the same success status.</p>"
        )
    return f"<div class='section'><h2>Asymmetries</h2>{body}</div>"


def _build_schema_warning(missing: dict) -> str:
    if not missing:
        return ""
    lines = []
    for k, where in sorted(missing.items()):
        if where == "A":
            lines.append(f"  - <b>{_html.escape(k)}</b>: missing from batch A; column dropped")
        elif where == "B":
            lines.append(f"  - <b>{_html.escape(k)}</b>: missing from batch B; column dropped")
        else:
            lines.append(f"  - <b>{_html.escape(k)}</b>: missing from both batches; column dropped")
    return (
        "<div class='section' style='background:#fff7d6;border:1px solid #d2b400'>"
        "<h2>Schema mismatch warning</h2>"
        f"<ul>{''.join(f'<li>{l}</li>' for l in lines)}</ul>"
        "</div>"
    )


def _build_plots_section(plot_files: dict, metrics: list[dict],
                         reforecast_section, out_dir: str) -> str:
    """Per-metric: GFS and ERA5 thumbnails side-by-side (bar above, scatter below)."""
    def metric_block(m: dict, reforecast: bool = False) -> str:
        forecasts = ["GFS"] if reforecast else ["GFS", "ERA5"]
        bar_cells = ""
        sc_cells  = ""
        for ft in forecasts:
            bar_key = (m["key"], ft, "bar")
            sc_key  = (m["key"], ft, "scatter")
            if bar_key in plot_files:
                bar_cells += (
                    f"<td><div class='plot-label'>{_html.escape(ft)}</div>"
                    f"<img src='{plot_files[bar_key]}' alt='{m['key']} {ft} bar'></td>"
                )
            if sc_key in plot_files:
                sc_cells += (
                    f"<td><div class='plot-label'>{_html.escape(ft)}</div>"
                    f"<img src='{plot_files[sc_key]}' alt='{m['key']} {ft} scatter'></td>"
                )
        return (
            f"<div class='metric-plot-card'>"
            f"<h3>{_html.escape(m['label'])}</h3>"
            f"<table class='plot-grid'><tr>{bar_cells}</tr><tr>{sc_cells}</tr></table>"
            f"</div>"
        )

    parts = [metric_block(m) for m in metrics]
    body = "".join(parts)
    return f"<div class='section'><h2>Per-metric plots</h2>{body}</div>"


# ── Main entry ────────────────────────────────────────────────────────────────


_SORT_JS = """\
<script>
(function () {
  var states = {};
  document.querySelectorAll('th[data-col]').forEach(function (th) {
    th.style.cursor = 'pointer';
    th.title = 'Click to sort';
    th.addEventListener('click', function () {
      var table = th.closest('table');
      var tid   = table.id || (table.id = 'tbl_' + Math.random().toString(36).slice(2));
      var col   = parseInt(th.dataset.col, 10);
      var st    = states[tid] || {col: -1, dir: 0};
      table.querySelectorAll('th[data-col] .si').forEach(function (s) { s.textContent = ''; });
      if (st.col === col) {
        st.dir = st.dir === 1 ? -1 : (st.dir === -1 ? 0 : 1);
      } else { st.col = col; st.dir = 1; }
      states[tid] = st;
      var icon = th.querySelector('.si');
      if (icon) icon.textContent = st.dir === 1 ? ' ↑' : (st.dir === -1 ? ' ↓' : '');
      var tbody = table.querySelector('tbody');
      if (!tbody) return;
      var rows  = Array.from(tbody.querySelectorAll('tr'));
      if (st.dir === 0) {
        rows.sort(function (a, b) {
          return parseInt(a.dataset.orig || 0, 10) - parseInt(b.dataset.orig || 0, 10);
        });
      } else {
        rows.sort(function (a, b) {
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


_CSS = f"""
  body {{ font-family: monospace; margin: 24px; background: #f5f5f5; color: #222; }}
  h1   {{ margin-bottom: 4px; }}
  h2   {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 4px; margin-top: 0; }}
  h3   {{ color: #2c3e50; margin: 6px 0; }}
  .section {{
    background: #fff; border: 1px solid #ccc; border-radius: 6px;
    padding: 16px 24px; margin-bottom: 28px; box-shadow: 0 2px 6px rgba(0,0,0,0.06);
  }}
  .batch-row  {{ display: flex; gap: 16px; margin-bottom: 24px; }}
  .batch-card {{ flex: 1; background: #fff; border: 1px solid #ccc; border-radius: 6px;
                 padding: 12px 18px; box-shadow: 0 2px 6px rgba(0,0,0,0.06); font-size: 0.9em; }}
  .legend  {{ display: flex; gap: 12px; flex-wrap: wrap; margin: 6px 0 18px 0;
              font-size: 0.85em; align-items: center; }}
  .legend span {{ padding: 3px 10px; border: 1px solid #bbb; border-radius: 3px; }}
  table   {{ border-collapse: collapse; width: 100%; margin-bottom: 6px; font-size: 0.78em; }}
  th, td  {{ border: 1px solid #bbb; padding: 4px 6px; text-align: center; white-space: nowrap; }}
  th      {{ background: #2c3e50; color: #fff; font-weight: bold; user-select: none; }}
  th.sub  {{ background: #34495e; }}
  .si     {{ font-size: 0.8em; }}
  tr:hover {{ filter: brightness(0.96); }}
  .metric-plot-card {{
    background: #fafafa; border: 1px solid #ddd; border-radius: 4px;
    padding: 10px 14px; margin-bottom: 14px;
  }}
  .plot-grid img {{ max-width: 100%; height: auto; display: block; }}
  .plot-grid td  {{ vertical-align: top; padding: 6px; border: none; }}
  .plot-label    {{ font-size: 0.8em; color: #666; margin-bottom: 2px; text-align: center; }}
"""


def write_compare_html(out_dir: str, batch_a: BatchData, batch_b: BatchData,
                       per_launch: list[dict], aggregates: list[dict],
                       reforecast_section: Optional[dict],
                       asymmetries: dict, missing: dict,
                       available_metrics: list[dict], tldr_lines: list[str],
                       plot_files: dict) -> str:
    leading = ["Launch", "Forecast", "Campaign"]
    per_launch_header = _build_metric_header(available_metrics, leading)
    per_launch_body   = _build_per_launch_body(per_launch, available_metrics)

    agg_leading = ["Group", "—", "N"]
    aggregate_header = _build_metric_header(available_metrics, agg_leading)
    aggregate_body   = _build_aggregate_body(aggregates, available_metrics)

    winloss_table = _build_winloss_table(aggregates, available_metrics)

    plots_section = _build_plots_section(plot_files, available_metrics,
                                         reforecast_section, out_dir)

    reforecast_html = ""
    if reforecast_section:
        ref_m = reforecast_section["metric"]
        ref_rows = reforecast_section["rows"]
        ref_agg  = reforecast_section["aggregate"]
        ref_header = _build_metric_header([ref_m], ["Launch", "Forecast", "Campaign"])
        ref_body   = _build_per_launch_body(ref_rows, [ref_m])
        ref_agg_header = _build_metric_header([ref_m], ["Group", "—", "N"])
        ref_agg_body   = _build_aggregate_body([ref_agg], [ref_m])
        # Plot
        ref_plot = ""
        bar_key = (ref_m["key"], "GFS", "bar")
        sc_key  = (ref_m["key"], "GFS", "scatter")
        if bar_key in plot_files or sc_key in plot_files:
            cells = ""
            if bar_key in plot_files:
                cells += f"<td><img src='{plot_files[bar_key]}' alt='reforecast bar'></td>"
            if sc_key in plot_files:
                cells += f"<td><img src='{plot_files[sc_key]}' alt='reforecast scatter'></td>"
            ref_plot = (
                f"<div class='metric-plot-card'><h3>{_html.escape(ref_m['label'])}</h3>"
                f"<table class='plot-grid'><tr>{cells}</tr></table></div>"
            )
        reforecast_html = (
            "<div class='section'><h2>Reforecast diagnostics (GFS-only)</h2>"
            "<p style='color:#555;font-size:0.85em'>Reforecast = sim re-run against truth altitude profile, "
            "isolating wind-forecast error from altitude-model error.</p>"
            f"<h3>Aggregate</h3><table id='tbl-ref-agg'>{ref_agg_header}{ref_agg_body}</table>"
            f"<h3>Per launch</h3><table id='tbl-ref-detail'>{ref_header}{ref_body}</table>"
            f"{ref_plot}"
            "</div>"
        )

    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>EarthSHAB Batch Comparison: {_html.escape(batch_a.batch_id)} vs {_html.escape(batch_b.batch_id)}</title>
  <style>{_CSS}</style>
</head>
<body>
  <h1>EarthSHAB Batch Comparison</h1>
  {_build_header_section(batch_a, batch_b)}
  {_build_schema_warning(missing)}

  <div class='legend'>
    <b>Cell colors (Δ = B − A):</b>
    <span style='background:{_GREEN}'>improved</span>
    <span style='background:{_RED}'>worsened</span>
    <span style='background:{_GRAY}'>unchanged (below noise floor)</span>
    &nbsp;&nbsp;<span style='color:#666'>Click any column header to sort ↑↓</span>
  </div>

  {_build_tldr_section(tldr_lines)}
  {_build_asymmetries_section(asymmetries)}

  <div class='section'>
    <h2>Aggregates</h2>
    <table id='tbl-aggregates'>{aggregate_header}{aggregate_body}</table>
  </div>

  <div class='section'>
    <h2>Win / Unchanged / Loss counts</h2>
    <p style='font-size:0.85em;color:#666'>
      Per-launch deltas classified by metric. Green count = launches that improved,
      red count = launches that worsened, middle = within noise floor.
    </p>
    {winloss_table}
  </div>

  <div class='section'>
    <h2>Per-launch deltas</h2>
    <table id='tbl-per-launch'>{per_launch_header}{per_launch_body}</table>
  </div>

  {plots_section}
  {reforecast_html}

  {_SORT_JS}
</body>
</html>"""

    path = os.path.join(out_dir, "compare.html")
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    return path
