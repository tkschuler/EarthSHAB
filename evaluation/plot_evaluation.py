import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import fluids


def plot_comparison(time_sim, el_sim, v_sim, T_atm_sim, df_truth,
                    sim_phases=None, truth_phases=None,
                    t_sunset_sim=None, t_sunset_truth=None):
    """Evaluation comparison figure.

    Layout (3-row gridspec):
      Row 0: Altitude (full width) — coloured phase boxes + float level lines.
             Boxes are derived from float_mask (slow motion only), so transition
             curves at the top are left unshaded — a visible gap appears between
             the ascent box, the float box, and the descent box.
      Row 1: Ascent velocity | Float vertical motion | Descent velocity
             Each subplot is zoomed to the corresponding phase window.
      Row 2: Temperature (2/3 width) | Pressure (1/3 width)
    """
    sns.set_style("darkgrid")
    plt.rcParams.update({'font.size': 11})

    el    = np.array(el_sim)
    v     = np.array(v_sim)
    T_atm = np.array(T_atm_sim)

    # Unify sim timestamps to pd.Timestamp for consistent comparisons
    t_sim_pd = pd.DatetimeIndex([pd.Timestamp(t) for t in time_sim])

    # Subsample simulation to ~1 point per minute for plotting speed
    step    = 60
    t_plot  = t_sim_pd[::step]
    el_plot = el[::step]
    v_plot  = v[::step]
    T_plot  = T_atm[::step]

    # Truth timestamps as a DatetimeIndex (used throughout)
    if df_truth is not None and len(df_truth) > 0:
        t_truth_all = pd.DatetimeIndex([pd.Timestamp(t) for t in df_truth['time']])
    else:
        t_truth_all = pd.DatetimeIndex([])

    # ── Compute all phase boundaries upfront ─────────────────────────────────
    if sim_phases is not None:
        i_enter = sim_phases['i_enter']
        i_exit  = sim_phases['i_exit']
        asc_idx = np.where(sim_phases['ascent_mask'])[0]
        flt_idx = np.where(sim_phases['float_mask'])[0]
        des_idx = np.where(sim_phases['descent_mask'])[0]

        t_asc_start = t_sim_pd[0]                                              # always from the start
        t_asc_end   = t_sim_pd[asc_idx[-1]] if asc_idx.size else t_sim_pd[i_enter]
        # float_mask uses BOTH altitude and velocity criteria, so its boundaries
        # sit inside the curved transition regions → creates a visible gap
        t_flt_start = t_sim_pd[flt_idx[0]]  if flt_idx.size else t_sim_pd[i_enter]
        t_flt_end   = t_sim_pd[flt_idx[-1]] if flt_idx.size else t_sim_pd[i_exit]
        t_des_start = t_sim_pd[des_idx[0]]  if des_idx.size else t_sim_pd[i_exit]
        t_des_end   = t_sim_pd[des_idx[-1]] if des_idx.size else t_sim_pd[-1]
    else:
        n = len(t_sim_pd)
        t_asc_start, t_asc_end = t_sim_pd[0],          t_sim_pd[n // 3]
        t_flt_start, t_flt_end = t_sim_pd[n // 3],     t_sim_pd[2 * n // 3]
        t_des_start, t_des_end = t_sim_pd[2 * n // 3], t_sim_pd[-1]

    # Truth phase boundaries (all three phases, indexed to df_truth)
    if truth_phases is not None and len(t_truth_all) > 0:
        asc_t_idx = np.where(truth_phases['ascent_mask'])[0]
        flt_t_idx = np.where(truth_phases['float_mask'])[0]
        des_t_idx = np.where(truth_phases['descent_mask'])[0]
        i_enter_t = truth_phases['i_enter']
        i_exit_t  = truth_phases['i_exit']
        t_asc_start_t = t_truth_all[0]                                                  # always from the start
        t_asc_end_t   = t_truth_all[asc_t_idx[-1]] if asc_t_idx.size else t_truth_all[i_enter_t]
        t_flt_start_t = t_truth_all[flt_t_idx[0]]  if flt_t_idx.size else t_truth_all[i_enter_t]
        t_flt_end_t   = t_truth_all[flt_t_idx[-1]] if flt_t_idx.size else t_truth_all[i_exit_t]
        t_des_start_t = t_truth_all[des_t_idx[0]]  if des_t_idx.size else t_truth_all[i_exit_t]
        t_des_end_t   = t_truth_all[des_t_idx[-1]] if des_t_idx.size else t_truth_all[-1]
    else:
        t_asc_start_t = t_asc_start
        t_asc_end_t   = t_asc_end
        t_flt_start_t = t_flt_start
        t_flt_end_t   = t_flt_end
        t_des_start_t = t_des_start
        t_des_end_t   = t_des_end

    # ── Sunset timestamps and elapsed helpers ────────────────────────────────
    t_ss_sim   = pd.Timestamp(t_sunset_sim)   if t_sunset_sim   is not None else None
    t_ss_truth = pd.Timestamp(t_sunset_truth) if t_sunset_truth is not None else None

    _sim_sunset_kw   = dict(color='royalblue',  linewidth=1.5, linestyle=':', zorder=5)
    _truth_sunset_kw = dict(color='darkorange', linewidth=1.5, linestyle=':', zorder=5)

    t0_sim   = t_sim_pd[0]
    t0_truth = t_truth_all[0] if len(t_truth_all) > 0 else t0_sim

    # Launch-elapsed minutes for temp/pressure plots
    sim_elapsed_launch   = (t_plot - t0_sim).total_seconds() / 60
    truth_elapsed_launch = ((t_truth_all - t0_truth).total_seconds() / 60
                            if len(t_truth_all) > 0 else np.array([]))

    ss_sim_elapsed   = (t_ss_sim   - t0_sim).total_seconds()   / 60 if t_ss_sim   is not None else None
    ss_truth_elapsed = (t_ss_truth - t0_truth).total_seconds() / 60 if t_ss_truth is not None else None

    # ── Layout ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 12))
    gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.46, wspace=0.35)

    ax_alt  = fig.add_subplot(gs[0, :])    # full-width altitude
    ax_asc  = fig.add_subplot(gs[1, 0])    # ascent velocity window
    ax_flt  = fig.add_subplot(gs[1, 1])    # float velocity window
    ax_des  = fig.add_subplot(gs[1, 2])    # descent velocity window
    ax_temp = fig.add_subplot(gs[2, 0:2])  # temperature
    ax_pres = fig.add_subplot(gs[2, 2])    # pressure

    sim_line_kw   = dict(color='royalblue',  label='Simulation')
    truth_line_kw = dict(color='darkorange', linestyle='--', marker='.', markersize=4,
                         label='Ground Truth (APRS)')

    # ── Altitude plot ─────────────────────────────────────────────────────────
    ax_alt.plot(t_plot, el_plot, **sim_line_kw)
    ax_alt.plot(df_truth['time'], df_truth['altitude'], **truth_line_kw)
    ax_alt.set_xlabel('Time (MST)')
    ax_alt.set_ylabel('Altitude (m)')
    ax_alt.set_title('Altitude')

    # ── 5 phase boxes: sim ascent, truth ascent, float, sim descent, truth descent ──
    # Sim boxes: solid fill.  Truth boxes: hatched, same color family.
    _span_kw = dict(zorder=0)
    ax_alt.axvspan(t_asc_start,   t_asc_end,   alpha=0.18, color='limegreen',
                   label='Sim ascent window',   **_span_kw)
    ax_alt.axvspan(t_asc_start_t, t_asc_end_t, alpha=0.12, color='limegreen',
                   hatch='//', label='Truth ascent window', **_span_kw)
    ax_alt.axvspan(t_flt_start,   t_flt_end,   alpha=0.13, color='mediumpurple',
                   label='Float window',        **_span_kw)
    ax_alt.axvspan(t_des_start,   t_des_end,   alpha=0.18, color='tomato',
                   label='Sim descent window',  **_span_kw)
    ax_alt.axvspan(t_des_start_t, t_des_end_t, alpha=0.12, color='tomato',
                   hatch='//', label='Truth descent window', **_span_kw)

    # Boundary vertical lines — one per phase edge (sim + truth)
    _bnd_kw = dict(linewidth=1.5, alpha=0.65, zorder=2)
    for t_b, clr in [
        (t_asc_start,   'limegreen'),    (t_asc_end,   'limegreen'),
        (t_asc_start_t, 'limegreen'),    (t_asc_end_t, 'limegreen'),
        (t_flt_start,   'mediumpurple'), (t_flt_end,   'mediumpurple'),
        (t_des_start,   'tomato'),       (t_des_end,   'tomato'),
        (t_des_start_t, 'tomato'),       (t_des_end_t, 'tomato'),
    ]:
        ax_alt.axvline(t_b, color=clr, **_bnd_kw)

    # Float altitude estimate — bounded box (mean line + ±std bracket)
    def _draw_float_alt_box(ax, fa_mean, fa_std, t0, t1, color, label):
        """Draw mean line and ±std bracket bounded to [t0, t1]."""
        if math.isnan(fa_mean):
            return
        # Mean line spanning the float window
        ax.plot([t0, t1], [fa_mean, fa_mean], color=color, linestyle='--',
                linewidth=1.8, alpha=0.85, label=label, zorder=4)
        if not math.isnan(fa_std) and fa_std > 0:
            lo, hi = fa_mean - fa_std, fa_mean + fa_std
            bracket_kw = dict(color=color, linestyle=':', linewidth=1.2, alpha=0.75, zorder=4)
            # Top and bottom horizontal rails
            ax.plot([t0, t1], [hi, hi], **bracket_kw)
            ax.plot([t0, t1], [lo, lo], **bracket_kw)
            # Vertical side connectors (only span mean±std, not down to zero)
            ax.plot([t0, t0], [lo, hi], **bracket_kw)
            ax.plot([t1, t1], [lo, hi], **bracket_kw)

    if sim_phases is not None:
        fa_mean = sim_phases.get('float_alt_mean', math.nan)
        fa_std  = sim_phases.get('float_alt_std',  math.nan)
        lbl = (f'Sim float: {fa_mean:.0f} ± {fa_std:.0f} m'
               if not math.isnan(fa_mean) and not math.isnan(fa_std)
               else f'Sim float: {fa_mean:.0f} m' if not math.isnan(fa_mean) else None)
        if lbl:
            _draw_float_alt_box(ax_alt, fa_mean, fa_std,
                                t_flt_start, t_flt_end, 'royalblue', lbl)

    if truth_phases is not None:
        fa_mean_t = truth_phases.get('float_alt_mean', math.nan)
        fa_std_t  = truth_phases.get('float_alt_std',  math.nan)
        lbl = (f'Truth float: {fa_mean_t:.0f} ± {fa_std_t:.0f} m'
               if not math.isnan(fa_mean_t) and not math.isnan(fa_std_t)
               else f'Truth float: {fa_mean_t:.0f} m' if not math.isnan(fa_mean_t) else None)
        if lbl:
            _draw_float_alt_box(ax_alt, fa_mean_t, fa_std_t,
                                t_flt_start_t, t_flt_end_t, 'darkorange', lbl)

    # Sunsets on the altitude plot (absolute time axis)
    for t_ss, kw, lbl in [
        (t_ss_sim,   _sim_sunset_kw,   'Sim Sunset'),
        (t_ss_truth, _truth_sunset_kw, 'Truth Sunset'),
    ]:
        if t_ss is not None:
            ax_alt.axvline(t_ss, label=lbl, **kw)
    ax_alt.legend(loc='lower left', ncol=3, fontsize=9)

    # ── Velocity windows (elapsed time, phase-relative) ───────────────────────
    def _vel_window(ax, sim_ts, sim_te, truth_ts, truth_te, title):
        """Compare sim vs truth velocity on a shared elapsed-time axis.

        Sim elapsed starts at sim_ts; truth elapsed starts at truth_ts.
        This aligns the two trajectories even when their phase windows
        do not overlap in wall-clock time.
        """
        sim_ts, sim_te     = pd.Timestamp(sim_ts),   pd.Timestamp(sim_te)
        truth_ts, truth_te = pd.Timestamp(truth_ts), pd.Timestamp(truth_te)

        # Sim
        sim_mask = (t_plot >= sim_ts) & (t_plot <= sim_te)
        if sim_mask.any():
            sim_e = (t_plot[sim_mask] - sim_ts).total_seconds() / 60
            ax.plot(sim_e, v_plot[sim_mask], color='royalblue', label='Simulation')

        # Truth
        if len(t_truth_all) > 0:
            truth_mask = (t_truth_all >= truth_ts) & (t_truth_all <= truth_te)
            if truth_mask.any():
                truth_e = (t_truth_all[truth_mask] - truth_ts).total_seconds() / 60
                ax.plot(truth_e, df_truth['v_truth'].values[truth_mask],
                        color='darkorange', linestyle='--', marker='.', markersize=4,
                        label='Ground Truth')

        ax.axhline(0, color='gray', linewidth=0.8, linestyle=':')

        # Sunsets in phase-relative elapsed time
        if t_ss_sim is not None and sim_ts <= t_ss_sim <= sim_te:
            ax.axvline((t_ss_sim - sim_ts).total_seconds() / 60,
                       label='Sim Sunset', **_sim_sunset_kw)
        if t_ss_truth is not None and truth_ts <= t_ss_truth <= truth_te:
            ax.axvline((t_ss_truth - truth_ts).total_seconds() / 60,
                       label='Truth Sunset', **_truth_sunset_kw)

        ax.set_xlabel('Elapsed Time (min)')
        ax.set_ylabel('Vertical Velocity (m/s)')
        ax.set_title(title)
        ax.legend(fontsize=8)

    _vel_window(ax_asc, t_asc_start, t_asc_end, t_asc_start_t, t_asc_end_t,  'Ascent Rate')
    _vel_window(ax_flt, t_flt_start, t_flt_end, t_flt_start_t, t_flt_end_t,  'Float Vertical Motion')
    _vel_window(ax_des, t_des_start, t_des_end, t_des_start_t, t_des_end_t,  'Descent Rate')

    # ── Temperature (elapsed from each trajectory's launch) ───────────────────
    ax_temp.plot(sim_elapsed_launch, T_plot, **sim_line_kw)
    if len(t_truth_all) > 0:
        valid_temp = df_truth['temp_k'].notna().values
        if valid_temp.any():
            ax_temp.scatter(truth_elapsed_launch[valid_temp],
                            df_truth['temp_k'].values[valid_temp],
                            color='darkorange', s=12, zorder=3, label='Ground Truth (APRS)')
    if ss_sim_elapsed is not None:
        ax_temp.axvline(ss_sim_elapsed,   label='Sim Sunset',   **_sim_sunset_kw)
    if ss_truth_elapsed is not None:
        ax_temp.axvline(ss_truth_elapsed, label='Truth Sunset', **_truth_sunset_kw)
    ax_temp.set_xlabel('Elapsed Time (min)')
    ax_temp.set_ylabel('Temperature (K)')
    ax_temp.set_title('Temperature')
    ax_temp.legend()

    # ── Pressure (elapsed from each trajectory's launch) ──────────────────────
    p_sim_kpa = np.array([fluids.atmosphere.ATMOSPHERE_1976(e).P / 1000 for e in el_plot])
    ax_pres.plot(sim_elapsed_launch, p_sim_kpa, **sim_line_kw)
    if len(t_truth_all) > 0:
        valid_pres = df_truth['pressure_pa'].notna().values
        if valid_pres.any():
            ax_pres.scatter(truth_elapsed_launch[valid_pres],
                            df_truth['pressure_pa'].values[valid_pres] / 1000,
                            color='darkorange', s=12, zorder=3, label='Ground Truth (APRS)')
    if ss_sim_elapsed is not None:
        ax_pres.axvline(ss_sim_elapsed,   label='Sim Sunset',   **_sim_sunset_kw)
    if ss_truth_elapsed is not None:
        ax_pres.axvline(ss_truth_elapsed, label='Truth Sunset', **_truth_sunset_kw)
    ax_pres.set_xlabel('Elapsed Time (min)')
    ax_pres.set_ylabel('Pressure (kPa)')
    ax_pres.set_title('Pressure')
    ax_pres.legend()

    fig.suptitle('Simulation vs Ground Truth', fontsize=14, y=1.00)
    plt.tight_layout()
    return fig
