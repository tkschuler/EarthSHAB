import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import seaborn as sns
import fluids


def plot_comparison(time_sim, el_sim, v_sim, T_atm_sim, df_truth,
                    sim_phases=None, truth_phases=None):
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

    # ── Compute all phase boundaries upfront ─────────────────────────────────
    if sim_phases is not None:
        i_enter = sim_phases['i_enter']
        i_exit  = sim_phases['i_exit']
        asc_idx = np.where(sim_phases['ascent_mask'])[0]
        flt_idx = np.where(sim_phases['float_mask'])[0]
        des_idx = np.where(sim_phases['descent_mask'])[0]

        t_asc_start = t_sim_pd[asc_idx[0]]  if asc_idx.size else t_sim_pd[0]
        t_asc_end   = t_sim_pd[asc_idx[-1]] if asc_idx.size else t_sim_pd[i_enter]
        # float_mask uses BOTH altitude and velocity criteria, so its boundaries
        # sit inside the curved transition regions → creates a visible gap
        t_flt_start = t_sim_pd[flt_idx[0]]  if flt_idx.size else t_sim_pd[i_enter]
        t_flt_end   = t_sim_pd[flt_idx[-1]] if flt_idx.size else t_sim_pd[i_exit]
        t_des_start = t_sim_pd[des_idx[0]]  if des_idx.size else t_sim_pd[i_exit]
        t_des_end   = t_sim_pd[des_idx[-1]] if des_idx.size else t_sim_pd[-1]
    else:
        n = len(t_sim_pd)
        t_asc_start, t_asc_end = t_sim_pd[0],     t_sim_pd[n // 3]
        t_flt_start, t_flt_end = t_sim_pd[n // 3], t_sim_pd[2 * n // 3]
        t_des_start, t_des_end = t_sim_pd[2 * n // 3], t_sim_pd[-1]

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

    # Phase boxes — gaps between ascent/float/descent are the transition curves
    ax_alt.axvspan(t_asc_start, t_asc_end,  alpha=0.13, color='limegreen',
                   zorder=0, label='Ascent window')
    ax_alt.axvspan(t_flt_start, t_flt_end,  alpha=0.13, color='mediumpurple',
                   zorder=0, label='Float window')
    ax_alt.axvspan(t_des_start, t_des_end,  alpha=0.13, color='tomato',
                   zorder=0, label='Descent window')

    # Float altitude estimate lines
    if sim_phases is not None:
        fa_mean = sim_phases.get('float_alt_mean', math.nan)
        fa_std  = sim_phases.get('float_alt_std',  math.nan)
        if not math.isnan(fa_mean):
            lbl = (f'Sim float: {fa_mean:.0f} ± {fa_std:.0f} m'
                   if not math.isnan(fa_std) else f'Sim float: {fa_mean:.0f} m')
            ax_alt.axhline(fa_mean, color='royalblue', linestyle='--',
                           linewidth=1.8, alpha=0.85, label=lbl)

    if truth_phases is not None:
        fa_mean_t = truth_phases.get('float_alt_mean', math.nan)
        fa_std_t  = truth_phases.get('float_alt_std',  math.nan)
        if not math.isnan(fa_mean_t):
            lbl = (f'Truth float: {fa_mean_t:.0f} ± {fa_std_t:.0f} m'
                   if not math.isnan(fa_std_t) else f'Truth float: {fa_mean_t:.0f} m')
            ax_alt.axhline(fa_mean_t, color='darkorange', linestyle='--',
                           linewidth=1.8, alpha=0.85, label=lbl)

    ax_alt.legend(loc='lower left', ncol=3, fontsize=9)

    # ── Velocity windows ─────────────────────────────────────────────────────
    def _vel_window(ax, t_start, t_end, title):
        """Plot sim and APRS velocity zoomed to a given time window."""
        t_s = pd.Timestamp(t_start)
        t_e = pd.Timestamp(t_end)

        sim_mask  = (t_plot >= t_s) & (t_plot <= t_e)
        aprs_mask = (df_truth['time'] >= t_s) & (df_truth['time'] <= t_e)

        ax.plot(t_plot[sim_mask], v_plot[sim_mask], color='royalblue', label='Simulation')
        ax.plot(df_truth['time'][aprs_mask], df_truth['v_truth'][aprs_mask],
                color='darkorange', linestyle='--', marker='.', markersize=4,
                label='Ground Truth')
        ax.axhline(0, color='gray', linewidth=0.8, linestyle=':')
        ax.set_xlabel('Time (MST)')
        ax.set_ylabel('Vertical Velocity (m/s)')
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.setp(ax.get_xticklabels(), rotation=20, ha='right')

    # Each velocity subplot matches exactly the corresponding altitude box
    _vel_window(ax_asc, t_asc_start, t_asc_end,  'Ascent Rate')
    _vel_window(ax_flt, t_flt_start, t_flt_end,  'Float Vertical Motion')
    _vel_window(ax_des, t_des_start, t_des_end,  'Descent Rate')

    # ── Temperature ───────────────────────────────────────────────────────────
    ax_temp.plot(t_plot, T_plot, **sim_line_kw)
    valid_temp = df_truth['temp_k'].notna()
    ax_temp.scatter(df_truth['time'][valid_temp], df_truth['temp_k'][valid_temp],
                    color='darkorange', s=12, zorder=3, label='Ground Truth (APRS)')
    ax_temp.set_xlabel('Time (MST)')
    ax_temp.set_ylabel('Temperature (K)')
    ax_temp.set_title('Temperature')
    ax_temp.legend()

    # ── Pressure ─────────────────────────────────────────────────────────────
    p_sim_kpa = np.array([fluids.atmosphere.ATMOSPHERE_1976(e).P / 1000 for e in el_plot])
    ax_pres.plot(t_plot, p_sim_kpa, **sim_line_kw)
    valid_pres = df_truth['pressure_pa'].notna()
    ax_pres.scatter(df_truth['time'][valid_pres], df_truth['pressure_pa'][valid_pres] / 1000,
                    color='darkorange', s=12, zorder=3, label='Ground Truth (APRS)')
    ax_pres.set_xlabel('Time (MST)')
    ax_pres.set_ylabel('Pressure (kPa)')
    ax_pres.set_title('Pressure')
    ax_pres.legend()

    fig.suptitle('Simulation vs Ground Truth', fontsize=14, y=1.00)
    plt.tight_layout()
