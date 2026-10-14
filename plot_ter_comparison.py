#!/usr/bin/env python3
"""Plot TER flux comparison: truth / prior / optimized on a regular lat-lon grid.

Layout  (2 rows × 3 panels):
  row 1 — absolute TER flux (vmin=0, shared colorbar):  Truth | Prior | Optimized
  row 2 — departure from truth (diverging, shared colorbar):
           stats box  |  Prior-Truth  |  Opt-Truth

Months to average are selected from the assimilation window [--t-start, --t-end).
CS data are interpolated to a regular 1°×1° grid with scipy griddata.

Usage:
  python3 plot_ter_comparison.py \\
      --truth    <osse_4dvar>/obs/truth_emissions_monthly_2016Q1.nc \\
      --state    <exp_dir>/4dvar_state_osse.npz \\
      --out      <plot_dir>/ter_comparison_iter001.png \\
      --t-start  2016-01-01 \\
      --t-end    2016-04-01 \\
      [--iter-label 001] \\
      [--prior-scale 1.3] \\
      [--res 1.0]
"""
import argparse
import os
import sys
from datetime import datetime

import numpy as np
import netCDF4 as nc
from netCDF4 import num2date
from scipy.interpolate import griddata

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import plot_4dvar as _p4d
from plot_4dvar import NF

# kgCO2/m²/s → gC/m²/day
_TO_GC_DAY = 1e3 * (12.0 / 44.0) * 86400.0

MONTH_LABELS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']


def _window_month_indices(truth_file, t_start, t_end):
    """Return indices of months in truth_file that fall within [t_start, t_end)."""
    ds = nc.Dataset(truth_file)
    cftimes = num2date(ds['time'][:], ds['time'].units)
    ds.close()
    keep = [i for i, ct in enumerate(cftimes)
            if t_start <= datetime(ct.year, ct.month, ct.day) < t_end]
    return keep, [cftimes[i] for i in keep]


def _month_label(cftimes):
    """Human-readable label for a list of cftime objects."""
    names = [MONTH_LABELS[ct.month - 1] for ct in cftimes]
    year  = cftimes[0].year
    if len(names) == 1:
        return f'{names[0]} {year}'
    return f'{names[0]}–{names[-1]} {year} mean'


def cs_to_latlon(data_cs, lats_cs, lons_cs, res=1.0):
    lat_reg = np.arange(-90,  90 + res, res)
    lon_reg = np.arange(-180, 180 + res, res)
    lon2d, lat2d = np.meshgrid(lon_reg, lat_reg)
    pts  = np.column_stack([lons_cs.ravel(), lats_cs.ravel()])
    vals = data_cs.ravel()
    grid = griddata(pts, vals, (lon2d, lat2d), method='linear')
    nans = np.isnan(grid)
    if nans.any():
        grid[nans] = griddata(pts, vals, (lon2d[nans], lat2d[nans]), method='nearest')
    return lat_reg, lon_reg, grid


def _draw_panel(ax, lat, lon, field, vmin, vmax, cmap, title):
    lon2d, lat2d = np.meshgrid(lon, lat)
    pm = ax.pcolormesh(lon2d, lat2d, field,
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       shading='auto', rasterized=True)
    _p4d._load_coast()
    for pts in _p4d._COAST_LINES:
        ax.plot(pts[:, 0], pts[:, 1], color='0.25', linewidth=0.4, zorder=2)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-90,   91, 30))
    ax.tick_params(labelsize=7)
    ax.grid(True, alpha=0.25, linewidth=0.4)
    ax.set_title(title, fontsize=10, pad=4)
    return pm


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--truth',       required=True,
                    help='truth_emissions_monthly_*.nc from semi-obs forward run')
    ap.add_argument('--state',       required=True,
                    help='4dvar_state_osse.npz (contains x_best / it_best)')
    ap.add_argument('--out',         required=True, help='Output PNG path')
    ap.add_argument('--t-start',     required=True,
                    help='Assimilation window start date  YYYY-MM-DD')
    ap.add_argument('--t-end',       required=True,
                    help='Assimilation window end date    YYYY-MM-DD (exclusive)')
    ap.add_argument('--iter-label',  default=None,
                    help='Zero-padded iteration number for display (e.g. 001)')
    ap.add_argument('--prior-scale', type=float, default=1.3)
    ap.add_argument('--res',         type=float, default=1.0,
                    help='Output grid resolution in degrees')
    args = ap.parse_args()

    t_start = datetime.strptime(args.t_start, '%Y-%m-%d')
    t_end   = datetime.strptime(args.t_end,   '%Y-%m-%d')

    # ---- select months within the assimilation window ----
    keep_idx, keep_times = _window_month_indices(args.truth, t_start, t_end)
    if not keep_idx:
        raise SystemExit(
            f'No months in {args.truth} fall within [{args.t_start}, {args.t_end})')
    month_str = _month_label(keep_times)

    # ---- load truth emissions ----
    ds = nc.Dataset(args.truth)
    all_months = np.array(ds['EmisCO2_NetTerrExch'][:])  # (nmon, 6, 24, 24) kgCO2/m²/s
    lats_cs    = np.array(ds['lats'][:])
    lons_cs    = ((np.array(ds['lons'][:]) + 180.0) % 360.0) - 180.0
    ds.close()

    ter_raw      = all_months[keep_idx].mean(axis=0)   # (6, 24, 24)
    cs_res       = ter_raw.shape[-1]
    ter_truth_cs = ter_raw * _TO_GC_DAY                # gC/m²/day

    # ---- load σ_best ----
    if not os.path.exists(args.state):
        print(f'No state file yet ({args.state}); skipping TER comparison plot.')
        return
    st       = np.load(args.state, allow_pickle=True)
    sigma_cs = st['x_best'].reshape(NF, cs_res, cs_res)
    it_best  = int(st['it_best'])
    iter_tag = args.iter_label if args.iter_label else str(it_best)

    # ---- derive flux fields ----
    ter_prior_cs = ter_truth_cs * args.prior_scale
    ter_opt_cs   = ter_truth_cs * args.prior_scale * sigma_cs

    # ---- interpolate to regular grid ----
    lat, lon, truth_ll = cs_to_latlon(ter_truth_cs, lats_cs, lons_cs, args.res)
    _,   _,   prior_ll = cs_to_latlon(ter_prior_cs, lats_cs, lons_cs, args.res)
    _,   _,   opt_ll   = cs_to_latlon(ter_opt_cs,   lats_cs, lons_cs, args.res)

    diff_prior_ll = prior_ll - truth_ll
    diff_opt_ll   = opt_ll   - truth_ll

    # ---- colour limits ----
    abs_max = np.ceil(np.nanmax([truth_ll, prior_ll, opt_ll]) * 10) / 10
    dif_max = np.ceil(np.nanmax(np.abs([diff_prior_ll, diff_opt_ll])) * 10) / 10

    target_sigma = 1.0 / args.prior_scale

    # ---- layout: 2 rows × 3 data cols + 1 narrow cbar col ----
    fig = plt.figure(figsize=(18, 9))
    gs = fig.add_gridspec(2, 4,
                          width_ratios=[1, 1, 1, 0.04],
                          hspace=0.32, wspace=0.08,
                          left=0.04, right=0.96, top=0.90, bottom=0.05)

    ax_truth      = fig.add_subplot(gs[0, 0])
    ax_prior      = fig.add_subplot(gs[0, 1])
    ax_opt        = fig.add_subplot(gs[0, 2])
    cax1          = fig.add_subplot(gs[0, 3])
    ax_box        = fig.add_subplot(gs[1, 0])
    ax_diff_prior = fig.add_subplot(gs[1, 1])
    ax_diff_opt   = fig.add_subplot(gs[1, 2])
    cax2          = fig.add_subplot(gs[1, 3])

    # ---- row 1: absolute fluxes ----
    pm1 = _draw_panel(ax_truth, lat, lon, truth_ll,
                      0, abs_max, 'YlOrRd', f'Truth TER — {month_str}')
    _draw_panel(ax_prior, lat, lon, prior_ll,
                0, abs_max, 'YlOrRd',
                f'Prior TER (×{args.prior_scale:.1f}) — {month_str}')
    _draw_panel(ax_opt, lat, lon, opt_ll,
                0, abs_max, 'YlOrRd',
                f'Optimized TER (iter {iter_tag}) — {month_str}')

    cb1 = fig.colorbar(pm1, cax=cax1)
    cb1.set_label('gC m⁻² day⁻¹', fontsize=9)
    cb1.ax.tick_params(labelsize=8)

    # ---- row 2 left: stats box ----
    ax_box.set_axis_off()
    stats_text = (
        f'{"":10s}{"min":>8s}  {"mean":>8s}  {"max":>8s}\n'
        f'{"─" * 42}\n'
        f'{"Truth":10s}{np.nanmin(truth_ll):>8.3f}  {np.nanmean(truth_ll):>8.3f}  {np.nanmax(truth_ll):>8.3f}\n'
        f'{"Prior":10s}{np.nanmin(prior_ll):>8.3f}  {np.nanmean(prior_ll):>8.3f}  {np.nanmax(prior_ll):>8.3f}\n'
        f'{"Optimized":10s}{np.nanmin(opt_ll):>8.3f}  {np.nanmean(opt_ll):>8.3f}  {np.nanmax(opt_ll):>8.3f}\n'
        f'\n'
        f'Units: gC m⁻² day⁻¹'
    )
    ax_box.text(0.5, 0.5, stats_text,
                transform=ax_box.transAxes,
                ha='center', va='center', fontsize=10,
                fontfamily='monospace',
                bbox=dict(boxstyle='round,pad=0.7', fc='#f7f7f7', ec='0.6', lw=1.2))

    # ---- row 2: difference panels ----
    pm2 = _draw_panel(ax_diff_prior, lat, lon, diff_prior_ll,
                      -dif_max, dif_max, 'RdBu_r',
                      f'Prior − Truth — {month_str}')
    _draw_panel(ax_diff_opt, lat, lon, diff_opt_ll,
                -dif_max, dif_max, 'RdBu_r',
                f'Optimized − Truth — {month_str}')

    cb2 = fig.colorbar(pm2, cax=cax2)
    cb2.set_label('gC m⁻² day⁻¹', fontsize=9)
    cb2.ax.tick_params(labelsize=8)

    fig.suptitle(
        f'TER Flux Comparison — prior scale ×{args.prior_scale:.1f}'
        f'  |  {month_str}',
        fontsize=13, fontweight='bold', y=0.96)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {args.out}')


if __name__ == '__main__':
    main()
