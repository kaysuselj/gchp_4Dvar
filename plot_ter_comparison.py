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


def cell_areas(coord_file, im):
    """Cell areas (unit sphere) of a C{im} grid from its corner lat/lon."""
    d = nc.Dataset(coord_file)
    clat = np.radians(np.array(d['corner_lats'][:], dtype=np.float64)).reshape(NF, im + 1, im + 1)
    clon = np.radians(np.array(d['corner_lons'][:], dtype=np.float64)).reshape(NF, im + 1, im + 1)
    d.close()
    v = np.stack([np.cos(clat) * np.cos(clon), np.cos(clat) * np.sin(clon), np.sin(clat)], -1)

    def tri(a, b, c):   # spherical triangle area (Van Oosterom & Strackee)
        num = np.abs(np.einsum('...i,...i', a, np.cross(b, c)))
        den = 1 + np.einsum('...i,...i', a, b) + np.einsum('...i,...i', b, c) + np.einsum('...i,...i', c, a)
        return 2 * np.arctan2(num, den)
    p00, p10, p11, p01 = v[:, :-1, :-1], v[:, 1:, :-1], v[:, 1:, 1:], v[:, :-1, 1:]
    return tri(p00, p10, p11) + tri(p00, p11, p01)


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
                       shading='auto', rasterized=True, **_p4d._tr())
    if _p4d._is_geo(ax):                                   # Robinson (plot_4dvar._PROJ)
        ax.coastlines(resolution='110m', linewidth=0.4, color='0.25', zorder=2)
        ax.set_global()
        gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='0.5', alpha=0.3,
                          xlocs=range(-180, 181, 60), ylocs=range(-90, 91, 30))
        gl.top_labels = gl.right_labels = False
        gl.xlabel_style = gl.ylabel_style = {'size': 7}
    else:
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
    ap.add_argument('--prior-scale', type=float, default=1.3,
                    help='constant prior error E (used if no --prior-err-file)')
    ap.add_argument('--prior-err-file', default=None,
                    help='prior error field E (make_prior_err.py, CS stacked): '
                         'prior TER = truth TER * E, target sigma = 1/E')
    ap.add_argument('--prior-err-dir', default=None,
                    help='monthly control: E per month in <dir>/<YYYY>/<MM>.nc4 (make_prior_err.py --out-dir)')
    ap.add_argument('--coord-file', default=('/nobackupp27/ksuselj1/gchp_14.5.3_adjoint_surfaceF/'
                    'tests_final/P2_advection/P2a_edge_corner_8ms/results/coord.nc4'),
                    help='C-grid file with corner_lats/corner_lons (cell areas for the metric)')
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
    # annual states keep the best iterate (x_best, it_best); monthly states (lbfgsb_step_monthly.py) only
    # the current one (x_next, iteration; 4dvar_state_best.npz holds the best one in the same keys)
    x_best   = st['x_best'] if 'x_best' in st.files else st['x_next']
    it_best  = int(st['it_best']) if 'it_best' in st.files else int(st['iteration'])
    iter_tag = args.iter_label if args.iter_label else str(it_best)
    ncell    = NF * cs_res * cs_res
    if x_best.size != ncell:                                   # monthly control: one figure per month
        nmon = x_best.size // ncell
        if nmon * ncell != x_best.size or nmon != len(keep_idx):
            raise SystemExit(f'state has {x_best.size} values: not C{cs_res} x {len(keep_idx)} months')
        sig_m = x_best.reshape(nmon, NF, cs_res, cs_res)
        tru_m = all_months[keep_idx] * _TO_GC_DAY
        e_m, rows = [], []
        for m, tm in enumerate(keep_times):
            if args.prior_err_dir:
                f = os.path.join(args.prior_err_dir, f'{tm.year}', f'{tm.month:02d}.nc4')
                with nc.Dataset(f) as de:
                    e_m.append(np.array(de['Prior_Err'][0]).reshape(NF, cs_res, cs_res))
            else:
                e_m.append(np.full((NF, cs_res, cs_res), args.prior_scale))
        base, ext = os.path.splitext(args.out)
        for m, tm in enumerate(keep_times):
            lab = tm.strftime('%b %Y')
            r = _figure(args, tru_m[m], e_m[m], sig_m[m], lats_cs, lons_cs, cs_res, lab, it_best, iter_tag,
                        f'monthly, E {lab} mean {e_m[m].mean():.3g}', f'{base}_{tm:%Y%m}{ext}')
            rows.append((lab,) + r)
        # whole window, weights over all months (TER x area of each month)
        w = np.abs(tru_m) * cell_areas(args.coord_file, cs_res)[None]
        w = w / w.sum()
        ep, eo = np.array(e_m), np.array(e_m) * sig_m
        for lab, pb, pr, ob, orr in rows:
            print(f'  {lab}: prior bias {pb:+.4f} rms {pr:.4f} | optimized bias {ob:+.4f} rms {orr:.4f}')
        print(f'TER x area weighted flux ratio error (all months): prior bias {float((w * (ep - 1)).sum()):+.4f} '
              f'rms {float(np.sqrt((w * (ep - 1) ** 2).sum())):.4f}  |  optimized (iter {it_best}) bias '
              f'{float((w * (eo - 1)).sum()):+.4f} rms {float(np.sqrt((w * (eo - 1) ** 2).sum())):.4f}')
        return
    sigma_cs = x_best.reshape(NF, cs_res, cs_res)

    # ---- prior error E (constant or field) ----
    if args.prior_err_file:
        de = nc.Dataset(args.prior_err_file)
        e_cs = np.array(de['Prior_Err'][0]).reshape(NF, cs_res, cs_res)
        e_label = getattr(de, 'recipe', os.path.basename(args.prior_err_file))
        # keep only the field definition (drop --out / --coord-file / --cs-res paths of the recipe)
        _toks, _keep, _skip = e_label.split(), [], False
        for _t in _toks:
            if _skip:
                _skip = False; continue
            if _t in ('--out', '--coord-file', '--cs-res'):
                _skip = True; continue
            _keep.append(_t)
        e_label = ' '.join(_keep).replace('--type ', '').replace('--', '') or e_label
        de.close()
    else:
        e_cs = np.full((NF, cs_res, cs_res), args.prior_scale)
        e_label = f'constant {args.prior_scale:g}'

    _figure(args, ter_truth_cs, e_cs, sigma_cs, lats_cs, lons_cs, cs_res, month_str, it_best, iter_tag, e_label,
            args.out, monthly=False)


def _figure(args, ter_truth_cs, e_cs, sigma_cs, lats_cs, lons_cs, cs_res, month_str, it_best, iter_tag, e_label, out,
            monthly=True):
    """truth / prior / optimized TER maps, differences and the weighted flux error for one period"""
    # ---- derive flux fields ----
    ter_prior_cs = ter_truth_cs * e_cs
    ter_opt_cs   = ter_truth_cs * e_cs * sigma_cs

    # ---- metric: TER x area weighted error of the flux ratio to truth (sigma*E - 1) ----
    w = np.abs(ter_truth_cs) * cell_areas(args.coord_file, cs_res)
    w = w / w.sum()
    def _werr(r):
        return float((w * (r - 1.0)).sum()), float(np.sqrt((w * (r - 1.0) ** 2).sum()))
    prior_bias, prior_rms = _werr(e_cs)
    opt_bias, opt_rms = _werr(e_cs * sigma_cs)
    if not monthly:
        print(f'TER x area weighted flux ratio error: prior bias {prior_bias:+.4f} rms {prior_rms:.4f}  |  '
              f'optimized (iter {it_best}) bias {opt_bias:+.4f} rms {opt_rms:.4f}')

    # ---- interpolate to regular grid ----
    lat, lon, truth_ll = cs_to_latlon(ter_truth_cs, lats_cs, lons_cs, args.res)
    _,   _,   prior_ll = cs_to_latlon(ter_prior_cs, lats_cs, lons_cs, args.res)
    _,   _,   opt_ll   = cs_to_latlon(ter_opt_cs,   lats_cs, lons_cs, args.res)

    diff_prior_ll = prior_ll - truth_ll
    diff_opt_ll   = opt_ll   - truth_ll

    # ---- colour limits ----
    abs_max = np.ceil(np.nanmax([truth_ll, prior_ll, opt_ll]) * 10) / 10
    dif_max = np.ceil(np.nanmax(np.abs([diff_prior_ll, diff_opt_ll])) * 10) / 10


    # ---- layout: 2 rows × 3 data cols + 1 narrow cbar col ----
    fig = plt.figure(figsize=(18, 9))
    gs = fig.add_gridspec(2, 4,
                          width_ratios=[1, 1, 1, 0.04],
                          hspace=0.32, wspace=0.08,
                          left=0.04, right=0.96, top=0.90, bottom=0.05)

    ax_truth      = fig.add_subplot(gs[0, 0], **_p4d._map_kw())
    ax_prior      = fig.add_subplot(gs[0, 1], **_p4d._map_kw())
    ax_opt        = fig.add_subplot(gs[0, 2], **_p4d._map_kw())
    cax1          = fig.add_subplot(gs[0, 3])
    ax_box        = fig.add_subplot(gs[1, 0])
    ax_diff_prior = fig.add_subplot(gs[1, 1], **_p4d._map_kw())
    ax_diff_opt   = fig.add_subplot(gs[1, 2], **_p4d._map_kw())
    cax2          = fig.add_subplot(gs[1, 3])

    # ---- row 1: absolute fluxes ----
    pm1 = _draw_panel(ax_truth, lat, lon, truth_ll,
                      0, abs_max, 'YlOrRd', f'Truth TER — {month_str}')
    _draw_panel(ax_prior, lat, lon, prior_ll,
                0, abs_max, 'YlOrRd',
                f'Prior TER (truth × E) — {month_str}')
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
        f'Units: gC m⁻² day⁻¹\n\n'
        f'TER×area weighted (σE − 1):\n'
        f'{"Prior":10s} bias {prior_bias:+.3f}  rms {prior_rms:.3f}\n'
        f'{"Optimized":10s} bias {opt_bias:+.3f}  rms {opt_rms:.3f}'
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
        f'TER Flux Comparison — prior error E: {e_label}'
        f'  |  {month_str}',
        fontsize=13, fontweight='bold', y=0.96)

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')
    return prior_bias, prior_rms, opt_bias, opt_rms

if __name__ == '__main__':
    main()
