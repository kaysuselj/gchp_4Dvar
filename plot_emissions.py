#!/usr/bin/env python3
"""Compare truth vs optimized EmisCO2_Total on the native C24 cubed-sphere.

Produces per-iteration PNG files saved to --plot-dir:
  emiss_truth_jan_iter001.png   — truth map for each month
  emiss_opt_jan_iter001.png     — optimized map for each month
  emiss_diff_jan_iter001.png    — (optimized - truth) for each month
  emiss_global_mean_iter001.png — bar chart of spatial-mean flux per month

Usage (called from 4dvar_optimizer.osse.*.run after archive_emissions):
  python3 plot_emissions.py
      --truth    obs/truth_emissions_monthly_2016Q1.nc
      --opt      4dvar_output_osse/emissions_monthly_iter001.nc4
      --plot-dir plots_osse/20261001_123456
      --iter-label 001
      [--coord-file <any GEOSChem GCHP nc4 with lats/lons>]
      [--search-dir <dir to scan for lats/lons if --coord-file absent>]
"""
import argparse
import glob
import os
import sys

import numpy as np

# Import CS helpers from plot_4dvar (SCRIPT_DIR is on PYTHONPATH in the run script)
try:
    from plot_4dvar import load_cs_coords, plot_map_cs, to_faces, NF
except ImportError as e:
    sys.exit(f"Cannot import plot_4dvar: {e}")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

MONTH_NAMES = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# kgCO2/m2/s  →  gC/m2/day
_KG_CO2_TO_GC_DAY = 1e3 * (12.0 / 44.0) * 86400.0


def _load_emiss(path, cs_res):
    """Return EmisCO2_Total as (ntime, 6, cs_res, cs_res) in gC/m2/day.

    EmisCO2_Total shape from GCHP HISTORY: (time, lev, nf, Ydim, Xdim).
    Sum over lev to get total column-integrated flux (surface + aviation).
    """
    ds = xr.open_dataset(path)
    da = ds['EmisCO2_Total']
    if 'lev' in da.dims:
        da = da.sum('lev')
    # Now (time, nf, Ydim, Xdim)
    arr = da.values * _KG_CO2_TO_GC_DAY
    times = pd.DatetimeIndex(ds['time'].values)
    return arr, times


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--truth',       required=True,
                    help='Monthly-mean truth EmisCO2_Total nc file')
    ap.add_argument('--opt',         required=True,
                    help='Monthly-mean optimized EmisCO2_Total nc file for this iter')
    ap.add_argument('--plot-dir',    required=True)
    ap.add_argument('--iter-label',  required=True, help='Zero-padded iteration number')
    ap.add_argument('--coord-file',  default=None,
                    help='Any GEOSChem GCHP nc4 carrying "lats"/"lons" data variables')
    ap.add_argument('--search-dir',  action='append', default=[],
                    help='Directory to scan for a GCHP file with lats/lons (repeatable)')
    ap.add_argument('--cs-res',      type=int, default=24)
    args = ap.parse_args()

    # ---- existence checks ----
    if not os.path.exists(args.truth):
        print(f"Truth emissions not found: {args.truth} — skipping")
        return 0
    if not os.path.exists(args.opt):
        print(f"Optimized emissions not found: {args.opt} — skipping")
        return 0

    os.makedirs(args.plot_dir, exist_ok=True)
    cs_res = args.cs_res

    # ---- CS lat/lon grid ----
    try:
        lats, lons = load_cs_coords(cs_res,
                                    coord_file=args.coord_file,
                                    search_dirs=args.search_dir)
    except RuntimeError as e:
        print(f"WARNING: {e}\n  — emissions maps skipped; global-mean plot will still be made")
        lats = lons = None

    # ---- load fields ----
    truth_arr, truth_times = _load_emiss(args.truth, cs_res)
    opt_arr,   opt_times   = _load_emiss(args.opt,   cs_res)

    truth_months = truth_times.month
    opt_months   = opt_times.month
    common_months = sorted(set(truth_months) & set(opt_months))

    if not common_months:
        print("No overlapping months between truth and optimized files — nothing to plot")
        return 0

    # ---- per-month CS maps ----
    if lats is not None:
        for m in common_months:
            t_idx = int(np.where(truth_months == m)[0][0])
            o_idx = int(np.where(opt_months   == m)[0][0])
            mname = MONTH_NAMES[m - 1]
            il    = args.iter_label

            t_field = to_faces(truth_arr[t_idx], cs_res)
            o_field = to_faces(opt_arr[o_idx],   cs_res)
            diff    = o_field - t_field

            vabs = float(max(np.abs(t_field).max(), np.abs(o_field).max(), 1e-12))
            dabs = float(max(np.abs(diff).max(), 1e-12))
            cbar = 'gC m⁻² day⁻¹'

            plot_map_cs(t_field, lats, lons,
                        f'Truth EmisCO2_Total — {mname} 2016',
                        os.path.join(args.plot_dir, f'emiss_truth_{mname.lower()}_iter{il}.png'),
                        vmin=-vabs, vmax=vabs, cbar_label=cbar)
            plot_map_cs(o_field, lats, lons,
                        f'Optimized EmisCO2_Total — {mname} 2016  (iter {il})',
                        os.path.join(args.plot_dir, f'emiss_opt_{mname.lower()}_iter{il}.png'),
                        vmin=-vabs, vmax=vabs, cbar_label=cbar)
            plot_map_cs(diff, lats, lons,
                        f'Optimized − Truth EmisCO2_Total — {mname} 2016  (iter {il})',
                        os.path.join(args.plot_dir, f'emiss_diff_{mname.lower()}_iter{il}.png'),
                        vmin=-dabs, vmax=dabs, cbar_label=cbar)

    # ---- global spatial-mean bar chart ----
    truth_means = [float(truth_arr[np.where(truth_months == m)[0][0]].mean())
                   for m in common_months]
    opt_means   = [float(opt_arr[np.where(opt_months == m)[0][0]].mean())
                   for m in common_months]
    mnames = [MONTH_NAMES[m - 1] for m in common_months]

    x = np.arange(len(common_months))
    fig, ax = plt.subplots(figsize=(max(5, len(common_months) * 2 + 1), 4))
    ax.bar(x - 0.2, truth_means, width=0.4, label='Truth',
           color='steelblue', edgecolor='k', linewidth=0.5)
    ax.bar(x + 0.2, opt_means,   width=0.4, label=f'Optimized (iter {args.iter_label})',
           color='darkorange', edgecolor='k', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(mnames)
    ax.set_ylabel('Global mean EmisCO2_Total (gC m⁻² day⁻¹)')
    ax.set_title(f'Global mean CO₂ flux: truth vs optimized  (iter {args.iter_label})')
    ax.axhline(0, color='k', linewidth=0.6, linestyle='--')
    ax.legend()
    ax.grid(axis='y', alpha=0.4)
    plt.tight_layout()
    path = os.path.join(args.plot_dir, f'emiss_global_mean_iter{args.iter_label}.png')
    plt.savefig(path, dpi=120, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {path}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
