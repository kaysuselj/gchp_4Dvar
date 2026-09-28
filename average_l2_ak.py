#!/usr/bin/env python3
"""
average_l2_ak.py

Compute the mean OCO-2 averaging kernel, pressure levels, CO2 a priori, and
observation uncertainty over all ORCHIDEE-ECCO2 L2# soundings in a given time
window.  The output file is used by make_semi_obs_xco2_simplified.py and
co2_adjoint_forcing_osse.py (--simplified-ak-file) to run the 4D-Var without
per-sounding L2# AKs.

Usage:
    python3 average_l2_ak.py \\
        --t-start 2016-01-01 \\
        --t-end   2016-04-01 \\
        --out     mean_ak_Jan_Mar_2016.nc

Output variables (dimension lev_sat = 20):
    mean_ak         mean averaging kernel               [dimensionless]
    mean_prs        mean satellite pressure levels       [hPa]
    mean_co2_apr    mean CO2 a priori profile            [mol/mol]
    mean_xco2_apr   mean XCO2 a priori (column scalar)  [mol/mol]
    mean_xco2_std   mean XCO2 uncertainty               [ppm]
"""

import argparse
import calendar
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

# ORCHIDEE-ECCO2 L2# data location
BASE_FOLD_OBS = '/nobackupp17/jliu7/OCO2/L2#/OCO2-B10/PSEUDO/TRENDY/ORCHIDEE-ECCO2/'

# Add the script's own directory to sys.path so convert_satelite_tracks_osse is found
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from convert_satelite_tracks_osse import read_oco_daily_osse


def average_l2_ak(t_start, t_end, base_fold, out_path):
    cur = t_start
    acc_ak    = []   # list of (n_soundings, 20) arrays
    acc_prs   = []
    acc_co2apr = []
    acc_xco2apr = []
    acc_xco2std = []

    while cur < t_end:
        year  = cur.year
        month = cur.month
        n_days = calendar.monthrange(year, month)[1]
        n_ok = 0
        for day in range(1, n_days + 1):
            date_str = f'{year}{month:02d}{day:02d}'
            try:
                ds = read_oco_daily_osse(date_str, base_fold)
            except Exception as e:
                print(f'  Warning: {date_str}: {e}')
                continue
            if ds is None:
                continue

            # Filter to [t_start, t_end)
            times = pd.DatetimeIndex(pd.to_datetime(ds['time'].values))
            mask  = (times >= t_start) & (times < t_end)
            if not np.any(mask):
                continue
            ds = ds.isel(time=np.where(mask)[0])

            prs  = ds['pressure'].values                 # (n, 20)
            ak   = ds['xCO2-averagingKernel'].values     # (n, 20)
            co2a = ds['CO2-apriori'].values              # (n, 20)
            xa   = ds['xCO2-apriori'].values             # (n,)
            std  = ds['xCO2-uncertainty'].values         # (n,)

            # Skip soundings where pressure is all zero (bad retrieval)
            valid = np.any(prs > 0, axis=1)
            if not np.any(valid):
                continue

            acc_ak.append(ak[valid])
            acc_prs.append(prs[valid])
            acc_co2apr.append(co2a[valid])
            acc_xco2apr.append(xa[valid])
            acc_xco2std.append(std[valid])
            n_ok += np.sum(valid)

        print(f'{year}-{month:02d}: {n_ok} soundings')
        cur += pd.DateOffset(months=1)

    if not acc_ak:
        print('ERROR: no soundings found', file=sys.stderr)
        sys.exit(1)

    ak_all     = np.vstack(acc_ak)       # (N, 20)
    prs_all    = np.vstack(acc_prs)
    co2apr_all = np.vstack(acc_co2apr)
    xco2apr_all = np.concatenate(acc_xco2apr)
    xco2std_all = np.concatenate(acc_xco2std)
    N = len(ak_all)
    print(f'Total soundings: {N}')

    # For pressure: many soundings have some zero-padded levels at the top.
    # Compute the mean only over levels where pressure > 0; levels that are
    # always zero will have mean zero (above model top — not used in interp).
    with np.errstate(invalid='ignore'):
        prs_mean = np.where(prs_all > 0, prs_all, np.nan).mean(axis=0)
        prs_mean = np.nan_to_num(prs_mean, nan=0.0)

    ds_out = xr.Dataset({
        'mean_ak': (
            ['lev_sat'],
            ak_all.mean(axis=0).astype(np.float64),
            {'long_name': 'Mean OCO-2 averaging kernel (over all soundings in window)',
             'units': '1'}),
        'mean_prs': (
            ['lev_sat'],
            prs_mean.astype(np.float64),
            {'long_name': 'Mean satellite pressure levels (over valid soundings)',
             'units': 'hPa',
             'comment': 'Zero entries indicate levels above the model top (unused)'}),
        'mean_co2_apr': (
            ['lev_sat'],
            co2apr_all.mean(axis=0).astype(np.float64),
            {'long_name': 'Mean CO2 a priori profile',
             'units': 'mol mol-1'}),
        'mean_xco2_apr': (
            [],
            float(xco2apr_all.mean()),
            {'long_name': 'Mean XCO2 a priori (column scalar)',
             'units': 'mol mol-1'}),
        'mean_xco2_std': (
            [],
            float(xco2std_all.mean()) * 1e6,  # mol/mol → ppm
            {'long_name': 'Mean XCO2 observation uncertainty',
             'units': 'ppm'}),
        'n_soundings': (
            [],
            N,
            {'long_name': 'Number of soundings used to compute the mean'}),
    })
    ds_out.attrs['description'] = (
        f'Mean ORCHIDEE-ECCO2 OCO-2 L2# AK and a priori over '
        f'{t_start.date()} to {t_end.date()}. '
        f'Created by average_l2_ak.py.')
    ds_out.attrs['t_start'] = str(t_start.date())
    ds_out.attrs['t_end']   = str(t_end.date())

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    ds_out.to_netcdf(out_path)
    print(f'Wrote {out_path}')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--t-start', default='2016-01-01',
                    help='Start date inclusive YYYY-MM-DD')
    ap.add_argument('--t-end',   default='2016-04-01',
                    help='End date exclusive YYYY-MM-DD')
    ap.add_argument('--l2-dir',  default=BASE_FOLD_OBS,
                    help='ORCHIDEE-ECCO2 root directory')
    ap.add_argument('--out',     required=True,
                    help='Output netCDF file path')
    args = ap.parse_args()

    average_l2_ak(
        t_start   = pd.Timestamp(args.t_start),
        t_end     = pd.Timestamp(args.t_end),
        base_fold = args.l2_dir,
        out_path  = args.out,
    )


if __name__ == '__main__':
    main()
