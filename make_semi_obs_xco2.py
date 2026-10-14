#!/usr/bin/env python3
"""
make_semi_obs_xco2.py

Apply OCO-2 L2# averaging kernel to GCHP sat_track output to produce
truth xCO2 semi-observations for the OSSE.

For each L2# obs point, find the matching GCHP sat_track profile
(matched by time/lat/lon), interpolate GCHP CO2 to the L2# pressure
grid, then apply the column-averaging kernel:

    xco2_truth = xco2_apr + sum( AK * (co2_gchp_interp - co2_apr) )

Output netCDF contains xco2_truth alongside xco2_apriori, xco2_uncertainty,
and xco2_l2 (the original ORCHIDEE-ECCO2 pseudo-obs, kept for diagnostics).

Usage:
    python3 make_semi_obs_xco2.py \\
        --sat-track-dir  <path/to/forward_run/OutputDir/> \\
        --l2-dir         <path/to/L2#/ORCHIDEE-ECCO2/> \\
        --t-start        2016-01-01 \\
        --t-end          2016-04-01 \\
        --out            <path/to/semi_obs_xco2.nc>
"""

import argparse
import calendar
import glob
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
from scipy.interpolate import interp1d

L2_SUBFOLDERS = ['OCO2b10_29Aug2020_LG', 'OCO2b10_29Aug2020_LN']

TOL_TIME_S = 60    # max time mismatch (seconds)
TOL_LAT    = 0.1   # max lat mismatch (degrees)
TOL_LON    = 0.1   # max lon mismatch (degrees)


def load_sat_track_month(sat_track_dir, year, month):
    """Load sat_track file for the given calendar month.

    GCHP names the file with the *previous* month's timestamp, so
    January 2016 data lives in GEOSChem.sat_track.20151201_*.nc4.
    """
    prev = pd.Timestamp(year=year, month=month, day=1) - pd.DateOffset(months=1)
    pattern = os.path.join(sat_track_dir,
                           f'GEOSChem.sat_track.{prev.year:04d}{prev.month:02d}*.nc4')
    files = glob.glob(pattern)
    if not files:
        print(f'  WARNING: no sat_track file for {year}-{month:02d} '
              f'(pattern: {pattern})')
        return None
    ds = xr.open_dataset(files[0])
    co2  = ds['SpeciesConcVV_CO2'].values          # (lev, ntime)
    prs  = ds['Met_PMIDDRY'].values                # (lev, ntime)  hPa
    lats = ds['latitude'].values                   # (ntime,)
    lons = ds['longitude'].values % 360            # (ntime,)
    times = pd.DatetimeIndex(
        pd.to_datetime(ds['time'].values).round('s'))
    return dict(co2=co2, prs=prs, lats=lats, lons=lons, times=times)


def load_l2_day(l2_dir, year, month, day):
    """Load L2# ORCHIDEE-ECCO2 obs for one day (LG + LN), or None."""
    date_str  = f'{year:04d}/{month:02d}/{day:02d}'
    base_date = datetime(year, month, day)
    datasets  = []
    for subfold in L2_SUBFOLDERS:
        fpath = os.path.join(l2_dir, subfold, date_str + '.nc')
        if not os.path.exists(fpath):
            continue
        ds = xr.open_dataset(fpath)
        time_hrs = ds['time'].values
        time_dt  = np.array([base_date + timedelta(hours=float(h))
                              for h in time_hrs])
        ds = ds.assign(time=('nSamples', time_dt))
        ds['longitude'] = ds['longitude'] % 360
        datasets.append(ds)
    if not datasets:
        return None
    combined = xr.concat(datasets, dim='nSamples').sortby('time')
    return combined


def apply_ak(co2_gchp_vv, prs_gchp_hpa,
             prs_obs, ak, co2_apr_ppm, xco2_apr_ppm):
    """
    Interpolate GCHP CO2 profile to L2# pressure levels and apply AK.

    co2_gchp_vv  : (nlev,)  mol/mol
    prs_gchp_hpa : (nlev,)  hPa (can be surface-first or top-first)
    prs_obs      : (20,)    hPa (L2# levels, may contain zeros)
    ak           : (20,)    dimensionless averaging kernel
    co2_apr_ppm  : (20,)    ppm
    xco2_apr_ppm : float    ppm

    Returns truth xCO2 in ppm, or NaN on failure.
    """
    co2_ppm = co2_gchp_vv * 1e6          # mol/mol -> ppm

    # Ensure pressure is monotonically increasing (surface last is fine
    # for interp1d as long as we sort)
    order = np.argsort(prs_gchp_hpa)
    prs_sorted = prs_gchp_hpa[order]
    co2_sorted = co2_ppm[order]

    valid = prs_obs > 0
    if valid.sum() < 2:
        return np.nan

    try:
        f = interp1d(prs_sorted, co2_sorted,
                     bounds_error=False,
                     fill_value=(co2_sorted[0], co2_sorted[-1]))
        co2_interp = f(prs_obs)
    except Exception:
        return np.nan

    diff = co2_interp - co2_apr_ppm
    diff[~valid] = 0.0
    return float(xco2_apr_ppm) + float(np.sum(ak * diff))


def process_month(year, month, sat_track_dir, l2_dir):
    """Return list of matched-obs dicts for one calendar month."""
    st = load_sat_track_month(sat_track_dir, year, month)
    if st is None:
        return []

    num_days = calendar.monthrange(year, month)[1]
    results  = []

    for day in range(1, num_days + 1):
        l2 = load_l2_day(l2_dir, year, month, day)
        if l2 is None:
            continue

        t_arr   = l2['time'].values
        lat_arr = l2['latitude'].values
        lon_arr = l2['longitude'].values % 360
        n_obs   = len(t_arr)

        t_st_s  = st['times'].astype('int64') / 1e9   # epoch seconds

        for j in range(n_obs):
            t_obs   = pd.Timestamp(t_arr[j]).round('s')
            lat_obs = float(lat_arr[j])
            lon_obs = float(lon_arr[j])

            t_obs_s = t_obs.timestamp()
            dt_s    = np.abs(t_st_s - t_obs_s)
            mask    = (dt_s            <= TOL_TIME_S) & \
                      (np.abs(st['lats'] - lat_obs) <= TOL_LAT) & \
                      (np.abs(st['lons'] - lon_obs) <= TOL_LON)
            idx = np.where(mask)[0]
            if len(idx) == 0:
                continue
            k = idx[np.argmin(dt_s[idx])]

            prs_obs  = l2['pressure'].values[j]               # (20,) hPa
            ak       = l2['xCO2-averagingKernel'].values[j]   # (20,)
            co2_apr  = l2['CO2-apriori'].values[j]            # (20,) ppm
            xco2_apr = float(l2['xCO2-apriori'].values[j])    # ppm

            xco2_truth = apply_ak(
                st['co2'][:, k], st['prs'][:, k],
                prs_obs, ak, co2_apr, xco2_apr)

            if np.isnan(xco2_truth):
                continue

            results.append(dict(
                time             = t_obs.to_datetime64(),
                latitude         = lat_obs,
                longitude        = lon_obs,
                xco2_truth       = xco2_truth,
                xco2_apriori     = xco2_apr,
                xco2_uncertainty = float(l2['xCO2-uncertainty'].values[j]),
                xco2_l2          = float(l2['xCO2'].values[j]),
            ))

    print(f'  {year}-{month:02d}: {len(results)} matched obs')
    return results


def main():
    ap = argparse.ArgumentParser(
        description='Apply OCO-2 AK to GCHP sat_track to produce truth xCO2')
    ap.add_argument('--sat-track-dir', required=True,
                    help='Directory containing GEOSChem.sat_track.*.nc4')
    ap.add_argument('--l2-dir', required=True,
                    help='L2# ORCHIDEE-ECCO2 root directory')
    ap.add_argument('--t-start', default='2016-01-01',
                    help='Start date (inclusive), YYYY-MM-DD')
    ap.add_argument('--t-end',   default='2016-04-01',
                    help='End date (exclusive), YYYY-MM-DD')
    ap.add_argument('--out', required=True,
                    help='Output netCDF file path')
    args = ap.parse_args()

    t_start = pd.Timestamp(args.t_start)
    t_end   = pd.Timestamp(args.t_end)

    all_results = []
    cur = t_start
    while cur < t_end:
        print(f'Processing {cur.year}-{cur.month:02d}...')
        all_results += process_month(cur.year, cur.month,
                                     args.sat_track_dir, args.l2_dir)
        cur += pd.DateOffset(months=1)

    if not all_results:
        print('ERROR: no matched obs found', file=sys.stderr)
        sys.exit(1)

    times = np.array([r['time'] for r in all_results])
    order = np.argsort(times)
    times = times[order]

    def col(key, dtype=np.float32):
        return np.array([all_results[i][key] for i in order], dtype=dtype)

    ds_out = xr.Dataset({
        'xco2_truth':        ('time', col('xco2_truth')),
        'xco2_apriori':      ('time', col('xco2_apriori')),
        'xco2_uncertainty':  ('time', col('xco2_uncertainty')),
        'xco2_l2':           ('time', col('xco2_l2')),
        'latitude':          ('time', col('latitude')),
        'longitude':         ('time', col('longitude')),
    }, coords={'time': times})

    ds_out['xco2_truth'].attrs.update(
        units='ppm', long_name='Truth xCO2 from GCHP default run + OCO-2 AK')
    ds_out['xco2_l2'].attrs.update(
        units='ppm', long_name='ORCHIDEE-ECCO2 L2# pseudo-obs xCO2')
    ds_out['xco2_apriori'].attrs.update(
        units='ppm', long_name='OCO-2 L2# xCO2 a priori')
    ds_out['xco2_uncertainty'].attrs.update(
        units='ppm', long_name='OCO-2 L2# xCO2 observation uncertainty')
    ds_out.attrs['description'] = (
        'OSSE semi-observations: GCHP default CO2 + OCO-2 averaging kernel. '
        f'Period {args.t_start} to {args.t_end}.')

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    ds_out.to_netcdf(args.out)
    print(f'Wrote {len(times)} obs -> {args.out}')


if __name__ == '__main__':
    main()
