#!/usr/bin/env python3
"""
make_semi_obs_xco2_simplified.py

Apply a fixed (space-time mean) OCO-2 averaging kernel to GCHP sat_track output
to produce truth xCO2 semi-observations for the simplified OSSE.

Unlike make_semi_obs_xco2.py, this script does NOT require per-sounding L2# data.
The mean AK (from average_l2_ak.py) is applied to every sat_track profile regardless
of whether an ORCHIDEE-ECCO2 sounding exists at that position.

The AK variables are embedded in the output netCDF so the forcing script needs only
this one file (no separate --mean-ak-file at run time).

Usage:
    python3 make_semi_obs_xco2_simplified.py \\
        --sat-track-dir <path/to/OutputDir/> \\
        --mean-ak-file  mean_ak_Jan_Mar_2016.nc \\
        --t-start       2016-01-01 \\
        --t-end         2016-04-01 \\
        --out           <path/to/semi_obs_xco2_6h_simplified_2016Q1.nc>

Output variables:
    xco2_truth(time)        [ppm]   truth XCO2 with mean AK applied
    xco2_uncertainty(time)  [ppm]   mean uncertainty (same for all obs)
    latitude(time)          [deg]
    longitude(time)         [deg]
    time(time)
    mean_ak(lev_sat)        [1]     embedded mean averaging kernel
    mean_prs(lev_sat)       [hPa]   embedded mean pressure levels
    mean_co2_apr(lev_sat)   [mol/mol] embedded mean CO2 a priori
    mean_xco2_apr           [mol/mol] embedded mean XCO2 a priori scalar
"""

import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd
import xarray as xr


# MAPL writes sat_track files with the PREVIOUS month's timestamp.
# This map lets us find the right file for each data month.
def _sat_track_pattern(sat_track_dir, year, month):
    """Return glob pattern for the sat_track file containing data for year/month.

    MAPL writes two forms:
      canonical:  GEOSChem.sat_track.YYYYMMDD_HHMMz.nc4  (e.g. 20160101_0030z)
      short-form: GEOSChem.sat_track.YYYYMM_HHMMz.nc4    (e.g. 201601_0030z) -- duplicate

    Include '01' after the month to anchor to the canonical YYYYMMDD form and
    avoid matching both files, which would double-count the observations.
    """
    prev = pd.Timestamp(year=year, month=month, day=1) - pd.DateOffset(months=1)
    return os.path.join(sat_track_dir,
                        f'GEOSChem.sat_track.{prev.year:04d}{prev.month:02d}01*.nc4')


def load_mean_ak(ak_file):
    """Load mean AK from average_l2_ak.py output. Returns dict."""
    ds = xr.open_dataset(ak_file)
    out = {
        'mean_ak':       ds['mean_ak'].values.astype(np.float64),        # (20,)
        'mean_prs':      ds['mean_prs'].values.astype(np.float64),       # (20,) hPa
        'mean_co2_apr':  ds['mean_co2_apr'].values.astype(np.float64),   # (20,) mol/mol
        'mean_xco2_apr': float(ds['mean_xco2_apr'].values),              # mol/mol
        'mean_xco2_std': float(ds['mean_xco2_std'].values),              # ppm
    }
    ds.close()
    valid = out['mean_prs'] > 0
    if not np.any(valid):
        raise ValueError(f'mean_prs in {ak_file} has no positive values')
    print(f'Mean AK loaded from {ak_file}  '
          f'({valid.sum()} valid pressure levels, '
          f'mean_xco2_std={out["mean_xco2_std"]:.3f} ppm)')
    return out


def apply_mean_ak(co2_gchp_vv, prs_gchp_hpa, ak):
    """
    Interpolate GCHP CO2 to mean pressure levels and apply mean AK.

    co2_gchp_vv  : (nlev,)  mol/mol, surface-first (descending pressure OK)
    prs_gchp_hpa : (nlev,)  hPa
    ak           : dict from load_mean_ak

    Returns truth xCO2 in ppm, or NaN on failure.
    """
    mean_prs     = ak['mean_prs']
    mean_ak      = ak['mean_ak']
    mean_co2_apr = ak['mean_co2_apr']
    mean_xco2_apr = ak['mean_xco2_apr']

    valid = mean_prs > 0
    if valid.sum() < 2:
        return np.nan

    sort_idx = np.argsort(prs_gchp_hpa)
    try:
        co2_interp = np.interp(
            mean_prs[valid],
            prs_gchp_hpa[sort_idx],
            co2_gchp_vv[sort_idx],
        )
    except Exception:
        return np.nan

    diff = co2_interp - mean_co2_apr[valid]
    xco2 = mean_xco2_apr + float(np.sum(mean_ak[valid] * diff))
    return xco2 * 1e6   # mol/mol → ppm


def process_month(sat_track_dir, year, month, t_start, t_end, ak):
    """Process one calendar month of sat_track files. Returns list of dicts."""
    pattern = _sat_track_pattern(sat_track_dir, year, month)
    files   = sorted(glob.glob(pattern))
    if not files:
        print(f'  WARNING: no sat_track file for {year}-{month:02d} '
              f'(pattern: {os.path.basename(pattern)})')
        return []

    results = []
    n_nan = 0
    for fpath in files:
        ds = xr.open_dataset(fpath)
        times = pd.DatetimeIndex(pd.to_datetime(ds['time'].values).round('s'))
        mask  = (times >= t_start) & (times < t_end)
        if not np.any(mask):
            ds.close()
            continue

        idx_arr = np.where(mask)[0]
        co2_all = ds['SpeciesConcVV_CO2'].values   # (lev, time)
        prs_all = ds['Met_PMIDDRY'].values          # (lev, time)
        lat_all = ds['latitude'].values             # (time,)
        lon_all = ds['longitude'].values            # (time,)
        ds.close()

        for i in idx_arr:
            xco2 = apply_mean_ak(co2_all[:, i], prs_all[:, i], ak)
            if np.isnan(xco2):
                n_nan += 1
                continue
            results.append(dict(
                time      = times[i].to_datetime64(),
                latitude  = float(lat_all[i]),
                longitude = float(lon_all[i]),
                xco2_truth       = float(xco2),
                xco2_uncertainty = ak['mean_xco2_std'],
            ))

    print(f'  {year}-{month:02d}: {len(results)} obs'
          + (f'  ({n_nan} NaN skipped)' if n_nan else ''))
    return results


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--sat-track-dir', required=True,
                    help='Directory containing GEOSChem.sat_track.*.nc4 from truth run')
    ap.add_argument('--mean-ak-file',  required=True,
                    help='mean AK netCDF from average_l2_ak.py')
    ap.add_argument('--t-start', default='2016-01-01',
                    help='Start date inclusive YYYY-MM-DD')
    ap.add_argument('--t-end',   default='2016-04-01',
                    help='End date exclusive YYYY-MM-DD')
    ap.add_argument('--out',     required=True,
                    help='Output netCDF path')
    args = ap.parse_args()

    t_start = pd.Timestamp(args.t_start)
    t_end   = pd.Timestamp(args.t_end)

    ak = load_mean_ak(args.mean_ak_file)

    all_results = []
    cur = t_start
    while cur < t_end:
        print(f'Processing {cur.year}-{cur.month:02d}...')
        all_results += process_month(
            args.sat_track_dir, cur.year, cur.month, t_start, t_end, ak)
        cur += pd.DateOffset(months=1)

    if not all_results:
        print('ERROR: no obs produced — check sat_track_dir and t_start/t_end',
              file=sys.stderr)
        sys.exit(1)

    times = np.array([r['time'] for r in all_results])
    order = np.argsort(times)
    times = times[order]

    def col(key, dtype=np.float32):
        return np.array([all_results[i][key] for i in order], dtype=dtype)

    # Load mean AK dataset for embedding
    ds_ak = xr.open_dataset(args.mean_ak_file)

    ds_out = xr.Dataset(
        {
            'xco2_truth': (
                ['time'], col('xco2_truth'),
                {'units': 'ppm',
                 'long_name': 'Truth XCO2 from GCHP default run + mean OCO-2 AK'}),
            'xco2_uncertainty': (
                ['time'], col('xco2_uncertainty'),
                {'units': 'ppm',
                 'long_name': 'Mean OCO-2 XCO2 observation uncertainty (fixed for all obs)'}),
            'latitude': (
                ['time'], col('latitude'),
                {'units': 'degrees_north'}),
            'longitude': (
                ['time'], col('longitude'),
                {'units': 'degrees_east'}),
            # Embed the AK so the forcing script needs only this file
            'mean_ak': (
                ['lev_sat'], ds_ak['mean_ak'].values.astype(np.float32),
                {'units': '1',
                 'long_name': 'Mean OCO-2 averaging kernel used to compute xco2_truth'}),
            'mean_prs': (
                ['lev_sat'], ds_ak['mean_prs'].values.astype(np.float32),
                {'units': 'hPa',
                 'long_name': 'Mean satellite pressure levels',
                 'comment': 'Zero entries indicate levels above model top (unused)'}),
            'mean_co2_apr': (
                ['lev_sat'], ds_ak['mean_co2_apr'].values.astype(np.float32),
                {'units': 'mol mol-1',
                 'long_name': 'Mean CO2 a priori profile'}),
            'mean_xco2_apr': (
                [], float(ds_ak['mean_xco2_apr'].values),
                {'units': 'mol mol-1',
                 'long_name': 'Mean XCO2 a priori (column scalar)'}),
        },
        coords={'time': times},
    )
    ds_ak.close()

    ds_out.attrs['description'] = (
        'OSSE simplified semi-observations: GCHP default CO2 + mean OCO-2 AK '
        f'(averaged over {t_start.date()} to {t_end.date()}). '
        'AK is fixed (not per-sounding) and embedded in this file.')
    ds_out.attrs['t_start']     = str(t_start.date())
    ds_out.attrs['t_end']       = str(t_end.date())
    ds_out.attrs['mean_ak_src'] = os.path.abspath(args.mean_ak_file)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    ds_out.to_netcdf(args.out)
    print(f'Wrote {len(times)} obs → {args.out}')


if __name__ == '__main__':
    main()
