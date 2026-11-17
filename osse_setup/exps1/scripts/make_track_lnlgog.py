#!/usr/bin/env python3
"""exps1 satellite track file with all three OCO-2 sounding types of J. Liu's B10 "50 km buddy" set: land nadir (LN),
land glint (LG) and ocean glint (OG). Same format as gchp_4Dvar/convert_satelite_tracks_osse.py (latitude, longitude
float32; time float64, hours since 1900-01-01), which reads LG + LN only. Soundings sorted by time.
Source per day: the pseudo set (PSEUDO/TRENDY/ORCHIDEE-ECCO2/OCO2b10_29Aug2020_*), which starts 2016-01-01; before that
the real set (OCO2-B10/OCO2b10_29Aug2020_*). Both have identical positions, times, averaging kernels, uncertainties and
priors (checked 2016-07-15 and 2017-03-15); only xCO2 differs, which the OSSE does not use (2026-11-10).
    make_track_lnlgog.py --out <track_file.nc> --start YYYY-MM-DD --end YYYY-MM-DD [--subdirs LN,LG,OG]
"""
import argparse, os
from datetime import datetime, timedelta
import numpy as np, xarray as xr

B10 = '/nobackupp17/jliu7/OCO2/L2#/OCO2-B10/'
PSEUDO = B10 + 'PSEUDO/TRENDY/ORCHIDEE-ECCO2/'
SUB = {'LN': 'OCO2b10_29Aug2020_LN', 'LG': 'OCO2b10_29Aug2020_LG', 'OG': 'OCO2b10_29Aug2020_OG'}


def day_file(sub, day):
    for root in (PSEUDO, B10):
        f = os.path.join(root, sub, f'{day:%Y/%m/%d}.nc')
        if os.path.exists(f):
            return f, ('pseudo' if root == PSEUDO else 'real')
    return None, None


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--start', default='2016-01-01')
    ap.add_argument('--end', default='2017-01-01', help='exclusive')
    ap.add_argument('--subdirs', default='LN,LG,OG')
    a = ap.parse_args()
    subs = [SUB[s] for s in a.subdirs.split(',')]
    t0, t1 = datetime.fromisoformat(a.start), datetime.fromisoformat(a.end)
    lat, lon, tim = [], [], []
    count = {s: 0 for s in subs}; src_days = {'pseudo': 0, 'real': 0}
    day = t0
    while day < t1:
        for s in subs:
            f, src = day_file(s, day)
            if f is None:
                continue
            ds = xr.open_dataset(f)
            n = ds['latitude'].size
            lat.append(ds['latitude'].values); lon.append(ds['longitude'].values % 360)
            tim.append(np.array([day + timedelta(hours=float(h)) for h in ds['time'].values], 'datetime64[ns]'))
            count[s] += n; src_days[src] += 1
        day += timedelta(days=1)
    lat, lon, tim = map(np.concatenate, (lat, lon, tim))
    o = np.argsort(tim, kind='stable')
    ds = xr.Dataset({'latitude': (['time'], lat[o].astype('float32')), 'longitude': (['time'], lon[o].astype('float32'))},
                    coords={'time': tim[o]})
    enc = {'latitude': {'dtype': 'float32', '_FillValue': np.float32(np.nan)},
           'longitude': {'dtype': 'float32', '_FillValue': np.float32(np.nan)},
           'time': {'dtype': 'float64', '_FillValue': np.float64(np.nan), 'units': 'hours since 1900-01-01',
                    'calendar': 'proleptic_gregorian'}}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    ds.to_netcdf(a.out, encoding=enc)
    print(f'{a.out}: {ds.time.size} soundings {a.start} .. {a.end} (' + ', '.join(f'{k[-2:]} {v}' for k, v in count.items())
          + f'); daily files from the pseudo set {src_days["pseudo"]}, from the real set {src_days["real"]}')
