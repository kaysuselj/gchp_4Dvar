#!/usr/bin/env python3
"""
Write the OSSE prior-error field E for the 4D-Var: prior TER = truth TER * E.

The truth run uses the unmodified TER; the 4D-Var forward AND adjoint multiply TER by
E (HEMCO scale factor 752, read from ../prior_err_co2.nc4 in the experiment work tree,
exactly like the control variable sigma = factor 750).  The optimizer should then find
sigma = 1/E (posterior flux sigma*E*TER = truth TER).

E is written on the native C{im} cubed sphere in the GMAO stacked layout
(lat = 6*im, lon = im) that MAPL reads without regridding -- the same layout as the
sigma file written by write_sigma.py --grid cs.

Field types (--type):
  constant   E = --value everywhere                          (e.g. 1.3)
  zonal      latitude bands, --bands="lat0:lat1:val,..."  (use "=": values may start with -)
  regional   boxes, --boxes="lat0:lat1:lon0:lon1:val,..."  (lon in -180..180 or 0..360;
             later boxes overwrite earlier ones; outside: --background)
  random     correlated Gaussian field: --value (mean) + --amplitude (std) with a
             Gaussian correlation in great-circle distance, --length-km, --seed;
             clipped to [--min, --max]

The recipe (command line) is stored in the file attributes.

usage examples:
  make_prior_err.py --out prior_err_co2.nc4 --type constant --value 1.3
  make_prior_err.py --out e.nc4 --type zonal --bands="-30:30:1.3,30:90:0.8"
  make_prior_err.py --out e.nc4 --type random --value 1.0 --amplitude 0.3 --length-km 1500 --seed 1
"""
import argparse
import os
import sys

import warnings

import numpy as np
import xarray as xr

warnings.filterwarnings("ignore", message="Duplicate dimension names")
NF = 6
R_EARTH_KM = 6371.0
DEFAULT_COORD = ('/nobackupp27/ksuselj1/gchp_14.5.3_adjoint_surfaceF/tests_final/'
                 'P2_advection/P2a_edge_corner_8ms/results/coord.nc4')


def read_coords(path, im):
    with xr.open_dataset(path, drop_variables=['anchor']) as ds:
        lats = np.asarray(ds['lats'].values, dtype=np.float64).reshape(NF, im, im)
        lons = np.asarray(ds['lons'].values, dtype=np.float64).reshape(NF, im, im)
    return lats, (lons + 180.0) % 360.0 - 180.0          # lon in [-180, 180)


def parse_list(spec, n):
    out = []
    for item in spec.split(','):
        v = [float(x) for x in item.split(':')]
        if len(v) != n:
            raise SystemExit(f'bad item "{item}": need {n} numbers separated by ":"')
        out.append(v)
    return out


def in_lon_range(lon, lo0, lo1):
    lo0 = (lo0 + 180.0) % 360.0 - 180.0
    lo1 = (lo1 + 180.0) % 360.0 - 180.0
    if lo0 <= lo1:
        return (lon >= lo0) & (lon < lo1)
    return (lon >= lo0) | (lon < lo1)                    # box across the date line


def random_field(lats, lons, mean, amp, length_km, seed):
    phi, lam = np.radians(lats.ravel()), np.radians(lons.ravel())
    xyz = np.stack([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)], 1)
    ang = np.arccos(np.clip(xyz @ xyz.T, -1.0, 1.0))      # great-circle angle
    d = ang * R_EARTH_KM
    k = np.exp(-0.5 * (d / length_km) ** 2)
    noise = np.random.default_rng(seed).standard_normal(len(phi))
    f = k @ noise
    f = (f - f.mean()) / f.std()
    return (mean + amp * f).reshape(lats.shape)


def write_cs(field, im, path, attrs, t_stamp='2016-01-01T00:00:00'):
    stacked = np.asarray(field, dtype=np.float64).reshape(NF * im, im)
    ds = xr.Dataset(
        {'Prior_Err': (['time', 'lat', 'lon'], stacked[None, :, :],
                       {'long_name': 'OSSE prior error factor E (prior TER = truth TER * E)',
                        'units': '1'})},
        coords={'time': ('time', [np.datetime64(t_stamp, 'ns')], {'long_name': 'time'}),
                'lat': ('lat', np.arange(NF * im, dtype=np.float64),
                        {'long_name': 'Latitude', 'units': 'degrees_north'}),
                'lon': ('lon', np.arange(im, dtype=np.float64),
                        {'long_name': 'Longitude', 'units': 'degrees_east'})},
        attrs=attrs)
    tmp = path + '.tmp.nc4'
    ds.to_netcdf(tmp, format='NETCDF4', encoding={'Prior_Err': {'_FillValue': None},
                                                  'lat': {'_FillValue': None},
                                                  'lon': {'_FillValue': None}})
    os.replace(tmp, path)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', default=None, help='single (annual) E file')
    p.add_argument('--out-dir', default=None,
                   help='monthly E: <out-dir>/<YYYY>/<MM>.nc4 for each month of the window (+1 boundary month, '
                        'a copy of the last, for the ExtData right bracket), E_m = field x month factor m')
    p.add_argument('--month-factors', default='',
                   help='with --out-dir: comma-separated factor per month (e.g. "1.3,1.0,0.8" with --type constant '
                        '--value 1); default 1 for every month of --n-months')
    p.add_argument('--year', type=int, default=2016)
    p.add_argument('--start-month', type=int, default=1)
    p.add_argument('--n-months', type=int, default=0)
    p.add_argument('--type', required=True, choices=['constant', 'zonal', 'regional', 'random'])
    p.add_argument('--value', type=float, default=1.0, help='constant value / mean of random')
    p.add_argument('--background', type=float, default=1.0, help='value outside bands/boxes')
    p.add_argument('--bands', default='')
    p.add_argument('--boxes', default='')
    p.add_argument('--amplitude', type=float, default=0.3)
    p.add_argument('--length-km', type=float, default=1500.0)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--min', type=float, default=0.4,
                   help='lower clip of E (default 0.4: target sigma 1/E <= 2.5 < upper bound 3)')
    p.add_argument('--max', type=float, default=3.0)
    p.add_argument('--cs-res', type=int, default=24)
    p.add_argument('--coord-file', default=DEFAULT_COORD,
                   help='GCHP C-grid file with lats/lons (6,im,im) cell centres')
    a = p.parse_args()

    lats, lons = read_coords(a.coord_file, a.cs_res)
    if a.type == 'constant':
        e = np.full(lats.shape, a.value)
    elif a.type == 'zonal':
        e = np.full(lats.shape, a.background)
        for la0, la1, v in parse_list(a.bands, 3):
            e[(lats >= la0) & (lats < la1)] = v
    elif a.type == 'regional':
        e = np.full(lats.shape, a.background)
        for la0, la1, lo0, lo1, v in parse_list(a.boxes, 5):
            e[(lats >= la0) & (lats < la1) & in_lon_range(lons, lo0, lo1)] = v
    else:
        e = random_field(lats, lons, a.value, a.amplitude, a.length_km, a.seed)
    e = np.clip(e, a.min, a.max)
    if not np.all(np.isfinite(e)) or e.min() <= 0:
        sys.exit('E must be finite and > 0')

    attrs = {'recipe': ' '.join(sys.argv[1:]), 'type': a.type,
             'target_sigma': 'sigma_true = 1/E', 'coord_file': a.coord_file}
    if a.out_dir:                                  # monthly E (MODE=monthly), same layout as write_sigma_monthly.py
        fac = [float(x) for x in a.month_factors.split(',')] if a.month_factors else [1.0] * a.n_months
        if a.n_months and len(fac) != a.n_months:
            sys.exit(f'--month-factors has {len(fac)} values, --n-months is {a.n_months}')
        base = a.year * 12 + (a.start_month - 1)
        for m in range(len(fac) + 1):              # + boundary month (copy of the last)
            tot = base + m; yy, mm = tot // 12, tot % 12 + 1
            em = np.clip(e * fac[min(m, len(fac) - 1)], a.min, a.max)
            path = os.path.join(a.out_dir, f'{yy}', f'{mm:02d}.nc4')
            os.makedirs(os.path.dirname(path), exist_ok=True)
            write_cs(em, a.cs_res, path, dict(attrs, month_factor=fac[min(m, len(fac) - 1)]), f'{yy}-{mm:02d}-01T00:00:00')
            print(f'{path}: E min={em.min():.4f} max={em.max():.4f} mean={em.mean():.4f}'
                  + ('   (boundary month = copy of the last)' if m == len(fac) else ''))
        return
    if not a.out:
        sys.exit('give --out or --out-dir')
    write_cs(e, a.cs_res, a.out, attrs)
    print(f'{a.out}: E type={a.type}  min={e.min():.4f} max={e.max():.4f} mean={e.mean():.4f}  '
          f'(target sigma = 1/E: {1/e.max():.4f} .. {1/e.min():.4f})')


if __name__ == '__main__':
    main()
