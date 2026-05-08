#!/usr/bin/env python3
"""
Write sigma_co2.nc4 for HEMCO from the current 4D-Var state.

First call (no state file): sigma = 1.0 everywhere (prior = unperturbed emissions).
Subsequent calls: sigma is read from the x_next field in the state file.

Usage:
    python write_sigma.py --sigma-file /path/sigma_co2.nc4 \
                         [--state-file /path/4dvar_state.npz] \
                         [--nlat 46] [--nlon 72] \
                         [--t-start 2019-01-01]
"""
import argparse
import numpy as np
import pandas as pd
import xarray as xr


def make_latlon_grid(nlat, nlon):
    lats = np.linspace(-90.0, 90.0, nlat)
    lons = np.linspace(-180.0, 180.0 - 360.0 / nlon, nlon)
    return lats, lons


def write_sigma_file(sigma_2d, lats, lons, path, t_start):
    ds = xr.Dataset(
        {'Sigma_CO2': (
            ['time', 'lat', 'lon'],
            sigma_2d[np.newaxis].astype(np.float32),
            {'long_name': 'CO2 surface flux scaling factor',
             'units':     '1',
             'comment':   '4D-Var control variable; multiplies EmisCO2_Total'},
        )},
        coords={
            'time': [t_start],
            'lat':  ('lat', lats.astype(np.float64),
                     {'long_name': 'Latitude',  'units': 'degrees_north'}),
            'lon':  ('lon', lons.astype(np.float64),
                     {'long_name': 'Longitude', 'units': 'degrees_east'}),
        },
    )
    ds.to_netcdf(path, encoding={
        'time':      {'units': 'hours since 1900-01-01 00:00:00',
                      'calendar': 'proleptic_gregorian', 'dtype': 'float64'},
        'Sigma_CO2': {'dtype': 'float32', 'zlib': True, 'complevel': 4},
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sigma-file',  required=True,
                        help='Output path for sigma_co2.nc4')
    parser.add_argument('--state-file',  default=None,
                        help='4D-Var state file (4dvar_state.npz); '
                             'if absent sigma=1 everywhere')
    parser.add_argument('--nlat',        type=int, default=46)
    parser.add_argument('--nlon',        type=int, default=72)
    parser.add_argument('--t-start',     default='2019-01-01',
                        help='Simulation start date (YYYY-MM-DD)')
    args = parser.parse_args()

    import os
    lats, lons = make_latlon_grid(args.nlat, args.nlon)

    if args.state_file and os.path.exists(args.state_file):
        state    = np.load(args.state_file)
        sigma_1d = state['x_next']
        it       = int(state['iteration'])
    else:
        sigma_1d = np.ones(args.nlat * args.nlon, dtype=np.float64)
        it       = 1

    sigma_2d = sigma_1d.reshape(args.nlat, args.nlon)
    t_start  = pd.Timestamp(args.t_start)

    write_sigma_file(sigma_2d, lats, lons, args.sigma_file, t_start)
    print(f'[iter {it}] sigma written: '
          f'min={sigma_2d.min():.4f}  max={sigma_2d.max():.4f}  '
          f'mean={sigma_2d.mean():.4f}')


if __name__ == '__main__':
    main()
