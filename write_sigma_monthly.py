#!/usr/bin/env python3
"""
Write monthly sigma files for HEMCO/ExtData from the current 4D-Var state.

Twelve single-record files are written, one per month:

    <sigma-dir>/<YYYY>/<MM>.nc4        e.g.  sigma_monthly/2016/01.nc4

This mirrors the OCEAN_ECCO_DARWIN layout (%y4/%m2.nc), which ExtData reads
correctly with refresh template F%y4-%m2-01T00:00:00.  A single file holding
twelve time records is NOT used: ExtData's handling of a multi-record file
under a monthly refresh template is untested in this configuration.

Control vector layout is (nmon, nlat, nlon) flattened C-order, so element
(m, j, i) lives at index m*nlat*nlon + j*nlon + i.  Month index 0 is January.

First call (no state file): sigma = 1.0 everywhere (prior = unperturbed TER).
Subsequent calls: sigma is read from the x_next field in the state file.

Sigma multiplies the TER_CLASS_CTEM base emission only (HEMCO scale factor 750),
not total CO2 emissions.

The number of monthly fields is --n-months (default 12, the full-year run).
For a shorter window (e.g. a 3-month shakedown) pass --n-months 3; the files
are written for --start-month, --start-month+1, ... with year rollover, so the
control-vector length is n-months * nlat * nlon.

One EXTRA "boundary" file is written for the month after the last control
month (start-month + n-months).  The forward/adjoint window ends at that month
boundary, and while ExtData integrates the final control month it needs that
file as the RIGHT time-bracket -- without it the forward aborts at the last
month boundary in ExtData (ExtDataGridCompMod.F90 right-bracket file search).
Sigma_CO2 uses an F-persisted refresh (no interpolation), so the boundary
file's values never enter the solution; it duplicates the last control month
and is NOT part of the control vector.

Usage:
    python write_sigma_monthly.py --sigma-dir /path/sigma_monthly \
                                 [--state-file /path/4dvar_state.npz] \
                                 [--nlat 46] [--nlon 72] \
                                 [--year 2016] [--n-months 12] [--start-month 1]

Test override (bypasses the state file entirely):
    --const-per-month 1.0,0.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0
"""
import argparse
import os
import numpy as np
import pandas as pd
import xarray as xr

NMON = 12          # default control-vector length (full-year run)


def make_latlon_grid(nlat, nlon):
    lats = np.linspace(-90.0, 90.0, nlat)
    lons = np.linspace(-180.0, 180.0 - 360.0 / nlon, nlon)
    return lats, lons


def write_month(sigma_2d, lats, lons, path, t_stamp):
    ds = xr.Dataset(
        {'Sigma_CO2': (
            ['time', 'lat', 'lon'],
            sigma_2d[np.newaxis].astype(np.float32),
            {'long_name': 'CO2 TER flux scaling factor',
             'units':     '1',
             'comment':   '4D-Var control variable; multiplies EmisCO2_NetTerrExch'},
        )},
        coords={
            'time': [t_stamp],
            'lat':  ('lat', lats.astype(np.float64),
                     {'long_name': 'Latitude',  'units': 'degrees_north'}),
            'lon':  ('lon', lons.astype(np.float64),
                     {'long_name': 'Longitude', 'units': 'degrees_east'}),
        },
    )
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ds.to_netcdf(path, encoding={
        'time':      {'units': 'hours since 1900-01-01 00:00:00',
                      'calendar': 'proleptic_gregorian', 'dtype': 'float64'},
        'Sigma_CO2': {'dtype': 'float32', 'zlib': True, 'complevel': 4},
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sigma-dir',  required=True)
    parser.add_argument('--state-file', default=None)
    parser.add_argument('--nlat',       type=int, default=46)
    parser.add_argument('--nlon',       type=int, default=72)
    parser.add_argument('--year',       type=int, default=2016)
    parser.add_argument('--n-months',   type=int, default=NMON,
                        help='Number of monthly fields / control-vector months '
                             '(default 12, the full-year run).')
    parser.add_argument('--start-month', type=int, default=1,
                        help='Calendar month (1-12) of the first field.')
    parser.add_argument('--const-per-month', default=None,
                        help='Comma-separated n-months values; bypasses state '
                             'file. For plumbing tests only.')
    args = parser.parse_args()

    nlat, nlon = args.nlat, args.nlon
    nmon = args.n_months
    n = nmon * nlat * nlon
    lats, lons = make_latlon_grid(nlat, nlon)

    if args.const_per_month:
        vals = [float(v) for v in args.const_per_month.split(',')]
        if len(vals) != nmon:
            raise SystemExit(f'--const-per-month needs {nmon} values, got {len(vals)}')
        sigma = np.empty((nmon, nlat, nlon))
        for m, v in enumerate(vals):
            sigma[m] = v
        it = 0
        print(f'TEST MODE: constant sigma per month = {vals}')
    elif args.state_file and os.path.exists(args.state_file):
        state = np.load(args.state_file)
        x = state['x_next']
        if x.size != n:
            raise SystemExit(f'state x_next has {x.size} elements, expected {n} '
                             f'({nmon}x{nlat}x{nlon})')
        sigma = x.reshape(nmon, nlat, nlon)
        it = int(state['iteration'])
    else:
        sigma = np.ones((nmon, nlat, nlon))
        it = 1

    # Month index m maps to calendar (start_month + m), rolling over into the
    # next year(s) as needed, so a non-January or year-spanning window works.
    base = args.year * 12 + (args.start_month - 1)
    for m in range(nmon):
        tot = base + m
        yy, mm = tot // 12, tot % 12 + 1
        t_stamp = pd.Timestamp(f'{yy}-{mm:02d}-01')
        path = os.path.join(args.sigma_dir, f'{yy}', f'{mm:02d}.nc4')
        write_month(sigma[m], lats, lons, path, t_stamp)
        print(f'[iter {it}] {yy}-{mm:02d}: min={sigma[m].min():.4f} '
              f'max={sigma[m].max():.4f} mean={sigma[m].mean():.4f} -> {path}')

    # Boundary file at month (start_month + nmon): the window ends at this month
    # boundary and ExtData needs the file as the RIGHT time-bracket while it
    # integrates the final control month.  Sigma_CO2 is F-persisted (no
    # interpolation), so these values never enter the solution -- the file only
    # has to exist -- so duplicate the last control month.  Not a control var.
    tot = base + nmon
    yy, mm = tot // 12, tot % 12 + 1
    t_stamp = pd.Timestamp(f'{yy}-{mm:02d}-01')
    path = os.path.join(args.sigma_dir, f'{yy}', f'{mm:02d}.nc4')
    write_month(sigma[nmon - 1], lats, lons, path, t_stamp)
    print(f'[iter {it}] {yy}-{mm:02d}: boundary (=last control month) -> {path}')


if __name__ == '__main__':
    main()
