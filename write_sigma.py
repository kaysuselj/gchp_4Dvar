#!/usr/bin/env python3
"""
Write sigma_co2.nc4 for HEMCO from the current 4D-Var state.

First call (no state file): sigma = 1.0 everywhere (prior = unperturbed emissions).
Subsequent calls: sigma is read from the x_next field in the state file.

Two control-variable grids are supported:

  --grid latlon   sigma on a regular lat-lon grid (nlat x nlon); MAPL regrids
                  it to the model cubed sphere at read time.  n = nlat*nlon.
  --grid cs       sigma on the model's native C{im} cubed sphere, written in the
                  GMAO "stacked" layout (lat = 6*im, lon = im) that MAPL reads
                  with NO regrid.  n = 6*im*im.  This removes the lat-lon<->CS
                  interpolation mismatch in the gradient.

Usage:
    python write_sigma.py --sigma-file /path/sigma_co2.nc4 \
                         [--state-file /path/4dvar_state.npz] \
                         [--grid latlon|cs] \
                         [--nlat 46] [--nlon 72] \
                         [--cs-res 24] \
                         [--t-start 2019-01-01]
"""
import argparse
import os
import numpy as np
import pandas as pd
import xarray as xr

NF = 6   # cubed-sphere faces


# ---------------------------------------------------------------------------
# lat-lon
# ---------------------------------------------------------------------------

def make_latlon_grid(nlat, nlon):
    lats = np.linspace(-90.0, 90.0, nlat)
    lons = np.linspace(-180.0, 180.0 - 360.0 / nlon, nlon)
    return lats, lons


def write_sigma_latlon(sigma_2d, lats, lons, path, t_start):
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


# ---------------------------------------------------------------------------
# native cubed sphere (stacked layout: lat = 6*im, lon = im)
# ---------------------------------------------------------------------------

def write_sigma_cs(sigma_faces, im, path, t_stamp):
    """Write one single-record sigma file in the stacked CS layout.

    MAPL's GridManager identifies the cubed sphere from lat == 6*lon and reads
    the file with an identity regrid.  Matches adjoint_validation/cs_fd_tools.py
    (write_month_cs) so the sigma written here and the gradient read back use the
    same (6, im, im) C-order face convention.
    """
    stacked = np.asarray(sigma_faces, dtype=np.float64).reshape(NF * im, im)
    ds = xr.Dataset(
        {'Sigma_CO2': (['time', 'lat', 'lon'], stacked[None, :, :],
                       {'long_name': 'CO2 surface flux scaling factor (native CS)',
                        'units': '1',
                        'comment': '4D-Var control variable on the model C-grid'})},
        coords={
            'time': ('time', [t_stamp], {'long_name': 'time'}),
            # fake index coordinates: MAPL computes the geometry itself
            'lat': ('lat', np.arange(NF * im, dtype=np.float64),
                    {'long_name': 'Latitude', 'units': 'degrees_north'}),
            'lon': ('lon', np.arange(im, dtype=np.float64),
                    {'long_name': 'Longitude', 'units': 'degrees_east'}),
        })
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + '.tmp.nc4'
    ds.to_netcdf(tmp, format='NETCDF4',
                 encoding={'Sigma_CO2': {'_FillValue': None},
                           'lat': {'_FillValue': None},
                           'lon': {'_FillValue': None}})
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sigma-file',  required=True,
                        help='Output path for sigma_co2.nc4')
    parser.add_argument('--state-file',  default=None,
                        help='4D-Var state file (4dvar_state.npz); '
                             'if absent sigma=1 everywhere')
    parser.add_argument('--grid', choices=['latlon', 'cs'], default='latlon',
                        help='Control-variable grid (default latlon)')
    parser.add_argument('--nlat',        type=int, default=46)
    parser.add_argument('--nlon',        type=int, default=72)
    parser.add_argument('--cs-res',      type=int, default=24,
                        help='Cubed-sphere face size im (C{im}) for --grid cs')
    parser.add_argument('--t-start',     default='2019-01-01',
                        help='Simulation start date (YYYY-MM-DD)')
    args = parser.parse_args()

    t_start = pd.Timestamp(args.t_start)

    if args.grid == 'cs':
        im = args.cs_res
        n = NF * im * im
    else:
        n = args.nlat * args.nlon

    if args.state_file and os.path.exists(args.state_file):
        state    = np.load(args.state_file)
        sigma_1d = state['x_next']
        it       = int(state['iteration'])
        if sigma_1d.size != n:
            raise SystemExit(
                f'state x_next size {sigma_1d.size} != expected {n} for '
                f'--grid {args.grid}')
    else:
        sigma_1d = np.ones(n, dtype=np.float64)
        it       = 1

    if args.grid == 'cs':
        sigma_faces = sigma_1d.reshape(NF, im, im)
        write_sigma_cs(sigma_faces, im, args.sigma_file, t_start)
        s = sigma_faces
        print(f'[iter {it}] CS sigma written (C{im}): '
              f'min={s.min():.4f}  max={s.max():.4f}  mean={s.mean():.4f}')
    else:
        lats, lons = make_latlon_grid(args.nlat, args.nlon)
        sigma_2d = sigma_1d.reshape(args.nlat, args.nlon)
        write_sigma_latlon(sigma_2d, lats, lons, args.sigma_file, t_start)
        print(f'[iter {it}] lat-lon sigma written: '
              f'min={sigma_2d.min():.4f}  max={sigma_2d.max():.4f}  '
              f'mean={sigma_2d.mean():.4f}')


if __name__ == '__main__':
    main()
