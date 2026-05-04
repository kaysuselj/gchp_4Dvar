#!/usr/bin/env python3
"""
Compute xCO2 model–observation biases from a GCHP sat-track file.

Usage:
    python ./co2_biases.py <gchp_file> <output_file>

Arguments:
    gchp_file   : GCHP netCDF file (GEOSChem.sat_track.*) with a time dimension
    output_file : Output netCDF file; written variables: xco2_obs, xco2_hat,
                  lat, lon, time

Environment setup:
    python -m venv gchp-env
    source gchp-env/bin/activate
    pip install numpy pandas xarray netCDF4 metpy

    The following file must be in the same directory:
        co2_sat_compare_monthly.py
"""

import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr
from metpy.interpolate import interpolate_1d
from metpy.units import units

from co2_sat_compare_monthly import read_oco_monthly

FOLD_OBS = '/nobackupp17/jliu7/OCO2/BAKER/B11.2_L2#/'


def co2_biases(gchp_file, output_file):
    ds_gchp = xr.open_dataset(gchp_file)
    ds_gchp['time'] = pd.to_datetime(ds_gchp['time'].values).round('s')

    # Find all unique (year, month) combinations present in the GCHP file
    times = pd.DatetimeIndex(ds_gchp['time'].values)
    year_months = sorted(set(zip(times.year, times.month)))

    xco2_obs_list = []
    xco2_hat_list = []
    lat_list      = []
    lon_list      = []
    time_list     = []

    for year, month in year_months:
        print(f'Processing {year}-{month:02d}')

        # Select GCHP profiles for this month
        mask = (times.year == year) & (times.month == month)
        ds_mod = ds_gchp.isel(time=np.where(mask)[0])

        co2_mod = ds_mod['SpeciesConcVV_CO2'].transpose('time', 'lev')
        prs_mod = ds_mod['Met_PMIDDRY'].transpose('time', 'lev')

        # Read matching OCO-2 observations
        ds_obs = read_oco_monthly(year, month, FOLD_OBS)
        if ds_obs is None:
            print(f'  No OCO-2 data for {year}-{month:02d}, skipping')
            continue

        prs_mod = prs_mod * units.hPa
        co2_mod = co2_mod * units.K  # fake unit required by interpolate_1d
        prs_obs = ds_obs['pressure'] * units.hPa

        nobs, nlev = prs_obs.shape
        co2_interp = np.empty((nobs, nlev), dtype=np.float32)
        for i in range(nobs):
            co2_interp[i] = interpolate_1d(prs_obs[i], prs_mod[i], co2_mod[i])
            if np.isnan(co2_interp[i, 0]):
                co2_interp[i, 0] = co2_mod[i, 0]

        # Apply OCO-2 observation operator: x_hat = x_a + A (x_m - x_a)
        co2_pert = co2_interp - ds_obs['CO2-apriori']
        xco2_hat = ds_obs['xCO2-apriori'] + np.sum(
            ds_obs['xCO2-averagingKernel'] * co2_pert, axis=1
        )

        xco2_obs_list.append(ds_obs['xCO2'].values)
        xco2_hat_list.append(np.array(xco2_hat))
        lat_list.append(ds_obs['latitude'].values)
        lon_list.append(ds_obs['longitude'].values)
        time_list.append(ds_obs['time'].values)

    if not xco2_obs_list:
        print('No data processed — output file not written.')
        return

    xco2_obs_all = np.concatenate(xco2_obs_list) * 1e6  # ppm
    xco2_hat_all = np.concatenate(xco2_hat_list) * 1e6  # ppm
    lat_all      = np.concatenate(lat_list)
    lon_all      = np.concatenate(lon_list)
    time_all     = np.concatenate(time_list)

    ds_out = xr.Dataset(
        {
            'xco2_obs': ('time', xco2_obs_all, {'long_name': 'Observed xCO2', 'units': 'ppm'}),
            'xco2_hat': ('time', xco2_hat_all, {'long_name': 'Modelled xCO2 (obs operator applied)', 'units': 'ppm'}),
            'lat':      ('time', lat_all,       {'long_name': 'Latitude',  'units': 'degrees_north'}),
            'lon':      ('time', lon_all,       {'long_name': 'Longitude', 'units': 'degrees_east'}),
        },
        coords={'time': time_all},
    )

    encoding = {
        'time': {
            'units': 'hours since 1900-01-01 00:00:00',
            'calendar': 'proleptic_gregorian',
            'dtype': 'float64',
        }
    }
    ds_out.to_netcdf(output_file, encoding=encoding)
    print(f'Written {len(time_all)} profiles to {output_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gchp_file',   help='GCHP sat-track netCDF file')
    parser.add_argument('output_file', help='Output netCDF file')
    args = parser.parse_args()

    co2_biases(args.gchp_file, args.output_file)
