#!/usr/bin/env python3
"""
Compute 3D adjoint forcing from OCO-2 xCO2 observations, binned to GCHP
chemistry checkpoint times.  One netCDF file is written per checkpoint time
that contains at least one observation.  Output can be on a regular lat/lon
grid or regridded to a GEOS-5 gnomonic cubed-sphere grid.

Usage:
    # Regular lat/lon output (default)
    python ./co2_adjoint_forcing.py gchp_file output_dir \\
        --ts-chem 1200 \\
        --t-start 2015-01-01T00:00:00 --t-end 2015-01-31T23:59:59 \\
        --nlat 91 --nlon 144

    # Cubed-sphere output (C90)
    python ./co2_adjoint_forcing.py gchp_file output_dir \\
        --ts-chem 1200 \\
        --t-start 2015-01-01T00:00:00 --t-end 2015-01-31T23:59:59 \\
        --grid cubedsphere --cs-res 90 \\
        --nlat 181 --nlon 360          # intermediate accumulation resolution

Arguments:
    gchp_file    : GCHP sat-track netCDF file (GEOSChem.sat_track.*)
    output_dir   : output directory (created if absent); one file per checkpoint

    --ts-chem    : chemistry timestep [s] for checkpoint binning (required)
    --t-start    : start of assimilation window YYYY-MM-DDTHH:MM:SS (required)
    --t-end      : end of assimilation window (required)
    --grid       : latlon (default) or cubedsphere
    --cs-res     : cubed-sphere resolution N for C{N} grid (required if cubedsphere)
    --nlat       : intermediate lat/lon accumulation grid latitude points (default 181)
    --nlon       : intermediate lat/lon accumulation grid longitude points (default 360)

Environment setup:
    python -m venv gchp-env
    source gchp-env/bin/activate
    pip install numpy pandas xarray netCDF4 metpy xesmf

    The following file must be in the same directory:
        co2_sat_compare_monthly.py

Output files:
    CO2_adjoint_forcing_YYYYMMDD_HHMMz.nc4  (one per checkpoint with obs)

Output variables:
    latlon output     : (time, lev, lat, lon)
    cubedsphere output: (time, lev, nf, Ydim, Xdim)  nf=6 faces

    forcing : adjoint forcing dJ/d(xCO2_model) [ppm^-1]
    n_obs   : number of observations per cell per checkpoint [count]

Notes:
    - Accumulation always occurs on the intermediate lat/lon grid (--nlat/--nlon).
      For cubed sphere, this is then regridded with xesmf bilinear interpolation.
    - Level 1 = surface (highest pressure), level LLPAR = TOA.
    - Unit conversion v/v -> kg/box must be applied in Fortran before adding
      to State_Chm%%SpeciesAdj.
    - The GEOS-5 gnomonic cubed sphere face orientation used here follows the
      Putman & Lin (2007) convention. Verify face 0 covers lon~[-45,45] at the
      equator against your GCHP output.
    - OCO-2 data is stored by month so read_oco_monthly always reads all days
      in the month; observations outside [t_start, t_end] are discarded in the
      per-observation loop.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr
from metpy.interpolate import interpolate_1d
from metpy.units import units

from co2_sat_compare_monthly import read_oco_monthly

FOLD_OBS = '/nobackupp17/jliu7/OCO2/BAKER/B11.2_L2#/'


# ---------------------------------------------------------------------------
# Vertical interpolation weight matrix  (mirrors GET_INTMAP in oco2_Baker.f)
# ---------------------------------------------------------------------------

def get_intmap(pres_model, pres_sat):
    """
    Linear interpolation weight matrix W (nmodel x nsat).

    Implements GET_INTMAP from oco2_Baker.f.

    Convention:
        Forward:  x_sat[j]    = sum_i W[i,j] * x_model[i]   (= W.T @ x_model)
        Adjoint:  x_model_adj = W @ x_sat_adj

    Parameters
    ----------
    pres_model : (nmodel,) decreasing pressure [hPa], surface first
    pres_sat   : (nsat,)   satellite pressure levels [hPa]
    """
    nmodel = len(pres_model)
    nsat   = len(pres_sat)
    W = np.zeros((nmodel, nsat), dtype=np.float64)

    for j in range(nsat):
        p = float(pres_sat[j])
        if p > pres_model[0]:       # below surface: pin to surface level
            W[0, j] = 1.0
            continue
        found = False
        for i in range(nmodel - 1):
            hi  = pres_model[i]
            low = pres_model[i + 1]
            if p <= hi and p > low:
                diff        = hi - low
                W[i,     j] = (p  - low) / diff
                W[i + 1, j] = (hi - p  ) / diff
                found = True
                break
        if not found:               # above model top: pin to top level
            W[nmodel - 1, j] = 1.0

    return W


# ---------------------------------------------------------------------------
# Per-observation adjoint forcing profile
# ---------------------------------------------------------------------------

def obs_forcing_profile(pres_model, pres_sat, xAK, xco2_hat, xco2_obs, xco2_std):
    """
    Adjoint forcing on model levels for one observation (mirrors
    CALC_GOSAT_CO2_FORCE in oco2_Baker.f).

        diff         = xco2_hat - xco2_obs              [ppm]
        force_scalar = diff / sigma^2                    [ppm^-1]
        force_sat    = xAK  * force_scalar               [ppm^-1, nsat]
        force_model  = W    @ force_sat                  [ppm^-1, nmodel]
    """
    W            = get_intmap(np.asarray(pres_model, dtype=np.float64),
                              np.asarray(pres_sat,   dtype=np.float64))
    diff         = float(xco2_hat) - float(xco2_obs)
    force_scalar = diff / float(xco2_std) ** 2
    force_sat    = np.asarray(xAK, dtype=np.float64) * force_scalar
    force_model  = W @ force_sat
    return force_model, diff


# ---------------------------------------------------------------------------
# Checkpoint time grid helpers
# ---------------------------------------------------------------------------

def make_checkpoint_grid(t_start, t_end, ts_chem_s):
    times, t = [], t_start
    while t <= t_end:
        times.append(t)
        t += pd.Timedelta(seconds=ts_chem_s)
    return pd.DatetimeIndex(times)


def nearest_checkpoint(obs_time, checkpoints, ts_chem_s):
    half    = pd.Timedelta(seconds=ts_chem_s / 2)
    deltas  = np.abs(checkpoints - obs_time)
    idx     = deltas.argmin()
    return idx if deltas[idx] <= half else -1


# ---------------------------------------------------------------------------
# GEOS-5 gnomonic cubed-sphere grid  (Putman & Lin 2007 convention)
# ---------------------------------------------------------------------------

def cubedsphere_cell_centers(cs_res):
    """
    Compute cell-centre lat/lon for a C{cs_res} GEOS-5 gnomonic cubed sphere.

    Returns
    -------
    lats, lons : each (6, cs_res, cs_res) in degrees
        Face ordering (0-indexed):
            0 → equatorial, centred at lon=0
            1 → equatorial, centred at lon=90E
            2 → equatorial, centred at lon=180
            3 → equatorial, centred at lon=270E
            4 → North pole face
            5 → South pole face
    """
    n     = cs_res
    alpha = np.linspace(-np.pi / 4, np.pi / 4, n, endpoint=False) + np.pi / (4 * n)
    A, B  = np.meshgrid(alpha, alpha)   # A varies along x (lon), B along y (lat)

    lats = np.zeros((6, n, n))
    lons = np.zeros((6, n, n))

    def _to_latlon(px, py, pz):
        r   = np.sqrt(px**2 + py**2 + pz**2)
        lat = np.rad2deg(np.arcsin(np.clip(pz / r, -1, 1)))
        lon = np.rad2deg(np.arctan2(py, px))
        return lat, lon

    # Equatorial faces: face normal points in +x, +y, -x, -y respectively.
    # Local s=A (east–west gnomonic angle), t=B (north–south gnomonic angle).
    lats[0], lons[0] = _to_latlon( np.ones_like(A),  np.tan(A),  np.tan(B))
    lats[1], lons[1] = _to_latlon(-np.tan(A),         np.ones_like(A),  np.tan(B))
    lats[2], lons[2] = _to_latlon(-np.ones_like(A), -np.tan(A),  np.tan(B))
    lats[3], lons[3] = _to_latlon( np.tan(A),        -np.ones_like(A),  np.tan(B))

    # Polar faces
    lats[4], lons[4] = _to_latlon(-np.tan(B),  np.tan(A),  np.ones_like(A))   # North
    lats[5], lons[5] = _to_latlon( np.tan(B),  np.tan(A), -np.ones_like(A))   # South

    return lats, lons


# ---------------------------------------------------------------------------
# Regrid lat/lon dataset to cubed sphere via xesmf
# ---------------------------------------------------------------------------

def regrid_latlon_to_cubedsphere(ds_ll, cs_res):
    """
    Regrid an (lat, lon) xarray Dataset to a C{cs_res} cubed sphere.

    Parameters
    ----------
    ds_ll  : xr.Dataset with coords lat (nlat,) and lon (nlon,)
    cs_res : int, cubed-sphere resolution N

    Returns
    -------
    xr.Dataset with spatial dimensions (nf, Ydim, Xdim)
    """
    try:
        import xesmf as xe
    except ImportError:
        raise ImportError(
            'xesmf is required for cubed-sphere output.\n'
            'Install it with:  pip install xesmf'
        )

    cs_lat, cs_lon = cubedsphere_cell_centers(cs_res)  # (6, n, n)

    # xesmf target grid: provide 2-D lat/lon arrays with arbitrary shape
    # Flatten the 6 faces into a single (6*n, n) layout so xesmf sees a
    # 2-D structured target; we reshape back afterwards.
    n = cs_res
    target = xr.Dataset({
        'lat': (['y', 'x'], cs_lat.reshape(6 * n, n)),
        'lon': (['y', 'x'], cs_lon.reshape(6 * n, n)),
    })

    regridder = xe.Regridder(ds_ll, target, method='bilinear',
                             periodic=True, ignore_degenerate=True)
    ds_cs_flat = regridder(ds_ll)

    # Restore face dimension: (6*n, n) → (6, n, n)
    def _reshape_var(da):
        new_shape = da.shape[:-2] + (6, n, n)
        return da.values.reshape(new_shape)

    data_vars = {}
    for vname, da in ds_cs_flat.data_vars.items():
        dims_new = [d for d in da.dims if d not in ('y', 'x')] + ['nf', 'Ydim', 'Xdim']
        data_vars[vname] = (dims_new, _reshape_var(da), da.attrs)

    coords = {k: v for k, v in ds_cs_flat.coords.items() if k not in ('lat', 'lon', 'y', 'x')}
    coords['nf']   = np.arange(6)
    coords['Ydim'] = np.arange(n)
    coords['Xdim'] = np.arange(n)
    coords['cs_lat'] = (['nf', 'Ydim', 'Xdim'], cs_lat)
    coords['cs_lon'] = (['nf', 'Ydim', 'Xdim'], cs_lon)

    return xr.Dataset(data_vars, coords=coords)


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def _accumulate_forcing(gchp_file, t_start, t_end, ts_chem_s, nlat, nlon):
    """
    Load GCHP sat-track file, apply obs operator, and accumulate forcing
    on an intermediate regular lat/lon grid.

    Returns
    -------
    active_times   : pd.DatetimeIndex  checkpoint times with observations
    active_forcing : ndarray (n_ckpt, nlev, nlat, nlon)
    active_n_obs   : ndarray (n_ckpt, nlat, nlon)
    levs           : model level coordinate array
    lats, lons     : 1-D coordinate arrays for the accumulation grid
    """
    ds_gchp = xr.open_dataset(gchp_file)
    ds_gchp['time'] = pd.to_datetime(ds_gchp['time'].values).round('s')

    nlev = ds_gchp.sizes['lev']
    levs = ds_gchp['lev'].values

    lats = np.linspace(-90,  90,  nlat)
    lons = np.linspace(-180, 180, nlon, endpoint=False)
    dlat = lats[1] - lats[0]
    dlon = lons[1] - lons[0]

    checkpoints = make_checkpoint_grid(t_start, t_end, ts_chem_s)
    n_ckpt      = len(checkpoints)
    forcing     = np.zeros((n_ckpt, nlev, nlat, nlon), dtype=np.float64)
    n_obs       = np.zeros((n_ckpt, nlat, nlon),        dtype=np.int32)

    times_gchp  = pd.DatetimeIndex(ds_gchp['time'].values)
    year_months = sorted(set(zip(times_gchp.year, times_gchp.month)))
    total       = 0

    for year, month in year_months:
        print(f'Processing {year}-{month:02d}')

        mask_gchp = (times_gchp.year == year) & (times_gchp.month == month)
        ds_mod    = ds_gchp.isel(time=np.where(mask_gchp)[0])
        co2_mod   = ds_mod['SpeciesConcVV_CO2'].values   # (nobs, nlev) v/v, no units
        prs_mod   = ds_mod['Met_PMIDDRY'].values          # (nobs, nlev) hPa, no units

        ds_obs = read_oco_monthly(year, month, FOLD_OBS)
        if ds_obs is None:
            print('  No OCO-2 data, skipping')
            continue

        prs_mod_u = prs_mod * units.hPa
        co2_mod_u = co2_mod * units.K   # fake unit required by interpolate_1d
        prs_obs_u = ds_obs['pressure'].values * units.hPa
        nobs, nlev_sat = prs_obs_u.shape

        co2_interp = np.empty((nobs, nlev_sat), dtype=np.float32)
        for i in range(nobs):
            result = interpolate_1d(prs_obs_u[i], prs_mod_u[i], co2_mod_u[i])
            co2_interp[i] = result.magnitude
            if np.isnan(co2_interp[i, 0]):
                co2_interp[i, 0] = co2_mod[i, 0]   # plain array, no units

        co2_pert     = co2_interp - ds_obs['CO2-apriori'].values
        xco2_hat_all = (ds_obs['xCO2-apriori'].values
                        + np.sum(ds_obs['xCO2-averagingKernel'].values * co2_pert, axis=1))
        xco2_obs_all = ds_obs['xCO2'].values

        if 'xCO2-uncertainty' in ds_obs:
            xco2_std_all = ds_obs['xCO2-uncertainty'].values
        elif 'xstd' in ds_obs:
            xco2_std_all = ds_obs['xstd'].values
        else:
            print('  WARNING: no uncertainty variable found, using 1 ppm')
            xco2_std_all = np.ones(nobs)

        obs_times    = pd.DatetimeIndex(ds_obs['time'].values)
        obs_lats     = ds_obs['latitude'].values
        obs_lons     = ds_obs['longitude'].values
        xAK_all      = ds_obs['xCO2-averagingKernel'].values
        xco2_hat_ppm = xco2_hat_all * 1e6
        xco2_obs_ppm = xco2_obs_all * 1e6
        xco2_std_ppm = xco2_std_all * 1e6

        for i in range(nobs):
            t_obs = obs_times[i]
            if t_obs < t_start or t_obs > t_end:
                continue
            ckpt_idx = nearest_checkpoint(t_obs, checkpoints, ts_chem_s)
            if ckpt_idx < 0:
                continue

            ilat = int(np.round((obs_lats[i] - lats[0]) / dlat))
            ilon = int(np.round(((obs_lons[i] + 180) % 360 - 180 - lons[0]) / dlon))
            ilat = np.clip(ilat, 0, nlat - 1)
            ilon = np.clip(ilon, 0, nlon - 1)

            force_model, _ = obs_forcing_profile(
                prs_mod[i], ds_obs['pressure'].values[i],
                xAK_all[i], xco2_hat_ppm[i], xco2_obs_ppm[i], xco2_std_ppm[i],
            )
            forcing[ckpt_idx, :, ilat, ilon] += force_model
            n_obs[ckpt_idx, ilat, ilon]       += 1
            total += 1

    print(f'Total observations matched: {total}')

    has_obs        = n_obs.sum(axis=(1, 2)) > 0
    active_times   = checkpoints[has_obs]
    active_forcing = forcing[has_obs]
    active_n_obs   = n_obs[has_obs]
    return active_times, active_forcing, active_n_obs, levs, lats, lons


# ---------------------------------------------------------------------------
# Write helpers  (one file per checkpoint)
# ---------------------------------------------------------------------------

FORCING_ATTRS = {
    'long_name': 'Adjoint forcing dJ/d(xCO2_model)',
    'units': 'ppm^-1',
    'comment': ('H^T S^{-1}(Hx-y) summed over obs in each checkpoint window. '
                'Unit conversion v/v->kg/box must be applied in Fortran.'),
}
N_OBS_ATTRS = {'long_name': 'Number of observations accumulated', 'units': '1'}
TIME_ENCODING = {'units': 'hours since 1900-01-01 00:00:00',
                 'calendar': 'proleptic_gregorian', 'dtype': 'float64'}
FORCING_ENCODING = {'dtype': 'float32', 'zlib': True, 'complevel': 4}


def _checkpoint_filename(output_dir, t):
    tstr = pd.Timestamp(t).strftime('%Y%m%d_%H%Mz')
    return os.path.join(output_dir, f'CO2_adjoint_forcing_{tstr}.nc4')


def write_latlon(output_path, t, forcing_t, n_obs_t, levs, lats, lons):
    """Write a single checkpoint lat/lon forcing file."""
    ds = xr.Dataset(
        {
            'forcing': (['time', 'lev', 'lat', 'lon'], forcing_t[np.newaxis],  FORCING_ATTRS),
            'n_obs':   (['time', 'lat', 'lon'],         n_obs_t[np.newaxis],    N_OBS_ATTRS),
        },
        coords={'time': [t], 'lev': levs, 'lat': lats, 'lon': lons},
    )
    ds['lat'].attrs = {'long_name': 'Latitude',  'units': 'degrees_north'}
    ds['lon'].attrs = {'long_name': 'Longitude', 'units': 'degrees_east'}
    ds['lev'].attrs = {'long_name': 'Model level (1=surface, LLPAR=TOA)'}
    ds.to_netcdf(output_path,
                 encoding={'time': TIME_ENCODING, 'forcing': FORCING_ENCODING})
    print(f'  Written {output_path}')


def write_cubedsphere(output_path, t, forcing_t, n_obs_t, levs, lats, lons, cs_res):
    """Write a single checkpoint cubed-sphere forcing file."""
    ds_ll = xr.Dataset(
        {
            'forcing': (['time', 'lev', 'lat', 'lon'],
                        forcing_t[np.newaxis].astype(np.float32), FORCING_ATTRS),
            'n_obs':   (['time', 'lat', 'lon'],
                        n_obs_t[np.newaxis].astype(np.float32),   N_OBS_ATTRS),
        },
        coords={'time': [t], 'lev': levs, 'lat': lats, 'lon': lons},
    )
    ds_ll['lat'].attrs = {'units': 'degrees_north'}
    ds_ll['lon'].attrs = {'units': 'degrees_east'}

    ds_cs = regrid_latlon_to_cubedsphere(ds_ll, cs_res)
    ds_cs['lev'].attrs = {'long_name': 'Model level (1=surface, LLPAR=TOA)'}
    ds_cs['cs_lat'].attrs = {'long_name': 'Latitude of cubed-sphere cell centre',
                             'units': 'degrees_north'}
    ds_cs['cs_lon'].attrs = {'long_name': 'Longitude of cubed-sphere cell centre',
                             'units': 'degrees_east'}
    ds_cs.to_netcdf(output_path,
                    encoding={'time': TIME_ENCODING, 'forcing': FORCING_ENCODING})
    print(f'  Written {output_path}')


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def co2_adjoint_forcing(gchp_file, output_dir,
                        ts_chem_s, t_start, t_end,
                        nlat, nlon,
                        grid='latlon', cs_res=None):

    os.makedirs(output_dir, exist_ok=True)

    active_times, active_forcing, active_n_obs, levs, lats, lons = \
        _accumulate_forcing(gchp_file, t_start, t_end, ts_chem_s, nlat, nlon)

    if not len(active_times):
        print('No observations matched any checkpoint — output not written.')
        return

    print(f'Writing {len(active_times)} checkpoint file(s) to {output_dir}/')
    for i, t in enumerate(active_times):
        fpath = _checkpoint_filename(output_dir, t)
        if grid == 'latlon':
            write_latlon(fpath, t, active_forcing[i], active_n_obs[i], levs, lats, lons)
        else:
            write_cubedsphere(fpath, t, active_forcing[i], active_n_obs[i],
                              levs, lats, lons, cs_res)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('gchp_file',   help='GCHP sat-track netCDF file')
    parser.add_argument('output_dir',  help='Output directory (created if absent)')
    parser.add_argument('--ts-chem', type=int, required=True, dest='ts_chem_s',
                        help='Chemistry timestep in seconds')
    parser.add_argument('--t-start', required=True,
                        help='Assimilation window start YYYY-MM-DDTHH:MM:SS')
    parser.add_argument('--t-end',   required=True,
                        help='Assimilation window end   YYYY-MM-DDTHH:MM:SS')
    parser.add_argument('--grid', choices=['latlon', 'cubedsphere'], default='latlon',
                        help='Output grid (default: latlon)')
    parser.add_argument('--cs-res', type=int, default=None,
                        help='Cubed-sphere resolution N for C{N} (required if --grid cubedsphere)')
    parser.add_argument('--nlat', type=int, default=181,
                        help='Intermediate accumulation grid latitude points (default 181 → 1°)')
    parser.add_argument('--nlon', type=int, default=360,
                        help='Intermediate accumulation grid longitude points (default 360 → 1°)')

    args = parser.parse_args()

    if args.grid == 'cubedsphere' and args.cs_res is None:
        parser.error('--cs-res is required when --grid cubedsphere')

    co2_adjoint_forcing(
        gchp_file   = args.gchp_file,
        output_dir  = args.output_dir,
        ts_chem_s   = args.ts_chem_s,
        t_start     = pd.Timestamp(args.t_start),
        t_end       = pd.Timestamp(args.t_end),
        nlat        = args.nlat,
        nlon        = args.nlon,
        grid        = args.grid,
        cs_res      = args.cs_res,
    )
