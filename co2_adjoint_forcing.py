#!/usr/bin/env python3
"""
Compute 3D adjoint forcing from OCO-2 xCO2 observations, binned to GCHP
chemistry checkpoint times.  One netCDF file is written per checkpoint,
including empty files (n_obs=0) so Fortran can distinguish 'no obs' from
'file missing'.  Forcing is stored per-observation (sparse); Fortran does
the cell assignment using State_Grid%{X,Y}{Min,Max} cell boundaries.

Usage:
    python ./co2_adjoint_forcing.py gchp_file output_dir \\
        --ts-chem 1200 \\
        --t-start 2015-01-01T00:00:00 --t-end 2015-01-31T23:59:59

Arguments:
    gchp_file    : GCHP sat-track netCDF file (GEOSChem.sat_track.*)
    output_dir   : output directory (created if absent); one file per checkpoint

    --ts-chem    : chemistry timestep [s] for checkpoint binning (required)
    --t-start    : start of assimilation window YYYY-MM-DDTHH:MM:SS (required)
    --t-end      : end of assimilation window (required)

Environment setup:
    python -m venv gchp-env
    source gchp-env/bin/activate
    pip install numpy pandas xarray netCDF4 metpy

    The following file must be in the same directory:
        co2_sat_compare_monthly.py

Output files:
    CO2_adjoint_forcing_YYYYMMDD_HHMMz.nc4  (one per checkpoint, always)

Output variables (sparse per-obs format):
    lat_obs  : (obs,)      observation latitudes  [degrees_north]
    lon_obs  : (obs,)      observation longitudes [-180, 180) [degrees_east]
    forcing  : (obs, lev)  dJ/d(CO2_mmr) [1/(kg_CO2/kg_dry_air)]

Notes:
    - Level 1 = surface (highest pressure), level LLPAR = TOA.
    - Fortran reads forcing as (nlev, n_obs) in column-major (matches Python
      (obs, lev) row-major layout in netCDF).
    - OCO-2 data is stored by month; observations outside [t_start, t_end]
      are discarded in the per-observation loop.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr

# Molecular weights [g/mol] for unit conversion ppm^-1 → 1/(kg_CO2/kg_dry_air)
M_AIR = 28.97
M_CO2 = 44.01
# SpeciesAdj [1/(kg_CO2/kg_dry)] = forcing [ppm^-1] * (M_AIR/M_CO2) * 1e6
PPM_TO_MMR_ADJ = (M_AIR / M_CO2) * 1e6   # ≈ 6.585e5

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
    except ImportError as _e:
        raise ImportError(
            f'xesmf import failed: {_e}\n'
            'Try:  pip install xesmf esmpy\n'
            'On HPC systems the ESMF C library must also be available.\n'
            'If "module load esmf" is available, run that before activating the venv.'
        ) from _e

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

def _accumulate_forcing(gchp_file, t_start, t_end, ts_chem_s):
    """
    Load GCHP sat-track file, apply obs operator, and collect per-obs forcing.

    Each GCHP profile is matched to the nearest OCO-2 observation by time
    (within OBS_MATCH_TOL).  Forcing is accumulated per observation, with
    longitude normalized to [-180, 180).  No intermediate lat/lon grid is used;
    Fortran does the nearest-cell assignment directly on the cubed sphere.

    Returns
    -------
    checkpoints  : pd.DatetimeIndex  all checkpoint times in the window
    obs_by_ckpt  : dict {ckpt_idx: [(lat, lon, force_profile), ...]}
    levs         : model level coordinate array
    """
    OBS_MATCH_TOL = pd.Timedelta(seconds=60)

    ds_gchp = xr.open_dataset(gchp_file)
    ds_gchp['time'] = pd.to_datetime(ds_gchp['time'].values).round('s')

    nlev = ds_gchp.sizes['lev']
    levs = ds_gchp['lev'].values

    checkpoints = make_checkpoint_grid(t_start, t_end, ts_chem_s)
    n_ckpt      = len(checkpoints)
    obs_by_ckpt = {i: [] for i in range(n_ckpt)}
    J_total     = 0.0

    # Restrict GCHP profiles to the assimilation window
    times_gchp  = pd.DatetimeIndex(ds_gchp['time'].values)
    window_mask = (times_gchp >= t_start) & (times_gchp <= t_end)
    ds_window   = ds_gchp.isel(time=np.where(window_mask)[0])
    times_win   = pd.DatetimeIndex(ds_window['time'].values)

    if len(times_win) == 0:
        print('No GCHP profiles within the assimilation window.')
        return checkpoints, obs_by_ckpt, levs, J_total

    # Only read OCO-2 months that are actually needed
    year_months = sorted(set(zip(times_win.year, times_win.month)))
    total = 0

    for year, month in year_months:
        print(f'Processing {year}-{month:02d}')

        # GCHP profiles for this month (already within window)
        mask_month   = (times_win.year == year) & (times_win.month == month)
        ds_mod       = ds_window.isel(time=np.where(mask_month)[0])
        gchp_times_m = pd.DatetimeIndex(ds_mod['time'].values)
        co2_mod      = ds_mod['SpeciesConcVV_CO2'].transpose('time', 'lev').values  # (n_gchp, nlev)
        prs_mod      = ds_mod['Met_PMIDDRY'].transpose('time', 'lev').values         # (n_gchp, nlev)

        ds_obs = read_oco_monthly(year, month, FOLD_OBS)
        if ds_obs is None:
            print('  No OCO-2 data, skipping')
            continue

        obs_times    = pd.DatetimeIndex(ds_obs['time'].values)
        obs_lats     = ds_obs['latitude'].values
        obs_lons     = ds_obs['longitude'].values
        prs_obs_all  = ds_obs['pressure'].values           # (n_obs, nlev_sat)
        xAK_all      = ds_obs['xCO2-averagingKernel'].values
        co2_apr_all  = ds_obs['CO2-apriori'].values
        xco2_apr_all = ds_obs['xCO2-apriori'].values
        xco2_obs_all = ds_obs['xCO2'].values

        if 'xCO2-uncertainty' in ds_obs:
            xco2_std_all = ds_obs['xCO2-uncertainty'].values
        elif 'xstd' in ds_obs:
            xco2_std_all = ds_obs['xstd'].values
        else:
            print('  WARNING: no uncertainty variable found, using 1 ppm')
            xco2_std_all = np.ones(len(obs_times))

        # Loop over GCHP profiles; match each to an OCO-2 obs by time
        for j in range(len(gchp_times_m)):
            t_gchp  = gchp_times_m[j]
            deltas  = np.abs(obs_times - t_gchp)
            idx_obs = int(deltas.argmin())
            if deltas[idx_obs] > OBS_MATCH_TOL:
                continue

            prs_mod_j = prs_mod[j]            # (nlev,)
            co2_mod_j = co2_mod[j]            # (nlev,)
            prs_obs_j = prs_obs_all[idx_obs]  # (nlev_sat,)

            # Interpolate model CO2 to satellite pressure levels.
            # np.interp requires ascending x; model levels run surface→TOA
            # (descending pressure), so flip before interpolating.
            # Out-of-bounds satellite levels (above model top or below surface)
            # receive the boundary model value — no NaN, no warning.
            sort_idx   = np.argsort(prs_mod_j)          # ascending pressure indices
            co2_interp = np.interp(prs_obs_j,
                                   prs_mod_j[sort_idx],
                                   co2_mod_j[sort_idx])

            # Apply OCO-2 observation operator
            co2_pert     = co2_interp - co2_apr_all[idx_obs]
            xco2_hat     = xco2_apr_all[idx_obs] + np.sum(xAK_all[idx_obs] * co2_pert)
            xco2_hat_ppm = float(xco2_hat)                 * 1e6
            xco2_obs_ppm = float(xco2_obs_all[idx_obs])    * 1e6
            xco2_std_ppm = float(xco2_std_all[idx_obs])    * 1e6

            ckpt_idx = nearest_checkpoint(t_gchp, checkpoints, ts_chem_s)
            if ckpt_idx < 0:
                continue

            force_model, diff = obs_forcing_profile(
                prs_mod_j, prs_obs_j,
                xAK_all[idx_obs], xco2_hat_ppm, xco2_obs_ppm, xco2_std_ppm,
            )
            J_total += 0.5 * (diff / xco2_std_ppm) ** 2
            # Normalize lon to [-180, 180) and store per-obs entry
            lon_norm = float((float(obs_lons[idx_obs]) + 180) % 360 - 180)
            obs_by_ckpt[ckpt_idx].append((
                float(obs_lats[idx_obs]),
                lon_norm,
                force_model * PPM_TO_MMR_ADJ,   # ppm^-1 → 1/(kg_CO2/kg_dry_air)
            ))
            total += 1

    n_with_obs = sum(1 for obs_list in obs_by_ckpt.values() if obs_list)
    print(f'Total observations matched: {total}  '
          f'({n_with_obs}/{n_ckpt} checkpoints have at least one obs)')
    print(f'Cost function J = {J_total:.6e}')

    return checkpoints, obs_by_ckpt, levs, J_total


# ---------------------------------------------------------------------------
# Write helpers  (one file per checkpoint + optional consolidated diagnostics)
# ---------------------------------------------------------------------------

FORCING_ATTRS = {
    'long_name': 'Adjoint forcing dJ/d(CO2_mass_mixing_ratio)',
    'units': '1/(kg_CO2/kg_dry_air)',
    'comment': ('H^T S^{-1}(Hx-y) * (M_air/M_CO2) * 1e6, summed over obs in '
                'each checkpoint window. Ready to load directly into '
                'State_Chm%SpeciesAdj [J/(kg_CO2/kg_dry_air)].'),
}
N_OBS_ATTRS = {'long_name': 'Number of observations accumulated', 'units': '1'}
TIME_ENCODING = {'units': 'hours since 1900-01-01 00:00:00',
                 'calendar': 'proleptic_gregorian', 'dtype': 'float64'}
FORCING_ENCODING = {'dtype': 'float32', 'zlib': True, 'complevel': 4}


def _checkpoint_filename(output_dir, t):
    tstr = pd.Timestamp(t).strftime('%Y%m%d_%H%Mz')
    return os.path.join(output_dir, f'CO2_adjoint_forcing_{tstr}.nc4')


def write_sparse(output_path, t, obs_lats, obs_lons, forcing_k, levs):
    """Write a per-obs sparse forcing file for one checkpoint.

    Always written, even when n_obs == 0 (empty obs dimension), so Fortran
    can distinguish 'no obs this checkpoint' from 'file missing'.

    forcing_k : (n_obs, nlev) float32
        Dim order ('obs', 'lev') → Fortran reads as forcing_f(nlev, n_obs).
    """
    n_obs = len(obs_lats)
    nlev  = len(levs)
    forcing_arr = np.asarray(forcing_k, dtype=np.float32).reshape(n_obs, nlev)
    ds = xr.Dataset(
        {
            'lat_obs': (['obs'], np.asarray(obs_lats, dtype=np.float32),
                        {'long_name': 'Observation latitude',  'units': 'degrees_north'}),
            'lon_obs': (['obs'], np.asarray(obs_lons, dtype=np.float32),
                        {'long_name': 'Observation longitude', 'units': 'degrees_east',
                         'comment': 'normalized to [-180, 180)'}),
            'forcing': (['obs', 'lev'], forcing_arr, FORCING_ATTRS),
        },
        coords={'time': [t], 'lev': levs},
    )
    ds['lev'].attrs = {'long_name': 'Model level (1=surface, LLPAR=TOA)'}
    # zlib compression requires at least one element; skip for empty files
    forcing_enc = FORCING_ENCODING if n_obs > 0 else {'dtype': 'float32'}
    ds.to_netcdf(output_path, encoding={'time': TIME_ENCODING, 'forcing': forcing_enc})
    print(f'  Written: {output_path}  ({n_obs} obs)')


# ---------------------------------------------------------------------------
# Consolidated diagnostics file  (optional; one per co2_adjoint_forcing call)
# ---------------------------------------------------------------------------

def write_consolidated(output_path, obs_by_ckpt, levs):
    """Write all observations from all checkpoints into a single netCDF file.

    Concatenates lat_obs, lon_obs, and forcing across every checkpoint that
    has at least one observation.  Empty checkpoints are skipped.

    The file is intended for post-processing only (plots, diagnostics); it is
    NOT read by GCHP.  Disable with save_diagnostics=False if the observation
    count is very large.

    forcing : (n_total_obs, nlev)  — same units as the per-checkpoint files.
    """
    all_lats, all_lons, all_forcing = [], [], []
    for obs_list in obs_by_ckpt.values():
        for lat, lon, force in obs_list:
            all_lats.append(lat)
            all_lons.append(lon)
            all_forcing.append(force)

    if not all_lats:
        print('  No observations — consolidated forcing file not written.')
        return

    n_obs = len(all_lats)
    ds = xr.Dataset(
        {
            'lat_obs': (['obs'], np.array(all_lats, dtype=np.float32),
                        {'long_name': 'Observation latitude',  'units': 'degrees_north'}),
            'lon_obs': (['obs'], np.array(all_lons, dtype=np.float32),
                        {'long_name': 'Observation longitude', 'units': 'degrees_east',
                         'comment': 'normalized to [-180, 180)'}),
            'forcing': (['obs', 'lev'],
                        np.array(all_forcing, dtype=np.float32),
                        FORCING_ATTRS),
        },
        coords={'lev': levs},
    )
    ds['lev'].attrs = {'long_name': 'Model level (1=surface, LLPAR=TOA)'}
    ds.to_netcdf(output_path, encoding={'forcing': FORCING_ENCODING})
    print(f'  Consolidated forcing written: {output_path}  ({n_obs} obs total)')


# ---------------------------------------------------------------------------
# Daily mean forcing profile  (always written; small file)
# ---------------------------------------------------------------------------

def write_daily_forcing(output_path, obs_by_ckpt, levs, checkpoints):
    """Write one mean forcing profile (nlev,) per calendar day.

    Groups all observations by the calendar day of their checkpoint time,
    then computes the mean forcing profile across all obs within that day.
    The result is a 2-D array mean_forcing(n_days, nlev).

    This file is always written regardless of save_diagnostics — it is small
    (n_days × nlev × 4 bytes) and useful for comparing how the adjoint signal
    evolves across iterations.
    """
    daily = {}  # date → list of (nlev,) arrays
    for ckpt_idx, obs_list in obs_by_ckpt.items():
        if not obs_list:
            continue
        t_day = pd.Timestamp(checkpoints[ckpt_idx]).normalize()
        if t_day not in daily:
            daily[t_day] = []
        for _, _, force in obs_list:
            daily[t_day].append(force)   # force is (nlev,)

    if not daily:
        print('  No observations — daily forcing file not written.')
        return

    days     = sorted(daily.keys())
    profiles = np.array([np.array(daily[d]).mean(axis=0) for d in days],
                        dtype=np.float32)   # (n_days, nlev)
    n_obs_day = np.array([len(daily[d]) for d in days], dtype=np.int32)
    dates     = np.array([np.datetime64(d.date(), 'D') for d in days])

    ds = xr.Dataset(
        {
            'mean_forcing': (
                ['date', 'lev'], profiles,
                {'long_name': 'Daily mean adjoint forcing profile',
                 'units': '1/(kg_CO2/kg_dry_air)',
                 'comment': 'Mean of forcing(obs, lev) over all obs in the calendar day'},
            ),
            'n_obs': (
                ['date'], n_obs_day,
                {'long_name': 'Number of observations in this calendar day'},
            ),
        },
        coords={'date': dates, 'lev': levs},
    )
    ds['lev'].attrs  = {'long_name': 'Model level (1=surface, LLPAR=TOA)'}
    ds['date'].attrs = {'long_name': 'Calendar day (UTC)'}
    ds.to_netcdf(output_path,
                 encoding={'mean_forcing': {'dtype': 'float32', 'zlib': True,
                                            'complevel': 4}})
    print(f'  Daily forcing written: {output_path}  '
          f'({len(days)} day(s), {n_obs_day.sum()} obs total)')


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def co2_adjoint_forcing(gchp_file, output_dir, ts_chem_s, t_start, t_end,
                        save_diagnostics=True):

    os.makedirs(output_dir, exist_ok=True)

    checkpoints, obs_by_ckpt, levs, J = \
        _accumulate_forcing(gchp_file, t_start, t_end, ts_chem_s)

    j_path = os.path.join(output_dir, 'J_value.txt')
    with open(j_path, 'w') as fh:
        fh.write(f'{J:.10e}\n')
    print(f'J written to {j_path}')

    nlev = len(levs)
    print(f'Writing {len(checkpoints)} checkpoint file(s) to {output_dir}/ '
          f'(zero-obs files included for all checkpoints)')
    for i, t in enumerate(checkpoints):
        fpath     = _checkpoint_filename(output_dir, t)
        obs_list  = obs_by_ckpt[i]
        if obs_list:
            obs_lats_k = np.array([d[0] for d in obs_list], dtype=np.float32)
            obs_lons_k = np.array([d[1] for d in obs_list], dtype=np.float32)
            forcing_k  = np.array([d[2] for d in obs_list], dtype=np.float32)
        else:
            obs_lats_k = np.empty(0, dtype=np.float32)
            obs_lons_k = np.empty(0, dtype=np.float32)
            forcing_k  = np.empty((0, nlev), dtype=np.float32)
        write_sparse(fpath, t, obs_lats_k, obs_lons_k, forcing_k, levs)

    if save_diagnostics:
        diag_path = os.path.join(output_dir, 'forcing_all_obs.nc4')
        write_consolidated(diag_path, obs_by_ckpt, levs)

    daily_path = os.path.join(output_dir, 'daily_forcing_stats.nc4')
    write_daily_forcing(daily_path, obs_by_ckpt, levs, checkpoints)

    return J


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
    parser.add_argument('--no-save-diagnostics', dest='save_diagnostics',
                        action='store_false', default=True,
                        help='Skip writing forcing_all_obs.nc4 (use for large obs counts)')
    args = parser.parse_args()

    co2_adjoint_forcing(
        gchp_file        = args.gchp_file,
        output_dir       = args.output_dir,
        ts_chem_s        = args.ts_chem_s,
        t_start          = pd.Timestamp(args.t_start),
        t_end            = pd.Timestamp(args.t_end),
        save_diagnostics = args.save_diagnostics,
    )
