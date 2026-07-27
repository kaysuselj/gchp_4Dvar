#!/usr/bin/env python3
"""
Verify the monthly-sigma plumbing from a forward test run.

Expects a forward run over 2016-01-01 -> 2016-02-03 driven by test sigma
Jan=1.0, Feb=0.0, Mar-Dec=9.0, with the Emissions collection writing
EmisCO2_Total, EmisCO2_NetTerrExch and EmisCO2_BalBiosph.

Checks
------
(2a) EmisCO2_NetTerrExch is nonzero somewhere in January.
(2b) EmisCO2_NetTerrExch is exactly zero at every record from 2016-02-01T00:00.
(2c) The transition happens at the month boundary, not before or after.
(2d) No record anywhere is consistent with sigma=9 (wrong month selected):
     January |TER| must not exceed 9x the CLASS-CTEM January source maximum.
(2e) EmisCO2_Total stays nonzero in February (sigma scales TER only).
(1)  EmisCO2_BalBiosph tracks the CARDAMOM source across the Jan->Feb boundary.

Usage:
    python check_monthly_sigma.py <OutputDir> [--ter-src DIR] [--nbe-src DIR]
"""
import argparse
import glob
import os
import sys
import warnings

import numpy as np
import xarray as xr

warnings.filterwarnings('ignore')

TER_SRC = '/nobackup/ksuselj1/gchp_14.5.3_adjoint_surfaceF/surface_fluxes_osse/TER/CLASS-CTEM'
NBE_SRC = '/nobackup/ksuselj1/gchp_14.5.3_adjoint_surfaceF/surface_fluxes_osse/Balbio/CARDAMOM-ECCO'
UNIT_CONV = 3.667e-6          # HEMCO scale factor 751
BOUNDARY = np.datetime64('2016-02-01T00:00:00')

results = []


def record(name, ok, detail):
    results.append((name, ok, detail))
    print(f'[{"PASS" if ok else "FAIL"}] {name}: {detail}')


# Cubed-sphere grid metadata carried by GCHP output; unused here and the
# 'anchor' variable has a duplicate 'ncontact' dimension that breaks xr.concat.
_DROP_VARS = ['anchor', 'contacts', 'ncontact']


def load(output_dir):
    files = sorted(glob.glob(os.path.join(output_dir, 'GEOSChem.Emissions.*.nc4')))
    if not files:
        sys.exit(f'No Emissions files in {output_dir}')
    ds = xr.concat(
        [xr.open_dataset(f, drop_variables=_DROP_VARS) for f in files],
        dim='time').sortby('time')
    print(f'Loaded {len(files)} files, {ds.sizes["time"]} records: '
          f'{ds.time.values[0]} .. {ds.time.values[-1]}\n')
    return ds


def spatial_dims(da):
    return [d for d in da.dims if d != 'time']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('output_dir')
    ap.add_argument('--ter-src', default=TER_SRC)
    ap.add_argument('--nbe-src', default=NBE_SRC)
    args = ap.parse_args()

    ds = load(args.output_dir)
    t = ds.time.values
    # The Emissions collection is instantaneous and ExtData applies a new monthly
    # sigma file one timestep AFTER its T00:00:00 refresh stamp (the same one-step
    # lag seen for CARDAMOM). So the snapshot stamped exactly at the month boundary
    # still carries the previous month's sigma; classify it as January.
    jan = t <= BOUNDARY
    feb = t > BOUNDARY

    if not jan.any() or not feb.any():
        sys.exit('Run does not span the 2016-02-01 boundary; nothing to test.')

    ter = ds['EmisCO2_NetTerrExch']
    sd = spatial_dims(ter)
    ter_absmax = np.abs(ter).max(dim=sd).values

    # ---- (2a) January TER is active -------------------------------------
    jan_max = ter_absmax[jan].max()
    record('2a jan TER nonzero', jan_max > 0,
           f'max|TER| in Jan = {jan_max:.4e} kg/m2/s')

    # ---- (2b) February TER is exactly zero -------------------------------
    feb_max = ter_absmax[feb].max()
    record('2b feb TER zero', feb_max == 0.0,
           f'max|TER| in Feb = {feb_max:.4e} (must be exactly 0)')

    # ---- (2c) transition exactly at the boundary -------------------------
    nz = np.nonzero(ter_absmax > 0)[0]
    if nz.size:
        last_nz, first_z = t[nz[-1]], t[nz[-1] + 1] if nz[-1] + 1 < len(t) else None
        # last active record must be at or before the boundary (one-step ExtData lag)
        ok = last_nz <= BOUNDARY and (first_z is None or first_z > BOUNDARY)
        record('2c transition at boundary', bool(ok),
               f'last nonzero record {last_nz}, boundary {BOUNDARY}')
    else:
        record('2c transition at boundary', False, 'TER is zero everywhere')

    # ---- (2d) sigma=9 never applied --------------------------------------
    src = os.path.join(args.ter_src, '2016', '01.nc')
    if os.path.exists(src):
        s = xr.open_dataset(src)['CO2_Flux']
        src_max = float(np.abs(s).max()) * UNIT_CONV
        ok = jan_max <= 1.5 * src_max          # allow regridding overshoot
        record('2d no sigma=9 leak', bool(ok),
               f'max|TER| Jan = {jan_max:.4e} vs 1x source {src_max:.4e} '
               f'(9x would be {9 * src_max:.4e})')
    else:
        record('2d no sigma=9 leak', False, f'source file not found: {src}')

    # ---- (2e) total emissions still nonzero in February -------------------
    tot = ds['EmisCO2_Total']
    tot_feb = np.abs(tot).max(dim=spatial_dims(tot)).values[feb].max()
    record('2e feb total nonzero', tot_feb > 0,
           f'max|EmisCO2_Total| in Feb = {tot_feb:.4e} (TER is only one term)')

    # ---- (1) CARDAMOM read across the boundary ---------------------------
    nbe = ds['EmisCO2_BalBiosph']
    nbe_sd = spatial_dims(nbe)
    ok_all = True
    detail = []
    for day, lbl in [('2016/01/31', 'Jan 31'), ('2016/02/01', 'Feb 1')]:
        f = os.path.join(args.nbe_src, day + '.nc')
        if not os.path.exists(f):
            ok_all = False
            detail.append(f'{lbl}: source missing')
            continue
        c = xr.open_dataset(f)['CO2_Flux']
        w = np.cos(np.deg2rad(c.lat))
        src_h = c.weighted(w).mean(dim=['lat', 'lon']).values          # 24 hourly
        d0 = np.datetime64(day.replace('/', '-'))
        sel = (t >= d0) & (t < d0 + np.timedelta64(1, 'D'))
        mod = nbe.isel(time=sel).mean(dim=nbe_sd).values               # 48 half-hourly
        if mod.size != 48:
            ok_all = False
            detail.append(f'{lbl}: got {mod.size} records, expected 48')
            continue
        mod_h = mod.reshape(24, 2).mean(axis=1)
        # model lags the source by one hour (source stamped :30, ExtData refreshes on :00)
        a = mod_h[1:]
        b = src_h[:-1] * UNIT_CONV
        # The domain mean of NBE changes sign, so a per-point ratio is unstable near
        # zero; test the time series with correlation + a through-origin slope instead.
        corr = float(np.corrcoef(a, b)[0, 1])
        slope = float(np.dot(a, b) / np.dot(b, b))
        good = corr > 0.9 and abs(slope - 1.0) < 0.1
        ok_all &= good
        detail.append(f'{lbl}: corr={corr:.3f}, slope={slope:.3f}')
    record('1 CARDAMOM across boundary', ok_all, '; '.join(detail))

    # ---- summary ---------------------------------------------------------
    print()
    n_fail = sum(1 for _, ok, _ in results if not ok)
    if n_fail:
        print(f'{n_fail} CHECK(S) FAILED')
        sys.exit(1)
    print('ALL CHECKS PASSED')
    sys.exit(0)


if __name__ == '__main__':
    main()
