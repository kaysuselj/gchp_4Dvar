#!/usr/bin/env python3
"""Offline replacement for the MAPL `track_file` History sampler (absent from MAPL 2.69 / GCHP 14.8.0).

Reads every-step instantaneous History output of the forward run (one 3-D/2-D field set per
model step) and the satellite track file, and writes GEOSChem.sat_track.* files with the same
names, variables, dimensions and time convention as the 14.5.3 MAPL sampler, so that
co2_adjoint_forcing_osse.py and the 4D-Var scripts work unchanged.

Sampling rules of the 14.5.3 sampler (MAPL_HistoryTrajectoryMod.F90, get_current_interval;
MAPL_LocstreamRegridder default ESMF_REGRIDMETHOD_NEAREST_STOD):
  * a sounding at time t with T_{k-1} < t <= T_k is sampled from the model state at T_k
    (History is called at the end of each step); soundings at or before the first output time
    go to the first output time;
  * horizontal: value of the grid cell whose centre is nearest on the sphere (no interpolation);
  * fields are copied as stored (no vertical flip): GEOS-Chem state fields are surface first,
    HEMCO emission fields top-of-atmosphere first, exactly as in the MAPL output.
Output: one file per calendar month of soundings, named like MAPL did (one month EARLIER than the
data, e.g. January soundings -> GEOSChem.sat_track.20151201_0030z.nc4) with
time units 'days since <that name's date> 00:30:00'.

Optional --method bilinear is reserved for the planned bilinear option (not implemented yet).

  make_sat_track_offline.py --track track_file_osse_2016.nc --inst-glob 'OutputDir/GEOSChem.Inst30.*.nc4' \
      --out-dir OutputDir [--t-start 2016-01-01 --t-end 2016-04-01]
"""
import argparse, glob, os
from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc

FIELDS3 = ['EmisCO2_Total', 'SpeciesConcVV_CO2', 'Met_AD', 'Met_BXHEIGHT', 'Met_DELP', 'Met_PMID', 'Met_PMIDDRY']
FIELDS2 = ['Met_PSC2WET', 'Met_PSC2DRY']


def unitvec(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def month_back(y, m):
    return (y - 1, 12) if m == 1 else (y, m - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--track', required=True)
    ap.add_argument('--inst-glob', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--t-start'); ap.add_argument('--t-end')
    ap.add_argument('--method', default='nearest', choices=['nearest'])
    a = ap.parse_args()

    tr = nc.Dataset(a.track)
    tt = nc.num2date(tr['time'][:], tr['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    tlat = np.asarray(tr['latitude'][:], 'f8'); tlon = np.asarray(tr['longitude'][:], 'f8')
    tt = np.array(tt)

    files = sorted(glob.glob(a.inst_glob))
    assert files, f'no files for {a.inst_glob}'
    # index of every output record: (time, file, record)
    recs = []
    for f in files:
        d = nc.Dataset(f)
        ts = nc.num2date(d['time'][:], d['time'].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        for k, t in enumerate(np.atleast_1d(ts)):
            recs.append((t, f, k))
        d.close()
    recs.sort(key=lambda r: r[0])
    T = np.array([r[0] for r in recs])
    t0, t1 = T[0], T[-1]
    if a.t_start: t0 = max(t0, datetime.fromisoformat(a.t_start))
    if a.t_end: t1 = min(t1, datetime.fromisoformat(a.t_end))

    # grid: cell centres from the first file
    d0 = nc.Dataset(files[0])
    glat = np.asarray(d0['lats'][:], 'f8').ravel(); glon = np.asarray(d0['lons'][:], 'f8').ravel()
    shape = d0['SpeciesConcVV_CO2'].shape[1:]        # (lev, nf, Ydim, Xdim)
    nlev = shape[0]
    f3 = [v for v in FIELDS3 if v in d0.variables]    # write what the run provides
    f2 = [v for v in FIELDS2 if v in d0.variables]    # (the forcing/obs scripts need SpeciesConcVV_CO2, Met_PMIDDRY)
    assert 'SpeciesConcVV_CO2' in f3 and 'Met_PMIDDRY' in f3, 'Inst output must contain SpeciesConcVV_CO2 and Met_PMIDDRY'
    d0.close()
    G = unitvec(glat, glon)

    # soundings inside the output period; step assignment T_{k-1} < t <= T_k
    sel = np.where((tt > T[0] - timedelta(seconds=1)) & (tt <= T[-1]))[0]
    sel = sel[(tt[sel] >= t0) & (tt[sel] <= t1 + timedelta(seconds=1))] if (a.t_start or a.t_end) else sel
    kstep = np.searchsorted(T, tt[sel], side='left')   # first T_k >= t
    kstep = np.clip(kstep, 0, len(T) - 1)
    cell = np.argmax(unitvec(tlat[sel], tlon[sel]) @ G.T, axis=1)  # nearest centre (max cosine)

    n = len(sel)
    out = {v: np.full((nlev, n), np.nan, 'f4') for v in f3}
    out.update({v: np.full(n, np.nan, 'f4') for v in f2})
    order = np.argsort(kstep, kind='stable')
    cur = None
    for j in order:
        _, f, k = recs[kstep[j]]
        if cur != (f, k):
            ds = nc.Dataset(f); cur = (f, k)
            arr3 = {v: np.asarray(ds[v][k], 'f4').reshape(nlev, -1) for v in f3}
            arr2 = {v: np.asarray(ds[v][k], 'f4').ravel() for v in f2}
            ds.close()
        for v in f3: out[v][:, j] = arr3[v][:, cell[j]]
        for v in f2: out[v][j] = arr2[v][cell[j]]

    # write one file per calendar month of the soundings, MAPL naming (one month earlier)
    os.makedirs(a.out_dir, exist_ok=True)
    months = sorted({(t.year, t.month) for t in tt[sel]})
    for (y, m) in months:
        idx = np.array([i for i in range(n) if (tt[sel][i].year, tt[sel][i].month) == (y, m)])
        yy, mm = month_back(y, m)
        ref = datetime(yy, mm, 1, 0, 30)
        path = os.path.join(a.out_dir, f'GEOSChem.sat_track.{yy:04d}{mm:02d}01_0030z.nc4')
        with nc.Dataset(path, 'w') as o:
            o.createDimension('lev', nlev); o.createDimension('time', len(idx))
            v = o.createVariable('lev', 'f8', ('lev',)); v[:] = np.arange(1, nlev + 1); v.units = 'layer'
            v = o.createVariable('time', 'f4', ('time',)); v.units = f'days since {ref:%Y-%m-%d %H:%M:%S}'
            v[:] = [(tt[sel][i] - ref).total_seconds() / 86400.0 for i in idx]
            v = o.createVariable('longitude', 'f8', ('time',)); v[:] = tlon[sel][idx]; v.units = 'degrees_east'
            v = o.createVariable('latitude', 'f8', ('time',)); v[:] = tlat[sel][idx]; v.units = 'degrees_east'
            for name in f3:
                v = o.createVariable(name, 'f4', ('lev', 'time')); v[:] = out[name][:, idx]
            for name in f2:
                v = o.createVariable(name, 'f4', ('time',)); v[:] = out[name][idx]
            o.description = 'offline satellite-track sampling of every-step instantaneous History output (make_sat_track_offline.py)'
        print(f'{path}: {len(idx)} soundings')


if __name__ == '__main__':
    main()
