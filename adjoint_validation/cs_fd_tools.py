#!/usr/bin/env python3
"""
Native cubed-sphere adjoint-vs-forward gradient test (test_adjoint_vs_forward).

Purpose: validate the adjoint gradient against finite differences with NO
regridding anywhere in the loop.  The lat-lon FD validation (fd_multicell)
cannot separate a chain bug from the fact that the control lives on lat-lon
(model feels sigma_cs = R sigma_ll via ExtData) while the reported gradient
is S g_cs via the HISTORY regrid, with S != R^T.  Here BOTH operators are
removed:

  * sigma is written directly on the model's C{im} cubed-sphere grid in the
    GMAO "stacked" layout (lat = 6*im, lon = im).  MAPL's grid manager
    detects lat == 6*lon as a cubed-sphere file and ExtData reads it with an
    identity regrid — the model feels exactly the sigma in the file.
  * the adjoint HISTORY writes SurfaceFluxAdj_CO2 on the native grid
    (no grid_label), so the gradient is read back cell-for-cell.

Face convention: internal arrays are (nf=6, j, i); the stacked file row
f*im + j is face f row j (the layout MAPL itself uses for restarts).  The
P0 plumbing check validates this mapping empirically: a sigma field with a
UNIQUE value at every cell must reproduce itself exactly in
EmisCO2_NetTerrExch_pert / EmisCO2_NetTerrExch_base.

Subcommands (orchestrated by test_adjoint_vs_forward.run):
  write-sigma   write monthly sigma files on the CS grid (base, single-cell,
                global, or unique-per-cell P0 pattern); perturbations go in
                the FIRST month of the window only
  check-p0      verify pert == sigma * base elementwise on
                EmisCO2_NetTerrExch between the two 1-day P0 forwards
  pick-cells    after the base adjoint: read the native-CS gradient, pick
                land cells spread over |g|, emit the run list (cells.txt)
                and the state npz
  analyze       parabola-fit the exact g per pattern from the +/-eps runs
                (J is exactly quadratic in sigma), compare with the adjoint
                gradient, classify, write report + scatter + map plots
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

NF = 6
LAND_THR = 1e-12     # |A(t_start)| above this marks a TER-active (land) cell

# CS mosaic grid metadata carried by GCHP output; 'anchor' has a duplicate
# 'ncontact' dimension that breaks xarray operations.
_DROP_VARS = ['anchor', 'contacts', 'ncontact']


# ---------------------------------------------------------------------------
#  layout helpers
# ---------------------------------------------------------------------------

def to_faces(arr, im):
    """Normalize a 2-D CS record to (6, im, im) from either layout."""
    a = np.asarray(arr)
    if a.shape == (NF, im, im):
        return a
    if a.shape == (NF * im, im):
        return a.reshape(NF, im, im)
    raise ValueError(f'cannot interpret shape {a.shape} as C{im} cubed sphere')


def cell_tag(kind, f, j, i, eps):
    base = {'cell': f'f{f}j{j:02d}i{i:02d}', 'global': 'global'}[kind]
    return f'{base}_e{eps:+.4f}'


def expected_g(g, kind, f, j, i):
    if kind == 'cell':
        return float(g[f, j, i])
    if kind == 'global':
        return float(g.sum())
    raise ValueError(kind)


# ---------------------------------------------------------------------------
#  write-sigma
# ---------------------------------------------------------------------------

def write_month_cs(sigma_faces, im, path, t_stamp):
    """Write one single-record sigma file in the stacked CS layout."""
    stacked = np.asarray(sigma_faces, dtype=np.float64).reshape(NF * im, im)
    ds = xr.Dataset(
        {'Sigma_CO2': (['time', 'lat', 'lon'], stacked[None, :, :],
                       {'long_name': 'CO2 TER flux scaling factor (native CS)',
                        'units': '1'})},
        coords={
            'time': ('time', [t_stamp],
                     {'long_name': 'time'}),
            # fake index coordinates: MAPL identifies the cubed sphere from
            # lat == 6*lon and computes the geometry itself
            'lat': ('lat', np.arange(NF * im, dtype=np.float64),
                    {'long_name': 'Latitude', 'units': 'degrees_north'}),
            'lon': ('lon', np.arange(im, dtype=np.float64),
                    {'long_name': 'Longitude', 'units': 'degrees_east'}),
        })
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp.nc4'
    ds.to_netcdf(tmp, format='NETCDF4',
                 encoding={'Sigma_CO2': {'_FillValue': None},
                           'lat': {'_FillValue': None},
                           'lon': {'_FillValue': None}})
    os.replace(tmp, path)


def cmd_write_sigma(a):
    im = a.im
    sigma = np.full((a.n_months + 1, NF, im, im), a.base, dtype=np.float64)
    if a.kind == 'cell':
        f, j, i = (int(v) for v in a.cell.split(','))
        sigma[0, f, j, i] += a.eps
        tag = f'+ eps={a.eps:+.4f} at cell (f={f}, j={j}, i={i}), month 1'
    elif a.kind == 'global':
        sigma[0] += a.eps
        tag = f'+ eps={a.eps:+.4f} at ALL cells, month 1'
    elif a.kind == 'unique':
        ncell = NF * im * im
        idx = np.arange(1, ncell + 1, dtype=np.float64).reshape(NF, im, im)
        sigma[0] += a.unique_amp * idx / ncell
        tag = (f'+ UNIQUE per-cell pattern (amp={a.unique_amp}), month 1 '
               '(P0 face-mapping probe)')
    elif a.kind == 'none':
        tag = '(base only)'
    else:
        raise SystemExit(f'unknown --kind {a.kind}')

    # one extra trailing month so the refresh at t_end never misses a file
    for m in range(a.n_months + 1):
        mon = a.start_month + m
        yr, mon = a.year + (mon - 1) // 12, (mon - 1) % 12 + 1
        path = os.path.join(a.sigma_dir, f'{yr}', f'{mon:02d}.nc4')
        write_month_cs(sigma[m], im, path, pd.Timestamp(f'{yr}-{mon:02d}-01'))
    print(f'CS sigma written to {a.sigma_dir}: constant {a.base} {tag} '
          f'({a.n_months + 1} files from {a.year}-{a.start_month:02d})')


# ---------------------------------------------------------------------------
#  check-p0
# ---------------------------------------------------------------------------

def _load_emis(d, im):
    files = sorted(glob.glob(os.path.join(d, 'GEOSChem.Emissions.*.nc4')))
    if not files:
        raise SystemExit(f'ERROR: no Emissions files in {d}')
    ds = xr.concat([xr.open_dataset(f, drop_variables=_DROP_VARS)
                    for f in files], dim='time').sortby('time')
    da = ds['EmisCO2_NetTerrExch']
    recs = np.stack([to_faces(da.isel(time=k).values, im)
                     for k in range(da.sizes['time'])])
    return pd.DatetimeIndex(ds['time'].values).round('s'), recs


def cmd_check_p0(a):
    im = a.im
    tb, base = _load_emis(a.base_dir, im)
    tp, pert = _load_emis(a.pert_dir, im)
    if len(tb) != len(tp) or (tb != tp).any():
        raise SystemExit(f'ERROR: time axes differ between {a.base_dir} '
                         f'({len(tb)} recs) and {a.pert_dir} ({len(tp)} recs)')
    with xr.open_dataset(a.sigma_file) as ds:
        sigma = to_faces(ds['Sigma_CO2'].isel(time=0).values, im)

    # With >= 2 records, skip the first: a mid-run ExtData refresh applies one
    # timestep after its stamp, so a record straddling a switch-on would be
    # mixed.  A SINGLE record is fine as-is: sigma is constant over the run
    # and ExtData loads it at initialization, so every step already carries it
    # (verified empirically 2026-07-14: single daily-averaged record obeys
    # pert == sigma*base to 1.1e-7 at all 3456 cells).
    skipped = len(tb) >= 2
    if skipped:
        base, pert = base[1:], pert[1:]
    base, pert = base.astype(np.float64), pert.astype(np.float64)

    ref = np.abs(sigma[None] * base)
    scale = ref.max()
    if scale <= 0:
        raise SystemExit('ERROR: EmisCO2_NetTerrExch is identically zero — '
                         'nothing to check')
    err = np.abs(pert - sigma[None] * base) / scale
    moved = np.abs(pert - base).max() / max(np.abs(base).max(), 1e-300)

    print(f'check-p0: {base.shape[0]} records compared'
          + (' (first record skipped for the ExtData 1-step lag)'
             if skipped else ' (single record: sigma constant from init)'))
    print(f'  max |pert - sigma*base| / max|sigma*base| = {err.max():.3e} '
          f'(tol {a.tol:.1e})')
    print(f'  max |pert - base| / max|base|             = {moved:.3e} '
          '(must be >> 0: perturbation felt)')
    if err.max() < a.tol and moved > 1e-4:
        print('P0 PASS: ExtData read the CS sigma file natively (no regrid) '
              'and the face mapping is correct.')
        return
    if moved <= 1e-4:
        raise SystemExit('P0 FAIL: perturbed run does not differ from base — '
                         'sigma file not applied (check ExtData/HEMCO wiring)')
    kmax = np.unravel_index(np.argmax(err.max(axis=0)), err.shape[1:])
    raise SystemExit(f'P0 FAIL: pert != sigma*base (max rel err {err.max():.3e} '
                     f'at (f,j,i)={kmax}) — regrid in the loop or wrong '
                     'face/row mapping')


# ---------------------------------------------------------------------------
#  gradient reading (native CS)
# ---------------------------------------------------------------------------

def read_boundary_accumulators_cs(adj_output_dir, boundaries, im):
    """
    Read SurfaceFluxAdj_CO2 at each boundary time, exactly, on the native
    CS grid.  Mirrors lbfgsb_step_monthly.read_boundary_accumulators
    including the backwards-HISTORY quirk: the record at t_end is never
    written; A(t_end) = 0 by construction, so it is zero-filled.  A missing
    INTERIOR boundary is a hard error.

    Returns A : (len(boundaries), 6, im, im) float64
    """
    pattern = os.path.join(adj_output_dir, 'GEOSChem.Adjoint.*.nc4')
    files = sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError(f'No adjoint output files matching {pattern}')

    targets = [np.datetime64(pd.Timestamp(t)) for t in boundaries]
    A = np.full((len(targets), NF, im, im), np.nan, dtype=np.float64)
    src = [None] * len(targets)

    for f in files:
        with xr.open_dataset(f, drop_variables=_DROP_VARS) as ds:
            times = pd.DatetimeIndex(ds['time'].values).round('s')
            for k, t in enumerate(targets):
                hits = np.where(times == t)[0]
                if hits.size == 0:
                    continue
                if src[k] is not None:
                    print(f'  NOTE: boundary {t} found again in '
                          f'{os.path.basename(f)} (first: {src[k]}); '
                          'keeping the first read')
                    continue
                da = ds['SurfaceFluxAdj_CO2'].isel(time=int(hits[0]))
                A[k] = to_faces(da.values, im).astype(np.float64)
                src[k] = os.path.basename(f)

    if src[-1] is None:
        A[-1] = 0.0
        src[-1] = '(not written by backwards HISTORY; A(t_end)=0 by construction)'
        print(f'  NOTE: no adjoint record at t_end={targets[-1]} — '
              'expected for backwards HISTORY; using A(t_end)=0')

    missing = [str(targets[k]) for k in range(len(targets)) if src[k] is None]
    if missing:
        raise RuntimeError('No adjoint record at month boundary(ies): '
                           + ', '.join(missing))
    for k, t in enumerate(targets):
        print(f'  A({t}) <- {src[k]}')
    return A


def read_cs_coords(paths, im):
    """lats/lons cell centres (6, im, im) from the first file that has them."""
    for p in paths:
        if not p or not os.path.exists(p):
            continue
        with xr.open_dataset(p, drop_variables=_DROP_VARS) as ds:
            if 'lats' in ds and 'lons' in ds:
                lats = to_faces(np.squeeze(ds['lats'].values), im)
                lons = to_faces(np.squeeze(ds['lons'].values), im)
                return lats, lons
    raise RuntimeError('no file with native-CS "lats"/"lons" coordinate '
                       f'variables among: {paths}')


# ---------------------------------------------------------------------------
#  pick-cells
# ---------------------------------------------------------------------------

def cmd_pick_cells(a):
    im = a.im
    t_start, t_end = pd.Timestamp(a.t_start), pd.Timestamp(a.t_end)

    # Whole-month windows: use month-start boundaries so each month's
    # accumulated gradient can be read separately.  Sub-monthly windows
    # (e.g. a 5-day test) only have two meaningful timestamps: t_start
    # (where the backwards HISTORY writes the total gradient) and t_end
    # (never written; A=0 by construction).
    ms_bounds = pd.date_range(t_start, t_end, freq='MS')
    if len(ms_bounds) >= 2 and ms_bounds[0] == t_start and ms_bounds[-1] == t_end:
        bounds = ms_bounds
    else:
        bounds = pd.DatetimeIndex([t_start, t_end])
        print(f'Sub-monthly window {t_start} -> {t_end}: '
              'using [t_start, t_end] as accumulator boundaries')

    print(f'Reading SurfaceFluxAdj_CO2 at {len(bounds)} boundaries (native CS)')
    A = read_boundary_accumulators_cs(a.adj_output, bounds, im)
    g = A[0] - A[1]          # month-1 (or sub-monthly) gradient

    adj_files = sorted(glob.glob(os.path.join(a.adj_output,
                                              'GEOSChem.Adjoint.*.nc4')))
    lats, lons = read_cs_coords(adj_files[:1] + [a.coord_file], im)

    land = np.abs(A[0]) > LAND_THR

    # Optionally restrict to cells that are at least `face_margin` cells away
    # from every face edge on all four sides (0-indexed j,i in [m, im-1-m]).
    # Used to test whether face-boundary halo errors cause the adjoint mismatch.
    interior = np.ones((6, im, im), dtype=bool)
    m = getattr(a, 'face_margin', 0)
    if m > 0:
        mask = np.zeros((6, im, im), dtype=bool)
        mask[:, m:im - m, m:im - m] = True
        interior = mask
        n_int = int(mask.sum())
        print(f'face_margin={m}: restricting to {n_int} interior cells '
              f'(j,i in [{m},{im-1-m}]); '
              f'{6*im*im - n_int} edge-zone cells excluded')

    e = a.eps

    if getattr(a, 'fixed_cells', None):
        # Use cell indices from an existing cells.txt so two runs with different
        # physics test the adjoint at exactly the same grid points.
        fixed = read_cells(a.fixed_cells)
        seen, cell_order = set(), []
        for kind, f, j, i, _, _ in fixed:
            if kind == 'cell' and (f, j, i) not in seen:
                cell_order.append((int(f), int(j), int(i)))
                seen.add((f, j, i))
        print(f'Using {len(cell_order)} fixed cells from {a.fixed_cells}')
        entries = [('global', 0, 0, 0, +e), ('global', 0, 0, 0, -e)]
        for n, (f, j, i) in enumerate(cell_order):
            entries.append(('cell', f, j, i, +e))
            entries.append(('cell', f, j, i, -e))
            if n == 0:
                entries.append(('cell', f, j, i, +2 * e))
    else:
        cand = land & interior & (np.abs(g) >= a.gmin)
        ff, jj, ii = np.where(cand)
        if ff.size == 0:
            raise SystemExit(f'ERROR: no land cell with |g| >= {a.gmin}; '
                             f'max |g| on land = {np.abs(g)[land].max():.3e}')
        mags = np.abs(g[ff, jj, ii])

        # log-spaced |g| targets, greedily matched subject to a minimum
        # geographic separation (relaxed if necessary) — same idea as the
        # lat-lon multicell picker
        chosen, order = [], []
        targets = np.geomspace(mags.max(), max(a.gmin, mags.min()), a.n)
        for tgt in targets:
            by_tgt = np.argsort(np.abs(np.log(mags) - np.log(tgt)))
            pick = None
            for relax in (1.0, 0.5, 0.0):
                for k in by_tgt:
                    if k in order:
                        continue
                    dmin = min((np.hypot(lats[ff[k], jj[k], ii[k]] - lats[c],
                                         min(abs(lons[ff[k], jj[k], ii[k]] - lons[c]),
                                             360 - abs(lons[ff[k], jj[k], ii[k]]
                                                       - lons[c])))
                                for c in chosen), default=np.inf)
                    if dmin >= a.min_sep * relax:
                        pick = k
                        break
                if pick is not None:
                    break
            if pick is None:
                continue
            order.append(pick)
            chosen.append((ff[pick], jj[pick], ii[pick]))

        entries = [('global', 0, 0, 0, +e), ('global', 0, 0, 0, -e)]
        for n, k in enumerate(order):
            f, j, i = int(ff[k]), int(jj[k]), int(ii[k])
            entries.append(('cell', f, j, i, +e))
            entries.append(('cell', f, j, i, -e))
            if n == 0:                       # quadraticity probe at the top cell
                entries.append(('cell', f, j, i, +2 * e))

    lines = ['# native-CS FD run list: kind f j i eps tag']
    print(f'{"tag":24s} {"lat":>7s} {"lon":>8s} {"g_adj(expected)":>16s}')
    for kind, f, j, i, ee in entries:
        tag = cell_tag(kind, f, j, i, ee)
        lines.append(f'{kind} {f} {j} {i} {ee:+.4f} {tag}')
        gx = expected_g(g, kind, f, j, i)
        loc = ('   all', '    all') if kind == 'global' else \
              (f'{lats[f, j, i]:7.2f}', f'{lons[f, j, i]:8.2f}')
        print(f'{tag:24s} {loc[0]:>7s} {loc[1]:>8s} {gx:+16.4e}')

    with open(a.out, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    np.savez(a.state_out, g=g, A=A, lats=lats, lons=lons, im=np.int64(im),
             boundaries=np.array([str(b) for b in bounds]))
    print(f'{len(entries)} GCHP forwards -> {a.out}; state -> {a.state_out}')


# ---------------------------------------------------------------------------
#  analyze
# ---------------------------------------------------------------------------

def read_cells(path):
    entries = []
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        kind, f, j, i, e, tag = ln.split()
        entries.append((kind, int(f), int(j), int(i), float(e), tag))
    return entries


def fit_parabola(J0, Js):
    """J(e) = J0 + e*g + (e^2/2)*H through (0, J0) and the sampled points."""
    if len(Js) < 2:
        return None
    ee = np.array(sorted(Js))
    dJ = np.array([Js[e] - J0 for e in ee])
    M = np.column_stack([ee, 0.5 * ee ** 2])
    (g, H), *_ = np.linalg.lstsq(M, dJ, rcond=None)
    resid = dJ - M @ np.array([g, H])
    return float(g), float(H), float(np.abs(resid).max())


def cmd_analyze(a):
    st = np.load(a.state)
    g, lats, lons = st['g'], st['lats'], st['lons']
    J0 = float(open(a.j_base).read().strip())

    cells = {}
    for kind, f, j, i, e, tag in read_cells(a.cells):
        jf = os.path.join(a.jdir, f'J_{tag}.txt')
        if os.path.exists(jf):
            cells.setdefault((kind, f, j, i), {})[e] = \
                float(open(jf).read().strip())

    lines = ['Native cubed-sphere adjoint-vs-forward gradient validation',
             '=' * 60,
             'sigma control AND gradient on the native C24 grid — '
             'NO regrid operator anywhere in the loop.',
             f'J_base = {J0:.10e}   base sigma = {a.base}',
             '']
    rows = []
    for (kind, f, j, i), Js in sorted(cells.items()):
        g_adj = expected_g(g, kind, f, j, i)
        one_sided = {e: (Jp - J0) / e for e, Jp in Js.items()}
        fit = fit_parabola(J0, Js)
        g_fit, H = (fit[0], fit[1]) if fit else (None, None)
        rows.append((kind, f, j, i, g_adj, g_fit, H, one_sided))

        loc = 'GLOBAL uniform' if kind == 'global' else \
              f'({lats[f, j, i]:7.2f}N, {lons[f, j, i]:8.2f}E)'
        lines.append(f'{kind} (f={f}, j={j:2d}, i={i:2d}) {loc}')
        lbl = 'g_adj (A1-A2)' if kind == 'cell' else 'g_adj (summed)'
        lines.append(f'  {lbl:22s} = {g_adj:+.6e}')
        for e in sorted(one_sided):
            lines.append(f'  one-sided FD (eps={e:+.4f}) = {one_sided[e]:+.6e}')
        if fit:
            rel = abs(g_fit - g_adj) / max(abs(g_fit), 1e-30)
            lines.append(f'  FIT g (exact quadratic) = {g_fit:+.6e}   '
                         f'rel.err(g_adj) = {rel:.2%}   '
                         f'ratio g_adj/g_fit = {g_adj / g_fit:+.3f}')
            lines.append(f'  FIT H (curvature)       = {H:+.6e}')
            if len(Js) > 2:
                dJmax = max(abs(np.array(list(Js.values())) - J0).max(), 1e-30)
                ok = fit[2] < 1e-3 * dJmax
                lines.append(f'  quadraticity: max|resid| = {fit[2]:.3e} '
                             f'({"OK" if ok else "NON-QUADRATIC — check"})')
        lines.append('')

    def stats(kinds):
        have = [(r[5], r[4]) for r in rows if r[5] is not None and r[0] in kinds]
        if not have:
            return None
        fd = np.array([h[0] for h in have])
        ga = np.array([h[1] for h in have])
        return {'fd': fd, 'ga': ga,
                'rel': np.abs(fd - ga) / np.maximum(np.abs(fd), 1e-30),
                'ratio': ga / fd,
                'slope': float((ga * fd).sum() / (fd * fd).sum()),
                'corr': float(np.corrcoef(fd, ga)[0, 1]) if len(fd) > 1 else 1.0}

    sc = stats({'cell'})
    gl = stats({'global'})
    if sc:
        r = sc['ratio']
        lines += [f"Single cells (fitted g): n={len(sc['fd'])}   through-origin "
                  f"slope g_adj = {sc['slope']:+.4f} * g_fit   corr = {sc['corr']:+.4f}",
                  f'  ratio g_adj/g_fit: mean = {r.mean():+.4f}   std = {r.std():.4f}'
                  f'   range = [{r.min():+.3f}, {r.max():+.3f}]']
    if gl:
        lines += [f"global ratio g_adj/g_fit = {gl['ratio'][0]:+.4f}"]
    lines += ['']

    verdict, rc = 'INCOMPLETE: no pattern has two eps values yet.', 5
    if sc:
        r = sc['ratio']
        cells_ok = bool(np.all(sc['rel'] < a.tol))
        global_ok = bool(gl and np.all(gl['rel'] < a.tol))
        cells_const = r.std() / max(abs(r.mean()), 1e-30) < 0.15
        if cells_ok and (gl is None or global_ok):
            verdict = (f'PASS: on the NATIVE grid the adjoint gradient matches '
                       f"exact-FD within {a.tol:.0%} at all {len(sc['fd'])} "
                       'cells' + (' and globally' if global_ok else '') +
                       '. The adjoint chain is correct; any lat-lon FD '
                       'mismatch is the ExtData/HISTORY interpolation '
                       '(S != R^T) — move the control vector to the CS grid '
                       'or smooth increments (correlated B).')
            rc = 0
        elif cells_const:
            verdict = (f'(b) SYSTEMATIC factor: g_adj = {r.mean():+.3f} x FD '
                       'consistently across native-grid cells — with regrids '
                       'ruled out this is a REAL constant-factor bug in the '
                       'adjoint chain (units/dt/convention).')
            rc = 3
        else:
            verdict = ('(c) SPATIALLY VARYING mismatch on the native grid — '
                       'regrids ruled out: approximate adjoint physics '
                       '(advection limiter / convection) or chaotic noise. '
                       f"corr = {sc['corr']:+.3f}, slope = {sc['slope']:+.3f}."
                       + (' Global probe DOES match, so aggregate gradients '
                          'remain usable.' if global_ok else ''))
            rc = 4
    lines += ['VERDICT: ' + verdict]

    report = '\n'.join(lines)
    print(report)
    if a.report:
        with open(a.report, 'w') as fh:
            fh.write(report + '\n')
        print(f'Report written -> {a.report}')
    if a.plot_dir:
        os.makedirs(a.plot_dir, exist_ok=True)
        _plots(a.plot_dir, g, lats, lons, rows)
    sys.exit(rc)


def _plots(plot_dir, g, lats, lons, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # --- scatter: individual cells only, global excluded ---
    fig, ax = plt.subplots(figsize=(8, 8))
    for kind, f, j, i, g_adj, g_fit, H, one_sided in rows:
        if kind != 'cell':
            continue
        for e, g1 in one_sided.items():
            ax.plot(g1, g_adj, 'x', color='0.6', ms=6, zorder=2)
        if g_fit is None:
            continue
        ax.plot(g_fit, g_adj, marker='o', color='tab:red', ms=9,
                zorder=3, ls='none')
        ax.annotate(f'({f},{j},{i})', (g_fit, g_adj),
                    textcoords='offset points', xytext=(6, 4), fontsize=8)
    lo = min(ax.get_xlim()[0], ax.get_ylim()[0])
    hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([lo, hi], [lo, hi], 'k--', lw=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel('finite-difference dJ/dsigma (fitted exact; x = one-sided)')
    ax.set_ylabel('adjoint g = A(t1) - A(t2), native CS cell')
    ax.set_title('Native-CS adjoint gradient vs FD truth (no regrid) — '
                 'individual cells only, 1:1 dashed')
    ax.grid(alpha=0.3)
    fig.savefig(os.path.join(plot_dir, 'fd_scatter_cs.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # --- spatial map with clearly labelled FD cells ---
    try:
        import cartopy.crs as ccrs
        fig = plt.figure(figsize=(11, 5.5))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.coastlines(lw=0.6)
    except Exception:
        fig, ax = plt.subplots(figsize=(11, 5.5))
    vmax = np.abs(g).max()
    pm = ax.scatter(lons.ravel(), lats.ravel(), c=g.ravel(), s=14,
                    cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    cell_rows = [(kind, f, j, i, g_adj, g_fit, H, one_sided)
                 for (kind, f, j, i, g_adj, g_fit, H, one_sided) in rows
                 if kind == 'cell']
    for n, (kind, f, j, i, g_adj, g_fit, *_) in enumerate(cell_rows, start=1):
        lat, lon = float(lats[f, j, i]), float(lons[f, j, i])
        # black halo then yellow fill so the marker pops on any background
        ax.plot(lon, lat, marker='*', color='black', ms=18, zorder=4, ls='none')
        ax.plot(lon, lat, marker='*', color='yellow', ms=14, zorder=5, ls='none')
        ratio_str = (f'{g_adj/g_fit:+.2f}x' if g_fit else '')
        ax.annotate(f'#{n} ({f},{j},{i})\n{ratio_str}',
                    (lon, lat), textcoords='offset points',
                    xytext=(8, 4), fontsize=7, color='black',
                    bbox=dict(boxstyle='round,pad=0.2', fc='white',
                              alpha=0.7, ec='none'),
                    zorder=6)
    plt.colorbar(pm, ax=ax, shrink=0.8, label='g = dJ/dsigma (native CS cell)')
    ax.set_title('Month-1 adjoint sensitivity on the native C24 grid '
                 'with FD test cells (★ = perturbed cell, label = g_adj/g_FD ratio)')
    fig.savefig(os.path.join(plot_dir, 'g_map_cs.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)
    print(f'Plots written -> {plot_dir}/fd_scatter_cs.png, g_map_cs.png')


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    w = sub.add_parser('write-sigma')
    w.add_argument('--sigma-dir', required=True)
    w.add_argument('--base', type=float, default=1.0)
    w.add_argument('--year', type=int, default=2016)
    w.add_argument('--start-month', type=int, default=1)
    w.add_argument('--n-months', type=int, default=1)
    w.add_argument('--im', type=int, default=24)
    w.add_argument('--kind', default='none',
                   choices=['none', 'cell', 'global', 'unique'])
    w.add_argument('--cell', default=None, help='"f,j,i" for --kind cell')
    w.add_argument('--eps', type=float, default=0.0)
    w.add_argument('--unique-amp', type=float, default=0.5)
    w.set_defaults(func=cmd_write_sigma)

    c = sub.add_parser('check-p0')
    c.add_argument('--base-dir', required=True)
    c.add_argument('--pert-dir', required=True)
    c.add_argument('--sigma-file', required=True)
    c.add_argument('--im', type=int, default=24)
    c.add_argument('--tol', type=float, default=3e-5)
    c.set_defaults(func=cmd_check_p0)

    k = sub.add_parser('pick-cells')
    k.add_argument('--adj-output', required=True)
    k.add_argument('--t-start', required=True)
    k.add_argument('--t-end', required=True)
    k.add_argument('--im', type=int, default=24)
    k.add_argument('--n', type=int, default=5)
    k.add_argument('--gmin', type=float, default=30.0)
    k.add_argument('--min-sep', type=float, default=20.0)
    k.add_argument('--eps', type=float, default=0.05)
    k.add_argument('--face-margin', type=int, default=0,
                   help='exclude cells within this many cells of any face edge '
                        '(0=all cells; 4 tests interior-only hypothesis)')
    k.add_argument('--coord-file', default=None,
                   help='fallback file with native-CS lats/lons (e.g. a P0 '
                        'Emissions file) if the adjoint output has none')
    k.add_argument('--fixed-cells', default=None,
                   help='cells.txt from another run: use its (f,j,i) indices '
                        'instead of picking dynamically, so two runs with '
                        'different physics test the same grid points')
    k.add_argument('--out', required=True)
    k.add_argument('--state-out', required=True)
    k.set_defaults(func=cmd_pick_cells)

    z = sub.add_parser('analyze')
    z.add_argument('--state', required=True)
    z.add_argument('--cells', required=True)
    z.add_argument('--j-base', required=True)
    z.add_argument('--jdir', required=True)
    z.add_argument('--base', type=float, default=1.0)
    z.add_argument('--tol', type=float, default=0.05)
    z.add_argument('--report', default=None)
    z.add_argument('--plot-dir', default=None)
    z.set_defaults(func=cmd_analyze)

    a = p.parse_args()
    a.func(a)


if __name__ == '__main__':
    main()
