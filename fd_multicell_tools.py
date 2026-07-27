#!/usr/bin/env python3
"""
Multi-cell, two-batch finite-difference validation of the monthly adjoint
gradient (follow-up to fd_gradient_tools.py after the 2026-07-14 one-sided
FD FAIL, job 24849168: g_adj/FD = 2.91 at one cell).

Key insight driving the design: CO2 is a passive linear tracer, so J_obs is
EXACTLY quadratic in sigma:  J(sigma + e*p) = J0 + e*g_p + (e^2/2) H_p  for a
perturbation pattern p.  Therefore TWO one-sided FDs at different eps solve
the parabola exactly:
    FD(e) = (J(e) - J0)/e = g + (e/2) H
    => g and H from any two distinct eps values (no symmetric pair needed).
The 24849168 mismatch had exactly the sign of a curvature bias (FD less
negative than g_adj), so the two-batch fit discriminates between
  (a) adjoint correct, one-sided FD curvature-biased,
  (b) constant systematic factor (units/dt/convention bug),
  (c) spatially varying error (approximate advection/convection adjoints),
  (d) interpolation/footprint artifact (see below).

Perturbation patterns:
  cell     single lat-lon cell (regrid-SENSITIVE: the control lives on
           lat-lon, the model feels sigma_cs = R sigma_ll via ExtData, and
           the reported gradient is S g_cs via the HISTORY lat-lon regrid;
           S != R^T, and a one-cell spike is the worst case)
  block3   3x3 block around the original cell (mostly regrid-immune)
  global   uniform over all cells (fully regrid-immune: constants regrid
           exactly) — separates an interpolation artifact from a true bug.

Batches: every pattern is run at eps1 (+1% default) AND eps2 (-5% default);
the scatter plot overlays both batches plus the parabola-fitted exact g.

Subcommands (orchestrated by fd_validate_multicell.run):
  pick-cells        select patterns/cells, emit cells.txt run list
  write-sigma-cell  write the 12 monthly sigma files for one entry
  analyze           parabola-fit g per pattern, compare with g_diff,
                    classify (a)-(d), write report + scatter + map plots
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd

from write_sigma_monthly import make_latlon_grid, write_month, NMON

LAND_THR = 1e-12          # |A(t_start)| above this marks a TER-active (land) cell


def cell_tag(j, i, eps, kind='cell'):
    base = {'cell': f'j{j:02d}_i{i:02d}', 'block3': f'b3j{j:02d}_i{i:02d}',
            'global': 'global'}[kind]
    return f'{base}_e{eps:+.4f}'


def block3_cells(j, i, nlat, nlon):
    """(j, i) pairs of the 3x3 block centred on (j, i), lon wrapped."""
    return [(jj, (i + di) % nlon)
            for jj in range(max(j - 1, 0), min(j + 2, nlat))
            for di in (-1, 0, 1)]


def expected_g(g_diff, kind, j, i):
    """Adjoint prediction of dJ/d(eps) for the given perturbation pattern."""
    if kind == 'cell':
        return float(g_diff[j, i])
    if kind == 'block3':
        return float(sum(g_diff[jj, ii]
                         for jj, ii in block3_cells(j, i, *g_diff.shape)))
    if kind == 'global':
        return float(g_diff.sum())
    raise ValueError(kind)


# ---------------------------------------------------------------------------

def cmd_pick_cells(a):
    st     = np.load(a.state)
    g_diff = st['g_diff']
    A1     = st['A'][0]
    j0, i0 = int(st['jlat']), int(st['ilon'])
    lats, lons = make_latlon_grid(*g_diff.shape)

    land = np.abs(A1) > LAND_THR
    cand = land & (np.abs(g_diff) >= a.gmin)
    cand[j0, i0] = False                       # original cell handled separately
    jj, ii = np.where(cand)
    mags   = np.abs(g_diff[jj, ii])

    # log-spaced |g| targets between gmin and the max candidate magnitude,
    # greedily matched to the nearest-|g| candidate that keeps >= min-sep
    # degrees from all previously chosen cells (geographic spread).
    chosen  = [(j0, i0)]
    targets = np.geomspace(mags.max(), a.gmin, a.n)
    order   = []
    for tgt in targets:
        by_tgt = np.argsort(np.abs(np.log(mags) - np.log(tgt)))
        pick   = None
        for relax in (1.0, 0.5, 0.0):          # relax separation if needed
            for k in by_tgt:
                if (jj[k], ii[k]) in chosen or k in order:
                    continue
                dmin = min(np.hypot(lats[jj[k]] - lats[cj],
                                    min(abs(lons[ii[k]] - lons[ci]),
                                        360 - abs(lons[ii[k]] - lons[ci])))
                           for cj, ci in chosen)
                if dmin >= a.min_sep * relax:
                    pick = k
                    break
            if pick is not None:
                break
        if pick is None:
            continue
        order.append(pick)
        chosen.append((jj[pick], ii[pick]))

    e1, e2 = a.eps, a.eps2
    lines = ['# multi-cell FD run list: kind j i eps tag',
             f'# seed_tag={cell_tag(j0, i0, e1)}   (J from fd_validation/J_pert.txt)']
    entries = [('cell',   j0, i0, e1),         # batch 1, seeded from 24849168
               ('cell',   j0, i0, e2),         # batch 2 -> exact g at the cell
               ('cell',   j0, i0, 2 * e1),     # 4-point parabola / quadraticity
               ('global', 0,  0,  e1),         # regrid-immune, decisive early
               ('global', 0,  0,  e2),
               ('block3', j0, i0, e1),
               ('block3', j0, i0, e2)]
    for k in order:
        entries.append(('cell', int(jj[k]), int(ii[k]), e1))
        entries.append(('cell', int(jj[k]), int(ii[k]), e2))

    print(f'{"tag":26s} {"lat":>7s} {"lon":>8s} {"g_adj(expected)":>16s}')
    for kind, j, i, e in entries:
        tag = cell_tag(j, i, e, kind)
        lines.append(f'{kind} {j} {i} {e:+.4f} {tag}')
        g = expected_g(g_diff, kind, j, i)
        loc = ('   all  ', '   all  ') if kind == 'global' else \
              (f'{lats[j]:7.2f}', f'{lons[i]:8.2f}')
        print(f'{tag:26s} {loc[0]:>7s} {loc[1]:>8s} {g:+16.4e}')

    with open(a.out, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    n_runs = len(entries) - 1                  # first entry is seeded
    print(f'{len(entries)} entries ({n_runs} GCHP forwards needed) -> {a.out}')


# ---------------------------------------------------------------------------

def cmd_write_sigma_cell(a):
    lats, lons = make_latlon_grid(a.nlat, a.nlon)
    sigma = np.full((NMON, a.nlat, a.nlon), a.base, dtype=np.float64)
    j, i  = (int(v) for v in a.cell.split(','))
    if a.kind == 'cell':
        sigma[0, j, i] += a.eps
        where = f'cell (j={j}, i={i}) = ({lats[j]:.2f}N, {lons[i]:.2f}E)'
    elif a.kind == 'block3':
        for jj, ii in block3_cells(j, i, a.nlat, a.nlon):
            sigma[0, jj, ii] += a.eps
        where = f'3x3 block around (j={j}, i={i}) = ({lats[j]:.2f}N, {lons[i]:.2f}E)'
    elif a.kind == 'global':
        sigma[0, :, :] += a.eps
        where = 'ALL cells (global uniform)'
    else:
        raise SystemExit(f'unknown --kind {a.kind}')
    for m in range(NMON):
        path = os.path.join(a.sigma_dir, f'{a.year}', f'{m + 1:02d}.nc4')
        write_month(sigma[m], lats, lons, path,
                    pd.Timestamp(f'{a.year}-{m + 1:02d}-01'))
    print(f'FD sigma written to {a.sigma_dir}: constant {a.base} '
          f'+ eps={a.eps:+.4f} at month=1, {where}')


# ---------------------------------------------------------------------------

def read_cells(path):
    entries = []
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        kind, j, i, e, tag = ln.split()
        entries.append((kind, int(j), int(i), float(e), tag))
    return entries


def fit_parabola(J0, Js):
    """
    Fit J(e) = J0 + e*g + (e^2/2)*H through (0, J0) and the sampled points.
    Exact for 2 eps values; overdetermined (quadraticity residual) for >2.
    Returns (g, H, max_residual) or None if fewer than 2 eps values.
    """
    if len(Js) < 2:
        return None
    ee = np.array(sorted(Js))
    dJ = np.array([Js[e] - J0 for e in ee])
    M  = np.column_stack([ee, 0.5 * ee ** 2])
    (g, H), *_ = np.linalg.lstsq(M, dJ, rcond=None)
    resid = dJ - M @ np.array([g, H])
    return float(g), float(H), float(np.abs(resid).max())


def cmd_analyze(a):
    st     = np.load(a.state)
    g_diff = st['g_diff']
    lats, lons = make_latlon_grid(*g_diff.shape)
    J0 = float(open(a.j_base).read().strip())

    # gather available J values per perturbation pattern: {(kind,j,i): {eps: J}}
    cells = {}
    for kind, j, i, e, tag in read_cells(a.cells):
        jf = os.path.join(a.jdir, f'J_{tag}.txt')
        if os.path.exists(jf):
            cells.setdefault((kind, j, i), {})[e] = float(open(jf).read().strip())

    lines = ['Multi-cell two-batch finite-difference gradient validation',
             '=' * 60,
             f'J_base = {J0:.10e}   base sigma = {a.base}',
             '']
    rows = []  # (kind, j, i, g_adj, g_fit or None, H or None, one_sided, Js)
    for (kind, j, i), Js in sorted(cells.items()):
        g_adj = expected_g(g_diff, kind, j, i)
        one_sided = {e: (Jp - J0) / e for e, Jp in Js.items()}
        fit = fit_parabola(J0, Js)
        g_fit, H = (fit[0], fit[1]) if fit else (None, None)
        rows.append((kind, j, i, g_adj, g_fit, H, one_sided, Js))

        loc = 'GLOBAL uniform' if kind == 'global' else \
              f'({lats[j]:7.2f}N, {lons[i]:8.2f}E)'
        lines.append(f'{kind} (j={j:2d}, i={i:2d}) {loc}')
        lbl = 'g_adj (A1-A2)' if kind == 'cell' else 'g_adj (A1-A2, summed)'
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
                lines.append(f'  quadraticity: max|resid| = {fit[2]:.3e} '
                             f'({"OK" if fit[2] < 1e-3 * max(abs(np.array(list(Js.values())) - J0).max(), 1e-30) else "NON-QUADRATIC — check"})')
        lines.append('')

    def stats(kinds):
        have = [(r[4], r[3]) for r in rows if r[4] is not None and r[0] in kinds]
        if not have:
            return None
        fd = np.array([h[0] for h in have]); ga = np.array([h[1] for h in have])
        return {'fd': fd, 'ga': ga,
                'rel':   np.abs(fd - ga) / np.maximum(np.abs(fd), 1e-30),
                'ratio': ga / fd,
                'slope': float((ga * fd).sum() / (fd * fd).sum()),
                'corr':  float(np.corrcoef(fd, ga)[0, 1]) if len(fd) > 1 else 1.0}

    sc  = stats({'cell'})                      # single cells: regrid-sensitive
    im  = stats({'block3', 'global'})          # block/global: regrid-immune
    verdict = 'INCOMPLETE: no pattern has two eps values yet.'
    rc = 5
    if sc:
        r = sc['ratio']
        lines += [f"Single cells (fitted g): n={len(sc['fd'])}   through-origin "
                  f"slope g_adj = {sc['slope']:+.4f} * g_fit   corr = {sc['corr']:+.4f}",
                  f'  ratio g_adj/g_fit: mean = {r.mean():+.4f}   std = {r.std():.4f}'
                  f'   range = [{r.min():+.3f}, {r.max():+.3f}]']
    if im:
        for row in rows:
            if row[0] in ('block3', 'global') and row[4] is not None:
                lines += [f'  {row[0]:6s} ratio g_adj/g_fit = {row[3] / row[4]:+.4f}']
    lines += ['']
    if sc:
        r = sc['ratio']
        cells_ok    = bool(np.all(sc['rel'] < a.tol))
        cells_const = r.std() / max(abs(r.mean()), 1e-30) < 0.15
        im_ok = bool(im and np.all(im['rel'] < a.tol))
        if cells_ok:
            verdict = (f'(a) ADJOINT VALIDATED: fitted-exact FD matches g_diff '
                       f"within {a.tol:.0%} at all {len(sc['fd'])} cells"
                       + (' AND on block/global probes' if im_ok else '')
                       + '. The 2026-07-14 one-sided FAIL was curvature bias '
                       '(J exactly quadratic in sigma). Fix '
                       'fd_validate_gradient.run to fit the parabola (two eps).')
            rc = 0
        elif im_ok and not cells_ok:
            verdict = ('(d) INTERPOLATION/FOOTPRINT artifact: block/global '
                       '(regrid-immune) probes match FD but single cells do not '
                       '— the HISTORY lat-lon regrid S is not the transpose of '
                       'the ExtData regrid R. Aggregate/smoothed gradients are '
                       'trustworthy; single-cell values are not. Consider a '
                       'correlated-B / smoothing transform (already planned) '
                       'or computing g on the native CS grid.')
            rc = 6
        elif cells_const:
            verdict = (f'(b) SYSTEMATIC factor: g_adj = {r.mean():+.3f} x FD '
                       'consistently across cells'
                       + ('' if im_ok else ' (and block/global also off)')
                       + ' — hunt a constant-factor convention/units bug in '
                       'the adjoint chain.')
            rc = 3
        else:
            verdict = ('(c) SPATIALLY VARYING mismatch: no constant factor — '
                       'approximate adjoint physics (advection limiters / '
                       f"convection) or noise. corr = {sc['corr']:+.3f}, "
                       f"slope = {sc['slope']:+.3f}."
                       + (' Block/global probes DO match FD, so aggregate '
                          'gradients remain usable.' if im_ok else ''))
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
        _plots(a.plot_dir, g_diff, lats, lons, rows, a.eps, a.eps2)
    sys.exit(rc)


def _plots(plot_dir, g_diff, lats, lons, rows, eps1, eps2):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # -- scatter: adjoint vs FD, both batches overlaid ----------------------
    batch = {eps1: dict(marker='o', mfc='none', color='tab:orange',
                        label=f'batch 1: one-sided FD (eps={eps1:+.2%})'),
             eps2: dict(marker='^', mfc='none', color='tab:purple',
                        label=f'batch 2: one-sided FD (eps={eps2:+.2%})')}
    fitsty = {'cell':   dict(marker='o', color='tab:red',   label='fitted exact g (cell)'),
              'block3': dict(marker='s', color='tab:blue',  label='fitted exact g (3x3 block)'),
              'global': dict(marker='D', color='tab:green', label='fitted exact g (global)')}

    def draw(ax, sub, legend):
        seen = set()
        for kind, j, i, g_adj, g_fit, H, one_sided, Js in sub:
            for e, g1 in one_sided.items():
                bs = batch.get(e)
                if bs is None:                 # e.g. the +2*eps1 parabola point
                    ax.plot(g1, g_adj, 'x', color='0.6', ms=6, zorder=2)
                    continue
                key = ('b', e)
                ax.plot(g1, g_adj, bs['marker'], mfc=bs['mfc'], color=bs['color'],
                        ms=8, zorder=2,
                        label=bs['label'] if legend and key not in seen else None)
                seen.add(key)
            if g_fit is not None:
                fs = fitsty[kind]
                ax.plot(g_fit, g_adj, fs['marker'], color=fs['color'], ms=9,
                        zorder=3,
                        label=fs['label'] if legend and kind not in seen else None)
                seen.add(kind)
                name = kind if kind == 'global' else f'({j},{i})'
                ax.annotate(name, (g_fit, g_adj), textcoords='offset points',
                            xytext=(6, 4), fontsize=8)
        lo = min(ax.get_xlim()[0], ax.get_ylim()[0])
        hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lo, hi], [lo, hi], 'k--', lw=1, label='1:1' if legend else None)
        have = [(r[4], r[3]) for r in sub if r[4] is not None and r[0] == 'cell']
        if have:
            fd = np.array([h[0] for h in have]); ga = np.array([h[1] for h in have])
            s  = (ga * fd).sum() / (fd * fd).sum()
            ax.plot([lo, hi], [s * lo, s * hi], 'r-', lw=1,
                    label=f'through-origin fit (cells): slope {s:+.3f}'
                          if legend else None)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel('finite-difference dJ/dsigma')
        ax.set_ylabel('adjoint g_diff = A(t1) - A(t2)')
        ax.grid(alpha=0.3)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 7))
    draw(axL, rows, legend=True)
    axL.set_title('all patterns (incl. block / global)')
    draw(axR, [r for r in rows if r[0] == 'cell'], legend=False)
    axR.set_title('zoom: single cells only')
    axL.legend(fontsize=8, loc='best')
    fig.suptitle('Monthly-sigma adjoint gradient vs FD truth  '
                 f'(Jan, 2-month window; batches eps={eps1:+.2%} and {eps2:+.2%})')
    fig.savefig(os.path.join(plot_dir, 'fd_scatter.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)

    # -- map: g_diff with test cells marked ---------------------------------
    try:
        import cartopy.crs as ccrs
        fig = plt.figure(figsize=(11, 5.5))
        ax  = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.coastlines(lw=0.6)
    except Exception:
        fig, ax = plt.subplots(figsize=(11, 5.5))
    vmax = np.abs(g_diff).max()
    pm = ax.pcolormesh(lons, lats, g_diff, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    for kind, j, i, g_adj, g_fit, H, one_sided, Js in rows:
        if kind != 'cell':
            continue
        ax.plot(lons[i], lats[j], 'kx', ms=9, mew=2)
        ax.annotate(f'({j},{i})', (lons[i], lats[j]),
                    textcoords='offset points', xytext=(5, 5), fontsize=8)
    plt.colorbar(pm, ax=ax, shrink=0.8,
                 label='g_diff = A(t1) - A(t2)  [dJ/dsigma_Jan]')
    ax.set_title('January adjoint sensitivity dJ/dsigma (differencing convention) '
                 'with FD test cells')
    fig.savefig(os.path.join(plot_dir, 'g_diff_map.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)
    print(f'Plots written -> {plot_dir}/fd_scatter.png, g_diff_map.png')


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    c = sub.add_parser('pick-cells')
    c.add_argument('--state', required=True, help='fd_state.npz from pick-cell')
    c.add_argument('--n', type=int, default=4, help='extra cells beyond the original')
    c.add_argument('--gmin', type=float, default=30.0, help='min |g_diff|')
    c.add_argument('--min-sep', type=float, default=15.0, help='degrees between cells')
    c.add_argument('--eps',  type=float, default=0.01,  help='batch-1 eps')
    c.add_argument('--eps2', type=float, default=-0.05, help='batch-2 eps')
    c.add_argument('--out', required=True)
    c.set_defaults(func=cmd_pick_cells)

    w = sub.add_parser('write-sigma-cell')
    w.add_argument('--sigma-dir', required=True)
    w.add_argument('--base', type=float, default=1.0)
    w.add_argument('--year', type=int, default=2016)
    w.add_argument('--nlat', type=int, default=46)
    w.add_argument('--nlon', type=int, default=72)
    w.add_argument('--cell', required=True, help='"j,i" (ignored for --kind global)')
    w.add_argument('--eps', type=float, required=True)
    w.add_argument('--kind', default='cell', choices=['cell', 'block3', 'global'])
    w.set_defaults(func=cmd_write_sigma_cell)

    z = sub.add_parser('analyze')
    z.add_argument('--state',  required=True)
    z.add_argument('--cells',  required=True)
    z.add_argument('--j-base', required=True)
    z.add_argument('--jdir',   required=True)
    z.add_argument('--base', type=float, default=1.0)
    z.add_argument('--eps',  type=float, default=0.01,  help='batch-1 eps (plot label)')
    z.add_argument('--eps2', type=float, default=-0.05, help='batch-2 eps (plot label)')
    z.add_argument('--tol',  type=float, default=0.05)
    z.add_argument('--report',   default=None)
    z.add_argument('--plot-dir', default=None)
    z.set_defaults(func=cmd_analyze)

    a = p.parse_args()
    a.func(a)


if __name__ == '__main__':
    main()
