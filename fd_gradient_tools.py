#!/usr/bin/env python3
"""
Finite-difference validation of the monthly gradient convention
(MONTHLY_SIGMA_STATUS.md outstanding item (a); doc/monthly_sigma_gradient.tex
Section 7).

Three subcommands, orchestrated by fd_validate_gradient.run:

  write-sigma   Write 12 constant monthly sigma files (base value), optionally
                adding eps at the single test cell recorded in the pick-cell
                state file.

  pick-cell     After the BASE adjoint run: read SurfaceFluxAdj_CO2 at the
                window's month boundaries, form
                    g_diff = A(t_1) - A(t_2)     (differencing convention)
                    g_raw  = A(t_1)              (raw convention)
                and pick the test cell maximising min(|g_diff|, |A(t_2)|),
                i.e. a cell with both a strong January gradient AND a strong
                difference between the two conventions (A(t_2) = g_raw−g_diff
                is exactly the later-months contamination the raw convention
                would add).  Saves everything to an .npz state file.

  analyze       Compare the model truth  dJ/dsigma = (J_pert - J_base)/eps
                against g_diff and g_raw at the test cell.  PASS (exit 0) if
                FD matches g_diff within --tol AND g_raw is clearly worse;
                exit 3 on mismatch, exit 4 if the cell cannot discriminate.

The window is expected to be two months (t_start, t_start+1M, t_end
boundaries), but any span >= 2 months works; the perturbed month is always
the FIRST month of the window.
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd

from write_sigma_monthly import make_latlon_grid, write_month, NMON
from lbfgsb_step_monthly import read_boundary_accumulators


# ---------------------------------------------------------------------------

def cmd_write_sigma(a):
    lats, lons = make_latlon_grid(a.nlat, a.nlon)
    sigma = np.full((NMON, a.nlat, a.nlon), a.base, dtype=np.float64)
    tag = f'constant {a.base}'
    if a.perturb_state:
        st = np.load(a.perturb_state)
        m, j, i = int(st['month']), int(st['jlat']), int(st['ilon'])
        sigma[m, j, i] += a.eps
        tag += (f' + eps={a.eps} at month={m + 1}, '
                f'cell (j={j}, i={i}) = ({lats[j]:.2f}N, {lons[i]:.2f}E)')
    for m in range(NMON):
        path = os.path.join(a.sigma_dir, f'{a.year}', f'{m + 1:02d}.nc4')
        write_month(sigma[m], lats, lons, path,
                    pd.Timestamp(f'{a.year}-{m + 1:02d}-01'))
    print(f'FD sigma written to {a.sigma_dir}: {tag}')


# ---------------------------------------------------------------------------

def cmd_pick_cell(a):
    t_start, t_end = pd.Timestamp(a.t_start), pd.Timestamp(a.t_end)
    bounds = pd.date_range(t_start, t_end, freq='MS')
    if len(bounds) < 3 or bounds[0] != t_start or bounds[-1] != t_end:
        raise SystemExit(f'ERROR: window {t_start} -> {t_end} must span >= 2 '
                         'whole months (both ends on month boundaries)')
    print(f'Reading SurfaceFluxAdj_CO2 at {len(bounds)} boundaries')
    A = read_boundary_accumulators(a.adj_output, bounds, a.nlat, a.nlon)

    a_end, a_start = np.abs(A[-1]).max(), np.abs(A[0]).max()
    if a_end > 1e-6 * max(a_start, 1e-30):
        raise SystemExit(f'ERROR: A(t_end) is not ~0 '
                         f'(|A(t_end)|max={a_end:.3e}, |A(t_start)|max={a_start:.3e})')

    g_diff = A[0] - A[1]    # January gradient, differencing convention
    g_raw  = A[0]           # January gradient, raw convention
    contam = A[1]           # = g_raw - g_diff: later-months contamination

    if a.cell:
        j, i = (int(v) for v in a.cell.split(','))
        print(f'Using user-specified cell (j={j}, i={i})')
    else:
        score = np.minimum(np.abs(g_diff), np.abs(contam))
        flat  = np.argsort(score.ravel())[::-1][:5]
        print('Top candidate cells by min(|g_diff|, |A(t_2)|):')
        lats, lons = make_latlon_grid(a.nlat, a.nlon)
        for r, f in enumerate(flat):
            jj, ii = np.unravel_index(f, score.shape)
            print(f'  {r + 1}. (j={jj:2d}, i={ii:2d}) '
                  f'({lats[jj]:7.2f}N, {lons[ii]:8.2f}E)  '
                  f'g_diff={g_diff[jj, ii]:+.4e}  g_raw={g_raw[jj, ii]:+.4e}  '
                  f'A(t_2)={contam[jj, ii]:+.4e}')
        j, i = np.unravel_index(flat[0], score.shape)

    print(f'Selected cell: (j={j}, i={i})  '
          f'g_diff={g_diff[j, i]:+.6e}  g_raw={g_raw[j, i]:+.6e}  '
          f'A(t_2)={contam[j, i]:+.6e}')
    np.savez(a.out,
             month=np.int64(0), jlat=np.int64(j), ilon=np.int64(i),
             g_diff_cell=np.float64(g_diff[j, i]),
             g_raw_cell=np.float64(g_raw[j, i]),
             contam_cell=np.float64(contam[j, i]),
             g_diff=g_diff, g_raw=g_raw, A=A,
             boundaries=np.array([str(b) for b in bounds]))
    print(f'FD state saved -> {a.out}')


# ---------------------------------------------------------------------------

def cmd_analyze(a):
    st = np.load(a.state)
    J0 = float(open(a.j_base).read().strip())
    J1 = float(open(a.j_pert).read().strip())
    fd = (J1 - J0) / a.eps
    gd = float(st['g_diff_cell'])
    gr = float(st['g_raw_cell'])
    j, i = int(st['jlat']), int(st['ilon'])

    def relerr(g):
        return abs(fd - g) / max(abs(fd), 1e-30)

    lines = [
        'Finite-difference gradient validation',
        '=====================================',
        f'cell (j={j}, i={i}), month 1, eps={a.eps}, base sigma={a.base}',
        f'J_base = {J0:.10e}',
        f'J_pert = {J1:.10e}',
        f'dJ     = {J1 - J0:+.6e}',
        '',
        f'FD truth   dJ/dsigma = {fd:+.10e}',
        f'g_diff  (A1 - A2)    = {gd:+.10e}   rel.err = {relerr(gd):.3%}',
        f'g_raw   (A1)         = {gr:+.10e}   rel.err = {relerr(gr):.3%}',
        '',
    ]
    diff_ok  = relerr(gd) < a.tol
    raw_bad  = relerr(gr) > 2 * a.tol
    if diff_ok and raw_bad:
        lines.append(f'PASS: FD matches the DIFFERENCING convention within '
                     f'{a.tol:.1%} and clearly rejects the raw convention.')
        rc = 0
    elif diff_ok and not raw_bad:
        lines.append(f'INCONCLUSIVE: FD matches g_diff, but g_raw is also '
                     f'within {2 * a.tol:.1%} — the cell cannot discriminate '
                     '(|A(t_2)| too small there). Re-run pick-cell/analyze '
                     'with a different --cell.')
        rc = 4
    else:
        lines.append(f'FAIL: FD does not match the differencing convention '
                     f'(rel.err {relerr(gd):.3%} >= {a.tol:.1%}). '
                     'The gradient the optimizer uses would be wrong — do not '
                     'start the production run.')
        rc = 3

    report = '\n'.join(lines)
    print(report)
    if a.report:
        with open(a.report, 'w') as fh:
            fh.write(report + '\n')
        print(f'Report written -> {a.report}')
    sys.exit(rc)


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    w = sub.add_parser('write-sigma')
    w.add_argument('--sigma-dir', required=True)
    w.add_argument('--base', type=float, default=1.0)
    w.add_argument('--year', type=int, default=2016)
    w.add_argument('--nlat', type=int, default=46)
    w.add_argument('--nlon', type=int, default=72)
    w.add_argument('--perturb-state', default=None,
                   help='fd_state.npz from pick-cell; adds eps at its cell')
    w.add_argument('--eps', type=float, default=0.01)
    w.set_defaults(func=cmd_write_sigma)

    c = sub.add_parser('pick-cell')
    c.add_argument('--adj-output', required=True)
    c.add_argument('--t-start', required=True)
    c.add_argument('--t-end', required=True)
    c.add_argument('--nlat', type=int, default=46)
    c.add_argument('--nlon', type=int, default=72)
    c.add_argument('--cell', default=None, help='override: "j,i"')
    c.add_argument('--out', required=True)
    c.set_defaults(func=cmd_pick_cell)

    z = sub.add_parser('analyze')
    z.add_argument('--state',  required=True)
    z.add_argument('--j-base', required=True)
    z.add_argument('--j-pert', required=True)
    z.add_argument('--eps',  type=float, default=0.01)
    z.add_argument('--base', type=float, default=1.0)
    z.add_argument('--tol',  type=float, default=0.05)
    z.add_argument('--report', default=None)
    z.set_defaults(func=cmd_analyze)

    a = p.parse_args()
    a.func(a)


if __name__ == '__main__':
    main()
