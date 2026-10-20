#!/usr/bin/env python3
"""Compare the FD gradient of run_fd.pbs with the 4D-Var adjoint gradient at sigma = 1.

The adjoint gradient is g_prev of the experiment's iteration-1 state
(exps/<EXP>/4dvar_output_osse/4dvar_state_iter001.npz): at sigma = 1 the background
gradient is 0, so g_prev = dJ_obs/dsigma.  FD = (J+ - J-) / (2 eps).

usage: compare_fd.py <EXP> <cells file> [--eps 0.05]
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if len(args) != 2:
        sys.exit(__doc__)
    exp, cells = args
    eps = float(sys.argv[sys.argv.index('--eps') + 1]) if '--eps' in sys.argv else 0.05
    st = np.load(os.path.join(ROOT, 'exps', exp, '4dvar_output_osse', '4dvar_state_iter001.npz'))
    g = st['g_prev'].reshape(6, 24, 24)
    if abs(float(st['J_b_prev'])) > 1e-12:
        sys.exit('iteration-1 state has J_b != 0: not evaluated at sigma = 1')
    res = os.path.join(ROOT, 'fd_check', exp, 'results')
    rows = []
    for line in open(cells):
        t = line.split()
        if not t or t[0].startswith('#'):
            continue
        tag, f, j, i = t[0], int(t[1]), int(t[2]), int(t[3])
        pp, pm = (os.path.join(res, f'J_{tag}_{s}.txt') for s in 'pm')
        if not (os.path.exists(pp) and os.path.exists(pm)):
            rows.append((tag, f, j, i, g[f, j, i], np.nan))
            continue
        fd = (float(open(pp).read()) - float(open(pm).read())) / (2 * eps)
        rows.append((tag, f, j, i, g[f, j, i], fd))
    out = [f'FD check {exp}: adjoint dJ_obs/dsigma (iteration 1, sigma = 1) vs central FD, eps {eps}',
           f'J_obs(sigma=1) = {float(st["J_obs_prev"]):.6e}', '',
           f'{"tag":14s} {"f":>2s} {"j":>3s} {"i":>3s} {"adjoint":>13s} {"FD":>13s} {"adj/FD":>8s}']
    done = [r for r in rows if np.isfinite(r[5])]
    for tag, f, j, i, a, d in rows:
        ratio = a / d if np.isfinite(d) and d != 0 else np.nan
        out.append(f'{tag:14s} {f:2d} {j:3d} {i:3d} {a:13.5e} {d:13.5e} {ratio:8.4f}')
    if done:
        a = np.array([r[4] for r in done]); d = np.array([r[5] for r in done])
        slope = float((a * d).sum() / (d * d).sum())
        out += ['', f'cells done {len(done)}/{len(rows)}  slope adj vs FD {slope:.4f}  '
                    f'corr {np.corrcoef(a, d)[0, 1] if len(done) > 2 else np.nan:.4f}  '
                    f'weighted error sum|adj-FD|/sum|FD| {np.abs(a - d).sum() / np.abs(d).sum():.2%}']
    text = '\n'.join(out)
    print(text)
    with open(os.path.join(ROOT, 'fd_check', exp, 'report.txt'), 'w') as fh:
        fh.write(text + '\n')


if __name__ == '__main__':
    main()
