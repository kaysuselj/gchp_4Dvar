#!/usr/bin/env python3
"""Taylor test of the 4D-Var gradient (taylor_test.pbs, EXPERIMENT_MATRIX G1).

point   : state file with x_next = x0 + eps * d for the sigma writers (write_sigma[_monthly].py)
analyze : J_obs(x0 + eps d) - J_obs(x0) against eps g^T d for all eps directories of an experiment

x0, J_obs(x0) and g = dJ_obs/dsigma come from the iteration-1 state of the optimizer (x_prev, J_obs_prev,
g_prev; at sigma = 1 the background term and its gradient are 0, which is checked). Directions d:
  grad       d = -g / max|g|          (largest sigma change = eps)
  month:<k>  monthly control: d = -g in month k (1-based), 0 elsewhere, / max|.|
  rand:<s>   d = N(0,1) field (seed s) / max|.|
"""
import argparse, glob, os, re
import numpy as np


def load(state):
    s = np.load(state)
    x0, g, j0 = s['x_prev'], s['g_prev'], float(s['J_obs_prev'])
    if not np.allclose(x0, 1.0) or float(s['J_b_prev']) != 0.0:
        raise SystemExit(f'{state}: not an iteration-1 state (sigma != 1 or J_b != 0); the Taylor test uses g_obs = g')
    return x0, g, j0


def direction(g, dirn, nmon, monthly):
    if dirn == 'grad':
        d = -g.copy()
    elif dirn.startswith('month'):
        if not monthly:
            raise SystemExit('month:<k> needs a monthly control')
        k = int(re.sub(r'\D', '', dirn)) - 1
        d = np.zeros_like(g)
        blk = g.size // nmon
        d[k * blk:(k + 1) * blk] = -g[k * blk:(k + 1) * blk]
    elif dirn.startswith('rand'):
        d = np.random.default_rng(int(re.sub(r'\D', '', dirn) or 1)).standard_normal(g.size)
    else:
        raise SystemExit(f'unknown direction {dirn}')
    return d / np.abs(d).max()


def cmd_point(a):
    x0, g, _ = load(a.state)
    d = direction(g, a.dirn, a.nmon, a.monthly)
    x = x0 + a.eps * d
    np.savez(a.out, x_next=x, iteration=np.int64(0))
    print(f'x = x0 + {a.eps:g} d ({a.dirn}): sigma {x.min():.6f} .. {x.max():.6f}')


def cmd_analyze(a):
    x0, g, j0 = load(a.state)
    out = [f'Taylor test {a.exp}', f'J_obs(x0) = {j0:.6f}   (iteration-1 state {a.state})', '']
    rows = {}
    for w in sorted(glob.glob(os.path.join(a.exp, 'taylor', '*_*'))):
        jf = os.path.join(w, 'J_value.txt')
        if not os.path.isfile(jf):
            continue
        dirn, eps = os.path.basename(w).rsplit('_', 1)
        dirn = re.sub(r'^month(\d)', r'month:\1', re.sub(r'^rand(\d)', r'rand:\1', dirn))
        rows.setdefault(dirn, []).append((float(eps), float(open(jf).read().split()[0])))
    for dirn, lst in rows.items():
        d = direction(g, dirn, a.nmon, a.monthly)
        gd = float(g @ d)
        out.append(f'--- direction {dirn}:  g^T d = {gd:.6e}')
        out.append(f'{"eps":>10s} {"J_obs(x0+eps d)":>18s} {"dJ":>14s} {"eps g^T d":>14s} {"ratio":>10s} {"resid/eps^2":>14s}')
        for eps, j in sorted(lst, reverse=True):
            dj, lin = j - j0, eps * gd
            out.append(f'{eps:10.2e} {j:18.6f} {dj:14.6e} {lin:14.6e} {dj / lin:10.6f} {(dj - lin) / eps ** 2:14.6e}')
        out.append('')
    out.append('ratio = dJ / (eps g^T d) -> 1 for a correct gradient; resid/eps^2 constant = J_obs quadratic in sigma')
    print('\n'.join(out))


p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = p.add_subparsers(dest='cmd', required=True)
q = sub.add_parser('point'); q.add_argument('--state', required=True); q.add_argument('--dirn', default='grad')
q.add_argument('--eps', type=float, required=True); q.add_argument('--out', required=True)
q.add_argument('--nmon', type=int, default=3); q.add_argument('--monthly', type=int, default=0)
q = sub.add_parser('analyze'); q.add_argument('--exp', required=True); q.add_argument('--state', required=True)
q.add_argument('--nmon', type=int, default=3); q.add_argument('--monthly', type=int, default=0)
a = p.parse_args()
{'point': cmd_point, 'analyze': cmd_analyze}[a.cmd](a)
