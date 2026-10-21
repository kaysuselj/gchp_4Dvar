#!/usr/bin/env python3
"""Continue a 4D-Var run that stopped on the convergence test (2026-10-21: switch from ftol 0.02 to the
plateau rule). On a stop the step script delivers the BEST iterate as x_next; the L-BFGS state (anchor,
direction d, step alpha, s/y memory) of the accepted step is saved as well. This rebuilds the trial point the
step would have proposed, x_next = clip(x_anchor + alpha * d, sigma_min, sigma_max), so the run continues
exactly as if it had never stopped.
    resume_after_stop.py STATE_FILE --sigma-min 0.01 --sigma-max 10   (backup: STATE_FILE.stop_backup.npz)
"""
import argparse, shutil
import numpy as np
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('state'); p.add_argument('--sigma-min', type=float, required=True); p.add_argument('--sigma-max', type=float, required=True)
a = p.parse_args()
s = dict(np.load(a.state, allow_pickle=False))
it_done = int(s['iteration']) - 1
if int(s['it_best']) != it_done:
    raise SystemExit(f'last evaluated iteration {it_done} is not the best ({int(s["it_best"])}): stop was not on an accepted best step')
if not np.array_equal(s['x_anchor'], s['x_prev']):
    raise SystemExit('anchor != last evaluated point: the last step was not accepted')
if not np.array_equal(s['x_next'], s['x_best']):
    raise SystemExit('x_next is not the delivered best point: state was not written by a stop')
if float(s['alpha']) != 1.0 or int(s['n_backtrack']) != 0:
    raise SystemExit(f'unexpected step control alpha={float(s["alpha"])} n_backtrack={int(s["n_backtrack"])}')
x_raw = s['x_anchor'] + s['alpha'] * s['d_cur']
x_new = np.clip(x_raw, a.sigma_min, a.sigma_max)
gTd = float(np.dot(s['g_prev'], s['d_cur']))
if not gTd < 0:
    raise SystemExit(f'saved direction is not a descent direction (g.d = {gTd:.3e})')
shutil.copy2(a.state, a.state.replace('.npz', '') + '.stop_backup.npz')
s['x_next'] = x_new
np.savez(a.state, **s)
print(f'{a.state}: continue at iteration {it_done + 1}; trial = anchor + d (g.d {gTd:.3e}), max|dx| {np.abs(x_new - s["x_prev"]).max():.4e}, '
      f'sigma {x_new.min():.4f}-{x_new.max():.4f}, clipped {int(((x_raw < a.sigma_min) | (x_raw > a.sigma_max)).sum())}')
