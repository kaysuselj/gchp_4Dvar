#!/usr/bin/env python3
"""
One L-BFGS-B update step for the 4D-Var CO2 flux optimization.

Reads J_obs from the forcing directory and dJ/dsigma (SurfaceFluxAdj_CO2)
from the adjoint output, adds the background gradient, performs one
L-BFGS-B step with unit step size (no line search), and writes the updated
control vector to the state file.

L-BFGS-B state is persisted in --state-file between bash-loop iterations.
Bounds: sigma >= 0 (applied via projection after each step).

Exit codes:
    0  — step taken, continue optimization
    1  — converged (||grad||_inf < gtol  or  |dJ/J| < ftol)
    2  — error

Usage:
    python lbfgsb_step.py \
        --state-file   4dvar_state.npz \
        --forcing-dir  /path/forcing_files \
        --adj-output   /path/adjoint/OutputDir \
        --nlat 46 --nlon 72 \
        --sigma-b 0.5 --m 10 \
        --t-start 2019-01-01 \
        [--gtol 1e-5] [--ftol 1e-8]
"""
import sys
import os
import glob
import argparse
import numpy as np
import pandas as pd
import xarray as xr


# ---------------------------------------------------------------------------
# Gradient reading
# ---------------------------------------------------------------------------

def read_obs_gradient(adj_output_dir, t_start, nlat, nlon):
    """
    Read SurfaceFluxAdj_CO2 at the adjoint final time (= forward t_start).

    All GEOSChem.Adjoint.*.nc4 files are opened together; the time record
    nearest to t_start is selected. This handles the case where GCHP splits
    adjoint output across multiple daily files and none is named for t_start.

    The adjoint HISTORY.rc must be configured to output on a lat-lon grid
    (PC{nlon}x{nlat}-DC) so the array is already (nlat, nlon).
    """
    pattern = os.path.join(adj_output_dir, 'GEOSChem.Adjoint.*.nc4')
    files   = sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError(f'No adjoint output files matching {pattern}')

    t_target = pd.Timestamp(t_start)
    ds = (xr.open_dataset(files[0]) if len(files) == 1
          else xr.concat([xr.open_dataset(f) for f in files], dim='time'))
    grad = (ds['SurfaceFluxAdj_CO2']
            .sel(time=t_target, method='nearest')
            .values)   # (nlat, nlon) — already lat-lon from HISTORY.rc

    t_found = ds.time.sel(time=t_target, method='nearest').values
    print(f'  Gradient read from time {t_found} (target: {t_target})')

    if grad.shape != (nlat, nlon):
        raise ValueError(
            f'SurfaceFluxAdj_CO2 shape {grad.shape} != ({nlat}, {nlon}). '
            f'Check Adjoint.grid_label in HISTORY.rc (should be '
            f'PC{nlon}x{nlat}-DC).')
    return grad.ravel()


def read_J_obs(forcing_dir):
    j_path = os.path.join(forcing_dir, 'J_value.txt')
    return float(open(j_path).read().strip())


# ---------------------------------------------------------------------------
# L-BFGS-B two-loop recursion  (no line search; unit step)
# ---------------------------------------------------------------------------

def lbfgsb_direction(g, s_hist, y_hist, m_used):
    """
    Compute L-BFGS-B search direction d = -H_k g  via the two-loop recursion.
    s_hist, y_hist: arrays of shape (m, n); first m_used rows are valid
    (oldest first).
    Returns the descent direction d (same shape as g).
    """
    s = s_hist[:m_used]   # (m_used, n)
    y = y_hist[:m_used]

    q      = g.copy()
    alphas = np.zeros(m_used)
    rhos   = np.zeros(m_used)

    for i in range(m_used - 1, -1, -1):   # newest first
        rho_i      = 1.0 / np.dot(y[i], s[i])
        rhos[i]    = rho_i
        alphas[i]  = rho_i * np.dot(s[i], q)
        q          = q - alphas[i] * y[i]

    # Initial Hessian scaling: gamma = s_{k-1}^T y_{k-1} / y_{k-1}^T y_{k-1}
    if m_used > 0:
        gamma = np.dot(s[-1], y[-1]) / np.dot(y[-1], y[-1])
    else:
        # First iteration: scale so max per-element change ≈ 1
        gamma = 1.0 / max(np.abs(g).max(), 1e-8)

    r = gamma * q

    for i in range(m_used):               # oldest first
        beta = rhos[i] * np.dot(y[i], r)
        r    = r + s[i] * (alphas[i] - beta)

    return -r   # descent direction


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--state-file',  required=True)
    parser.add_argument('--forcing-dir', required=True)
    parser.add_argument('--adj-output',  required=True,
                        help='Adjoint OutputDir containing GEOSChem.Adjoint.*.nc4')
    parser.add_argument('--nlat',    type=int,   default=46)
    parser.add_argument('--nlon',    type=int,   default=72)
    parser.add_argument('--sigma-b', type=float, default=0.5)
    parser.add_argument('--m',       type=int,   default=10,
                        help='L-BFGS-B memory (number of vector pairs)')
    parser.add_argument('--t-start', default='2019-01-01')
    parser.add_argument('--gtol',    type=float, default=1e-5,
                        help='Convergence: max |grad| threshold')
    parser.add_argument('--ftol',    type=float, default=1e-8,
                        help='Convergence: relative |dJ/J| threshold')
    args = parser.parse_args()

    n = args.nlat * args.nlon
    m = args.m

    # ------------------------------------------------------------------
    # Load state (or initialise for iteration 1)
    # ------------------------------------------------------------------
    if os.path.exists(args.state_file):
        state   = np.load(args.state_file)
        x_cur   = state['x_next'].copy()      # sigma used in this iteration
        x_prev  = state['x_prev'].copy()      # sigma from iteration before
        g_prev  = state['g_prev'].copy()
        J_prev  = float(state['J_prev'])
        s_hist  = state['s_hist'].copy()      # (m, n)
        y_hist  = state['y_hist'].copy()
        m_used  = int(state['m_used'])
        it      = int(state['iteration'])
    else:
        x_cur   = np.ones(n)
        x_prev  = np.empty(n)
        g_prev  = np.empty(n)
        J_prev  = np.inf
        s_hist  = np.zeros((m, n))
        y_hist  = np.zeros((m, n))
        m_used  = 0
        it      = 1

    # ------------------------------------------------------------------
    # Read current J and gradient
    # ------------------------------------------------------------------
    try:
        J_obs = read_J_obs(args.forcing_dir)
        g_obs = read_obs_gradient(args.adj_output, args.t_start,
                                  args.nlat, args.nlon)
    except Exception as exc:
        print(f'ERROR reading J / gradient: {exc}', file=sys.stderr)
        sys.exit(2)

    # Background cost and gradient: J_b = 0.5 * ||x-1||^2 / sigma_b^2
    delta  = x_cur - 1.0
    J_b    = 0.5 * np.dot(delta, delta) / args.sigma_b ** 2
    g_b    = delta / args.sigma_b ** 2

    J_cur  = J_obs + J_b
    g_cur  = g_obs + g_b

    print(f'[iter {it}]  J_obs={J_obs:.6e}  J_b={J_b:.6e}  J={J_cur:.6e}')
    print(f'          |g_obs|_inf={np.abs(g_obs).max():.4e}  '
          f'|g_b|_inf={np.abs(g_b).max():.4e}  '
          f'|g|_inf={np.abs(g_cur).max():.4e}')

    # ------------------------------------------------------------------
    # Convergence check
    # ------------------------------------------------------------------
    g_inf   = np.abs(g_cur).max()
    dJ_rel  = abs(J_cur - J_prev) / max(abs(J_cur), 1.0)

    if g_inf < args.gtol:
        print(f'CONVERGED: ||grad||_inf={g_inf:.4e} < gtol={args.gtol}')
        sys.exit(1)
    if it > 1 and dJ_rel < args.ftol:
        print(f'CONVERGED: |dJ/J|={dJ_rel:.4e} < ftol={args.ftol}')
        sys.exit(1)

    # ------------------------------------------------------------------
    # Update L-BFGS-B history with the pair from the previous iteration
    # ------------------------------------------------------------------
    if it > 1:
        s_new = x_cur - x_prev
        y_new = g_cur - g_prev
        ys    = np.dot(y_new, s_new)
        if ys > 1e-10 * np.dot(y_new, y_new):   # curvature condition
            if m_used < m:
                s_hist[m_used] = s_new
                y_hist[m_used] = y_new
                m_used += 1
            else:
                # Shift oldest out (circular buffer: oldest is row 0)
                s_hist[:-1] = s_hist[1:]
                y_hist[:-1] = y_hist[1:]
                s_hist[-1]  = s_new
                y_hist[-1]  = y_new
        else:
            print(f'  Skipping history update: curvature condition failed '
                  f'(y^T s = {ys:.4e})')

    # ------------------------------------------------------------------
    # L-BFGS-B search direction and update
    # ------------------------------------------------------------------
    d      = lbfgsb_direction(g_cur, s_hist, y_hist, m_used)
    x_next = np.clip(x_cur + d, 0.0, None)   # unit step, enforce sigma >= 0

    dx_max = np.abs(x_next - x_cur).max()
    print(f'          max|dx|={dx_max:.4e}  '
          f'sigma: min={x_next.min():.4f}  max={x_next.max():.4f}  '
          f'mean={x_next.mean():.4f}')

    # ------------------------------------------------------------------
    # Save state for next iteration
    # ------------------------------------------------------------------
    np.savez(
        args.state_file,
        x_next    = x_next,
        x_prev    = x_cur,
        g_prev    = g_cur,
        J_prev    = np.float64(J_cur),
        s_hist    = s_hist,
        y_hist    = y_hist,
        m_used    = np.int64(m_used),
        iteration = np.int64(it + 1),
    )
    print(f'  State saved → {args.state_file}')

    sys.exit(0)   # continue


if __name__ == '__main__':
    main()
