#!/usr/bin/env python3
"""
One L-BFGS-B update step for the MONTHLY 4D-Var CO2 flux optimization.

The control vector is sigma(month, lat, lon) — 12 monthly TER scaling
fields, flattened C-order with month index 0 = January (must match
write_sigma_monthly.py).

Gradient convention (doc/monthly_sigma_gradient.tex):
    SurfaceFluxAdj_CO2 is a BACKWARD ACCUMULATOR,
        A(x,t) = sum_{s >= t} lambda(x,s) * F_prior(x,s) * dt,
    so the gradient for month m is the DIFFERENCE

        g_m = A(t_m) - A(t_{m+1}),        A(t_end) = 0,

    NOT the raw A(t_m) (raw values are nested partial sums).  This script
    reads SurfaceFluxAdj_CO2 at all 13 month boundaries t_1 .. t_12, t_end
    and differences them.  Runtime guards:
      - every boundary record must exist at its exact timestamp;
      - A(t_end) must be ~0 (it is exactly 0 by construction);
      - telescoping: sum_m g_m == A(t_start) to float roundoff.

Background term per month (same prior each month):
    J_b = 0.5 * sum_m ||sigma_m - 1||^2 / sigma_b^2
    g_b_m = (sigma_m - 1) / sigma_b^2

State-file schema is identical to lbfgsb_step.py (x_next, x_prev, g_prev,
J_prev, s_hist, y_hist, m_used, iteration) with n = 12*nlat*nlon, so
write_sigma_monthly.py consumes it directly.

Exit codes:
    0  — step taken, continue optimization
    1  — converged (||grad||_inf < gtol  or  |dJ/J| < ftol)
    2  — error

Usage:
    python lbfgsb_step_monthly.py \
        --state-file   4dvar_state_monthly.npz \
        --forcing-dir  /path/forcing_files \
        --adj-output   /path/adjoint/OutputDir \
        --nlat 46 --nlon 72 \
        --sigma-b 0.2 \
        --t-start 2016-01-01 --t-end 2017-01-01 \
        [--gtol 1e-5] [--ftol 1e-8]
"""
import sys
import os
import glob
import argparse
import numpy as np
import pandas as pd
import xarray as xr

NMON = 12          # default control-vector length (full-year run)


# ---------------------------------------------------------------------------
# Month boundaries
# ---------------------------------------------------------------------------

def month_boundaries(t_start, t_end, nmon=NMON):
    """The nmon+1 boundaries t_1 .. t_nmon, t_end of an nmon-month window.

    Both t_start and t_end must fall exactly on month starts and span
    exactly nmon months.
    """
    t_start, t_end = pd.Timestamp(t_start), pd.Timestamp(t_end)
    for name, t in (('t-start', t_start), ('t-end', t_end)):
        if not (t.day == 1 and t == t.normalize()):
            raise ValueError(f'{name}={t} is not a month boundary (day 1, 00:00)')
    bounds = pd.date_range(t_start, t_end, freq='MS')
    if len(bounds) != nmon + 1:
        raise ValueError(f'{t_start} -> {t_end} spans {len(bounds) - 1} months, '
                         f'expected {nmon}')
    return bounds


# ---------------------------------------------------------------------------
# Gradient reading
# ---------------------------------------------------------------------------

def read_boundary_accumulators(adj_output_dir, boundaries, nlat, nlon):
    """
    Read SurfaceFluxAdj_CO2 at each boundary time, exactly.

    Files are opened one at a time (a year of daily adjoint output is far
    too large to concatenate).  Each boundary must be found at its exact
    timestamp — a nearest-record fallback could silently read a
    half-timestep-off accumulator and bias every monthly gradient.

    The adjoint HISTORY.rc must output on a lat-lon grid (PC{nlon}x{nlat}-DC)
    so each record is already (nlat, nlon).

    Returns A : (len(boundaries), nlat, nlon) float64
    """
    pattern = os.path.join(adj_output_dir, 'GEOSChem.Adjoint.*.nc4')
    files   = sorted(glob.glob(pattern))
    if not files:
        raise RuntimeError(f'No adjoint output files matching {pattern}')

    targets = [np.datetime64(pd.Timestamp(t)) for t in boundaries]
    A       = np.full((len(targets), nlat, nlon), np.nan, dtype=np.float64)
    src     = [None] * len(targets)

    for f in files:
        with xr.open_dataset(f) as ds:
            times = pd.DatetimeIndex(ds['time'].values).round('s')
            for k, t in enumerate(targets):
                hits = np.where(times == t)[0]
                if hits.size == 0:
                    continue
                da = ds['SurfaceFluxAdj_CO2'].isel(time=int(hits[0]))
                if da.values.shape != (nlat, nlon):
                    raise ValueError(
                        f'SurfaceFluxAdj_CO2 shape {da.values.shape} != '
                        f'({nlat}, {nlon}) in {os.path.basename(f)}. Check '
                        f'Adjoint.grid_label in HISTORY.rc '
                        f'(should be PC{nlon}x{nlat}-DC).')
                if src[k] is not None:
                    print(f'  NOTE: boundary {t} found again in '
                          f'{os.path.basename(f)} (first: {src[k]}); '
                          'keeping the first read')
                    continue
                A[k]   = da.values.astype(np.float64)
                src[k] = os.path.basename(f)

    # MAPL HISTORY in backwards mode never writes the first backward day,
    # so the record at t_end (the last boundary) does not exist on disk.
    # A(t_end) = 0 by construction (backward accumulator initial condition),
    # so substitute zeros for it; any INTERIOR boundary missing is still a
    # hard error (truncated/misaligned adjoint output).
    if src[-1] is None:
        A[-1]   = 0.0
        src[-1] = '(not written by backwards HISTORY; A(t_end)=0 by construction)'
        print(f'  NOTE: no adjoint record at t_end={targets[-1]} — '
              'expected for backwards HISTORY; using A(t_end)=0')

    missing = [str(targets[k]) for k in range(len(targets)) if src[k] is None]
    if missing:
        raise RuntimeError(
            'No adjoint record at month boundary(ies): ' + ', '.join(missing)
            + f'. Adjoint output in {adj_output_dir} is incomplete or its '
              'HISTORY frequency does not hit month starts.')

    for k, t in enumerate(targets):
        print(f'  A({t}) read from {src[k]}  '
              f'|A|_max={np.abs(A[k]).max():.6e}')
    return A


def read_J_obs(forcing_dir):
    j_path = os.path.join(forcing_dir, 'J_value.txt')
    return float(open(j_path).read().strip())


def read_N_obs(forcing_dir):
    """Matched-obs count written by co2_adjoint_forcing_osse.py (None if absent)."""
    n_path = os.path.join(forcing_dir, 'N_obs.txt')
    if not os.path.exists(n_path):
        return None
    return int(open(n_path).read().strip())


def print_balance_diagnostics(J_obs, J_b, g_obs, g_b, n_obs, n_ctrl):
    """Chi-square and obs/background balance diagnostics.

    2*J_obs/N_obs ~ 1 when R is statistically consistent; the gradient ratio
    shows how strongly observations out-pull the background at the extremes
    (>> 1 means the background barely constrains the most-active cells).
    """
    if n_obs:
        print(f'          chi2: 2*J_obs/N_obs={2*J_obs/n_obs:.3f} (N_obs={n_obs})  '
              f'2*J_b/N_ctrl={2*J_b/n_ctrl:.4f} (N_ctrl={n_ctrl})')
    gb_inf = np.abs(g_b).max()
    gb_rms = np.sqrt((g_b**2).mean())
    if gb_inf < 1e-15:
        print('          balance: g_b = 0 (sigma at prior) — ratios undefined')
    else:
        print(f'          balance: |g_obs|inf/|g_b|inf='
              f'{np.abs(g_obs).max()/gb_inf:.1f}  '
              f'rms(g_obs)/rms(g_b)='
              f'{np.sqrt((g_obs**2).mean())/gb_rms:.1f}')


def monthly_obs_gradient(adj_output_dir, t_start, t_end, nlat, nlon, nmon=NMON):
    """g_obs as (nmon*nlat*nlon,) via boundary differencing, with guards."""
    boundaries = month_boundaries(t_start, t_end, nmon)
    A = read_boundary_accumulators(adj_output_dir, boundaries, nlat, nlon)

    # Guard 1: the accumulator is exactly 0 at t_end by construction.
    a_end   = np.abs(A[-1]).max()
    a_start = np.abs(A[0]).max()
    if a_end > 1e-6 * max(a_start, 1e-30):
        raise RuntimeError(
            f'A(t_end) is not ~0 (|A(t_end)|_max={a_end:.3e}, '
            f'|A(t_start)|_max={a_start:.3e}). SurfaceFluxAdj_CO2 is not the '
            'expected backward accumulator — wrong adjoint output?')

    g = A[:-1] - A[1:]          # (NMON, nlat, nlon):  g_m = A(t_m) - A(t_m+1)

    # Guard 2: telescoping — sum of monthly gradients must recover the
    # full-window gradient A(t_start).  Catches indexing/read errors.
    resid = np.abs(g.sum(axis=0) - A[0]).max()
    tol   = 1e-4 * max(a_start, 1e-30)
    if resid > tol:
        raise RuntimeError(
            f'Telescoping guard failed: |sum_m g_m - A(t_start)|_max='
            f'{resid:.3e} > {tol:.3e}')
    print(f'  Telescoping guard OK (residual {resid:.3e}); '
          f'|A(t_end)|_max={a_end:.3e}')

    for m in range(nmon):
        print(f'  g_{m + 1:02d}: |g|_inf={np.abs(g[m]).max():.4e}  '
              f'mean={g[m].mean():.4e}')
    return g.reshape(-1), A, boundaries


# ---------------------------------------------------------------------------
# Boundary diagnostics archive
# ---------------------------------------------------------------------------

def save_boundary_diagnostics(A, g_obs, boundaries, archive_dir, iteration, nlat, nlon):
    """Save compact NC4 with month-boundary accumulators and monthly obs gradients.

    A         : (nmon+1, nlat, nlon)  backward accumulator at each boundary
    g_obs     : (nmon*nlat*nlon,)     monthly obs gradient (flattened)
    boundaries: DatetimeIndex of nmon+1 boundary timestamps
    """
    nmon    = len(boundaries) - 1
    b_strs  = [t.strftime('%Y-%m-%d') for t in boundaries]
    m_strs  = b_strs[:-1]
    g_2d    = g_obs.reshape(nmon, nlat, nlon)

    ds = xr.Dataset(
        {
            'accumulator': xr.DataArray(
                A.astype(np.float32),
                dims=['boundary', 'lat', 'lon'],
                coords={'boundary': b_strs},
                attrs={'long_name': 'Backward accumulator A(t_m) at month boundaries',
                       'note': 'g_m = A[m] - A[m+1]'}
            ),
            'gradient_obs': xr.DataArray(
                g_2d.astype(np.float32),
                dims=['month', 'lat', 'lon'],
                coords={'month': m_strs},
                attrs={'long_name': 'Monthly obs gradient g_m = A(t_m) - A(t_{m+1})'}
            ),
        },
        attrs={'iteration': iteration}
    )
    out = os.path.join(archive_dir, f'adj_boundaries_iter{iteration:03d}.nc4')
    ds.to_netcdf(out)
    print(f'  Boundary diagnostics saved → {out}')


# ---------------------------------------------------------------------------
# L-BFGS-B two-loop recursion  (no line search; unit step)
# — identical to lbfgsb_step.py
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
    parser.add_argument('--n-months', type=int,  default=NMON,
                        help='Number of monthly control-vector fields (default '
                             '12); must equal the months spanned by t-start..t-end.')
    parser.add_argument('--sigma-b', type=float, default=0.2)
    parser.add_argument('--m',       type=int,   default=10,
                        help='L-BFGS-B memory (number of vector pairs)')
    parser.add_argument('--t-start', required=True,
                        help='Window start, must be a month boundary')
    parser.add_argument('--t-end',   required=True,
                        help='Window end, must be a month boundary, '
                             't-start + 12 months')
    parser.add_argument('--gtol',    type=float, default=1e-5,
                        help='Convergence: max |grad| threshold')
    parser.add_argument('--ftol',    type=float, default=1e-8,
                        help='Convergence: relative |dJ/J| threshold')
    parser.add_argument('--archive-dir', default=None,
                        help='If set, write adj_boundaries_iter{N}.nc4 (boundary '
                             'accumulators + monthly obs gradients) here each iteration')
    args = parser.parse_args()

    nmon = args.n_months
    n = nmon * args.nlat * args.nlon
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
        if x_cur.size != n:
            print(f'ERROR: state x_next has {x_cur.size} elements, expected '
                  f'{n} ({nmon}x{args.nlat}x{args.nlon}) — state file from a '
                  'different (non-monthly?) experiment', file=sys.stderr)
            sys.exit(2)
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
        g_obs, A_bnd, boundaries = monthly_obs_gradient(
            args.adj_output, args.t_start, args.t_end,
            args.nlat, args.nlon, nmon)
    except Exception as exc:
        print(f'ERROR reading J / gradient: {exc}', file=sys.stderr)
        sys.exit(2)

    # Background cost and gradient: J_b = 0.5 * ||x-1||^2 / sigma_b^2
    # (per month; x already stacks the 12 monthly fields)
    delta  = x_cur - 1.0
    J_b    = 0.5 * np.dot(delta, delta) / args.sigma_b ** 2
    g_b    = delta / args.sigma_b ** 2

    J_cur  = J_obs + J_b
    g_cur  = g_obs + g_b

    print(f'[iter {it}]  J_obs={J_obs:.6e}  J_b={J_b:.6e}  J={J_cur:.6e}')
    print(f'          |g_obs|_inf={np.abs(g_obs).max():.4e}  '
          f'|g_b|_inf={np.abs(g_b).max():.4e}  '
          f'|g|_inf={np.abs(g_cur).max():.4e}')
    print_balance_diagnostics(J_obs, J_b, g_obs, g_b,
                              read_N_obs(args.forcing_dir), n)

    if args.archive_dir:
        save_boundary_diagnostics(A_bnd, g_obs, boundaries, args.archive_dir,
                                  it, args.nlat, args.nlon)

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
    sig    = x_next.reshape(nmon, args.nlat, args.nlon)
    print(f'          max|dx|={dx_max:.4e}  '
          f'sigma: min={x_next.min():.4f}  max={x_next.max():.4f}  '
          f'mean={x_next.mean():.4f}')
    for mm in range(nmon):
        print(f'          sigma month {mm + 1:02d}: min={sig[mm].min():.4f}  '
              f'max={sig[mm].max():.4f}  mean={sig[mm].mean():.4f}')

    # ------------------------------------------------------------------
    # Save state for next iteration
    # ------------------------------------------------------------------
    np.savez(
        args.state_file,
        x_next    = x_next,
        x_prev    = x_cur,
        g_prev    = g_cur,
        J_prev    = np.float64(J_cur),
        J_obs_prev = np.float64(J_obs),
        J_b_prev   = np.float64(J_b),
        s_hist    = s_hist,
        y_hist    = y_hist,
        m_used    = np.int64(m_used),
        iteration = np.int64(it + 1),
        nlat      = np.int64(args.nlat),
        nlon      = np.int64(args.nlon),
        nmon      = np.int64(nmon),
        start_month = np.int64(pd.Timestamp(args.t_start).month),
    )
    print(f'  State saved → {args.state_file}')

    sys.exit(0)   # continue


if __name__ == '__main__':
    main()
