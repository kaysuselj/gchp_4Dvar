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

STEP CONTROL (added 2026-08-05, see project analysis of 3months_osse):
    The bare unit L-BFGS step diverged (J rose, sigma blew up to ~29) because a
    step that INCREASED J was accepted unconditionally and its (s,y) pair
    poisoned the Hessian.  This version adds a *deferred backtracking line
    search* + *keep-best-iterate*:
      - accept the just-evaluated iterate only if it decreased J (Armijo, c1);
      - on reject, keep the anchor + direction, halve the step, and DON'T touch
        the L-BFGS memory (each backtrack trial is the next outer forward+adjoint);
      - a non-descent L-BFGS direction (g.d >= 0) resets the memory to a scaled
        steepest-descent step, which is always descent;
      - the lowest-J iterate ever seen is tracked (x_best) and written to
        4dvar_state_best.npz, and delivered as the solution on convergence/stall.
    Set --max-backtrack 0 to recover the old pure-unit-step behaviour (keep-best
    still active, which is harmless).

State-file schema extends lbfgsb_step.py's (x_next, x_prev, g_prev, J_prev,
s_hist, y_hist, m_used, iteration, nlat, nlon, nmon, start_month) with the
anchor / step-control / best-iterate fields (x_anchor, g_anchor, J_anchor,
d_cur, alpha, n_backtrack, x_best, J_best, it_best).  Old-format states are
migrated on load (anchor := previous iterate).  write_sigma_monthly.py consumes
x_next exactly as before.

Exit codes:
    0  — step taken, continue optimization
    1  — converged (||grad||_inf < gtol  or  |dJ/J| < ftol) or line-search stall;
         the delivered x_next is the best iterate
    2  — error (handled/expected: bad state, missing inputs, ...)
    3  — unexpected crash (unhandled exception).  Distinct from 1 on purpose so
         a bug can never be misread as convergence by the driver.

Usage:
    python lbfgsb_step_monthly.py \
        --state-file   4dvar_state_monthly.npz \
        --forcing-dir  /path/forcing_files \
        --adj-output   /path/adjoint/OutputDir \
        --nlat 46 --nlon 72 \
        --sigma-b 0.2 \
        --t-start 2016-01-01 --t-end 2017-01-01 \
        [--gtol 1e-5] [--ftol 1e-8] [--c1 1e-4] [--max-backtrack 6]
"""
import sys
import os
import glob
import csv
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
# Per-iteration diagnostics history (machine-readable CSV)
# ---------------------------------------------------------------------------

HISTORY_FIELDS = [
    'iteration', 'action', 'reject', 'reject_reason', 'J', 'J_obs', 'J_b',
    'g_inf', 'g_2', 'chi2', 'dJ_rel', 'actual_dJ', 'pred_dJ', 'rho', 'gTd',
    'alpha', 'n_backtrack', 'n_reject_total', 'max_dx',
    'sigma_min', 'sigma_max', 'sigma_mean', 'sigma_prop_max',
    'n_at_bound', 'n_clipped', 'reset_memory', 'm_used', 'gamma',
    'curvature_ok', 'is_best', 'J_best', 'it_best',
]


def append_history(csv_path, row):
    """Append one diagnostics row (dict keyed by HISTORY_FIELDS) to csv_path,
    writing the header first if the file does not yet exist.  Best-effort: a
    diagnostics-write failure must never abort the optimization."""
    try:
        os.makedirs(os.path.dirname(csv_path), exist_ok=True)
        new = not os.path.exists(csv_path)
        with open(csv_path, 'a', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=HISTORY_FIELDS)
            if new:
                w.writeheader()
            w.writerow({k: row.get(k, '') for k in HISTORY_FIELDS})
        print(f'  History row appended → {csv_path}')
    except Exception as exc:
        print(f'  WARNING: could not write history CSV ({exc})')


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


def _gamma_scale(g, s_hist, y_hist, m_used):
    """Initial-Hessian scaling gamma used by lbfgsb_direction (for logging)."""
    if m_used > 0:
        return float(np.dot(s_hist[m_used - 1], y_hist[m_used - 1])
                     / np.dot(y_hist[m_used - 1], y_hist[m_used - 1]))
    return float(1.0 / max(np.abs(g).max(), 1e-8))


# ---------------------------------------------------------------------------
# Deferred backtracking line search + keep-best  (pure, unit-tested)
# ---------------------------------------------------------------------------

def lbfgs_backtrack_step(x_cur, g_cur, J_cur, st, m, c1, max_bt):
    """One projected L-BFGS step with a line search deferred across outer
    iterations.  PURE (no I/O) so it can be unit-tested on a toy quadratic.

    Parameters
    ----------
    x_cur, g_cur, J_cur : the iterate just evaluated (proposed last call).
    st : dict, the persisted optimizer state.  Mutated in place.  Keys:
        x_anchor,g_anchor,J_anchor : last ACCEPTED iterate (None on iter 1)
        d, alpha, n_bt             : current search direction / step / backtracks
        s_hist,y_hist,m_used       : L-BFGS memory
        it                         : iteration number (of x_cur)
        x_best,J_best,it_best      : best iterate seen
    m, c1, max_bt : memory size, Armijo constant, max backtracks along one dir.

    Returns (x_next, st, info).  x_next is None iff info['stop'] is set
    (line-search stalled — caller should deliver the best iterate and exit).
    """
    info = {'curvature_ok': '', 'reset_memory': False, 'is_best': False,
            'reject': False, 'reject_reason': '',
            'pred_dJ': float('nan'), 'actual_dJ': float('nan'),
            'rho': float('nan')}

    # ---- keep-best-iterate --------------------------------------------------
    if J_cur < st['J_best']:
        st['J_best'] = float(J_cur)
        st['x_best'] = x_cur.copy()
        st['it_best'] = int(st['it'])
        info['is_best'] = True

    first = st['x_anchor'] is None
    if first:
        accept = True
    else:
        # Judge the step ACTUALLY taken (anchor -> x_cur, post-clip) vs the anchor.
        step      = x_cur - st['x_anchor']
        pred_dJ   = float(np.dot(st['g_anchor'], step))     # < 0 for a descent step
        actual_dJ = float(J_cur - st['J_anchor'])
        info['pred_dJ'], info['actual_dJ'] = pred_dJ, actual_dJ
        info['rho'] = actual_dJ / pred_dJ if pred_dJ != 0.0 else float('nan')
        # Armijo sufficient decrease (pred_dJ<0): accept only a genuine decrease.
        accept = (actual_dJ < 0.0) and (actual_dJ <= c1 * pred_dJ)
        if not accept:
            info['reject'] = True
            if actual_dJ >= 0.0:
                info['reject_reason'] = 'J_rose'          # step increased the cost
            elif pred_dJ >= 0.0:
                info['reject_reason'] = 'nondescent_step'  # clip made the step non-descent
            else:
                info['reject_reason'] = 'insufficient_decrease'  # decreased < Armijo bound

    if accept:
        # form the L-BFGS pair from the accepted move (anchor_old -> x_cur)
        if not first:
            s_new = x_cur - st['x_anchor']
            y_new = g_cur - st['g_anchor']
            ys    = float(np.dot(y_new, s_new))
            if ys > 1e-10 * float(np.dot(y_new, y_new)):     # curvature condition
                if st['m_used'] < m:
                    st['s_hist'][st['m_used']] = s_new
                    st['y_hist'][st['m_used']] = y_new
                    st['m_used'] += 1
                else:
                    st['s_hist'][:-1] = st['s_hist'][1:]
                    st['y_hist'][:-1] = st['y_hist'][1:]
                    st['s_hist'][-1]  = s_new
                    st['y_hist'][-1]  = y_new
                info['curvature_ok'] = True
            else:
                info['curvature_ok'] = False
        # advance the anchor to the accepted iterate
        st['x_anchor'], st['g_anchor'], st['J_anchor'] = x_cur.copy(), g_cur.copy(), float(J_cur)
        # fresh L-BFGS direction from the new anchor
        d = lbfgsb_direction(g_cur, st['s_hist'], st['y_hist'], st['m_used'])
        gTd = float(np.dot(g_cur, d))
        if gTd >= 0.0:      # non-descent safeguard: drop memory -> scaled steepest descent
            st['m_used'] = 0
            st['s_hist'][:] = 0.0
            st['y_hist'][:] = 0.0
            d = lbfgsb_direction(g_cur, st['s_hist'], st['y_hist'], 0)
            gTd = float(np.dot(g_cur, d))
            info['reset_memory'] = True
        st['d'], st['alpha'], st['n_bt'] = d, 1.0, 0
        info['action'] = 'init' if first else 'accept'
        anchor_g = g_cur
    else:
        # reject the overshoot: keep anchor + direction, halve the step,
        # DO NOT update the L-BFGS memory (a bad step must not poison H)
        st['n_reject_total'] = int(st.get('n_reject_total', 0)) + 1
        st['n_bt'] += 1
        if st['n_bt'] > max_bt:
            info['action'] = 'stall'
            info['stop'] = ('stall',
                            f'no decrease after {max_bt} backtracks along the '
                            f'search direction')
            info['gTd'] = float(np.dot(st['g_anchor'], st['d']))
            return None, st, info
        st['alpha'] *= 0.5
        d = st['d']
        gTd = float(np.dot(st['g_anchor'], d))
        info['action'] = 'backtrack'
        anchor_g = st['g_anchor']

    # projected trial point along d from the current anchor
    x_raw  = st['x_anchor'] + st['alpha'] * st['d']
    x_next = np.clip(x_raw, 0.0, None)
    info['gTd']        = gTd
    info['alpha']      = st['alpha']
    info['n_backtrack'] = st['n_bt']
    info['n_clipped']  = int((x_raw < 0.0).sum())
    info['n_at_bound'] = int((x_next <= 0.0).sum())
    info['max_dx']     = float(np.abs(x_next - x_cur).max())
    info['gamma']      = _gamma_scale(anchor_g, st['s_hist'], st['y_hist'], st['m_used'])
    return x_next, st, info


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
                        help='Convergence: relative |dJ/J| threshold (accepted steps)')
    parser.add_argument('--c1',      type=float, default=1e-4,
                        help='Armijo sufficient-decrease constant for the '
                             'backtracking line search')
    parser.add_argument('--max-backtrack', type=int, default=6,
                        help='Max step halvings along one direction before '
                             'declaring a line-search stall (0 = old unit-step '
                             'behaviour, no backtracking)')
    parser.add_argument('--history-csv', default=None,
                        help='Per-iteration diagnostics CSV (default: '
                             '<state dir>/4dvar_output_osse/4dvar_history.csv)')
    parser.add_argument('--archive-dir', default=None,
                        help='If set, write adj_boundaries_iter{N}.nc4 (boundary '
                             'accumulators + monthly obs gradients) here each iteration')
    args = parser.parse_args()

    nmon = args.n_months
    n = nmon * args.nlat * args.nlon
    m = args.m

    history_csv = args.history_csv or os.path.join(
        os.path.dirname(os.path.abspath(args.state_file)),
        '4dvar_output_osse', '4dvar_history.csv')
    best_file = os.path.join(os.path.dirname(os.path.abspath(args.state_file)),
                             '4dvar_state_best.npz')

    # ------------------------------------------------------------------
    # Load state (or initialise for iteration 1).  New step-control / anchor /
    # best fields are read with defaults so OLD-format states migrate cleanly
    # (anchor := previous iterate; best starts tracking from here).
    # ------------------------------------------------------------------
    if os.path.exists(args.state_file):
        state   = np.load(args.state_file)
        keys    = set(state.files)
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
        if 'x_anchor' in keys:                # new-format state
            x_anchor = state['x_anchor'].copy()
            g_anchor = state['g_anchor'].copy()
            J_anchor = float(state['J_anchor'])
            d_cur    = state['d_cur'].copy()
            alpha    = float(state['alpha'])
            n_bt     = int(state['n_backtrack'])
            x_best   = state['x_best'].copy()
            J_best   = float(state['J_best'])
            it_best  = int(state['it_best'])
            n_reject_total = int(state['n_reject_total']) if 'n_reject_total' in keys else 0
        elif it > 1:                          # migrate old-format state
            print('  NOTE: migrating old-format state — anchor := previous iterate')
            x_anchor, g_anchor, J_anchor = x_prev, g_prev, J_prev
            d_cur = x_cur - x_prev            # the unit step that produced x_cur
            alpha, n_bt = 1.0, 0
            x_best, J_best, it_best = x_cur.copy(), float('inf'), 0
            n_reject_total = 0
        else:
            x_anchor = g_anchor = None; J_anchor = float('inf')
            d_cur = np.zeros(n); alpha, n_bt = 1.0, 0
            x_best, J_best, it_best = x_cur.copy(), float('inf'), 0
            n_reject_total = 0
    else:
        x_cur   = np.ones(n)
        x_prev  = np.empty(n)
        g_prev  = np.empty(n)
        J_prev  = np.inf
        s_hist  = np.zeros((m, n))
        y_hist  = np.zeros((m, n))
        m_used  = 0
        it      = 1
        x_anchor = g_anchor = None; J_anchor = float('inf')
        d_cur = np.zeros(n); alpha, n_bt = 1.0, 0
        x_best, J_best, it_best = x_cur.copy(), float('inf'), 0
        n_reject_total = 0

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
    n_obs = read_N_obs(args.forcing_dir)
    print_balance_diagnostics(J_obs, J_b, g_obs, g_b, n_obs, n)

    if args.archive_dir:
        save_boundary_diagnostics(A_bnd, g_obs, boundaries, args.archive_dir,
                                  it, args.nlat, args.nlon)

    g_inf = float(np.abs(g_cur).max())
    g_2   = float(np.sqrt(np.dot(g_cur, g_cur)))
    chi2  = (2.0 * J_obs / n_obs) if n_obs else float('nan')
    J_anchor_pre = J_anchor      # anchor J BEFORE this step (for ftol on accept)

    # ------------------------------------------------------------------
    # Deferred backtracking L-BFGS step + keep-best
    # ------------------------------------------------------------------
    st = dict(x_anchor=x_anchor, g_anchor=g_anchor, J_anchor=J_anchor,
              d=d_cur, alpha=alpha, n_bt=n_bt,
              s_hist=s_hist, y_hist=y_hist, m_used=m_used, it=it,
              x_best=x_best, J_best=J_best, it_best=it_best,
              n_reject_total=n_reject_total)
    x_next, st, info = lbfgs_backtrack_step(x_cur, g_cur, J_cur, st,
                                            m, args.c1, args.max_backtrack)

    # console step report
    if info.get('action') in ('accept', 'init'):
        rmsg = '  [memory reset -> steepest descent]' if info.get('reset_memory') else ''
        print(f'          step ACCEPT ({info["action"]}): rho={info["rho"]:.3f}  '
              f'g.d={info["gTd"]:.4e}{rmsg}')
    elif info.get('action') in ('backtrack', 'stall'):
        # NB: on the 'stall' action lbfgs_backtrack_step returns early, before
        # info['n_backtrack']/['alpha'] are populated — so read the backtrack
        # count from st['n_bt'] (always present) and .get() every info field to
        # keep this diagnostic line from ever crashing the step.
        print(f'          step REJECT [{info.get("reject_reason")}] '
              f'backtrack #{st["n_bt"]} (rho={info.get("rho", float("nan")):.3f}, '
              f'actual_dJ={info.get("actual_dJ", float("nan")):.3e}), '
              f'alpha->{info.get("alpha", float("nan")):.4g}, memory kept  '
              f'(rejects this run: {st["n_reject_total"]})')

    def _persist_best():
        """Write the best iterate to a standalone file for recovery / delivery."""
        try:
            np.savez(best_file, x_next=st['x_best'], J=np.float64(st['J_best']),
                     iteration=np.int64(st['it_best']),
                     nlat=np.int64(args.nlat), nlon=np.int64(args.nlon),
                     nmon=np.int64(nmon),
                     start_month=np.int64(pd.Timestamp(args.t_start).month))
        except Exception as exc:
            print(f'  WARNING: could not write best-state file ({exc})')

    if info.get('is_best'):
        _persist_best()
        print(f'  New best iterate: iter {st["it_best"]}  J={st["J_best"]:.6e} '
              f'→ {best_file}')

    # ------------------------------------------------------------------
    # Convergence / stall.  On any stop we DELIVER the best iterate (write it as
    # the final x_next) so the optimizer never returns a worse-than-best sigma.
    # ------------------------------------------------------------------
    def _save_state(x_deliver, converged):
        np.savez(
            args.state_file,
            x_next    = x_deliver,
            x_prev    = x_cur,
            g_prev    = g_cur,
            J_prev    = np.float64(J_cur),
            J_obs_prev = np.float64(J_obs),
            J_b_prev   = np.float64(J_b),
            s_hist    = st['s_hist'],
            y_hist    = st['y_hist'],
            m_used    = np.int64(st['m_used']),
            iteration = np.int64(it + 1),
            nlat      = np.int64(args.nlat),
            nlon      = np.int64(args.nlon),
            nmon      = np.int64(nmon),
            start_month = np.int64(pd.Timestamp(args.t_start).month),
            # step-control / anchor / best
            x_anchor  = (st['x_anchor'] if st['x_anchor'] is not None else x_cur),
            g_anchor  = (st['g_anchor'] if st['g_anchor'] is not None else g_cur),
            J_anchor  = np.float64(st['J_anchor']),
            d_cur     = st['d'],
            alpha     = np.float64(st['alpha']),
            n_backtrack = np.int64(st['n_bt']),
            x_best    = st['x_best'],
            J_best    = np.float64(st['J_best']),
            it_best   = np.int64(st['it_best']),
            n_reject_total = np.int64(st['n_reject_total']),
        )

    dJ_rel = float('nan')
    stop_reason = None
    if info.get('stop'):
        stop_reason = f'line-search stall ({info["stop"][1]})'
    elif g_inf < args.gtol:
        stop_reason = f'||grad||_inf={g_inf:.4e} < gtol={args.gtol}'
    elif info.get('action') == 'accept' and np.isfinite(J_anchor_pre):
        dJ_rel = abs(J_cur - J_anchor_pre) / max(abs(J_cur), 1.0)
        if dJ_rel < args.ftol:
            stop_reason = f'|dJ/J|={dJ_rel:.4e} < ftol={args.ftol}'

    # diagnostics CSV row (best-effort).  sigma_* describe the EVALUATED iterate
    # x_cur (the sigma this row's J/g correspond to); sigma_prop_max is the max of
    # the NEXT proposal x_next (blank when we stop and deliver best).
    sigma_prop_max = float(x_next.max()) if x_next is not None else ''
    append_history(history_csv, dict(
        iteration=it, action=info.get('action', ''), reject=info.get('reject'),
        reject_reason=info.get('reject_reason'), J=J_cur, J_obs=J_obs, J_b=J_b,
        g_inf=g_inf, g_2=g_2, chi2=chi2, dJ_rel=dJ_rel,
        actual_dJ=info.get('actual_dJ'), pred_dJ=info.get('pred_dJ'),
        rho=info.get('rho'), gTd=info.get('gTd'), alpha=info.get('alpha'),
        n_backtrack=info.get('n_backtrack'), n_reject_total=st['n_reject_total'],
        max_dx=info.get('max_dx'),
        sigma_min=float(x_cur.min()), sigma_max=float(x_cur.max()),
        sigma_mean=float(x_cur.mean()), sigma_prop_max=sigma_prop_max,
        n_at_bound=info.get('n_at_bound'), n_clipped=info.get('n_clipped'),
        reset_memory=info.get('reset_memory'), m_used=st['m_used'],
        gamma=info.get('gamma'), curvature_ok=info.get('curvature_ok'),
        is_best=info.get('is_best'), J_best=st['J_best'], it_best=st['it_best']))

    if stop_reason:
        print(f'CONVERGED/STOP: {stop_reason}')
        print(f'  Delivering BEST iterate (iter {st["it_best"]}, '
              f'J={st["J_best"]:.6e}) as the solution')
        _persist_best()
        _save_state(st['x_best'], converged=True)
        print(f'  State saved → {args.state_file}')
        sys.exit(1)

    # normal continue: deliver the trial point for the next forward
    dx_max = info.get('max_dx', float('nan'))
    signx = x_next.reshape(nmon, args.nlat, args.nlon)
    print(f'          max|dx|={dx_max:.4e}  '
          f'sigma: min={x_next.min():.4f}  max={x_next.max():.4f}  '
          f'mean={x_next.mean():.4f}  at_bound={info.get("n_at_bound")}  '
          f'clipped={info.get("n_clipped")}')
    for mm in range(nmon):
        print(f'          sigma month {mm + 1:02d}: min={signx[mm].min():.4f}  '
              f'max={signx[mm].max():.4f}  mean={signx[mm].mean():.4f}')

    _save_state(x_next, converged=False)
    print(f'  State saved → {args.state_file}')
    sys.exit(0)   # continue


if __name__ == '__main__':
    # Guard against exit-code collision: the wrapper reads exit code 1 as
    # "converged/stall", but an unhandled exception ALSO exits 1 — which would
    # silently masquerade as convergence.  sys.exit(0/1/2) raise SystemExit
    # (not an Exception subclass), so those intended codes pass through; only a
    # genuine crash is caught here and re-mapped to a distinct code 3.
    try:
        main()
    except Exception:
        import traceback
        traceback.print_exc()
        sys.stderr.write('UNEXPECTED ERROR in lbfgsb_step_monthly — '
                         'exiting 3 (this is a crash, NOT a convergence)\n')
        sys.exit(3)
