#!/usr/bin/env python3
"""
Tests of the correlated B (b_corr.py) and its use in lbfgsb_step_monthly.py --b-corr-file (2026-11-17).

1. Operator: L L^T = C (each field block), <L v, w> = <v, L^T w> (adjoint), v -> x -> v round trip, C_h(d) = exp(-d/L).
2. Gradient: the v-space gradient of J(v) = J_obs(1 + Sigma L v) + 1/2 v.v against central finite differences.
3. Toy 4D-Var through the REAL step script (its main(), i.e. the CLI, the state file and the line search), with a known
   linear observation operator: J_obs(x) = 1/2 (y - H x)^T R^-1 (y - H x). The script reads J_obs and the gradient
   from files that this test writes in place of the forward / adjoint runs (J_value.txt, N_obs.txt and the adjoint
   accumulators GEOSChem.Adjoint.*.nc4 at the month boundaries, as the real runs write them).
   a) diagonal B (no --b-corr-file): converges to the exact minimizer x* = 1 + B H^T (H B H^T + R)^-1 (y - H 1);
   b) correlated B: converges to the exact minimizer with B = Sigma C Sigma;
   c) regression: --b-corr-file with L = 0 km and tau = 0 (C = I) gives the diagonal-B solution;
   d) the correlated solution is closer to a smooth truth than the diagonal one (the point of the exercise).
Grid: a small cubed sphere C4 (6x4x4 = 96 cells) from the same equiangular construction as GCHP, 2 fields
(land / ocean halves), 3 months.
    python test_b_corr.py [--keep]
"""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import b_corr                                             # noqa: E402

IM, NMON, NF = 4, 3, 2
T0, T1 = '2016-01-01', '2016-04-01'
FIELDS = ['SurfaceFluxAdj_CO2', 'SurfaceFluxAdj2_CO2']


def cs_centres(im):
    """Cell centres of an equiangular gnomonic cubed sphere, GMAO face order (6, im, im) -> flattened."""
    a = np.linspace(-np.pi / 4, np.pi / 4, im + 1)
    c = 0.5 * (a[1:] + a[:-1])
    X, Y = np.meshgrid(np.tan(c), np.tan(c), indexing='xy')
    one = np.ones_like(X)
    faces = [(one, X, Y), (-X, one, Y), (-one, -X, Y), (X, -one, Y), (-Y, X, one), (Y, X, -one)]
    lon, lat = [], []
    for fx, fy, fz in faces:
        r = np.sqrt(fx ** 2 + fy ** 2 + fz ** 2)
        lat.append(np.degrees(np.arcsin(fz / r)))
        lon.append(np.degrees(np.arctan2(fy, fx)))
    return np.array(lon).reshape(-1), np.array(lat).reshape(-1)


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}{": " + detail if detail else ""}')
    if not ok:
        check.failed += 1


check.failed = 0


def write_adjoint(adj_dir, g_x, nmon, ncell):
    """Adjoint accumulators A(t_m) = sum_{s >= m} g_s at each month start (A(t_end) = 0 is implied), one file per
    boundary, (time, lat=6*im, lon=im) as the native-CS adjoint output, one variable per control field."""
    os.makedirs(adj_dir, exist_ok=True)
    for f in os.listdir(adj_dir):
        os.remove(os.path.join(adj_dir, f))
    g = g_x.reshape(NF, nmon, ncell)
    bounds = pd.date_range(T0, T1, freq='MS')
    for mi, t in enumerate(bounds[:-1]):
        data = {}
        for k, var in enumerate(FIELDS):
            A = g[k, mi:].sum(0).reshape(1, 6 * IM, IM)
            data[var] = (('time', 'lat', 'lon'), A)
        xr.Dataset(data, coords={'time': [t]}).to_netcdf(
            os.path.join(adj_dir, f'GEOSChem.Adjoint.{t:%Y%m%d}_0000z.nc4'))


def run_step(work, extra):
    cmd = [sys.executable, os.path.join(HERE, 'lbfgsb_step_monthly.py'),
           '--state-file', os.path.join(work, 'state.npz'), '--forcing-dir', os.path.join(work, 'forcing'),
           '--adj-output', os.path.join(work, 'adj'), '--grid', 'cs', '--cs-res', str(IM),
           '--fields', ','.join(FIELDS), '--n-months', str(NMON), '--t-start', T0, '--t-end', T1,
           '--ftol', '1e-12', '--gtol', '1e-9', '--plateau-n', '1', '--first-step-sigma-b', '1.0',
           '--history-csv', os.path.join(work, 'hist.csv')] + extra
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode not in (0, 1):
        print(r.stdout[-3000:], r.stderr[-3000:])
        raise RuntimeError(f'step script exit {r.returncode}')
    return r.returncode, r.stdout


def toy_4dvar(work, H, Rinv, y, sig_b, extra, max_it=80):
    """Drive the step script like the 4D-Var driver: evaluate J_obs and its gradient at the state's x_next."""
    for d in ('forcing', 'adj'):
        os.makedirs(os.path.join(work, d), exist_ok=True)
    for f in ('state.npz', 'hist.csv', '4dvar_state_best.npz'):
        if os.path.exists(os.path.join(work, f)):
            os.remove(os.path.join(work, f))
    np.savez(os.path.join(work, 'sigma_b.npz'), sigma_b=sig_b)
    n = H.shape[1]
    x = np.ones(n)
    for it in range(1, max_it + 1):
        r = y - H @ x
        J_obs = 0.5 * r @ Rinv @ r
        g_x = -(H.T @ (Rinv @ r))
        open(os.path.join(work, 'forcing', 'J_value.txt'), 'w').write(f'{J_obs:.17e}\n')
        open(os.path.join(work, 'forcing', 'N_obs.txt'), 'w').write(f'{len(y)}\n')
        write_adjoint(os.path.join(work, 'adj'), g_x, NMON, 6 * IM * IM)
        rc, out = run_step(work, ['--sigma-b-file', os.path.join(work, 'sigma_b.npz')] + extra)
        st = np.load(os.path.join(work, 'state.npz'))
        x = st['x_next']
        if rc == 1:
            return x, it, st
    return x, max_it, st


def main():
    keep = '--keep' in sys.argv
    rng = np.random.default_rng(7)
    lon, lat = cs_centres(IM)
    nc_ = lon.size
    n = NF * NMON * nc_

    # ---------------- 1. operator ----------------
    print('1. operator')
    land = lat > 0                                       # field masks: "land" north, "ocean" south (+ a few cells)
    ocean = (lat <= 0) | (np.abs(lon) < 30)
    tmp = tempfile.mkdtemp(prefix='bcorr_')
    opf = os.path.join(tmp, 'b.npz')
    b_corr.build(lon, lat, [land, ocean], ['TER', 'ocean'], [1500.0, 3000.0], 1.0, NMON, opf)
    B = b_corr.BCorr(opf)
    u = b_corr.unit_vectors(lon, lat)
    for k, (m, Lk) in enumerate(((land, 1500.0), (ocean, 3000.0))):
        i = np.flatnonzero(m)
        Ch = np.exp(-b_corr.great_circle_km(u[i]) / Lk)
        Ct = b_corr.corr_time(NMON, 1.0)
        LL = np.kron(B.Lt, B.Lh[k]) @ np.kron(B.Lt, B.Lh[k]).T
        err = np.abs(LL - np.kron(Ct, Ch)).max()
        check(f'L L^T = C_t (x) C_h, field {k}', err < 1e-10, f'max err {err:.1e}')
    d01 = b_corr.great_circle_km(u[:2])[0, 1]
    ang = np.degrees(np.arccos(np.clip(u[0] @ u[1], -1, 1)))
    check('great-circle distance', abs(d01 - np.radians(ang) * 6371.0) < 1e-6, f'{d01:.1f} km')
    v, w = rng.standard_normal(n), rng.standard_normal(n)
    lhs, rhs = B.apply_L(v) @ w, v @ B.apply_LT(w)
    check('adjoint <Lv,w> = <v,L^T w>', abs(lhs - rhs) < 1e-10 * abs(lhs), f'{lhs:.12e} vs {rhs:.12e}')
    # full B from the operator vs Sigma C Sigma; cells outside the mask uncorrelated
    Lfull = np.column_stack([B.apply_L(e) for e in np.eye(n)])
    C = Lfull @ Lfull.T
    blk = C.reshape(NF, NMON, nc_, NF, NMON, nc_)
    check('no cross-field correlation', np.abs(blk[0, :, :, 1]).max() < 1e-12)
    out0 = np.flatnonzero(~land)
    sub = blk[0, :, :, 0][:, out0][:, :, :, out0]          # field 0, (month, cell outside, month, cell outside)
    check('cells outside the mask: identity', np.allclose(sub.reshape(NMON * len(out0), -1), np.eye(NMON * len(out0))))
    ins = np.flatnonzero(land)
    check('cells outside the mask: uncorrelated with the inside', np.abs(blk[0, :, :, 0][:, out0][:, :, :, ins]).max() < 1e-12)
    check('diag(C) = 1 (unit variance)', np.allclose(np.diag(C), 1.0))
    check('L invertible (v -> x -> v)', np.linalg.cond(Lfull) < 1e8, f'cond {np.linalg.cond(Lfull):.2e}')

    # ---------------- 2. v-gradient vs finite differences ----------------
    print('2. v-space gradient vs central finite differences')
    sig = 0.2 + 0.3 * rng.random(n)
    m_obs = 150
    H = rng.standard_normal((m_obs, n)) * (rng.random((m_obs, n)) < 0.1)
    Rinv = np.diag(1.0 / (0.5 + rng.random(m_obs)) ** 2)
    y = H @ (1 + 0.3 * rng.standard_normal(n))

    def Jv(vv):
        xx = 1 + sig * B.apply_L(vv)
        r = y - H @ xx
        return 0.5 * r @ Rinv @ r + 0.5 * vv @ vv

    v0 = 0.3 * rng.standard_normal(n)
    x0 = 1 + sig * B.apply_L(v0)
    g_an = B.apply_LT(sig * (-(H.T @ (Rinv @ (y - H @ x0))))) + v0
    for trial in range(3):
        dvec = rng.standard_normal(n)
        eps = 1e-4
        fd = (Jv(v0 + eps * dvec) - Jv(v0 - eps * dvec)) / (2 * eps)
        an = g_an @ dvec
        check(f'directional derivative {trial}', abs(fd - an) < 1e-7 * max(1.0, abs(an)), f'FD {fd:.10e} adj {an:.10e}')

    # ---------------- 3. toy 4D-Var through the step script ----------------
    print('3. toy 4D-Var through lbfgsb_step_monthly.py')
    # smooth truth: a large-scale pattern per field and month, so correlations should help
    tru = np.ones((NF, NMON, nc_))
    for k in range(NF):
        for mm in range(NMON):
            tru[k, mm] = 1 + 0.25 * np.sin(np.radians(lat) * (1 + k)) + 0.1 * np.cos(np.radians(lon)) + 0.05 * mm
    tru[0][:, ~land] = 1.0
    tru[1][:, ~ocean] = 1.0
    x_true = tru.reshape(-1)
    # observations: sparse, each sees a few cells of a field in one month (a crude "column")
    m_obs = 60
    H = np.zeros((m_obs, n))
    for o in range(m_obs):
        cols = rng.choice(n, 4, replace=False)
        H[o, cols] = rng.random(4)
    sig_o = 0.02
    Rinv = np.eye(m_obs) / sig_o ** 2
    y = H @ x_true + sig_o * rng.standard_normal(m_obs)
    sig_b = np.full(n, 0.3)

    def exact(Bm):
        return 1 + Bm @ H.T @ np.linalg.solve(H @ Bm @ H.T + np.linalg.inv(Rinv), y - H @ np.ones(n))

    work = os.path.join(tmp, 'run')
    os.makedirs(work)
    xa, it_a, _ = toy_4dvar(work, H, Rinv, y, sig_b, [])
    x_ex_d = exact(np.diag(sig_b ** 2))
    e = np.abs(xa - x_ex_d).max()
    check(f'a) diagonal B: exact minimizer ({it_a} iterations)', e < 2e-3, f'max |x - x*| {e:.2e}')

    xb, it_b, stb = toy_4dvar(work, H, Rinv, y, sig_b, ['--b-corr-file', opf])
    Bc = np.diag(sig_b) @ C @ np.diag(sig_b)
    x_ex_c = exact(Bc)
    e = np.abs(xb - x_ex_c).max()
    check(f'b) correlated B: exact minimizer ({it_b} iterations)', e < 2e-3, f'max |x - x*| {e:.2e}')
    check('b) state x_next = 1 + sigma_b L v_next', np.allclose(stb['x_next'], 1 + sig_b * B.apply_L(stb['v_next'])))
    check('b) state x_best physical (= 1 + sigma_b L v_best)', np.allclose(stb['x_best'], 1 + sig_b * B.apply_L(stb['v_best'])))
    # first step of the correlated run: largest |dx| / sigma_b = 1 (--first-step-sigma-b 1.0)
    h = [l.split(',') for l in open(os.path.join(work, 'hist.csv')).read().splitlines()]
    col = h[0].index('max_dx')
    check('b) first step max |dx| = 1 sigma_b', abs(float(h[1][col]) / 0.3 - 1) < 1e-9, f'max |dx| {float(h[1][col]):.6f}')

    op0 = os.path.join(tmp, 'b0.npz')
    b_corr.build(lon, lat, [land, ocean], ['TER', 'ocean'], [0.0, 0.0], 0.0, NMON, op0)
    xc, it_c, _ = toy_4dvar(work, H, Rinv, y, sig_b, ['--b-corr-file', op0])
    e = np.abs(xc - xa).max()
    check(f'c) C = I gives the diagonal-B solution ({it_c} iterations)', e < 2e-3, f'max |x_I - x_diag| {e:.2e}')

    ed = np.sqrt(np.mean((x_ex_d - x_true) ** 2)); ec = np.sqrt(np.mean((x_ex_c - x_true) ** 2))
    e0 = np.sqrt(np.mean((1 - x_true) ** 2))
    check('d) correlated B closer to the smooth truth', ec < ed,
          f'RMSE prior {e0:.4f}, diagonal {ed:.4f}, correlated {ec:.4f}')

    if keep:
        print(f'  files kept in {tmp}')
    else:
        shutil.rmtree(tmp)
    print(f'{"ALL PASSED" if check.failed == 0 else f"{check.failed} FAILED"}')
    sys.exit(1 if check.failed else 0)


if __name__ == '__main__':
    main()
