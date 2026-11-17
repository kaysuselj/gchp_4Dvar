#!/usr/bin/env python3
"""
Correlated background-error covariance B for the monthly 4D-Var control (2026-11-17).

B = Sigma C Sigma, with Sigma = diag(sigma_b) (the diagonal background error used so far) and C a correlation matrix
that is separable in space and time and block-diagonal in the control fields (TER, ocean; no cross-field correlation):

    C_field = C_t (x) C_h,      C_h[i, j] = exp(-d_ij / L_field),      C_t[m, n] = exp(-|m - n| / tau)

d_ij is the great-circle distance between the cell centres (exponential correlations are positive definite on the
sphere with the great-circle distance). C_h acts only on the cells of the field's mask (cells with a non-zero prior
flux in any month of the window); the other cells are not correlated with anything (identity), so they stay inert as
before (zero gradient, sigma stays 1).

Change of variable (the optimizer works in v, J_b = 1/2 v^T v):

    x = 1 + Sigma L v,       L = L_t (x) L_h,       L_h = U_h Lambda_h^1/2  (C_h = U_h Lambda_h U_h^T), same for L_t
    grad_v J = v + L^T Sigma grad_x J_obs

so L L^T = C and Sigma L L^T Sigma = B. L is square (full rank up to round-off), so v has the size of x.

Control layout (as lbfgsb_step_monthly.py): x = (n_fields, n_months, n_cells) flattened C-order, n_cells = 6*im*im on
the native cubed sphere (GMAO stacked layout, the cell order of the HEMCO / adjoint output files).

Build (writes the operator; seconds at C24, the dense C_h is n_mask x n_mask):
    b_corr.py build --grid-file <HEMCO emissions file with lons/lats> --prior-file <prior emissions file> \\
        --fields EmisCO2_NetTerrExch,EmisCO2_Ocean --length-km 500,1000 --tau-months 1.0 --n-months 19 \\
        [--first-month-index 0] --out b_corr.npz
Defaults: TER 500 km, ocean 1000 km, tau = 1 month (adjacent months correlated exp(-1) = 0.37).
    b_corr.py info b_corr.npz
"""
import argparse
import sys
import numpy as np

R_EARTH_KM = 6371.0


def unit_vectors(lon_deg, lat_deg):
    lo, la = np.radians(lon_deg), np.radians(lat_deg)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def great_circle_km(u):
    """(k, 3) unit vectors -> (k, k) great-circle distances [km] (arctan2 form, accurate at small distances)."""
    cr = np.linalg.norm(np.cross(u[:, None, :], u[None, :, :]), axis=-1)
    dt = np.clip(u @ u.T, -1.0, 1.0)
    return R_EARTH_KM * np.arctan2(cr, dt)


def sqrt_corr(C, name, rtol=1e-10):
    """Square root L with L L^T = C (symmetric eigendecomposition); tiny negative eigenvalues (round-off) -> 0."""
    lam, U = np.linalg.eigh(C)
    lmin = lam.min()
    if lmin < -rtol * lam.max() * 1e3:
        raise ValueError(f'{name}: correlation matrix not positive semi-definite (min eigenvalue {lmin:.3e})')
    lam = np.clip(lam, 0.0, None)
    return U * np.sqrt(lam)[None, :], lam


def corr_time(nmon, tau):
    m = np.arange(nmon, dtype=np.float64)
    if tau <= 0:
        return np.eye(nmon)
    return np.exp(-np.abs(m[:, None] - m[None, :]) / tau)


class BCorr:
    """The operator x - 1 = Sigma L v (L only; Sigma is applied by the caller, it may change per run)."""

    def __init__(self, path):
        z = np.load(path, allow_pickle=False)
        self.path = path
        self.n_fields = int(z['n_fields'])
        self.n_months = int(z['n_months'])
        self.n_cells = int(z['n_cells'])
        self.n = self.n_fields * self.n_months * self.n_cells
        self.fields = [str(s) for s in z['fields']]
        self.length_km = np.asarray(z['length_km'], float)
        self.tau = float(z['tau_months'])
        self.Lt = np.asarray(z['Lt'], np.float64)                   # (n_months, n_months)
        self.idx = [np.asarray(z[f'idx_{k}'], np.int64) for k in range(self.n_fields)]
        self.Lh = [np.asarray(z[f'Lh_{k}'], np.float64) for k in range(self.n_fields)]

    def _blocks(self, a):
        return a.reshape(self.n_fields, self.n_months, self.n_cells)

    def apply_L(self, v):
        """L v. Inside each field's mask: L_t V L_h^T; outside: identity."""
        V = self._blocks(np.asarray(v, np.float64))
        out = V.copy()
        for k in range(self.n_fields):
            i = self.idx[k]
            out[k][:, i] = self.Lt @ V[k][:, i] @ self.Lh[k].T
        return out.reshape(-1)

    def apply_LT(self, w):
        """L^T w (the adjoint of apply_L)."""
        W = self._blocks(np.asarray(w, np.float64))
        out = W.copy()
        for k in range(self.n_fields):
            i = self.idx[k]
            out[k][:, i] = self.Lt.T @ W[k][:, i] @ self.Lh[k]
        return out.reshape(-1)

    def summary(self):
        return (f'B correlations {self.path}: {self.n_fields} field(s) {",".join(self.fields)}, {self.n_months} months, '
                f'{self.n_cells} cells; L = ' + ', '.join(f'{f} {l:g} km ({len(i)} cells)' for f, l, i in
                                                         zip(self.fields, self.length_km, self.idx))
                + f'; tau = {self.tau:g} month(s)')


def build(lons, lats, masks, fields, length_km, tau, n_months, out, meta=''):
    """lons, lats: (n_cells,) cell centres in control order; masks: list of (n_cells,) bool per field."""
    u = unit_vectors(lons, lats)
    Lt, lam_t = sqrt_corr(corr_time(n_months, tau), 'C_t')
    arrays = dict(n_fields=np.int64(len(fields)), n_months=np.int64(n_months), n_cells=np.int64(lons.size),
                  fields=np.array(fields), length_km=np.array(length_km, float), tau_months=np.float64(tau),
                  Lt=Lt, meta=np.str_(meta))
    for k, (m, L) in enumerate(zip(masks, length_km)):
        i = np.flatnonzero(m)
        d = great_circle_km(u[i])
        Ch = np.exp(-d / L) if L > 0 else np.eye(i.size)
        Lh, lam = sqrt_corr(Ch, f'C_h {fields[k]}')
        arrays[f'idx_{k}'] = i
        arrays[f'Lh_{k}'] = Lh
        n_eff = lam.sum() ** 2 / (lam ** 2).sum()            # participation ratio: effective number of spatial modes
        print(f'  {fields[k]}: {i.size} cells, L = {L:g} km, eigenvalues {lam.min():.2e}..{lam.max():.2e}, '
              f'effective number of modes {n_eff:.0f}')
    print(f'  time: {n_months} months, tau = {tau:g}, eigenvalues {lam_t.min():.3f}..{lam_t.max():.3f}')
    np.savez(out, **arrays)
    print(f'  written {out}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build')
    b.add_argument('--grid-file', required=True, help='netCDF with lons, lats (nf, Ydim, Xdim): the cell centres')
    b.add_argument('--prior-file', required=True, help='prior emissions (time, nf, Ydim, Xdim), for the field masks')
    b.add_argument('--fields', default='EmisCO2_NetTerrExch,EmisCO2_Ocean')
    b.add_argument('--length-km', default='500,1000', help='horizontal correlation length per field [km]')
    b.add_argument('--tau-months', type=float, default=1.0, help='temporal correlation time [months] (0 = none)')
    b.add_argument('--n-months', type=int, required=True)
    b.add_argument('--first-month-index', type=int, default=0, help='record of the prior file = first control month')
    b.add_argument('--out', required=True)
    i = sub.add_parser('info')
    i.add_argument('file')
    a = ap.parse_args()
    if a.cmd == 'info':
        print(BCorr(a.file).summary())
        return
    import netCDF4 as nc
    g = nc.Dataset(a.grid_file)
    lons = np.asarray(g['lons'][:], float).reshape(-1)
    lats = np.asarray(g['lats'][:], float).reshape(-1)
    fields = [f for f in a.fields.split(',') if f]
    Ls = [float(v) for v in a.length_km.split(',')]
    if len(Ls) != len(fields):
        sys.exit(f'--length-km needs one value per field ({len(fields)})')
    p = nc.Dataset(a.prior_file)
    masks = []
    for f in fields:
        v = np.asarray(p[f][a.first_month_index:a.first_month_index + a.n_months], float)
        if v.shape[0] != a.n_months:
            sys.exit(f'{a.prior_file}: {f} has {v.shape[0]} records from {a.first_month_index}, need {a.n_months}')
        masks.append(np.abs(v.reshape(a.n_months, -1)).max(0) > 0)
    build(lons, lats, masks, fields, Ls, a.tau_months, a.n_months, a.out,
          meta=f'grid {a.grid_file}; masks from {a.prior_file}')


if __name__ == '__main__':
    main()
