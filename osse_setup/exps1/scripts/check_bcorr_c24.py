"""Real-size check of the C24 correlated-B operator with the OSSE 1 iteration-1 gradient (archived by the optimizer)."""
import sys, time, numpy as np, netCDF4 as nc
sys.path.insert(0, "/nobackupp27/ksuselj1/gchp_14.8.0_adj/gchp_4Dvar")
from b_corr import BCorr
X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
B = BCorr(f"{X}/bcorr/b_corr_c24_L500_1000_tau1.npz"); print(B.summary())
n = B.n; rng = np.random.default_rng(1)
t0 = time.time(); v = rng.standard_normal(n); Lv = B.apply_L(v); w = rng.standard_normal(n); LTw = B.apply_LT(w)
print(f"apply_L + apply_LT on n = {n}: {time.time() - t0:.2f} s")
print(f"adjoint test: {Lv @ w:.12e} vs {v @ LTw:.12e}")
# sample variance of L v (should be ~1 inside the masks)
V = np.array([B.apply_L(rng.standard_normal(n)) for _ in range(200)])
var = V.var(0).reshape(B.n_fields, B.n_months, B.n_cells)
for k in range(B.n_fields):
    print(f"  field {k}: sample var of (L v) in mask {var[k][:, B.idx[k]].mean():.3f} (expect 1)")
# iteration-1 obs gradient from the state archive (g_anchor of iteration 1 = full gradient at the prior; J_b grad = 0)
st = np.load(f"{X}/OSSE1_C24_DL/4dvar_output_osse/4dvar_state_iter001.npz")
g = st["g_anchor"].astype(np.float64); assert g.size == n
sig = 0.2
gv = B.apply_LT(sig * g)
# diagonal first step (sigma_b-scaled): d_x = -c sig^2 g with max|d_x|/sig = 1
dx_diag = -(1.0 / np.abs(sig * g).max()) * sig ** 2 * g
# correlated first step: d_v = -c g_v with max |sig L d_v| = 1 sigma_b
Ld = B.apply_L(gv); dx_corr = -sig * Ld / np.abs(Ld).max()
for name, d in (("diagonal", dx_diag), ("correlated", dx_corr)):
    D = d.reshape(B.n_fields, B.n_months, B.n_cells)
    for k, f in enumerate(("TER", "ocean")):
        i = B.idx[k]
        print(f"  first step {name:10s} {f:5s}: max|dx| {np.abs(D[k][:, i]).max():.3f}, mean dx {D[k][:, i].mean():+.4f}, "
              f"rms {np.sqrt((D[k][:, i] ** 2).mean()):.4f}")
print(f"  J_obs decrease predicted by the first step (g.dx): diagonal {g @ dx_diag:.4e}, correlated {g @ dx_corr:.4e}")
