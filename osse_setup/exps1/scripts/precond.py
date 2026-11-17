import numpy as np
X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
for e, sbf in (("OSSE1_C24", None), ("OSSE2_C24", f"{X}/OSSE2_C24/sigma_b.npz")):
    st = np.load(f"{X}/{e}/4dvar_state_osse_monthly.npz")
    g = st["g_anchor"].reshape(2, -1)                         # gradient at the prior (iteration 1), J_b = 0
    sb = (np.load(sbf)["sigma_b"] if sbf else np.full(g.size, 0.2)).reshape(2, -1)
    print(e)
    for k, f in enumerate(("TER", "ocean")):
        nz = np.abs(g[k]) > 0
        print(f"  {f:5s}: |g| max {np.abs(g[k]).max():.3g}, median of nonzero {np.median(np.abs(g[k][nz])):.3g}; "
              f"sigma_b median {np.median(sb[k][nz]):.2f}")
    # current first step: dx = -0.571 (median sigma_b) * g / max|g|  (or 0.2 for OSSE1)
    fs = 0.2 if sbf is None else float(np.median(sb))
    dx_now = -fs * g / np.abs(g).max()
    # sigma_b-scaled (preconditioned) steepest descent: dx = -c * sb^2 * g, c so that max |dx / sb| = 1
    p = sb ** 2 * g
    c = 1.0 / np.abs(p / sb).max()
    dx_pc = -c * p
    for name, dx in (("now (scalar)", dx_now), ("sigma_b^2-scaled", dx_pc)):
        r = [np.abs(dx[k]).max() for k in (0, 1)]
        m = [np.abs(dx[k][np.abs(g[k]) > 0]).mean() for k in (0, 1)]
        print(f"  first step {name:17s}: max|dx| TER {r[0]:.3f} ocean {r[1]:.3f}; mean|dx| TER {m[0]:.4f} ocean {m[1]:.4f}")
