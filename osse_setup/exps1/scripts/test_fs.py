import sys, numpy as np
sys.path.insert(0, sys.argv[1])
import lbfgsb_step_monthly as L
rng = np.random.default_rng(1); n = 12
g = rng.standard_normal(n) * 100; sb = np.r_[np.full(6, 0.3), np.full(6, 5.0)]
z = np.zeros((5, n))
# constant sigma_b: sigma_b-scaled step with frac f == default step with first_step = f * sigma_b
d0 = L.lbfgsb_direction(g, z, z, 0, 0.2); d1 = L.lbfgsb_direction(g, z, z, 0, 1.0, np.full(n, 0.2))
print("constant sigma_b: identical to the default step:", np.allclose(d0, d1))
d = L.lbfgsb_direction(g, z, z, 0, 1.0, sb)
print("max |d|/sigma_b =", round(float(np.abs(d / sb).max()), 6), "| descent:", float(g @ d) < 0,
      "| direction proportional to -sigma_b^2 g:", np.allclose(d / (-(sb ** 2) * g), (d / (-(sb ** 2) * g))[0]))
s_hist = rng.standard_normal((5, n)); y_hist = s_hist * 2.0
print("with memory: sigma_b option ignored:", np.allclose(L.lbfgsb_direction(g, s_hist, y_hist, 3, 1.0, sb), L.lbfgsb_direction(g, s_hist, y_hist, 3, 1.0)))
st = dict(x_anchor=None, g_anchor=None, J_anchor=None, d=None, alpha=1.0, n_bt=0, s_hist=np.zeros((5, n)), y_hist=np.zeros((5, n)),
          m_used=0, it=1, x_best=np.ones(n), J_best=np.inf, it_best=0, n_reject_total=0)
xn, st, info = L.lbfgs_backtrack_step(np.ones(n), g, 10.0, st, 5, 1e-4, 6, np.r_[np.zeros(6), np.full(6, -4.0)], np.full(n, 4.0), 1.0, sb)
print("full step: x min/max", xn.min().round(3), xn.max().round(3), "| clipped", info["n_clipped"], "| gamma", info["gamma"])
