import numpy as np, netCDF4 as nc
X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
d = nc.Dataset(f"{X}/obs_c24_1510_1705/semi_obs_xco2_201510_201705.nc")
n = np.asarray(d["xco2_noise"][:], float); u = np.asarray(d["xco2_uncertainty"][:], float) * 1e6
z = n / u
N = z.size
print(f"N = {N}; J_obs(truth) = 0.5 sum (noise/sigma)^2 = {0.5 * (z ** 2).sum():.4e}  (expected N/2 = {N / 2:.4e}); mean z^2 {np.mean(z ** 2):.4f}")
# J_b at the truth: true scale factor = truth / prior per cell and month, sigma_b from the file; same bounds as the run
t = nc.Dataset(f"{X}/obs_c24_1510_1705/truth_emissions_monthly_201510_201705.nc")
p = nc.Dataset(f"{X}/OSSE2_C24_DL/4dvar_output_osse/emissions_monthly_iter001.nc4")
sb = np.load(f"{X}/OSSE2_C24_DL/sigma_b.npz")["sigma_b"].reshape(2, 19, 6, 24, 24)
tot = 0
for k, (v, lo, hi) in enumerate((("EmisCO2_NetTerrExch", 0, 4), ("EmisCO2_Ocean", -4, 4))):
    T = np.asarray(t[v][:], float); P = np.asarray(p[v][:], float)
    m = np.abs(P) > 0
    x = np.ones_like(P); x[m] = T[m] / P[m]
    xc = np.clip(x, lo, hi)
    jb = 0.5 * (((xc - 1) / sb[k]) ** 2)[m].sum()
    tot += jb
    print(f"{v}: cells with flux {m.sum()}, true SF p5/p50/p95 {np.percentile(x[m], [5, 50, 95]).round(2)}, "
          f"clipped {(x[m] != xc[m]).sum()}, J_b(truth) {jb:.4e}")
print(f"J_b(truth) total {tot:.4e};  J(truth) approx {0.5 * (z ** 2).sum() + tot:.4e}")
