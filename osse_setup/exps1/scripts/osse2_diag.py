import sys, numpy as np, netCDF4 as nc
X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
I = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/inputs/surface_fluxes_osse"
# (1) split the predicted decrease g.d of the rejected direction by field
st = np.load(f"{X}/OSSE2_C24/4dvar_state_osse_monthly.npz")
g, d = st["g_anchor"].reshape(2, -1), st["d_cur"].reshape(2, -1)
gd = [float(g[k] @ d[k]) for k in (0, 1)]
print(f"rejected direction (anchor = iteration 5): g.d TER {gd[0]:.4g}, ocean {gd[1]:.4g} (total {sum(gd):.4g}); "
      f"max|d| TER {np.abs(d[0]).max():.3g} ocean {np.abs(d[1]).max():.3g}; alpha now {float(st['alpha']):.4g}")
# (2) J_b at the truth scale factors (TER ORCHIDEE/CTEM, ocean ECCO/MOM6 on C24, clipped to the bounds)
sys.path.insert(0, "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_realcase/scripts")
from make_inputs import COORD, subcell_points
c = nc.Dataset(COORD)
lon, lat = subcell_points(c["corner_lons"][:], c["corner_lats"][:])
E = np.concatenate([[-90.0], np.arange(-88.0, 89.0, 4.0), [90.0]])
jj = np.clip(np.searchsorted(E, lat, side="right") - 1, 0, 45); ii = np.floor(((lon + 182.5) % 360.0) / 5.0).astype(int) % 72
w = np.cos(np.radians(lat))
to_c24 = lambda f: ((f[jj, ii] * w).sum(axis=(-1, -2)) / w.sum(axis=(-1, -2))).reshape(-1)
rd = lambda p, y, m: to_c24(np.nan_to_num(np.squeeze(nc.Dataset(f"{I}/{p}/{y}/{m:02d}.nc")["CO2_Flux"][:]).astype(float)))
months = [(2015, 10), (2015, 11), (2015, 12)] + [(2016, m) for m in range(1, 13)] + [(2017, m) for m in range(1, 5)]
sb = np.load(f"{X}/OSSE2_C24/sigma_b.npz")["sigma_b"].reshape(2, 19, -1)
xb = np.load(f"{X}/OSSE2_C24/4dvar_state_best.npz")["x_next"].reshape(2, 19, -1)
Jb_true = Jb_best = 0.0; lo, hi = (0.0, -4.0), (4.0, 4.0)
for k, (tr, pr) in enumerate((("TER/ORCHIDEE", "TER/CLASS-CTEM"), ("Ocean/ECCO-Darwin-MON-v05", "Ocean/MOM6-MON"))):
    jt = jbst = 0.0
    for m, (y, mo) in enumerate(months):
        ft, fp = rd(tr, y, mo), rd(pr, y, mo)
        with np.errstate(divide="ignore", invalid="ignore"):
            s = np.where(np.abs(fp) > 0, ft / fp, 1.0)
        s = np.clip(np.nan_to_num(s, nan=1.0, posinf=hi[k], neginf=lo[k]), lo[k], hi[k])
        jt += 0.5 * np.sum(((s - 1) / sb[k, m]) ** 2)
        jbst += 0.5 * np.sum(((xb[k, m] - 1) / sb[k, m]) ** 2)
    print(f"J_b {['TER', 'ocean'][k]:5s}: at the truth sigma (clipped to the bounds) {jt:.4g}; at the best iterate {jbst:.4g}")
    Jb_true += jt; Jb_best += jbst
print(f"J at the truth ~ J_obs(noise only) 174558 + J_b {Jb_true:.4g} = {174558 + Jb_true:.4g} (lower bound if clipping changes J_obs); "
      f"best iterate J = 521824 (J_b {Jb_best:.4g})")
