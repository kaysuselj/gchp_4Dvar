"""OSSE 2 TER, 2016 mean (prior - truth) vs (posterior - truth) on C24: same colour scale check, RMSE, fraction of cells
improved, and the error after averaging over larger areas (scale dependence)."""
import csv, numpy as np, netCDF4 as nc
X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
t = nc.Dataset(f"{X}/obs_c24_1510_1705/truth_emissions_monthly_201510_201705.nc")
lab = [f"{d:%Y-%m}" for d in nc.num2date(t["time"][:], t["time"].units, only_use_cftime_datetimes=False)]
S = slice(lab.index("2016-01"), lab.index("2016-12") + 1)
lons, lats = np.asarray(t["lons"][:], float).ravel(), np.asarray(t["lats"][:], float).ravel()
v = "EmisCO2_NetTerrExch"
tr = np.asarray(t[v][:], float)[S].mean(0).ravel()
pr = np.asarray(nc.Dataset(f"{X}/OSSE2_C24_DL/4dvar_output_osse/emissions_monthly_iter001.nc4")[v][:], float)[S].mean(0).ravel()
posts = {}
for e in ("OSSE2_C24", "OSSE2_C24_DL"):
    h = list(csv.DictReader(open(f"{X}/{e}/4dvar_output_osse/4dvar_history.csv")))
    it = [int(r["iteration"]) for r in h if r["action"] in ("init", "accept")][-1]
    posts[f"{e} it{it}"] = np.asarray(nc.Dataset(f"{X}/{e}/4dvar_output_osse/emissions_monthly_iter{it:03d}.nc4")[v][:], float)[S].mean(0).ravel()
m = (np.abs(tr) + np.abs(pr)) > 0
k = 86400 * 1000 * 12 / 44                                     # kg CO2 m-2 s-1 -> gC m-2 day-1
dp = (pr - tr) * k
vmax = np.nanpercentile(np.abs(dp[m & (np.abs(dp) > 0)]), 98)
print(f"colour limit (98th pct of |prior-truth|) = {vmax:.3f} gC m-2 day-1; cells with TER {m.sum()}")
u = np.stack([np.cos(np.radians(lats)) * np.cos(np.radians(lons)), np.cos(np.radians(lats)) * np.sin(np.radians(lons)), np.sin(np.radians(lats))], -1)[m]
D = 6371 * np.arccos(np.clip(u @ u.T, -1, 1))
def smooth(x, R):                                              # mean over all TER cells within R km
    W = (D <= R).astype(float); return W @ x / W.sum(1)
print(f"{'':22s} {'RMSE':>6s} {'|err|>lim':>9s} {'improved':>8s}   RMSE after averaging within  300   600  1000  2000 km")
for name, x in [("prior", pr)] + list(posts.items()):
    d = ((x - tr) * k)[m]
    imp = (np.abs(d) < np.abs(dp[m])).mean() if name != "prior" else float("nan")
    sm = [np.sqrt(np.mean(smooth(d, R) ** 2)) for R in (300, 600, 1000, 2000)]
    print(f"{name:22s} {np.sqrt(np.mean(d ** 2)):6.3f} {np.mean(np.abs(d) > vmax):9.3f} {imp:8.2f}   "
          + "                     " + " ".join(f"{s:5.3f}" for s in sm))
