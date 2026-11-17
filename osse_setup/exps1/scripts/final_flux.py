import sys, numpy as np, netCDF4 as nc
X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
I = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/inputs/surface_fluxes_osse"
sys.path.insert(0, "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_realcase/scripts")
from make_inputs import COORD, subcell_points
c = nc.Dataset(COORD)
lon, lat = subcell_points(c["corner_lons"][:], c["corner_lats"][:])
E = np.concatenate([[-90.0], np.arange(-88.0, 89.0, 4.0), [90.0]])
jj = np.clip(np.searchsorted(E, lat, side="right") - 1, 0, 45); ii = np.floor(((lon + 182.5) % 360.0) / 5.0).astype(int) % 72
w = np.cos(np.radians(lat))
to_c24 = lambda f: ((f[jj, ii] * w).sum(axis=(-1, -2)) / w.sum(axis=(-1, -2))).reshape(-1)
R = 6.371e6; a = w.sum(axis=(-1, -2)).reshape(-1); area = a / a.sum() * 4 * np.pi * R**2
months = [(2015, 10), (2015, 11), (2015, 12)] + [(2016, m) for m in range(1, 13)] + [(2017, m) for m in range(1, 5)]
pg = lambda f: float((f * area).sum() * 1e-6 * 30.4 * 86400 / 1e12)
rd = lambda p, y, m: to_c24(np.nan_to_num(np.squeeze(nc.Dataset(f"{I}/{p}/{y}/{m:02d}.nc")["CO2_Flux"][:]).astype(float)))
truth = {k: [rd(p, y, m) for y, m in months] for k, p in (("TER", "TER/ORCHIDEE"), ("ocean", "Ocean/ECCO-Darwin-MON-v05"))}
priors = {"OSSE1": ("TER/ORCHIDEE", "Ocean/ECCO-Darwin-MON-v05", 1.3), "OSSE2": ("TER/CLASS-CTEM", "Ocean/MOM6-MON", 1.0)}
cache = {}
for e in ("OSSE1_C24", "OSSE1_C24_DL", "OSSE2_C24", "OSSE2_C24_DL"):
    o = e.split("_")[0]; tp, op, E0 = priors[o]
    if o not in cache:
        cache[o] = {"TER": [rd(tp, y, m) * E0 for y, m in months], "ocean": [rd(op, y, m) * E0 for y, m in months]}
    st = np.load(f"{X}/{e}/4dvar_state_best.npz"); x = st["x_next"].reshape(2, 19, -1)
    out = []
    for k, f in enumerate(("TER", "ocean")):
        T = sum(pg(t) for t in truth[f]); P = sum(pg(p) for p in cache[o][f]); S = sum(pg(p * x[k, m]) for m, p in enumerate(cache[o][f]))
        e0 = sum(np.abs(p - t).sum() for p, t in zip(cache[o][f], truth[f])); e1 = sum(np.abs(p * x[k, m] - t).sum() for m, (p, t) in enumerate(zip(cache[o][f], truth[f])))
        rms0 = np.sqrt(np.mean(np.concatenate([(p - t) ** 2 for p, t in zip(cache[o][f], truth[f])])))
        rms1 = np.sqrt(np.mean(np.concatenate([(p * x[k, m] - t) ** 2 for m, (p, t) in enumerate(zip(cache[o][f], truth[f]))])))
        out.append(f"{f:5s} total err {P - T:+6.2f} -> {S - T:+6.2f} PgC, sum|err| {100 * (e1 / e0 - 1):+6.1f} %, rms {100 * (rms1 / rms0 - 1):+6.1f} %, sigma {x[k].min():.2f}..{x[k].max():.2f}")
    print(f"{e:13s} best it {int(st['iteration']):2d} J {float(st['J']):.4g} | " + " | ".join(out))
