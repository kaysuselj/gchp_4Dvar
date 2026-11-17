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
for e, tp, op, E0 in (("OSSE1_C24", "TER/ORCHIDEE", "Ocean/ECCO-Darwin-MON-v05", 1.3), ("OSSE2_C24", "TER/CLASS-CTEM", "Ocean/MOM6-MON", 1.0)):
    st = np.load(f"{X}/{e}/4dvar_state_best.npz")
    x = st["x_next"].reshape(2, 19, -1); it = int(st["iteration"])
    for k, (tt, pp, name) in enumerate((("TER/ORCHIDEE", tp, "TER"), ("Ocean/ECCO-Darwin-MON-v05", op, "ocean"))):
        T = P = S = 0.0; e0 = e1 = 0.0
        for m, (y, mo) in enumerate(months):
            ft, fp = rd(tt, y, mo), rd(pp, y, mo) * E0
            T += pg(ft); P += pg(fp); S += pg(fp * x[k, m])
            e0 += np.abs(fp - ft).sum(); e1 += np.abs(fp * x[k, m] - ft).sum()
        print(f"{e} best iterate (it {it}) {name:5s}: 19-month total truth {T:7.2f}, prior {P:7.2f} ({P - T:+6.2f}), now {S:7.2f} ({S - T:+6.2f}) PgC; "
              f"sum|flux - truth| {100 * (e1 / e0 - 1):+6.1f} % vs prior; sigma {x[k].min():.2f}..{x[k].max():.2f}")
