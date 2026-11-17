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
# C24 cell areas from the sub-cell sampling weights (relative), scaled to the Earth's area
R = 6.371e6; a = w.sum(axis=(-1, -2)).reshape(-1); area = a / a.sum() * 4 * np.pi * R**2
months = [(2015, 10), (2015, 11), (2015, 12)] + [(2016, m) for m in range(1, 13)] + [(2017, m) for m in range(1, 5)]
pg = lambda f: float((f * area).sum() * 1e-6 * 30.4 * 86400 / 1e12)  # kg C km-2 s-1 -> PgC per month
for e, tp, E0 in (("OSSE1_C24", "TER/ORCHIDEE", 1.3), ("OSSE2_C24", "TER/CLASS-CTEM", 1.0)):
    x = np.load(f"{X}/{e}/4dvar_state_osse_monthly.npz")["x_next"].reshape(2, 19, -1)
    T = P = S = 0.0
    for m, (y, mo) in enumerate(months):
        ft = to_c24(np.nan_to_num(np.squeeze(nc.Dataset(f"{I}/TER/ORCHIDEE/{y}/{mo:02d}.nc")["CO2_Flux"][:]).astype(float)))
        fp = to_c24(np.nan_to_num(np.squeeze(nc.Dataset(f"{I}/{tp}/{y}/{mo:02d}.nc")["CO2_Flux"][:]).astype(float))) * E0
        T += pg(ft); P += pg(fp); S += pg(fp * x[0, m])
    print(f"{e} TER total over 19 months: truth {T:.1f}, prior {P:.1f} ({P - T:+.1f}), after step 1 {S:.1f} ({S - T:+.1f}) PgC")
