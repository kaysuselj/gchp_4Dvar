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
pg = lambda f: (f * area).sum(axis=-1) * 1e-6 * 30.4 * 86400 / 1e12
st = np.load(f"{X}/OSSE2_C24/4dvar_state_osse_monthly.npz")
g = st["g_anchor"].reshape(2, 19, -1)
sb = np.load(f"{X}/OSSE2_C24/sigma_b.npz")["sigma_b"].reshape(2, 19, -1)
ft = np.array([to_c24(np.nan_to_num(np.squeeze(nc.Dataset(f"{I}/TER/ORCHIDEE/{y}/{m:02d}.nc")["CO2_Flux"][:]).astype(float))) for y, m in months])
fp = np.array([to_c24(np.nan_to_num(np.squeeze(nc.Dataset(f"{I}/TER/CLASS-CTEM/{y}/{m:02d}.nc")["CO2_Flux"][:]).astype(float))) for y, m in months])
T, P = pg(ft).sum(), pg(fp).sum()
def report(name, dx):
    S = pg(fp * (1 + dx[0])).sum()
    err0 = np.abs(fp - ft).sum(); err1 = np.abs(fp * (1 + dx[0]) - ft).sum()
    print(f"{name:40s} TER total - truth {S - T:+6.1f} PgC (prior {P - T:+.1f}); sum|TER - truth| {100 * (err1 / err0 - 1):+5.1f} %; "
          f"max|dx| TER {np.abs(dx[0]).max():.2f} ocean {np.abs(dx[1]).max():.2f}")
report("now: 0.571 x g / max|g|", -0.571 * g / np.abs(g).max())
p = sb ** 2 * g
for lab, scale in (("sigma_b^2 g, max|dx/sigma_b| = 1", 1 / np.abs(p / sb).max()),
                   ("sigma_b^2 g, max|dx/sigma_b| = 0.5", 0.5 / np.abs(p / sb).max())):
    report(lab, -scale * p)
for fs in (0.2, 0.17):
    report(f"scalar {fs} x g / max|g|", -fs * g / np.abs(g).max())
