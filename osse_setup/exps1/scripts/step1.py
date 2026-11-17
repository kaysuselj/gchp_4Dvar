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
months = [(2015, 10), (2015, 11), (2015, 12)] + [(2016, m) for m in range(1, 13)] + [(2017, m) for m in range(1, 5)]
for e, ter_t, ocn_t, ter_p, ocn_p, E0 in (("OSSE1_C24", "TER/ORCHIDEE", "Ocean/ECCO-Darwin-MON-v05", "TER/ORCHIDEE", "Ocean/ECCO-Darwin-MON-v05", 1.3),
                                         ("OSSE2_C24", "TER/ORCHIDEE", "Ocean/ECCO-Darwin-MON-v05", "TER/CLASS-CTEM", "Ocean/MOM6-MON", 1.0)):
    st = np.load(f"{X}/{e}/4dvar_state_osse_monthly.npz")
    x = st["x_next"].reshape(2, 19, -1)
    out = []
    for k, (tt, tp) in enumerate(((ter_t, ter_p), (ocn_t, ocn_p))):
        err0 = err1 = 0.0; n_dir = n_tot = 0
        for m, (y, mo) in enumerate(months):
            ft = to_c24(np.nan_to_num(np.squeeze(nc.Dataset(f"{I}/{tt}/{y}/{mo:02d}.nc")["CO2_Flux"][:]).astype(float)))
            fp = to_c24(np.nan_to_num(np.squeeze(nc.Dataset(f"{I}/{tp}/{y}/{mo:02d}.nc")["CO2_Flux"][:]).astype(float))) * E0
            s = x[k, m]
            prior_err, post_err = fp - ft, fp * s - ft
            err0 += np.abs(prior_err).sum(); err1 += np.abs(post_err).sum()
            moved = np.abs(s - 1) > 1e-4
            n_tot += moved.sum(); n_dir += (moved & (np.abs(post_err) < np.abs(prior_err))).sum()
        out.append(f"{['TER', 'ocean'][k]}: sigma {x[k].min():.3f}..{x[k].max():.3f} (mean {x[k].mean():.3f}); "
                   f"sum|flux - truth| prior -> step1 {100 * (err1 / err0 - 1):+.1f} %; cells moved toward truth {100 * n_dir / max(n_tot, 1):.0f} % of {n_tot}")
    print(e); [print("   ", s) for s in out]
