"""RMSE of the window-mean TER error vs averaging radius for the earlier OSSEs O1 (E = 1.3 everywhere) and O4 (regional
boxes), C24 vs C48, DL adjoint, Jan-Mar 2016. Truth TER = prior TER / E (prior = truth x E by construction, E =
prior_err_co2.nc4), posterior = HEMCO TER at the last accepted iterate."""
import csv, numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse"
GRID = {24: f"{B}/exps1/obs_c24_1510_1705/truth_emissions_monthly_201510_201705.nc",
        48: f"{B}/exps1/obs_c48_1510_1705/truth_emissions_monthly_201510_201705.nc"}
RADII = [0, 150, 300, 500, 750, 1000, 1500, 2000, 3000, 4000]
K = 86400 * 1000 * 12 / 44                                         # kg CO2 m-2 s-1 -> gC m-2 day-1


def unit(lo, la):
    lo, la = np.radians(lo), np.radians(la)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def tri(a, b, c):
    num = np.abs(np.einsum("...i,...i", a, np.cross(b, c)))
    den = 1 + np.einsum("...i,...i", a, b) + np.einsum("...i,...i", b, c) + np.einsum("...i,...i", c, a)
    return 2 * np.arctan2(num, den)


def grid(res):
    g = nc.Dataset(GRID[res])
    P = unit(np.asarray(g["corner_lons"][:], float), np.asarray(g["corner_lats"][:], float))
    area = (tri(P[:, :-1, :-1], P[:, 1:, :-1], P[:, 1:, 1:]) + tri(P[:, :-1, :-1], P[:, 1:, 1:], P[:, :-1, 1:])) * 6.371e6 ** 2
    return np.asarray(g["lons"][:], float).ravel(), np.asarray(g["lats"][:], float).ravel(), area.ravel()


rows, curves = [], {}
for osse, tag in (("O1", "O1_{r}_DL"), ("O4", "O4_regions_{r}_DL")):
    for res in (24, 48):
        e = tag.format(r=f"C{res}"); D = f"{B}/exps/{e}"
        lon, lat, area = grid(res)
        E = np.asarray(nc.Dataset(f"{D}/prior_err_co2.nc4")["Prior_Err"][:], float).reshape(-1)
        h = list(csv.DictReader(open(f"{D}/4dvar_output_osse/4dvar_history.csv")))
        it = [int(r["iteration"]) for r in h if r["action"] in ("init", "accept")][-1]
        rd = lambda i: np.asarray(nc.Dataset(f"{D}/4dvar_output_osse/emissions_monthly_iter{i:03d}.nc4")["EmisCO2_NetTerrExch"][:], float).mean(0).reshape(-1)
        pr, po = rd(1), rd(it)
        tr = pr / E
        m = np.abs(pr) > 0
        u = unit(lon[m], lat[m])
        dist = 6371.0 * np.arctan2(np.linalg.norm(np.cross(u[:, None], u[None]), axis=-1), np.clip(u @ u.T, -1, 1))
        w = area[m]
        out = {}
        for nm, x in (("prior", pr), ("posterior", po)):
            err = ((x - tr) * K)[m]
            s = []
            for R in RADII:
                W = (dist <= max(R, 1e-6)) * w[None, :]
                a = W @ err / W.sum(1)
                s.append(np.sqrt(np.sum(w * a ** 2) / w.sum()))
            out[nm] = np.array(s)
        curves[(osse, res)] = out
        ratio = out["posterior"] / out["prior"]
        cross = next((R for R, q in zip(RADII, ratio) if q < 1), None)
        rows.append(f"{e:18s} it {it:2d}  cells {m.sum():5d}  prior " + " ".join(f"{v:.3f}" for v in out["prior"])
                    + "\n" + " " * 42 + "post  " + " ".join(f"{v:.3f}" for v in out["posterior"])
                    + "\n" + " " * 42 + "ratio " + " ".join(f"{v:.2f}" for v in ratio))
print("window-mean TER error RMSE [gC m-2 day-1] vs averaging radius " + " ".join(str(R) for R in RADII) + " km")
print("\n".join(rows))
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
for ax, osse in zip(axs, ("O1", "O4")):
    for res, col in ((24, "#2a78d6"), (48, "#eb6834")):
        c = curves[(osse, res)]
        ax.plot(RADII, c["prior"], "--o", color=col, ms=3.5, lw=1.4, alpha=0.6, label=f"C{res} prior - truth")
        ax.plot(RADII, c["posterior"], "-o", color=col, ms=3.5, lw=1.8, label=f"C{res} posterior - truth")
    ax.set_ylim(bottom=0); ax.set_xlabel("averaging radius [km] (0 = single cell: C24 ~400 km, C48 ~200 km)", fontsize=8)
    ax.set_title(f"{osse} (Jan-Mar 2016, DL): RMSE of the window-mean TER error [gC m-2 day-1]", fontsize=9, loc="left")
    ax.grid(True, color="#e4e3df"); ax.legend(fontsize=7, frameon=False)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{B}/exps1/figures/fig_rmse_scale_O1_O4.png", dpi=150)
