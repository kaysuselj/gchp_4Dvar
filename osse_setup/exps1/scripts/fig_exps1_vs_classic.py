#!/usr/bin/env python3
"""exps1 (2026-11-16): our C24 OSSE 2 (GCHP, continuous and discrete advection adjoint) vs J. Liu's GEOS-Chem Classic
CMS-Flux OSSE LNLG2_TER2 (4x5), TER, 2016. Both start from the same prior TER (CLASS-CTEM) and have the same truth TER
(ORCHIDEE). Her posterior = EMS-sf-36.nc (evaluation 36 = iteration 31, the best J, 1.571e5) x CLASS-CTEM TER.
Everything is compared on her 4x5 grid: our C24 prior / posterior / truth (HEMCO diagnostics, PgC per cell per month) are
summed into the 4x5 cells with 12 x 12 sample points per C24 cell (each carries 1/144 of the cell total), the same
sampling as make_exps1_inputs.py.
Differences between the two set-ups (see README): her observations LN + LG only (ours LN + LG + OG), her window
2016-01 .. 2019-06 (ours 2015-10 .. 2017-04), her ocean prior = truth and not optimized (ours MOM6, optimized), her scale
factors unbounded (ours TER [0, 4]), 4x5 Classic vs C24 GCHP.
  fig_classic_monthly.png   2016 monthly global TER bias and 4x5 cell RMSE (prior, her posterior, our posteriors)
  fig_classic_regional.png  2016 TransCom land-region TER totals minus truth
  fig_classic_maps.png      2016 mean TER posterior minus truth on 4x5: hers and ours
  classic_compare.txt       the numbers
    fig_exps1_vs_classic.py <out dir>
"""
import csv, os, sys
import numpy as np, netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import shapefile

sys.path.insert(0, "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_realcase/scripts")
from make_inputs import subcell_points

X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
TRUTH = f"{X}/obs_c24_1510_1705/truth_emissions_monthly_201510_201705.nc"
JL = ("/nobackupp17/jliu7/ADJOINT/MEEMONG_NEW2/GHGF-CMS-16Dec2019/App_4x5/CARDAMOM-ECCO-biweekly/TRENDY/CLASS-CTEM/"
      "top-down/LNLG2_TER2/GDT-EMS/EMS-sf-36.nc")
V4 = "/nobackupp17/jliu7/INVENTORY/V4/TER2"
TC = "/nobackupp17/jliu7/ADJOINT/utility/data/transcom.regrid.land.ocn.nc"
COAST = os.path.expanduser("~/.local/share/cartopy/shapefiles/natural_earth/physical/ne_110m_coastline.shp")
RUNS = [("OSSE2_C24", "GCHP C24 continuous"), ("OSSE2_C24_DL", "GCHP C24 discrete")]
COL = {"GCHP C24 continuous": "#2a78d6", "GCHP C24 discrete": "#eb6834", "Classic 4x5 (J. Liu)": "#1baf7a"}
INK, INK2, GRID, PRIOR = "#0b0b0b", "#52514e", "#e4e3df", "#9a9993"
TCNAME = {1: "N. Am. Boreal", 2: "N. Am. Temperate", 3: "S. Am. Tropical", 4: "S. Am. Temperate", 5: "N. Africa",
          6: "S. Africa", 7: "Eurasia Boreal", 8: "Eurasia Temperate", 9: "Tropical Asia", 10: "Australia", 11: "Europe"}
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rep = []
KG2PGC = 12.0 / 44.0 / 1e12


def style(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    if grid:
        ax.grid(True, axis=grid, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


# ---- 4x5 grid (GEOS-Chem: half-size polar rows) ----
LAT_EDGES = np.concatenate([[-90.0], np.arange(-88.0, 89.0, 4.0), [90.0]])
R = 6.371e6
a45 = (R ** 2 * np.radians(5.0) * np.diff(np.sin(np.radians(LAT_EDGES))))[:, None] * np.ones((1, 72))     # m2
rep.append(f"4x5 total area {a45.sum():.4e} m2")
mdays = np.array([31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], float)

# ---- J. Liu: prior (CLASS-CTEM), truth (ORCHIDEE), posterior = SF x prior; kg C km-2 s-1 -> PgC per cell per month ----
to_pgc = lambda f, m: f * 1e-6 * a45 * mdays[m] * 86400 / 1e12
jl_prior = np.array([to_pgc(np.asarray(nc.Dataset(f"{V4}/CLASS-CTEM/2016/{m + 1:02d}.nc")["CO2_Flux"][:], float), m) for m in range(12)])
jl_truth = np.array([to_pgc(np.asarray(nc.Dataset(f"{V4}/ORCHIDEE/2016/{m + 1:02d}.nc")["CO2_Flux"][:], float), m) for m in range(12)])
sf = np.asarray(nc.Dataset(JL)["EMS-01"][:12], float)                  # mmscl 0..11 = 2016-01 .. 2016-12
jl_post = sf * jl_prior
rep.append(f"Classic SF 2016: min {sf.min():.3f} max {sf.max():.3f}; prior {jl_prior.sum():.2f}, truth {jl_truth.sum():.2f}, "
           f"posterior {jl_post.sum():.2f} PgC")

# ---- ours on C24 -> 4x5 ----
t = nc.Dataset(TRUTH)
clon, clat = np.asarray(t["corner_lons"][:], float), np.asarray(t["corner_lats"][:], float)
tt = nc.num2date(t["time"][:], t["time"].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
lab19 = [f"{d:%Y-%m}" for d in tt]
SEL = slice(lab19.index("2016-01"), lab19.index("2016-12") + 1)
md19 = np.array([((d.replace(month=d.month % 12 + 1, year=d.year + d.month // 12)) - d).days for d in tt], float)


def unit(lo, la):
    lo, la = np.radians(lo), np.radians(la)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def tri(a, b, c):
    num = np.abs(np.einsum("...i,...i", a, np.cross(b, c)))
    den = 1 + np.einsum("...i,...i", a, b) + np.einsum("...i,...i", b, c) + np.einsum("...i,...i", c, a)
    return 2 * np.arctan2(num, den)


P = unit(clon, clat)
p00, p10, p01, p11 = P[:, :-1, :-1], P[:, 1:, :-1], P[:, :-1, 1:], P[:, 1:, 1:]
area24 = (tri(p00, p10, p11) + tri(p00, p11, p01)) * R ** 2
slon, slat = subcell_points(clon, clat)                                 # (6, 24, 24, 12, 12)
jj = np.clip(np.searchsorted(LAT_EDGES, slat, side="right") - 1, 0, 45).reshape(-1)
ii = (np.floor(((slon + 182.5) % 360.0) / 5.0).astype(int) % 72).reshape(-1)
NS = slon.shape[-1] * slon.shape[-2]


def c24_to_45(v):
    """(12, 6, 24, 24) PgC per C24 cell -> (12, 46, 72) PgC per 4x5 cell"""
    o = np.zeros((v.shape[0], 46, 72))
    for m in range(v.shape[0]):
        w = np.repeat(v[m].reshape(-1), NS) / NS
        np.add.at(o[m], (jj, ii), w)
    return o


def read_ter(path):
    f = np.asarray(nc.Dataset(path)["EmisCO2_NetTerrExch"][:], float)
    return c24_to_45((f * area24[None] * (md19[:, None, None, None] * 86400) * KG2PGC)[SEL])


g_truth = read_ter(TRUTH)
ours = {}
for e, lab in RUNS:
    D = f"{X}/{e}"
    hist = list(csv.DictReader(open(f"{D}/4dvar_output_osse/4dvar_history.csv")))
    acc = [int(h["iteration"]) for h in hist if h["action"] in ("init", "accept")]
    final = os.path.exists(f"{D}/CONVERGED_monthly")
    lab2 = lab + ("" if final else f" (not final, iteration {int(hist[-1]['iteration'])})")
    ours[lab2] = read_ter(f"{D}/4dvar_output_osse/emissions_monthly_iter{acc[-1]:03d}.nc4")
g_prior = read_ter(f"{X}/OSSE2_C24_DL/4dvar_output_osse/emissions_monthly_iter001.nc4")
rep.append(f"GCHP 2016 on 4x5: prior {g_prior.sum():.2f}, truth {g_truth.sum():.2f} PgC")
# consistency of the two set-ups: same prior and truth files, so the 4x5 sums should agree
for nm, a, b in (("prior", g_prior, jl_prior), ("truth", g_truth, jl_truth)):
    m = (np.abs(a).sum(0) + np.abs(b).sum(0)) > 0
    rep.append(f"  {nm}: GCHP vs Classic 4x5 2016 total {a.sum():.2f} vs {b.sum():.2f}; cell correlation of the 2016 totals "
               f"{np.corrcoef(a.sum(0)[m], b.sum(0)[m])[0, 1]:.4f}")

# errors: each system against its own truth (identical prior and truth files; the grid mapping differs slightly)
LJ = "Classic 4x5 (J. Liu)"
cases = [("prior (CLASS-CTEM)", g_prior - g_truth, jl_prior - jl_truth)]
err = {LJ: jl_post - jl_truth}
err.update({k: v - g_truth for k, v in ours.items()})
mask = (np.abs(jl_truth).sum(0) + np.abs(jl_prior).sum(0)) > 0
w = (a45[mask] / a45[mask].sum())
per_m2 = 1e15 / a45[None] / mdays[:, None, None]                       # PgC per cell per month -> gC m-2 day-1
bias = lambda e: e.sum(axis=(1, 2))
rmse = lambda e: np.sqrt((((e * per_m2)[:, mask]) ** 2 * w).sum(1))
rmse_y = lambda e: np.sqrt((((e.sum(0) * 1e15 / a45 / 366)[mask]) ** 2 * w).sum())

rep.append("\n2016 TER, each against its own truth, on 4x5 (bias = global total minus truth; RMSE area-weighted over the "
           "cells with TER, gC m-2 day-1)")
rows = [("prior, GCHP mapping", cases[0][1]), ("prior, Classic grid", cases[0][2])] + list(err.items())
for nm, e in rows:
    rep.append(f"  {nm:48s} annual {e.sum():+6.2f} PgC; monthly |bias| mean {np.abs(bias(e)).mean():.3f} PgC; "
               f"monthly RMSE mean {rmse(e).mean():.3f}; annual-mean RMSE {rmse_y(e):.3f}")
    p0 = cases[0][2] if nm.startswith(("Classic", "prior, Classic")) else cases[0][1]
    rep.append(f"  {'':48s} monthly RMSE vs own prior {100 * (rmse(e).mean() / rmse(p0).mean() - 1):+.1f} %, "
               f"annual-mean RMSE vs own prior {100 * (rmse_y(e) / rmse_y(p0) - 1):+.1f} %")
    rep.append(f"  {'':48s} bias " + " ".join(f"{v:+.2f}" for v in bias(e)) + " | RMSE " + " ".join(f"{v:.3f}" for v in rmse(e)))

# ---- figure: monthly bias and RMSE ----
xm = np.arange(12)
fig, axs = plt.subplots(1, 2, figsize=(13, 4.2))
for ax, fn, ttl in ((axs[0], bias, "2016 TER monthly bias (global total minus truth) [PgC per month]"),
                    (axs[1], rmse, "2016 TER monthly 4x5 cell RMSE (area-weighted) [gC m-2 day-1]")):
    ax.plot(xm, fn(cases[0][2]), "--o", color=PRIOR, ms=3, lw=1.8, label="prior (CLASS-CTEM) - truth, Classic 4x5")
    if fn is rmse:      # the C24 -> 4x5 mapping smooths the field, so the GCHP prior RMSE on 4x5 is lower: plot it too
        ax.plot(xm, fn(cases[0][1]), ":s", color=PRIOR, ms=3, lw=1.8, label="prior - truth, GCHP C24 mapped to 4x5")
    for nm, e in err.items():
        ax.plot(xm, fn(e), "-o", color=COL[nm.split(" (not")[0]], ms=3, lw=1.8, label=f"posterior - truth, {nm}")
    if fn is bias:
        ax.axhline(0, color=INK, lw=0.8)
    else:
        ax.set_ylim(bottom=0)
    ax.set_xticks(xm); ax.set_xticklabels([f"{m + 1:02d}" for m in xm]); ax.set_xlabel("month of 2016", fontsize=8, color=INK2)
    ax.set_title(ttl, fontsize=9, color=INK, loc="left"); style(ax)
axs[0].legend(fontsize=7, frameon=False); axs[1].legend(fontsize=7, frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(f"{out}/fig_classic_monthly.png", dpi=150); plt.close(fig)

# ---- figure: TransCom land regions ----
tc = nc.Dataset(TC)
mland = np.asarray(tc["mask"][:]).squeeze().astype(int)
assert mland.shape == (46, 72)
ids = list(range(1, 12))
reg = lambda e: np.array([e.sum(0)[mland == r].sum() for r in ids])
fig, ax = plt.subplots(figsize=(13, 4.6))
xg = np.arange(len(ids)); wb = 0.2
ax.bar(xg - 1.5 * wb, reg(cases[0][2]), wb, color=PRIOR, label="prior (CLASS-CTEM)")
for j, (nm, e) in enumerate(err.items()):
    ax.bar(xg + (j - 0.5) * wb, reg(e), wb, color=COL[nm.split(" (not")[0]], label=f"posterior, {nm}")
ax.axhline(0, color=INK, lw=0.8)
ax.set_xticks(xg); ax.set_xticklabels([f"{r} {TCNAME[r]}" for r in ids], rotation=35, ha="right", fontsize=8)
ax.set_title("2016 TER, TransCom land regions: total minus truth [PgC]", fontsize=9, color=INK, loc="left")
style(ax); ax.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{out}/fig_classic_regional.png", dpi=150); plt.close(fig)
rep.append("\nTransCom land regions, 2016 total minus truth [PgC]: " + " ".join(f"{r:>6d}" for r in ids) + "  sum|err|")
for nm, e in [("prior (Classic grid)", cases[0][2]), ("prior (GCHP mapping)", cases[0][1])] + list(err.items()):
    v = reg(e); rep.append(f"  {nm:48s} " + " ".join(f"{x:+6.2f}" for x in v) + f"  {np.abs(v).sum():6.2f}")

# ---- figure: maps ----
fig, axs = plt.subplots(2, 2, figsize=(13, 7.6), gridspec_kw=dict(hspace=0.45, wspace=0.12))
panels = [("prior (CLASS-CTEM) - truth", cases[0][2])] + [(f"posterior - truth\n{nm}", e) for nm, e in err.items()]
vmax = np.percentile(np.abs((cases[0][2].sum(0) * 1e15 / a45 / 366)[mask]), 98)
for ax, (ttl, e) in zip(axs.ravel(), panels):
    v = np.ma.masked_where(~mask, e.sum(0) * 1e15 / a45 / 366)
    pc = ax.pcolormesh(np.arange(-182.5, 178, 5.0), LAT_EDGES, v, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="flat")
    if os.path.exists(COAST):
        for s in shapefile.Reader(COAST).shapes():
            p = np.array(s.points); ax.plot(p[:, 0], p[:, 1], color=INK2, lw=0.4)
    ax.set_facecolor("#efeeeb"); ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_aspect("equal")
    ax.set_xticks(range(-180, 181, 90)); ax.set_yticks(range(-90, 91, 45))
    ax.set_title(f"2016 mean TER: {ttl}  (RMSE {rmse_y(e):.3f})", fontsize=9, color=INK, loc="left"); style(ax, None)
cb = fig.colorbar(pc, ax=axs.ravel().tolist(), orientation="horizontal", shrink=0.5, pad=0.07, aspect=40)
cb.set_label("2016 mean TER difference on 4x5 [gC m-2 day-1]", fontsize=8, color=INK2); cb.ax.tick_params(labelsize=7)
fig.savefig(f"{out}/fig_classic_maps.png", dpi=150, bbox_inches="tight"); plt.close(fig)

open(f"{out}/classic_compare.txt", "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
