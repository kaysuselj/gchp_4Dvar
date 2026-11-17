#!/usr/bin/env python3
"""exps1 C24 OSSE results (2026-11-16): OSSE1 and OSSE2, continuous (RW, ADJ_METHOD 3) and discrete (DL, ADJ_METHOD 4)
advection adjoints. Fluxes are the ones the model applied (HEMCO diagnostics archived per iteration,
4dvar_output_osse/emissions_monthly_iterNNN.nc4, native C24): TER = EmisCO2_NetTerrExch (category 5), ocean =
EmisCO2_Ocean; truth = obs_c24_1510_1705/truth_emissions_monthly_201510_201705.nc; prior = iteration 1 (sigma = 1).
"Posterior" = the latest accepted iterate (= the best one). Runs still running are labelled "not final (iteration N)".
The fluxes are analyzed for 2016 only (January-December); Oct-Dec 2015 and Jan-Apr 2017 are spin-up / spin-down of the
19-month window. J and the cost-function figure are over the whole window (that is what was minimized).
Figures (in <out>):
  fig_cost.png            J (log) per iteration with the prior -> final drop, J_obs and J_b, and 2 J_obs / N_obs
  fig_flux_error.png      per iteration: 2016 global total error and cell RMSE (relative to the prior), TER and ocean
  fig_transcom_map.png    the 22 TransCom regions (J. Liu, transcom.regrid.land.ocn.nc, 4x5, mapped to C24 centres)
  fig_regional.png        2016 regional totals, prior and posterior minus truth, TransCom land (TER) / ocean regions
  fig_latbands.png        the same for latitude bands
  fig_monthly_global.png  monthly global TER and ocean: truth, prior, posterior
  fig_monthly_bias_rmse.png  per month of 2016: global bias and area-weighted cell RMSE, prior and posterior
  fig_maps_ter.png        2016 mean TER, prior and posterior minus truth (maps)
  results.txt             all numbers
    fig_exps1_c24_results.py <out dir>
"""
import csv, os, sys
import numpy as np, netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import shapefile

X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
TRUTH = f"{X}/obs_c24_1510_1705/truth_emissions_monthly_201510_201705.nc"
TC = "/nobackupp17/jliu7/ADJOINT/utility/data/transcom.regrid.land.ocn.nc"
COAST = os.path.expanduser("~/.local/share/cartopy/shapefiles/natural_earth/physical/ne_110m_coastline.shp")
RUNS = [("OSSE1_C24", "OSSE 1", "continuous"), ("OSSE1_C24_DL", "OSSE 1", "discrete"),
        ("OSSE2_C24", "OSSE 2", "continuous"), ("OSSE2_C24_DL", "OSSE 2", "discrete")]
COL = {"continuous": "#2a78d6", "discrete": "#eb6834"}
INK, INK2, GRID, PRIOR = "#0b0b0b", "#52514e", "#e4e3df", "#9a9993"
TCNAME = {1: "N. Am. Boreal", 2: "N. Am. Temperate", 3: "S. Am. Tropical", 4: "S. Am. Temperate", 5: "N. Africa",
          6: "S. Africa", 7: "Eurasia Boreal", 8: "Eurasia Temperate", 9: "Tropical Asia", 10: "Australia", 11: "Europe",
          12: "N. Pacific Temp.", 13: "W. Pacific Trop.", 14: "E. Pacific Trop.", 15: "S. Pacific Temp.",
          16: "Northern Ocean", 17: "N. Atlantic Temp.", 18: "Atlantic Trop.", 19: "S. Atlantic Temp.",
          20: "Southern Ocean", 21: "Indian Trop.", 22: "S. Indian Temp."}
BANDS = [(-90, -30), (-30, 0), (0, 30), (30, 60), (60, 90)]
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rep = []


def style(ax, grid="y"):
    from matplotlib.ticker import MaxNLocator
    if ax.get_xlabel().startswith("iteration") or grid == "both":
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    if grid:
        ax.grid(True, axis=grid, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


# ---- grid: C24 cell areas (spherical quadrilaterals from the corners) and centres ----
t = nc.Dataset(TRUTH)
clon, clat = np.asarray(t["corner_lons"][:], float), np.asarray(t["corner_lats"][:], float)
lons, lats = np.asarray(t["lons"][:], float), np.asarray(t["lats"][:], float)
tt = nc.num2date(t["time"][:], t["time"].units, only_use_cftime_datetimes=False, only_use_python_datetimes=True)
mdays19 = np.array([((d.replace(month=d.month % 12 + 1, year=d.year + d.month // 12)) - d).days for d in tt], float)
SEL = slice([f"{d:%Y-%m}" for d in tt].index("2016-01"), [f"{d:%Y-%m}" for d in tt].index("2016-12") + 1)
mdays = mdays19[SEL]
labels_m = [f"{d:%Y-%m}" for d in tt][SEL]
assert len(labels_m) == 12 and labels_m[0] == "2016-01" and labels_m[-1] == "2016-12"
NM = len(labels_m)


def unit(lo, la):
    lo, la = np.radians(lo), np.radians(la)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def tri(a, b, c):
    num = np.abs(np.einsum("...i,...i", a, np.cross(b, c)))
    den = 1 + np.einsum("...i,...i", a, b) + np.einsum("...i,...i", b, c) + np.einsum("...i,...i", c, a)
    return 2 * np.arctan2(num, den)


P = unit(clon, clat)
p00, p10, p01, p11 = P[:, :-1, :-1], P[:, 1:, :-1], P[:, :-1, 1:], P[:, 1:, 1:]
area = (tri(p00, p10, p11) + tri(p00, p11, p01)) * 6.371e6 ** 2          # (6, 24, 24) m2
rep.append(f"C24 total area {area.sum():.4e} m2 (Earth 5.1006e14)")
KG2PGC = 12.0 / 44.0 / 1e12


def totals(f):
    """(19, 6, 24, 24) kg CO2 m-2 s-1 -> PgC per cell per month, 2016 only (12, 6, 24, 24)"""
    return (f * area[None] * (mdays19[:, None, None, None] * 86400) * KG2PGC)[SEL]


def read_flux(path):
    d = nc.Dataset(path)
    return {"TER": totals(np.asarray(d["EmisCO2_NetTerrExch"][:], float)),
            "ocean": totals(np.asarray(d["EmisCO2_Ocean"][:], float))}


truth = read_flux(TRUTH)
rep.append("analysis period 2016-01..2016-12 (12 months)")
rep.append("truth 2016 totals: TER %.2f, ocean %.2f PgC" % (truth["TER"].sum(), truth["ocean"].sum()))

# ---- TransCom ids on C24 (nearest 4x5 cell of each C24 centre) ----
tc = nc.Dataset(TC)
mland, mocn = np.asarray(tc["mask"][:]).squeeze().astype(int), np.asarray(tc["mask_ocn"][:]).squeeze().astype(int)
la4 = np.asarray(tc["lat"][:], float); lo4 = np.asarray(tc["lon"][:], float)
jj = np.abs(lats[..., None] - la4).argmin(-1)
ii = np.abs(((lons[..., None] - lo4 + 180) % 360) - 180).argmin(-1)
rid = {"TER": mland[jj, ii], "ocean": mocn[jj, ii]}

# ---- runs: history, per-iteration fluxes ----
runs = {}
for e, osse, adj in RUNS:
    D = f"{X}/{e}"
    hist = list(csv.DictReader(open(f"{D}/4dvar_output_osse/4dvar_history.csv")))
    final = os.path.exists(f"{D}/CONVERGED_monthly")
    last = int(hist[-1]["iteration"])
    lab = f"{adj}" + ("" if final else f" (not final, iteration {last})")
    it_acc = [int(h["iteration"]) for h in hist if h["action"] in ("init", "accept")]
    it_post = it_acc[-1]
    err = {k: [] for k in ("TER", "ocean")}
    rms = {k: [] for k in ("TER", "ocean")}
    flux = {}
    for it in it_acc:
        f = f"{D}/4dvar_output_osse/emissions_monthly_iter{it:03d}.nc4"
        if it == 1 and not os.path.exists(f):
            # iteration 1 = the prior (sigma = 1); the RW runs' iteration-1 forward was sampled by hand (2026-11-11) and
            # not archived, so take the DL run's: same prior, same forward model (J identical to 10 digits)
            f = f"{X}/{e}_DL/4dvar_output_osse/emissions_monthly_iter001.nc4"
        fl = read_flux(f)
        flux[it] = fl
        for k in ("TER", "ocean"):
            err[k].append(fl[k].sum() - truth[k].sum())
            m = np.abs(truth[k]) + np.abs(fl[k]) > 0
            rms[k].append(np.sqrt(np.mean((fl[k] - truth[k])[m] ** 2)))
    runs[e] = dict(osse=osse, adj=adj, hist=hist, final=final, last=last, lab=lab, it_acc=it_acc, it_post=it_post,
                   err=err, rms=rms, prior=flux[1], post=flux[it_post])
    rep.append(f"{e}: {'final' if final else 'running'}, last iteration {last}, posterior = iteration {it_post}, "
               f"J {float(hist[0]['J']):.4g} -> {min(float(h['J']) for h in hist if h['action'] in ('init', 'accept')):.4g}")

# ---- cost function at the truth (reference; same as scripts/j_truth.py) ----
OBS = f"{X}/obs_c24_1510_1705/semi_obs_xco2_201510_201705.nc"


def j_truth(e):
    """J_obs, J_b (TER, ocean) and J at the true fluxes for run e (its _DL twin has the prior and sigma_b).
    Identical-twin OSSE: J_obs(truth) = the added noise. J_b(truth): SF_true = truth / prior per cell and month (cells
    with prior flux), clipped to the optimizer bounds (TER [0, 4], ocean [-4, 4]), so J_b(truth) is a lower bound."""
    d = nc.Dataset(OBS)
    z = np.asarray(d["xco2_noise"][:], float) / (np.asarray(d["xco2_uncertainty"][:], float) * 1e6)
    jt = dict(N=z.size, J_obs=0.5 * (z ** 2).sum())
    D = f"{X}/{e}_DL" if not e.endswith("_DL") else f"{X}/{e}"
    p = nc.Dataset(f"{D}/4dvar_output_osse/emissions_monthly_iter001.nc4")
    sbf = f"{D}/sigma_b.npz"
    sb = np.load(sbf)["sigma_b"].reshape(2, 19, 6, 24, 24) if os.path.exists(sbf) else np.full((2, 19, 6, 24, 24), 0.2)
    jt["sigma_b"] = sbf if os.path.exists(sbf) else "constant 0.2"
    jt["n"] = 0
    for k, (k2, v, lo, hi) in enumerate((("TER", "EmisCO2_NetTerrExch", 0, 4), ("ocean", "EmisCO2_Ocean", -4, 4))):
        T = np.asarray(t[v][:], float); Pf = np.asarray(p[v][:], float)
        m = np.abs(Pf) > 0
        x = np.ones_like(Pf); x[m] = T[m] / Pf[m]
        xc = np.clip(x, lo, hi)
        jt[f"J_b_{k2}"] = 0.5 * (((xc - 1) / sb[k]) ** 2)[m].sum()
        jt[f"clip_{k2}"] = int((x[m] != xc[m]).sum()); jt[f"SF50_{k2}"] = np.median(x[m])
        jt["n"] += int(m.sum())
    jt["J_b"] = jt["J_b_TER"] + jt["J_b_ocean"]
    jt["J"] = jt["J_obs"] + jt["J_b"]
    return jt


JT = {}
rep.append("\nCost function at the truth (reference, not a target; J_b(truth) is a lower bound: SF_true clipped to the bounds)")
for e0, osse in (("OSSE1_C24", "OSSE 1"), ("OSSE2_C24", "OSSE 2")):
    jt = JT[osse] = j_truth(e0)
    rep.append(f"  {osse}: J_obs(truth) {jt['J_obs']:.4e} (N {jt['N']}, N/2 {jt['N'] / 2:.4e}); J_b(truth) TER "
               f"{jt['J_b_TER']:.4e} (median SF_true {jt['SF50_TER']:.3f}, clipped {jt['clip_TER']}), ocean "
               f"{jt['J_b_ocean']:.4e} (median SF_true {jt['SF50_ocean']:.3f}, clipped {jt['clip_ocean']}), total "
               f"{jt['J_b']:.4e} (n {jt['n']} cell-months with flux, n/2 {jt['n'] / 2:.4e}; sigma_b {jt['sigma_b']}); "
               f"J(truth) {jt['J']:.4e}")

# ---- (i) cost function: J (with prior -> final drop), J_obs and J_b, 2 J_obs / N ----
fig, axs = plt.subplots(3, 2, figsize=(11, 10.5))
for c, osse in enumerate(("OSSE 1", "OSSE 2")):
    notes = []
    for e, r in runs.items():
        if r["osse"] != osse:
            continue
        hh = r["hist"]
        it = np.array([int(h["iteration"]) for h in hh]); J = np.array([float(h["J"]) for h in hh])
        Jo = np.array([float(h["J_obs"]) for h in hh]); Jb = np.array([float(h["J_b"]) for h in hh])
        chi = np.array([float(h["chi2"]) for h in hh]); acc = np.array([h["action"] in ("init", "accept") for h in hh])
        col = COL[r["adj"]]
        for ax, y in ((axs[0, c], J), (axs[2, c], chi)):
            ax.plot(it[acc], y[acc], "-o", color=col, ms=4, lw=1.8, label=r["lab"])
            ax.plot(it[~acc], y[~acc], "x", color=col, ms=6, mew=1.5)
        axs[1, c].plot(it[acc], Jo[acc], "-o", color=col, ms=4, lw=1.8, label=f"J_obs, {r['adj']}")
        pb = acc & (Jb > 0)                     # J_b = 0 at the prior (sigma = 1): not drawn on the log axis
        axs[1, c].plot(it[pb], Jb[pb], "--s", color=col, ms=3.5, lw=1.5, label=f"J_b, {r['adj']}")
        J0, Jf = J[acc][0], J[acc][-1]
        notes.append(f"{r['lab']}: J {J0:.3g} -> {Jf:.3g} (factor {J0 / Jf:.0f}, -{100 * (1 - Jf / J0):.1f} %); "
                     f"final J_obs {Jo[acc][-1]:.3g}, J_b {Jb[acc][-1]:.3g} ({100 * Jb[acc][-1] / Jf:.2g} % of J)")
        rep.append(f"  cost {e}: " + notes[-1])
    # the cost function at the truth (reference lines, drawn after the runs so they come last in the legends)
    jt = JT[osse]
    axs[0, c].axhline(jt["J"], color=INK2, lw=1, ls=":", label=f"J(truth) = {jt['J']:.3g} (reference, not a target)")
    axs[1, c].axhline(jt["J_obs"], color=INK2, lw=0.9, ls="-", label=f"J_obs(truth) = {jt['J_obs']:.3g} (noise)")
    axs[1, c].axhline(jt["J_b"], color=INK2, lw=0.9, ls="--", label=f"J_b(truth) = {jt['J_b']:.3g}")
    axs[0, c].text(0.98, 0.62, "\n".join(n.split(":")[0] + ":\n  " + n.split(": ", 1)[1].split(";")[0] for n in notes),
                   transform=axs[0, c].transAxes, ha="right", va="top", fontsize=7.5, color=INK)
    axs[1, c].text(0.98, 0.03, "final:\n" + "\n".join(n.split(":")[0] + (":\n  " if "not final" in n else ": ")
                                                     + n.split("final ", 1)[1] for n in notes),
                   transform=axs[1, c].transAxes, ha="right", va="bottom", fontsize=7.5, color=INK)
    for ax in axs[:, c]:
        ax.set_yscale("log")
    axs[2, c].axhline(1.0, color=INK2, lw=1, ls=":")
    axs[2, c].text(1, 1.05, "noise level (2 J_obs / N = 1)", fontsize=7, color=INK2, va="bottom")
    axs[0, c].set_title(f"{osse}: cost function J = J_obs + J_b", fontsize=10, color=INK, loc="left")
    axs[1, c].set_title(f"{osse}: observation (J_obs) and background (J_b) terms", fontsize=10, color=INK, loc="left")
    axs[2, c].set_title(f"{osse}: 2 J_obs / N_obs (N = 350 112 soundings)", fontsize=10, color=INK, loc="left")
    for ax in axs[:, c]:
        ax.set_xlabel("iteration", fontsize=8, color=INK2); style(ax, "both")
    axs[0, c].legend(fontsize=8, frameon=False, loc="upper right")
    lo, hi = axs[1, c].get_ylim(); axs[1, c].set_ylim(lo / 8, hi * 30)  # headroom: legend above, "final" text below the curves
    axs[1, c].legend(fontsize=7, frameon=False, loc="upper right", ncol=3)
axs[0, 0].plot([], [], "x", color=INK2, label="rejected trial point (step halved)")
axs[0, 0].legend(fontsize=8, frameon=False, loc="upper right")
fig.tight_layout(); fig.savefig(f"{out}/fig_cost.png", dpi=150); plt.close(fig)

# ---- (ii) flux error per iteration ----
fig, axs = plt.subplots(2, 4, figsize=(15, 6.5))
for c, osse in enumerate(("OSSE 1", "OSSE 2")):
    for k_i, k in enumerate(("TER", "ocean")):
        ax_e, ax_r = axs[0, 2 * c + k_i], axs[1, 2 * c + k_i]
        for e, r in runs.items():
            if r["osse"] != osse:
                continue
            ax_e.plot(r["it_acc"], r["err"][k], "-o", color=COL[r["adj"]], ms=3.5, lw=1.8, label=r["lab"])
            ax_r.plot(r["it_acc"], 100 * (np.array(r["rms"][k]) / r["rms"][k][0] - 1), "-o", color=COL[r["adj"]], ms=3.5, lw=1.8)
        ax_e.axhline(0, color=INK, lw=0.8)
        ax_r.axhline(0, color=INK, lw=0.8)
        ax_e.set_title(f"{osse} {k}: 2016 total minus truth [PgC]", fontsize=9, color=INK, loc="left")
        ax_r.set_title(f"{osse} {k}: cell RMSE change vs prior [%]", fontsize=9, color=INK, loc="left")
        for ax in (ax_e, ax_r):
            style(ax, "both"); ax.set_xlabel("iteration (accepted steps)", fontsize=8, color=INK2)
    axs[0, 2 * c].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{out}/fig_flux_error.png", dpi=150); plt.close(fig)

# ---- cell polygons (dateline-safe) for maps ----
polys, idx = [], []
for f in range(6):
    for y in range(24):
        for x in range(24):
            lo = np.array([clon[f, y, x], clon[f, y, x + 1], clon[f, y + 1, x + 1], clon[f, y + 1, x]])
            la = np.array([clat[f, y, x], clat[f, y, x + 1], clat[f, y + 1, x + 1], clat[f, y + 1, x]])
            lo = (lo + 180) % 360 - 180
            if lo.max() - lo.min() > 180:
                lo = np.where(lo < 0, lo + 360, lo)
                polys.append(np.c_[lo - 360, la]); idx.append((f, y, x))
            polys.append(np.c_[lo, la]); idx.append((f, y, x))
idx = np.array(idx)


def coast(ax):
    if os.path.exists(COAST):
        for s in shapefile.Reader(COAST).shapes():
            p = np.array(s.points); ax.plot(p[:, 0], p[:, 1], color=INK2, lw=0.4)


def cs_map(ax, v, cmap, vmin, vmax, title):
    vals = np.array([v[f, y, x] for f, y, x in idx])
    pc = PolyCollection(polys, array=np.ma.masked_invalid(vals), cmap=cmap, edgecolors="face", linewidths=0.2)
    pc.set_clim(vmin, vmax); pc.cmap.set_bad("#efeeeb")
    ax.add_collection(pc); coast(ax)
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_aspect("equal")
    ax.set_xticks(range(-180, 181, 90)); ax.set_yticks(range(-90, 91, 45))
    ax.set_title(title, fontsize=9, color=INK, loc="left"); style(ax, None)
    return pc


# ---- TransCom map ----
fig, axs = plt.subplots(1, 2, figsize=(13, 3.9))
cm_l = plt.get_cmap("tab20", 11); cm_o = plt.get_cmap("tab20b", 11)
for ax, k, cm, rng in ((axs[0], "TER", cm_l, (1, 11)), (axs[1], "ocean", cm_o, (12, 22))):
    v = np.where(rid[k] > 0, rid[k], np.nan).astype(float)
    cs_map(ax, v, cm, rng[0] - 0.5, rng[1] + 0.5, f"TransCom {'land' if k == 'TER' else 'ocean'} regions {rng[0]}-{rng[1]} (on the C24 grid)")
    for r in range(rng[0], rng[1] + 1):
        s = rid[k] == r
        if s.sum() == 0:
            continue
        # label at the cell of the region nearest its mean position (unit-vector mean)
        u = unit(lons[s], lats[s]).mean(0); u /= np.linalg.norm(u)
        k0 = np.argmax(unit(lons[s], lats[s]) @ u)
        ax.text(lons[s][k0] if lons[s][k0] <= 180 else lons[s][k0] - 360, lats[s][k0], str(r), fontsize=8, ha="center",
                va="center", color=INK, fontweight="bold")
fig.text(0.01, 0.01, "Land: " + "; ".join(f"{r} {TCNAME[r]}" for r in range(1, 12)) + "\nOcean: "
         + "; ".join(f"{r} {TCNAME[r]}" for r in range(12, 23)), fontsize=7, color=INK2, va="bottom")
fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig(f"{out}/fig_transcom_map.png", dpi=150); plt.close(fig)


# ---- (iii) regional comparison ----
def regional(flux, k, ids):
    return np.array([flux[k][:, rid[k] == r].sum() for r in ids])


def band(flux, k):
    return np.array([flux[k][:, (lats >= a) & (lats < b)].sum() for a, b in BANDS])


rep.append("\nRegional 2016 totals minus truth [PgC] (prior | posterior per run)")
for name in ("regional", "latbands"):
    fig, axs = plt.subplots(2, 2, figsize=(14, 7.5), sharey=False)
    for c, osse in enumerate(("OSSE 1", "OSSE 2")):
        rr = [r for r in runs.values() if r["osse"] == osse]
        for row, k in enumerate(("TER", "ocean")):
            ax = axs[row, c]
            if name == "regional":
                ids = list(range(1, 12)) if k == "TER" else list(range(12, 23))
                T, Pr = regional(truth, k, ids), regional(rr[0]["prior"], k, ids)
                labs = [f"{r} {TCNAME[r]}" for r in ids]
                post = [regional(r["post"], k, ids) for r in rr]
            else:
                T, Pr = band(truth, k), band(rr[0]["prior"], k)
                labs = [f"{a}..{b}" for a, b in BANDS]
                post = [band(r["post"], k) for r in rr]
            xg = np.arange(len(labs)); w = 0.27
            ax.bar(xg - w, Pr - T, w, color=PRIOR, label="prior")
            for j, (r, p) in enumerate(zip(rr, post)):
                ax.bar(xg + j * w, p - T, w, color=COL[r["adj"]], label=f"posterior, {r['lab']}")
            ax.axhline(0, color=INK, lw=0.8)
            ax.set_xticks(xg); ax.set_xticklabels(labs, rotation=40 if name == "regional" else 0, ha="right" if name == "regional" else "center", fontsize=7)
            what = "TransCom land regions" if (name == "regional" and k == "TER") else ("TransCom ocean regions" if name == "regional" else "latitude bands")
            ax.set_title(f"{osse} {k}, {what}: 2016 total minus truth [PgC]", fontsize=9, color=INK, loc="left")
            style(ax)
            if row == 0 and c == 0:
                ax.legend(fontsize=7, frameon=False)
            if row == 0 and c == 1:
                ax.legend(fontsize=7, frameon=False)
            # numbers: summed |regional error|, prior vs posterior
            line = f"  {osse} {k} {name}: sum|err| prior {np.abs(Pr - T).sum():.2f}; " + "; ".join(
                f"{r['adj']} {np.abs(p - T).sum():.2f}" for r, p in zip(rr, post))
            rep.append(line)
    fig.tight_layout(); fig.savefig(f"{out}/fig_{name}.png", dpi=150); plt.close(fig)

# ---- monthly global totals ----
fig, axs = plt.subplots(2, 2, figsize=(13, 6.5))
xm = np.arange(NM)
for c, osse in enumerate(("OSSE 1", "OSSE 2")):
    rr = [r for r in runs.values() if r["osse"] == osse]
    for row, k in enumerate(("TER", "ocean")):
        ax = axs[row, c]
        ax.plot(xm, truth[k].sum(axis=(1, 2, 3)), "-", color=INK, lw=2.2, label="truth")
        ax.plot(xm, rr[0]["prior"][k].sum(axis=(1, 2, 3)), "--", color=PRIOR, lw=1.8, label="prior")
        for r in rr:
            ax.plot(xm, r["post"][k].sum(axis=(1, 2, 3)), "-o", color=COL[r["adj"]], ms=3, lw=1.5, label=f"posterior, {r['lab']}")
        ax.set_xticks(xm[::2]); ax.set_xticklabels(labels_m[::2], fontsize=7)
        ax.set_title(f"{osse} {k}: global monthly total [PgC per month]", fontsize=9, color=INK, loc="left")
        style(ax)
        if row == 0:
            ax.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{out}/fig_monthly_global.png", dpi=150); plt.close(fig)

# ---- monthly bias and RMSE of the posterior (best) and the prior vs truth ----
per_m2 = 1e15 / area[None] / mdays[:, None, None, None]        # PgC per cell per month -> gC m-2 day-1
rep.append("\nPer month of 2016: global bias [PgC per month] and area-weighted cell RMSE [gC m-2 day-1], prior | posterior")
fig, axs = plt.subplots(4, 2, figsize=(13, 11), sharex=True)
for c, osse in enumerate(("OSSE 1", "OSSE 2")):
    rr = [r for r in runs.values() if r["osse"] == osse]
    for k_i, k in enumerate(("TER", "ocean")):
        ax_b, ax_r = axs[2 * k_i, c], axs[2 * k_i + 1, c]
        m = (np.abs(truth[k]).sum(0) + np.abs(rr[0]["prior"][k]).sum(0)) > 0      # cells with this flux
        w = area[m] / area[m].sum()
        bias = lambda f: (f[k] - truth[k]).sum(axis=(1, 2, 3))
        rmse = lambda f: np.sqrt(((((f[k] - truth[k]) * per_m2)[:, m]) ** 2 * w).sum(1))
        ax_b.plot(xm, bias(rr[0]["prior"]), "--o", color=PRIOR, ms=3, lw=1.8, label="prior - truth")
        ax_r.plot(xm, rmse(rr[0]["prior"]), "--o", color=PRIOR, ms=3, lw=1.8, label="prior - truth")
        rep.append(f"  {osse} {k} prior: bias " + " ".join(f"{v:+.3f}" for v in bias(rr[0]["prior"]))
                   + " | RMSE " + " ".join(f"{v:.3f}" for v in rmse(rr[0]["prior"])))
        for r in rr:
            ax_b.plot(xm, bias(r["post"]), "-o", color=COL[r["adj"]], ms=3, lw=1.8, label=f"posterior - truth, {r['lab']}")
            ax_r.plot(xm, rmse(r["post"]), "-o", color=COL[r["adj"]], ms=3, lw=1.8, label=f"posterior - truth, {r['lab']}")
            rep.append(f"  {osse} {k} {r['adj']:10s}: bias " + " ".join(f"{v:+.3f}" for v in bias(r["post"]))
                       + " | RMSE " + " ".join(f"{v:.3f}" for v in rmse(r["post"])))
        ax_b.axhline(0, color=INK, lw=0.8)
        ax_r.set_ylim(bottom=0)
        ax_b.set_title(f"{osse} {k}: monthly bias (global total minus truth) [PgC per month]", fontsize=9, color=INK, loc="left")
        ax_r.set_title(f"{osse} {k}: monthly cell RMSE (area-weighted) [gC m-2 day-1]", fontsize=9, color=INK, loc="left")
        for ax in (ax_b, ax_r):
            style(ax)
    axs[0, c].legend(fontsize=7, frameon=False)
axs[-1, 0].set_xticks(xm); axs[-1, 0].set_xticklabels([l[5:] for l in labels_m], fontsize=7)
for ax in axs[-1]:
    ax.set_xlabel("month of 2016", fontsize=8, color=INK2)
fig.tight_layout(); fig.savefig(f"{out}/fig_monthly_bias_rmse.png", dpi=150); plt.close(fig)

# ---- maps: 2016 mean TER error ----
mean_flux = lambda f: (f * per_m2).mean(0)
fig, axs = plt.subplots(2, 3, figsize=(15, 7.2))
for row, osse in enumerate(("OSSE 1", "OSSE 2")):
    rr = [r for r in runs.values() if r["osse"] == osse]
    dprior = mean_flux(rr[0]["prior"]["TER"]) - mean_flux(truth["TER"])
    vmax = np.nanpercentile(np.abs(dprior[np.abs(dprior) > 0]), 98)
    m = np.abs(truth["TER"]).sum(0) + np.abs(rr[0]["prior"]["TER"]).sum(0) == 0
    pc = cs_map(axs[row, 0], np.where(m, np.nan, dprior), "RdBu_r", -vmax, vmax, f"{osse} TER: prior - truth")
    for j, r in enumerate(rr):
        d = mean_flux(r["post"]["TER"]) - mean_flux(truth["TER"])
        cs_map(axs[row, j + 1], np.where(m, np.nan, d), "RdBu_r", -vmax, vmax, f"{osse} TER: posterior - truth\n{r['lab']}")
    cb = fig.colorbar(pc, ax=axs[row, :].tolist(), orientation="horizontal", shrink=0.4, pad=0.08, aspect=40)
    cb.set_label(f"{osse}: 2016 mean TER difference [gC m-2 day-1]", fontsize=8, color=INK2); cb.ax.tick_params(labelsize=7)
fig.savefig(f"{out}/fig_maps_ter.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ---- numbers ----
rep.append("\nPer run: 2016 global total minus truth [PgC] and cell RMSE change vs prior, at the posterior")
for e, r in runs.items():
    rep.append(f"  {e:13s} {r['lab']:38s} posterior it {r['it_post']:2d}: " + "; ".join(
        f"{k} {r['err'][k][0]:+.2f} -> {r['err'][k][-1]:+.2f} PgC, RMSE {100 * (r['rms'][k][-1] / r['rms'][k][0] - 1):+.1f} %" for k in ("TER", "ocean")))
open(f"{out}/results.txt", "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
