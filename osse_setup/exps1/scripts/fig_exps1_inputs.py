#!/usr/bin/env python3
"""exps1 (2026-11-09) input figures, C24, 2016:
  fig_obs_per_day.png    soundings per day, land nadir (LN), land glint (LG), ocean glint (OG), stacked
  fig_obs_per_month.png  soundings per month, the same split
  fig_obs_map_1x1.png    soundings per 1x1 degree cell over the year (all types), log colour scale
  fig_sigma_b_osse2.png  OSSE 2 background error sigma_b (1-sigma of the scale factor): TER (CLASS-CTEM Uncertainty) and
                         ocean (MOM6 Uncertainty), January and July 2016, on the C24 grid
  obs_counts.txt         the numbers behind the observation figures
Reads J. Liu's daily pseudo B10 files directly (the same soundings as exps1/obs_c24/track_file_lnlgog_2016.nc) and
exps1/OSSE2_C24/sigma_b.npz (or --sigma-b <file>).
    fig_exps1_inputs.py --out-dir <dir> [--sigma-b <sigma_b.npz>]
"""
import argparse, calendar, os
from datetime import date
import numpy as np, netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import shapefile

L2 = "/nobackupp17/jliu7/OCO2/L2#/OCO2-B10/PSEUDO/TRENDY/ORCHIDEE-ECCO2/"
SUB = (("LN", "OCO2b10_29Aug2020_LN", "#2a78d6"), ("LG", "OCO2b10_29Aug2020_LG", "#1baf7a"),
       ("OG", "OCO2b10_29Aug2020_OG", "#eb6834"))
NAME = {"LN": "land nadir", "LG": "land glint", "OG": "ocean glint"}
COORD = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/tests_final/P2_advection/P2a/results/coord.nc4"
COAST = os.path.expanduser("~/.local/share/cartopy/shapefiles/natural_earth/physical/ne_110m_coastline.shp")
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

ap = argparse.ArgumentParser()
ap.add_argument("--out-dir", required=True)
ap.add_argument("--sigma-b", default="/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1/OSSE2_C24/sigma_b.npz")
ap.add_argument("--only-sigma-b", action="store_true", help="skip the observation figures")
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)


def style(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, axis=grid_axis, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def coast(ax):
    if os.path.exists(COAST):
        for s in shapefile.Reader(COAST).shapes():
            p = np.array(s.points)
            ax.plot(p[:, 0], p[:, 1], color=INK2, lw=0.4)


# ---- observations: per day, per month, per 1x1 cell ----
if not a.only_sigma_b:
    days = [date(2016, 1, 1).toordinal() + i for i in range(366)]
    per_day = {k: np.zeros(366, int) for k, _, _ in SUB}
    grid = np.zeros((180, 360), int)
    for i, o in enumerate(days):
        d = date.fromordinal(o)
        for k, sub, _ in SUB:
            f = f"{L2}/{sub}/{d:%Y/%m/%d}.nc"
            if not os.path.exists(f):
                continue
            g = nc.Dataset(f)
            lat, lon = np.asarray(g["latitude"][:], float), (np.asarray(g["longitude"][:], float) + 180) % 360 - 180
            per_day[k][i] = lat.size
            np.add.at(grid, (np.clip(np.floor(lat + 90).astype(int), 0, 179), np.clip(np.floor(lon + 180).astype(int), 0, 359)), 1)
    mon_idx = np.array([date.fromordinal(o).month - 1 for o in days])
    per_mon = {k: np.bincount(mon_idx, weights=v, minlength=12).astype(int) for k, v in per_day.items()}
    tot = sum(v.sum() for v in per_day.values())

    x = np.array([np.datetime64(date.fromordinal(o)) for o in days])
    fig, ax = plt.subplots(figsize=(10, 3.6))
    bottom = np.zeros(366)
    for k, _, col in SUB:
        ax.bar(x, per_day[k], bottom=bottom, width=1.0, color=col, label=f"{NAME[k]} ({k}), {per_day[k].sum():,}", lw=0)
        bottom += per_day[k]
    ax.set_ylabel("soundings per day", color=INK2, fontsize=9)
    ax.set_title(f"Pseudo OCO-2 B10 soundings per day, 2016 ({tot:,} in total)", fontsize=10, color=INK, loc="left")
    ax.set_xlim(x[0], x[-1])
    style(ax)
    ax.legend(fontsize=8, frameon=False, loc="upper left", ncol=3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out_dir, "fig_obs_per_day.png"), dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 3.6))
    xm = np.arange(12); w = 0.27
    for j, (k, _, col) in enumerate(SUB):
        ax.bar(xm + (j - 1) * w, per_mon[k], width=w - 0.02, color=col, label=f"{NAME[k]} ({k})")
    ax.set_xticks(xm); ax.set_xticklabels([calendar.month_abbr[m + 1] for m in xm])
    ax.set_ylabel("soundings per month", color=INK2, fontsize=9)
    ax.set_title("Pseudo OCO-2 B10 soundings per month, 2016", fontsize=10, color=INK, loc="left")
    style(ax)
    ax.legend(fontsize=8, frameon=False, loc="upper left", ncol=3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out_dir, "fig_obs_per_month.png"), dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    m = np.ma.masked_equal(grid, 0)
    pc = ax.pcolormesh(np.arange(-180, 181), np.arange(-90, 91), m, cmap="viridis", norm=LogNorm(1, max(2, grid.max())), shading="flat")
    coast(ax)
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_aspect("equal")
    ax.set_xticks(range(-180, 181, 60)); ax.set_yticks(range(-90, 91, 30))
    ax.set_title(f"Soundings per 1° x 1° cell, 2016, LN + LG + OG ({(grid > 0).sum():,} of 64,800 cells observed)",
                 fontsize=10, color=INK, loc="left")
    style(ax, grid_axis="both")
    cb = fig.colorbar(pc, ax=ax, shrink=0.8, pad=0.02); cb.set_label("soundings per cell (log scale)", fontsize=8, color=INK2)
    cb.ax.tick_params(labelsize=7, colors=INK2)
    fig.tight_layout(); fig.savefig(os.path.join(a.out_dir, "fig_obs_map_1x1.png"), dpi=150); plt.close(fig)

    with open(os.path.join(a.out_dir, "obs_counts.txt"), "w") as o:
        o.write("# exps1 pseudo OCO-2 B10 (ORCHIDEE-ECCO2) soundings, 2016 (fig_exps1_inputs.py)\n")
        o.write(f"total {tot}\n" + "".join(f"{k} ({NAME[k]}) {v.sum()}\n" for k, v in per_day.items()))
        o.write("\nmonth " + " ".join(f"{k:>7s}" for k, _, _ in SUB) + "   total\n")
        for i in range(12):
            o.write(f"{calendar.month_abbr[i + 1]:5s} " + " ".join(f"{per_mon[k][i]:7d}" for k, _, _ in SUB)
                    + f" {sum(per_mon[k][i] for k, _, _ in SUB):7d}\n")
        s = sum(per_day.values())
        o.write(f"\nper day: min {s.min()} median {int(np.median(s))} max {s.max()}; days with no soundings: {int((s == 0).sum())}\n")
        o.write(f"1x1 cells with soundings: {(grid > 0).sum()} of 64800; max per cell {grid.max()}; median of observed cells {int(np.median(grid[grid > 0]))}\n")
    print(open(os.path.join(a.out_dir, "obs_counts.txt")).read())

# ---- OSSE 2 sigma_b on the C24 grid ----
if os.path.exists(a.sigma_b):
    import sys
    from matplotlib.collections import PolyCollection
    sys.path.insert(0, "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_realcase/scripts")
    from make_inputs import subcell_points
    I = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/inputs/surface_fluxes_osse"
    sb = np.load(a.sigma_b)["sigma_b"].reshape(2, 12, 6, 24, 24)
    c = nc.Dataset(COORD)
    clon, clat = np.asarray(c["corner_lons"][:], float), np.asarray(c["corner_lats"][:], float)
    # the same 4x5 -> C24 sampling as make_exps1_inputs.py, to mark the cells where the flux is 0 (no gradient there)
    slon, slat = subcell_points(c["corner_lons"][:], c["corner_lats"][:])
    LAT_EDGES = np.concatenate([[-90.0], np.arange(-88.0, 89.0, 4.0), [90.0]])
    jj = np.clip(np.searchsorted(LAT_EDGES, slat, side="right") - 1, 0, 45)
    ii = np.floor(((slon + 182.5) % 360.0) / 5.0).astype(int) % 72
    has_flux = lambda f: (np.abs(f[jj, ii]).max(axis=(-1, -2)) > 0).reshape(6, 24, 24)
    # cell polygons; a cell crossing the dateline is drawn twice (shifted by +-360) so no gaps appear
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
    fig, axs = plt.subplots(2, 2, figsize=(12, 6.6))
    panels = ((0, "TER", "TER/CLASS-CTEM", 3.0), (1, "ocean", "Ocean/MOM6-MON", 12.5))
    for r, (fld, lab, src, vmax) in enumerate(panels):
        for cc, (mi, mn) in enumerate(((0, "Jan"), (6, "Jul"))):
            ax = axs[r, cc]
            flux = np.squeeze(nc.Dataset(f"{I}/{src}/2016/{mi + 1:02d}.nc")["CO2_Flux"][:]).astype(float)
            mask = has_flux(np.nan_to_num(flux))
            v = sb[fld, mi]
            val = np.array([v[f, y, x] if mask[f, y, x] else np.nan for f, y, x in idx])
            pcol = PolyCollection(polys, array=np.ma.masked_invalid(val), cmap="viridis", edgecolors="face", linewidths=0.2)
            pcol.set_clim(0, vmax); pcol.cmap.set_bad("#d9d8d4")
            ax.add_collection(pcol)
            coast(ax)
            ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_aspect("equal")
            ax.set_xticks(range(-180, 181, 90)); ax.set_yticks(range(-90, 91, 45))
            vm = v[mask]
            ax.set_title(f"{lab}, {mn} 2016: median {np.median(vm):.2f}, max {vm.max():.2f}", fontsize=9, color=INK, loc="left")
            style(ax, grid_axis="both")
        cb = fig.colorbar(pcol, ax=axs[r, :].tolist(), shrink=0.85, pad=0.02)
        cb.set_label(f"σ_b {lab}", fontsize=8, color=INK2); cb.ax.tick_params(labelsize=7, colors=INK2)
    fig.suptitle("OSSE 2 background error σ_b (1σ of the scale factor), C24 cells; TER from the CLASS-CTEM Uncertainty, "
                 "ocean from the MOM6 Uncertainty. Grey: no flux of that type (no gradient, σ stays 1)",
                 fontsize=9, color=INK, x=0.01, ha="left")
    fig.savefig(os.path.join(a.out_dir, "fig_sigma_b_osse2.png"), dpi=150); plt.close(fig)
    print("sigma_b figure written")
else:
    print(f"no {a.sigma_b}: sigma_b figure skipped")
