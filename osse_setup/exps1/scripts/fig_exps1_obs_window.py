#!/usr/bin/env python3
"""exps1 (2026-11-16) observation figures for the whole assimilation window 2015-10-01 .. 2017-05-01 (19 months):
  fig_obs_per_month_window.png  soundings per month, land nadir (LN), land glint (LG), ocean glint (OG); the analysis year
                                2016 is shaded
  fig_obs_density_window.png    soundings per 10^4 km2 over the window on 1x1 degree cells (log colour scale), and the
                                zonal-mean density per sounding type
  obs_counts_window.txt         the numbers behind the figures
Reads the same daily files as make_track_lnlgog.py (pseudo set from 2016-01-01, the real B10 set before; identical
positions) and checks the total against the track file the OSSEs use.
    fig_exps1_obs_window.py <out dir>
"""
import calendar, os, sys
from datetime import date, timedelta
import numpy as np, netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import shapefile

sys.path.insert(0, "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1/scripts")
from make_track_lnlgog import day_file, SUB as SUBDIR

TRACK = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1/obs_c24_1510_1705/track_file_lnlgog_201510_201705.nc"
COAST = os.path.expanduser("~/.local/share/cartopy/shapefiles/natural_earth/physical/ne_110m_coastline.shp")
TYPES = (("LN", "land nadir", "#2a78d6"), ("LG", "land glint", "#1baf7a"), ("OG", "ocean glint", "#eb6834"))
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
D0, D1 = date(2015, 10, 1), date(2017, 5, 1)                 # D1 exclusive
out = sys.argv[1]
os.makedirs(out, exist_ok=True)


def style(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, axis=grid_axis, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


months = []
d = D0
while d < D1:
    months.append(d); d = date(d.year + d.month // 12, d.month % 12 + 1, 1)
NM = len(months)
per_mon = {k: np.zeros(NM, int) for k, _, _ in TYPES}
grid = {k: np.zeros((180, 360), int) for k, _, _ in TYPES}
d = D0
while d < D1:
    im = (d.year - D0.year) * 12 + d.month - D0.month
    for k, _, _ in TYPES:
        f, _ = day_file(SUBDIR[k], d)
        if f is None:
            continue
        g = nc.Dataset(f)
        lat, lon = np.asarray(g["latitude"][:], float), (np.asarray(g["longitude"][:], float) + 180) % 360 - 180
        g.close()
        per_mon[k][im] += lat.size
        np.add.at(grid[k], (np.clip(np.floor(lat + 90).astype(int), 0, 179), np.clip(np.floor(lon + 180).astype(int), 0, 359)), 1)
    d += timedelta(days=1)
tot = sum(v.sum() for v in per_mon.values())
n_track = nc.Dataset(TRACK).dimensions["time"].size
assert tot == n_track, (tot, n_track)

# 1x1 cell areas [10^4 km2]
R = 6371.0
lat_e = np.radians(np.arange(-90, 91))
cell_area = (R ** 2 * np.radians(1.0) * (np.sin(lat_e[1:]) - np.sin(lat_e[:-1])))[:, None] * np.ones((1, 360)) / 1e4
allg = sum(grid.values())
dens = allg / cell_area

# ---- per month ----
fig, ax = plt.subplots(figsize=(11, 3.9))
xm = np.arange(NM)
i16 = [i for i, m in enumerate(months) if m.year == 2016]
ax.axvspan(i16[0] - 0.5, i16[-1] + 0.5, color="#f3f2ee", zorder=0)
bottom = np.zeros(NM)
for k, name, col in TYPES:
    ax.bar(xm, per_mon[k], bottom=bottom, width=0.8, color=col, label=f"{name} ({k}), {per_mon[k].sum():,}", lw=0,
           edgecolor="white", zorder=2)
    bottom += per_mon[k]
ax.set_xticks(xm); ax.set_xticklabels([f"{calendar.month_abbr[m.month]}\n{m.year}" if (m.month == 1 or i == 0) else calendar.month_abbr[m.month]
                                       for i, m in enumerate(months)], fontsize=7)
ax.set_xlim(-0.6, NM - 0.4)
ax.set_ylabel("soundings per month", color=INK2, fontsize=9)
ax.set_title(f"Pseudo OCO-2 B10 soundings per month, 2015-10 .. 2017-04 ({tot:,} in total, all assimilated; analysis year 2016 shaded)",
             fontsize=10, color=INK, loc="left")
style(ax)
ax.set_ylim(0, bottom.max() * 1.18)
ax.legend(fontsize=8, frameon=False, loc="upper left", ncol=3)
fig.tight_layout(); fig.savefig(f"{out}/fig_obs_per_month_window.png", dpi=150); plt.close(fig)

# ---- density map + zonal mean per type ----
fig, (ax, az) = plt.subplots(1, 2, figsize=(13.5, 5), gridspec_kw=dict(width_ratios=[4.2, 1]))
m = np.ma.masked_equal(dens, 0)
pc = ax.pcolormesh(np.arange(-180, 181), np.arange(-90, 91), m, cmap="viridis", norm=LogNorm(max(m.min(), 0.5), m.max()), shading="flat")
if os.path.exists(COAST):
    for s in shapefile.Reader(COAST).shapes():
        p = np.array(s.points); ax.plot(p[:, 0], p[:, 1], color="#ffffff", lw=0.4)
        ax.plot(p[:, 0], p[:, 1], color=INK2, lw=0.25)
ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_aspect("equal")
ax.set_xticks(range(-180, 181, 60)); ax.set_yticks(range(-90, 91, 30))
ax.set_title(f"Soundings per 10⁴ km², 2015-10 .. 2017-04, LN + LG + OG, 1° x 1° cells ({(allg > 0).sum():,} of 64,800 observed)",
             fontsize=10, color=INK, loc="left")
style(ax, grid_axis="both")
cb = fig.colorbar(pc, ax=ax, shrink=0.8, pad=0.02); cb.set_label("soundings per 10⁴ km² (log scale)", fontsize=8, color=INK2)
cb.ax.tick_params(labelsize=7, colors=INK2)
latc = np.arange(-89.5, 90)
band_area = cell_area.sum(1)
for k, name, col in TYPES:
    az.plot(grid[k].sum(1) / band_area, latc, color=col, lw=1.8, label=k)
az.plot(allg.sum(1) / band_area, latc, color=INK, lw=1.2, ls="--", label="all")
az.set_ylim(-90, 90); az.set_yticks(range(-90, 91, 30))
az.set_xlabel("soundings per 10⁴ km²\n(zonal mean, all longitudes)", fontsize=8, color=INK2)
az.set_title("zonal mean", fontsize=10, color=INK, loc="left")
style(az, grid_axis="both"); az.legend(fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(f"{out}/fig_obs_density_window.png", dpi=150); plt.close(fig)

with open(f"{out}/obs_counts_window.txt", "w") as o:
    o.write(f"# exps1 soundings, 2015-10-01 .. 2017-05-01 (fig_exps1_obs_window.py); total {tot} = track file {n_track}\n")
    o.write("".join(f"{k} ({n}) {per_mon[k].sum()}\n" for k, n, _ in TYPES))
    o.write("\nmonth    " + " ".join(f"{k:>7s}" for k, _, _ in TYPES) + "   total\n")
    for i, mo in enumerate(months):
        o.write(f"{mo:%Y-%m}  " + " ".join(f"{per_mon[k][i]:7d}" for k, _, _ in TYPES) + f" {sum(per_mon[k][i] for k, _, _ in TYPES):7d}\n")
    s16 = sum(per_mon[k][i16].sum() for k, _, _ in TYPES)
    o.write(f"\n2016 (analysis year): {s16}; spin-up 2015-10..12: {sum(per_mon[k][:i16[0]].sum() for k, _, _ in TYPES)}; "
            f"2017-01..04: {sum(per_mon[k][i16[-1] + 1:].sum() for k, _, _ in TYPES)}\n")
    o.write(f"1x1 cells with soundings: {(allg > 0).sum()} of 64800; density soundings per 1e4 km2: median of observed cells "
            f"{np.median(dens[allg > 0]):.2f}, 99th percentile {np.percentile(dens[allg > 0], 99):.2f}, max {dens.max():.2f}\n")
    o.write("global mean density (all cells) %.3f per 1e4 km2\n" % (allg.sum() / cell_area.sum()))
print(open(f"{out}/obs_counts_window.txt").read())
