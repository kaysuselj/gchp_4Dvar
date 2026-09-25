#!/usr/bin/env python3
"""
Two new diagnostic plots for osse-exp_constant_offset_v1:

  1. sigma_diff.png
       Difference between the final optimised scaling factor (sigma_final)
       and the theoretical target (sigma_theoretical = 1/1.3 over land,
       0 over ocean), plotted on the native C24 cubed-sphere grid.

  2. flux_diff.png
       Difference between the final optimised TER flux
           flux_opt  = TER * 1.3 * sigma_final
       and the theoretical TER flux (the "truth" the OSSE was built from)
           flux_theo = TER * 1.3 * (1/1.3) * land_mask  =  TER * land_mask
       where TER is the CLASS-CTEM monthly mean averaged over the
       assimilation window (Jan–Mar 2016).

Data flow
---------
  sigma_final   : 4dvar_state_osse.npz  → x_prev (3456,) = (6×24×24)
  CS coords     : adjoint/OutputDir/GEOSChem.Adjoint.20160401.nc4  → lats/lons (6,24,24)
  Land mask     : Cartopy 110-m land polygons (shapely point-in-polygon for each CS cell)
  TER fluxes    : surface_fluxes_osse/TER/CLASS-CTEM/2016/{01,02,03}.nc (kg m-2 s-1)
                  averaged over the 3-month window then interpolated to C24 centres
                  using scipy.interpolate.RegularGridInterpolator (bilinear)
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
import cartopy.feature as cf
import shapely.geometry as sgeom

# ── Paths ──────────────────────────────────────────────────────────────────
HERE       = os.path.dirname(os.path.abspath(__file__))
STATE_FILE = os.path.join(HERE, '4dvar_state_osse.npz')
COORD_FILE = os.path.join(HERE, 'adjoint', 'OutputDir',
                          'GEOSChem.Adjoint.20160401.nc4')
TER_DIR    = ('/nobackup/ksuselj1/gchp_14.5.3_adjoint_surfaceF/'
              'surface_fluxes_osse/TER/CLASS-CTEM')
PLOT_DIR   = os.path.join(HERE, 'plots_osse', '20260923_092236')

# Assimilation window: 2016-01 through 2016-03 (inclusive end month)
TER_MONTHS = [('2016', '01'), ('2016', '02'), ('2016', '03')]

SCALE_FACTOR = 1.3      # TER_OSSE_INFLATE in HEMCO (scale factor 752)
SIGMA_THEO   = 1.0 / SCALE_FACTOR   # = 1/1.3 over land: makes SCALE*sigma_theo=1

# TER file is in kgC/km²/s (confirmed by global total ~83 PgC/yr).
# HEMCO applies: TER × 3.667e-6 (CO2_UNIT_CONV=44/12×1e-6) × 1.3 × sigma
# to get kg CO2/m²/s.  For the carbon flux in gC/m²/yr we convert directly
# from kgC/km²/s:  × (1/1e6 km²/m²) × (1000 g/kg) × (3.156e7 s/yr)
KGC_KM2_S_TO_GC_M2_YR = 1e-6 * 1e3 * 3.156e7   # ≈ 31,560

NF, IM = 6, 24   # C24 cubed-sphere

os.makedirs(PLOT_DIR, exist_ok=True)

# ── Utility: shared Natural Earth coastlines ─────────────────────────────
_COAST_LINES = None

def _load_coast():
    global _COAST_LINES
    if _COAST_LINES is not None:
        return
    import shapefile as _shp
    candidates = [
        os.path.expanduser('~/.local/share/cartopy/shapefiles/'
                           'natural_earth/physical/ne_110m_coastline.shp'),
        os.path.expanduser('~/.local/share/cartopy/shapefiles/'
                           'natural_earth/physical/ne_50m_coastline.shp'),
    ]
    for path in candidates:
        if os.path.exists(path):
            try:
                sf = _shp.Reader(path)
                _COAST_LINES = [np.array(s.points) for s in sf.shapes()]
                return
            except Exception:
                pass

def _add_coast(ax):
    _load_coast()
    if _COAST_LINES:
        for pts in _COAST_LINES:
            ax.plot(pts[:, 0], pts[:, 1], color='0.3', linewidth=0.5, zorder=3)

def _set_ticks(ax):
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-90, 91, 30))
    ax.set_xticklabels(
        [f'{x}°{"W" if x<0 else ("E" if x>0 else "")}' for x in range(-180,181,60)],
        fontsize=8)
    ax.set_yticklabels(
        [f'{y}°{"S" if y<0 else ("N" if y>0 else "")}' for y in range(-90,91,30)],
        fontsize=8)


# ═══════════════════════════════════════════════════════════════════════════
# 1.  Load CS24 grid coordinates
# ═══════════════════════════════════════════════════════════════════════════
print('Loading CS24 coordinates…')
with xr.open_dataset(COORD_FILE) as ds:
    cs_lats = np.squeeze(ds['lats'].values).astype(float)   # (6,24,24)
    cs_lons = np.squeeze(ds['lons'].values).astype(float)   # (6,24,24) in 0–360
    cs_lons = ((cs_lons + 180.0) % 360.0) - 180.0          # wrap to -180..180

print(f'  CS lats: {cs_lats.min():.1f} … {cs_lats.max():.1f}')
print(f'  CS lons: {cs_lons.min():.1f} … {cs_lons.max():.1f}')


# ═══════════════════════════════════════════════════════════════════════════
# 2.  Build land mask on C24 grid (Cartopy 110-m land polygons, point-in-poly)
# ═══════════════════════════════════════════════════════════════════════════
print('Building C24 land mask (point-in-polygon, 110-m Natural Earth)…')

land_geoms = list(cf.NaturalEarthFeature(
    'physical', 'land', '110m').geometries())
from shapely.ops import unary_union
land_union  = unary_union(land_geoms)   # one MultiPolygon for fast lookup

flat_lats = cs_lats.ravel()
flat_lons = cs_lons.ravel()

land_mask = np.zeros(len(flat_lats), dtype=bool)
for k, (lat, lon) in enumerate(zip(flat_lats, flat_lons)):
    land_mask[k] = land_union.contains(sgeom.Point(lon, lat))

land_mask = land_mask.reshape(NF, IM, IM)
n_land = land_mask.sum()
print(f'  Land cells: {n_land} / {NF*IM*IM}  ({100*n_land/(NF*IM*IM):.1f} %)')


# ═══════════════════════════════════════════════════════════════════════════
# 3.  Load final sigma (C24)
# ═══════════════════════════════════════════════════════════════════════════
print('Loading final sigma from state file…')
sv         = np.load(STATE_FILE)
sigma_final = sv['x_prev'].reshape(NF, IM, IM)   # final completed step
print(f'  sigma_final: min={sigma_final.min():.4f}  max={sigma_final.max():.4f}  '
      f'mean={sigma_final.mean():.4f}')


# ═══════════════════════════════════════════════════════════════════════════
# 4.  Theoretical sigma (1/1.3 land, 0 ocean)
# ═══════════════════════════════════════════════════════════════════════════
sigma_theo       = np.zeros((NF, IM, IM))
sigma_theo[land_mask] = SIGMA_THEO
print(f'  sigma_theo : land={SIGMA_THEO:.4f} ({1/SCALE_FACTOR:.4f}), ocean=0')


# ═══════════════════════════════════════════════════════════════════════════
# 5.  Load TER fluxes and time-average over the assimilation window
# ═══════════════════════════════════════════════════════════════════════════
print(f'Loading TER flux for {len(TER_MONTHS)} months ({TER_MONTHS[0][1]}/{TER_MONTHS[0][0]}'
      f' – {TER_MONTHS[-1][1]}/{TER_MONTHS[-1][0]})…')

ter_accum = None
ter_lat   = None
ter_lon   = None

for (yr, mo) in TER_MONTHS:
    path = os.path.join(TER_DIR, yr, f'{mo}.nc')
    with xr.open_dataset(path) as ds:
        flux = ds['CO2_Flux'].values.squeeze().astype(float)   # (lat, lon)
        if ter_lat is None:
            ter_lat = ds['lat'].values.astype(float)
            ter_lon = ds['lon'].values.astype(float)
    ter_accum = flux if ter_accum is None else ter_accum + flux

ter_mean = ter_accum / len(TER_MONTHS)   # time-mean TER [kg m-2 s-1]
print(f'  TER lat-lon grid : {len(ter_lat)} lat × {len(ter_lon)} lon')
print(f'  TER mean: min={ter_mean.min():.3e}  max={ter_mean.max():.3e} kg m-2 s-1')


# ═══════════════════════════════════════════════════════════════════════════
# 6.  Interpolate TER onto C24 cell centres (bilinear)
# ═══════════════════════════════════════════════════════════════════════════
print('Interpolating TER onto C24 land cells (nearest-neighbour)…')
# Nearest-neighbour avoids averaging land and ocean source cells at coastlines,
# which bilinear interpolation would do (underestimating coastal TER values).
interp_fn = RegularGridInterpolator(
    (ter_lat, ter_lon),
    ter_mean,
    method='nearest',
    bounds_error=False,
    fill_value=None)

query_pts = np.column_stack([flat_lats, flat_lons])
ter_cs    = interp_fn(query_pts).reshape(NF, IM, IM)   # [kgC km-2 s-1]
print(f'  TER on C24: min={ter_cs.min():.3e}  max={ter_cs.max():.3e} kgC km-2 s-1')


# ═══════════════════════════════════════════════════════════════════════════
# 7.  Compute differences
# ═══════════════════════════════════════════════════════════════════════════

# --- Plot 1: sigma_final - sigma_theoretical ---
sigma_diff = sigma_final - sigma_theo   # positive = overestimate

# --- Plot 2: flux fields ---
# Prior flux:        TER * 1.3 * sigma_prior (sigma_prior = 1 everywhere)
# Optimised flux:    TER * 1.3 * sigma_final
# Theoretical flux:  TER * land_mask  (= TER * 1.3 * (1/1.3) * land_mask)
sigma_prior = np.ones((NF, IM, IM))
flux_prior = ter_cs * SCALE_FACTOR * sigma_prior        # [kgC km-2 s-1]
flux_opt  = ter_cs * SCALE_FACTOR * sigma_final         # [kgC km-2 s-1]
flux_theo = ter_cs * land_mask.astype(float)            # TER × 1.3 × (1/1.3) = TER on land
flux_diff_opt   = flux_opt   - flux_theo                # [kgC km-2 s-1]
flux_diff_prior = flux_prior - flux_theo                # [kgC km-2 s-1]

C = KGC_KM2_S_TO_GC_M2_YR
flux_prior_gC     = flux_prior     * C
flux_opt_gC       = flux_opt       * C
flux_theo_gC      = flux_theo      * C
flux_diff_opt_gC  = flux_diff_opt  * C
flux_diff_prior_gC = flux_diff_prior * C
print(f'  Prior  diff [gC m-2 yr-1]: min={flux_diff_prior_gC.min():.2f}  max={flux_diff_prior_gC.max():.2f}')
print(f'  Optim  diff [gC m-2 yr-1]: min={flux_diff_opt_gC.min():.2f}  max={flux_diff_opt_gC.max():.2f}')


# ═══════════════════════════════════════════════════════════════════════════
# 8.  Scatter-plot helper for CS fields
# ═══════════════════════════════════════════════════════════════════════════

def cs_scatter(ax, data, vmin, vmax, cmap='RdBu_r', land_only=False):
    mask = land_mask.ravel() if land_only else np.ones(len(land_mask.ravel()), dtype=bool)
    return ax.scatter(
        cs_lons.ravel()[mask], cs_lats.ravel()[mask],
        c=data.ravel()[mask], s=16, marker='s',
        cmap=cmap, vmin=vmin, vmax=vmax, zorder=1)

def _finish_ax(ax, title):
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    _set_ticks(ax)
    _add_coast(ax)
    ax.grid(True, alpha=0.3, linewidth=0.5)
    ax.set_title(title, fontsize=11)


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1: sigma_final - sigma_theoretical  (+ panels for each field)
# ═══════════════════════════════════════════════════════════════════════════
print('\nPlotting sigma difference…')

fig, axes = plt.subplots(1, 3, figsize=(19, 4.5))
fig.suptitle('Flux scaling factor σ  on C24 cubed-sphere  '
             '(after 11 iterations)', fontsize=12)

# Panel A — sigma_final (land only)
im0 = cs_scatter(axes[0], sigma_final, vmin=0.0, vmax=2.0, land_only=True)
_finish_ax(axes[0], 'σ final  (iteration 11)')
plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04, label='σ')

# Panel B — sigma_theoretical (land only; ocean is 0 by definition so omit)
im1 = cs_scatter(axes[1], sigma_theo, vmin=0.0, vmax=2.0, land_only=True)
_finish_ax(axes[1], f'σ theoretical  (1/{SCALE_FACTOR} over land)')
plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04, label='σ')

# Panel C — difference (land only)
land_diff = sigma_diff[land_mask]
amax_d = max(np.abs(land_diff).max(), 0.05)
im2 = cs_scatter(axes[2], sigma_diff, vmin=-amax_d, vmax=amax_d, land_only=True)
_finish_ax(axes[2], 'σ_final − σ_theoretical  (land only)')
plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04, label='Δσ')

# Stats annotation on diff panel (land cells only)
axes[2].text(0.01, 0.02,
             f'land: min={land_diff.min():.3f}  max={land_diff.max():.3f}  '
             f'RMSE={np.sqrt((land_diff**2).mean()):.4f}',
             transform=axes[2].transAxes, fontsize=8, color='0.3', va='bottom')

plt.tight_layout()
out1 = os.path.join(PLOT_DIR, 'sigma_diff.png')
plt.savefig(out1, dpi=150, bbox_inches='tight')
plt.close()
print(f'  Saved: {out1}')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2: 2-row flux comparison
#   Row 1: optimised | theoretical | Δ (optimised − theoretical)
#   Row 2: prior     | theoretical | Δ (prior − theoretical)
# Shared colour scale within column type (flux cols share one scale,
# difference cols share one scale).
# ═══════════════════════════════════════════════════════════════════════════
print('\nPlotting TER flux difference…')

# Shared flux colour scale (positive-only sequential)
vmax_f = max(np.nanmax(flux_opt_gC),
             np.nanmax(flux_prior_gC),
             np.nanmax(flux_theo_gC), 1.0)
vmin_f = 0.0

# Fixed difference colour scale as requested
diff_vmin, diff_vmax = -2000, 2000

fig, axes = plt.subplots(2, 3, figsize=(19, 10),
                         gridspec_kw={'bottom': 0.12})
fig.suptitle('TER flux  =  TER × 1.3 × σ   (Jan–Mar 2016 mean,  C24)  '
             '[gC m⁻² yr⁻¹]', fontsize=13)

# ── Row 1: optimised ────────────────────────────────────────────────────
im_opt = cs_scatter(axes[0, 0], flux_opt_gC,
                    vmin=vmin_f, vmax=vmax_f, cmap='YlOrRd', land_only=True)
_finish_ax(axes[0, 0], 'Optimised flux  (TER × 1.3 × σ_final)')

im_th0 = cs_scatter(axes[0, 1], flux_theo_gC,
                    vmin=vmin_f, vmax=vmax_f, cmap='YlOrRd', land_only=True)
_finish_ax(axes[0, 1], 'Theoretical flux  (TER over land)')

im_d0 = cs_scatter(axes[0, 2], flux_diff_opt_gC,
                   vmin=diff_vmin, vmax=diff_vmax, cmap='RdBu_r', land_only=True)
_finish_ax(axes[0, 2], 'Δ  (optimised − theoretical)')
ld0 = flux_diff_opt_gC[land_mask]
axes[0, 2].text(0.01, 0.02,
                f'land: min={ld0.min():.0f}  max={ld0.max():.0f}  '
                f'RMSE={np.sqrt((ld0**2).mean()):.1f}  gC m⁻² yr⁻¹',
                transform=axes[0, 2].transAxes,
                fontsize=8, color='0.3', va='bottom')

# ── Row 2: prior ────────────────────────────────────────────────────────
im_pr = cs_scatter(axes[1, 0], flux_prior_gC,
                   vmin=vmin_f, vmax=vmax_f, cmap='YlOrRd', land_only=True)
_finish_ax(axes[1, 0], 'Prior flux  (TER × 1.3,  σ_prior = 1)')

im_th1 = cs_scatter(axes[1, 1], flux_theo_gC,
                    vmin=vmin_f, vmax=vmax_f, cmap='YlOrRd', land_only=True)
_finish_ax(axes[1, 1], 'Theoretical flux  (TER over land)')

im_d1 = cs_scatter(axes[1, 2], flux_diff_prior_gC,
                   vmin=diff_vmin, vmax=diff_vmax, cmap='RdBu_r', land_only=True)
_finish_ax(axes[1, 2], 'Δ  (prior − theoretical)  ≡ 0.3 × TER')
ld1 = flux_diff_prior_gC[land_mask]
axes[1, 2].text(0.01, 0.02,
                f'land: min={ld1.min():.0f}  max={ld1.max():.0f}  '
                f'RMSE={np.sqrt((ld1**2).mean()):.1f}  gC m⁻² yr⁻¹',
                transform=axes[1, 2].transAxes,
                fontsize=8, color='0.3', va='bottom')

# ── Colorbars below the plots ────────────────────────────────────────────
# Flux colorbar: spans columns 0–1
ax_cbar_f = fig.add_axes([0.10, 0.04, 0.50, 0.025])
fig.colorbar(im_opt, cax=ax_cbar_f, orientation='horizontal',
             label='gC m⁻² yr⁻¹')

# Difference colorbar: spans column 2
ax_cbar_d = fig.add_axes([0.68, 0.04, 0.24, 0.025])
fig.colorbar(im_d0, cax=ax_cbar_d, orientation='horizontal',
             label='gC m⁻² yr⁻¹')

out2 = os.path.join(PLOT_DIR, 'flux_diff.png')
plt.savefig(out2, dpi=150, bbox_inches='tight')
plt.close()
print(f'  Saved: {out2}')

print('\nDone.')
