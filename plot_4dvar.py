#!/usr/bin/env python3
"""
Plot main results from a 4D-Var state file (4dvar_state.npz).

Produces:
  sigma_map.png       - current scaling factor field
  sigma_departure.png - sigma - 1 (departure from prior)
  gradient_map.png    - last gradient dJ/dsigma
  sigma_step.png      - last step (sigma_next - sigma_prev)
  convergence.png     - J and ||grad||_inf vs iteration
                        (requires --output-dir with per-iter snapshots)

Usage:
  python plot_4dvar.py --state-file 4dvar_state.npz \\
                       --output-dir 4dvar_output/   \\
                       --plot-dir   plots/          \\
                       [--nlat 46] [--nlon 72]
  sigma_animation.gif  - sigma evolving across all iterations
                         (requires --output-dir with per-iter snapshots)
  obs_locations.png    - observation locations coloured by iteration
  forcing_profile.png  - mean ± std adjoint forcing vs model level,
                         iteration 1 vs last  (requires forcing_iter_NNN.nc4
                         files in --output-dir, written when SAVE_DIAG=true)
  forcing_latlon.png   - vertically averaged forcing gridded on lat-lon,
                         iteration 1 vs last
"""

import argparse
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

try:
    import xarray as xr
    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False
    print('WARNING: xarray not found — forcing diagnostic plots will be skipped')

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print('WARNING: cartopy not found — land mask will not be drawn')


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def make_grid(nlat, nlon):
    lats = np.linspace(-90.0,  90.0,               nlat)
    lons = np.linspace(-180.0, 180.0 - 360.0/nlon, nlon)
    return lats, lons


def _add_land(ax):
    """Add land fill and coastlines to a cartopy GeoAxes."""
    ax.add_feature(cfeature.LAND,      facecolor='0.85', zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5,    zorder=2)
    ax.add_feature(cfeature.BORDERS,   linewidth=0.3,    zorder=2,
                   linestyle=':', edgecolor='0.4')


def _set_ticks(ax):
    ax.set_xticks(range(-180, 181, 60), crs=ccrs.PlateCarree())
    ax.set_yticks(range(-90,   91, 30), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())
    ax.tick_params(labelsize=8)


def plot_map(data, lats, lons, title, path,
             cmap='RdBu_r', vmin=None, vmax=None, cbar_label=''):
    """Save a single lat-lon map as a PNG."""
    if vmin is None and vmax is None:
        amax = max(np.abs(data).max(), 1e-12)
        vmin, vmax = -amax, amax

    lon2d, lat2d = np.meshgrid(lons, lats)

    if HAS_CARTOPY:
        proj = ccrs.PlateCarree()
        _, ax = plt.subplots(figsize=(12, 5),
                              subplot_kw={'projection': proj})
        _add_land(ax)
        im = ax.pcolormesh(lon2d, lat2d, data,
                           cmap=cmap, vmin=vmin, vmax=vmax,
                           shading='auto', transform=proj, zorder=1)
        ax.set_extent([-180, 180, -90, 90], crs=proj)
        _set_ticks(ax)
        ax.gridlines(alpha=0.3, linewidth=0.5, draw_labels=False)
    else:
        _, ax = plt.subplots(figsize=(12, 5))
        im = ax.pcolormesh(lon2d, lat2d, data,
                           cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_xticks(range(-180, 181, 60))
        ax.set_yticks(range(-90,   91, 30))
        ax.grid(True, alpha=0.3, linewidth=0.5)

    plt.colorbar(im, ax=ax, label=cbar_label, fraction=0.046, pad=0.04)
    ax.set_title(title, fontsize=12)

    stats = (f'min={data.min():.3f}  max={data.max():.3f}  '
             f'mean={data.mean():.3f}  |inf|={np.abs(data).max():.3e}')
    ax.text(0.01, 0.02, stats, transform=ax.transAxes,
            fontsize=8, color='0.3', va='bottom')

    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {path}')


# -----------------------------------------------------------------------
# Forcing diagnostic helpers
# -----------------------------------------------------------------------

def _load_forcing_file(path):
    """Return (lat, lon, forcing) from a forcing_iter_NNN.nc4 file,
    or None if the file cannot be opened."""
    if not HAS_XARRAY:
        return None
    try:
        ds  = xr.open_dataset(path)
        lat = ds['lat_obs'].values.astype(float)
        lon = ds['lon_obs'].values.astype(float)
        f   = ds['forcing'].values.astype(float)   # (n_obs, nlev)
        ds.close()
        return lat, lon, f
    except Exception as exc:
        print(f'  WARNING: could not read {path}: {exc}')
        return None


def plot_obs_time_series(forcing_dir, plot_dir):
    """Bar chart: number of observations per hour from per-checkpoint forcing files.

    Reads every CO2_adjoint_forcing_YYYYMMDD_HHMMz.nc4 in forcing_dir, parses
    the timestamp from the filename, and counts len(lat_obs) in each file.
    The obs count is the same every iteration, so any iteration's forcing dir
    works.
    """
    import re
    import matplotlib.dates as mdates

    pattern = os.path.join(forcing_dir, 'CO2_adjoint_forcing_*.nc4')
    files   = sorted(glob.glob(pattern))
    if not files:
        print(f'  No CO2_adjoint_forcing_*.nc4 files in {forcing_dir} — '
              'obs time series skipped.')
        return

    re_ts = re.compile(r'CO2_adjoint_forcing_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})z\.nc4$')
    times, counts = [], []
    for fpath in files:
        m = re_ts.search(os.path.basename(fpath))
        if not m:
            continue
        yr, mo, dy, hh, mm = (int(x) for x in m.groups())
        t = np.datetime64(f'{yr:04d}-{mo:02d}-{dy:02d}T{hh:02d}:{mm:02d}', 'ns')
        try:
            ds = xr.open_dataset(fpath)
            n  = ds.sizes.get('obs', 0)
            ds.close()
        except Exception:
            n = 0
        times.append(t)
        counts.append(n)

    if not times:
        print('  Could not parse any checkpoint files — obs time series skipped.')
        return

    times  = np.array(times,  dtype='datetime64[ns]')
    counts = np.array(counts, dtype=int)

    # Bin into hourly buckets (sum checkpoints within the same hour)
    hours_trunc = times.astype('datetime64[h]')
    unique_h, inv = np.unique(hours_trunc, return_inverse=True)
    hourly_counts = np.bincount(inv, weights=counts).astype(int)

    # Fill gaps so the bar chart shows zeros for quiet hours
    if len(unique_h) > 1:
        all_h = np.arange(unique_h[0], unique_h[-1] + np.timedelta64(1, 'h'),
                          np.timedelta64(1, 'h'), dtype='datetime64[h]')
        all_c = np.zeros(len(all_h), dtype=int)
        idx = np.searchsorted(all_h, unique_h)
        all_c[idx] = hourly_counts
    else:
        all_h, all_c = unique_h, hourly_counts

    dt_hours = all_h.astype('datetime64[ms]').astype(object)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(dt_hours, all_c, width=1 / 24, align='center',
           color='steelblue', edgecolor='white', linewidth=0.3)
    ax.set_ylabel('Number of observations')
    ax.set_xlabel('Time (UTC)')
    ax.set_title('Hourly observation count over assimilation window')
    ax.grid(True, axis='y', alpha=0.4, linewidth=0.5)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d\n%HZ'))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=max(1, len(all_h) // 12)))
    ax.text(0.01, 0.97, f'Total: {int(all_c.sum()):,} obs  |  {len(files)} checkpoints',
            transform=ax.transAxes, fontsize=9, va='top', color='0.3')
    plt.tight_layout()
    path = os.path.join(plot_dir, 'obs_time_series.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {path}')


def plot_obs_locations(lat, lon, n_iter, plot_dir):
    """Scatter plot of observation locations coloured by obs index (time proxy)."""
    title = f'Observation locations  (iteration {n_iter},  N={len(lat):,})'
    path  = os.path.join(plot_dir, 'obs_locations.png')
    c     = np.arange(len(lat))

    if HAS_CARTOPY:
        proj = ccrs.PlateCarree()
        fig, ax = plt.subplots(figsize=(12, 5), subplot_kw={'projection': proj})
        _add_land(ax)
        sc = ax.scatter(lon, lat, c=c, s=3, cmap='rainbow',
                        transform=proj, zorder=3, alpha=0.7)
        ax.set_extent([-180, 180, -90, 90], crs=proj)
        _set_ticks(ax)
        ax.gridlines(alpha=0.3, linewidth=0.5, draw_labels=False)
    else:
        fig, ax = plt.subplots(figsize=(12, 5))
        sc = ax.scatter(lon, lat, c=c, s=3, cmap='rainbow', alpha=0.7)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_xticks(range(-180, 181, 60))
        ax.set_yticks(range(-90,   91, 30))
        ax.grid(True, alpha=0.3, linewidth=0.5)

    plt.colorbar(sc, ax=ax, label='Observation index (proxy for time in window)',
                 fraction=0.046, pad=0.04)
    ax.set_title(title, fontsize=12)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {path}')


def plot_forcing_profile(data_first, data_last, iter_first, iter_last, plot_dir):
    """Mean ± std of adjoint forcing vs model level, iteration 1 and last overlaid.

    Blue = iteration 1,  red = last iteration.  Shading is ±1 std, alpha=0.20.
    data_first / data_last : (lat, lon, forcing, time) tuples or None.
    forcing shape : (n_obs, nlev),  level index 0 = surface, nlev-1 = TOA.
    """
    fig, ax = plt.subplots(figsize=(6, 8))
    ax.set_title('Adjoint forcing profile  ∂J/∂CO₂  per model level', fontsize=11)

    styles = [(data_first, iter_first, 'steelblue', 'Iteration 1'),
              (data_last,  iter_last,  'crimson',   f'Iteration {iter_last}')]

    for data, it, color, label in styles:
        if data is None:
            continue
        _, _, f = data                  # f: (n_obs, nlev)
        nlev   = f.shape[1]
        levs   = np.arange(1, nlev + 1)
        mean_f = f.mean(axis=0)
        std_f  = f.std(axis=0)
        n_obs  = f.shape[0]

        ax.plot(mean_f, levs, color=color, linewidth=1.8,
                label=f'{label}  (N={n_obs:,})')
        ax.fill_betweenx(levs, mean_f - std_f, mean_f + std_f,
                         alpha=0.20, color=color)

    ax.axvline(0, color='k', linewidth=0.7, linestyle='--')
    ax.set_xlabel('∂J/∂CO₂  [1/(kg CO₂/kg dry air)]', fontsize=9)
    ax.set_ylabel('Model level  (1 = surface)', fontsize=9)
    ax.invert_yaxis()
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    path = os.path.join(plot_dir, 'forcing_profile.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {path}')


def plot_forcing_latlon(data_first, data_last, iter_first, iter_last,
                        lats, lons, plot_dir):
    """Vertically averaged adjoint forcing gridded on the sigma lat-lon grid.

    Three panels: iteration 1 | last iteration | difference (last − first).
    Panels 1 and 2 share the same symmetric colour scale.
    The difference panel uses a separate symmetric scale centred on zero.
    """
    nlat, nlon = len(lats), len(lons)
    lat_edges  = np.concatenate([[lats[0]  - (lats[1]  - lats[0])  / 2],
                                 (lats[:-1] + lats[1:]) / 2,
                                 [lats[-1] + (lats[-1] - lats[-2]) / 2]])
    lon_edges  = np.concatenate([[lons[0]  - (lons[1]  - lons[0])  / 2],
                                 (lons[:-1] + lons[1:]) / 2,
                                 [lons[-1] + (lons[-1] - lons[-2]) / 2]])

    def _grid(data):
        if data is None:
            return None, 0
        lat_o, lon_o, f = data
        vals = f.mean(axis=1)
        sums, _, _ = np.histogram2d(lat_o, lon_o, bins=[lat_edges, lon_edges],
                                    weights=vals)
        cnts, _, _ = np.histogram2d(lat_o, lon_o, bins=[lat_edges, lon_edges])
        with np.errstate(invalid='ignore'):
            grid = np.where(cnts > 0, sums / cnts, np.nan)
        return grid, int(cnts.sum())

    grid1, n1 = _grid(data_first)
    grid2, n2 = _grid(data_last)

    if grid1 is None and grid2 is None:
        print('  No forcing data available for lat-lon plot — skipping.')
        return

    # Shared scale for panels 1 & 2
    valid = [g for g in [grid1, grid2] if g is not None]
    amax  = max(np.nanmax(np.abs(g)) for g in valid)
    amax  = max(amax, 1e-30)

    # Difference panel (last − first) and its independent scale
    if grid1 is not None and grid2 is not None:
        diff  = grid2 - grid1
        dmax  = np.nanmax(np.abs(diff))
        dmax  = max(dmax, 1e-30)
    else:
        diff, dmax = None, 1.0

    lon2d, lat2d = np.meshgrid(lons, lats)

    proj_kw = {'projection': ccrs.PlateCarree()} if HAS_CARTOPY else {}
    fig, axes = plt.subplots(1, 3, figsize=(19, 5), subplot_kw=proj_kw)
    fig.suptitle(
        'Vertically averaged adjoint forcing  ∂J/∂CO₂  [1/(kg CO₂/kg dry air)]',
        fontsize=11)

    def _draw(ax, grid, vmin, vmax, cmap, title, n_obs):
        if grid is None:
            ax.set_title(f'{title}\n(no data)', fontsize=10)
            return
        if HAS_CARTOPY:
            _add_land(ax)
            im = ax.pcolormesh(lon2d, lat2d, grid, cmap=cmap,
                               vmin=vmin, vmax=vmax, shading='auto',
                               transform=ccrs.PlateCarree(), zorder=1)
            ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
            _set_ticks(ax)
            ax.gridlines(alpha=0.3, linewidth=0.5, draw_labels=False)
        else:
            im = ax.pcolormesh(lon2d, lat2d, grid, cmap=cmap,
                               vmin=vmin, vmax=vmax, shading='auto')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            ax.set_xticks(range(-180, 181, 60))
            ax.set_yticks(range(-90,   91, 30))
            ax.grid(True, alpha=0.3, linewidth=0.5)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        lbl = f'{title}  (N={n_obs:,})' if n_obs else title
        ax.set_title(lbl, fontsize=10)

    _draw(axes[0], grid1, -amax, amax, 'RdBu_r',
          f'Iteration {iter_first}', n1)
    _draw(axes[1], grid2, -amax, amax, 'RdBu_r',
          f'Iteration {iter_last}', n2)
    _draw(axes[2], diff,  -dmax, dmax, 'PiYG',
          f'Difference  (iter {iter_last} − iter {iter_first})', 0)

    plt.tight_layout()
    path = os.path.join(plot_dir, 'forcing_latlon.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {path}')


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--state-file',  required=True,
                        help='Path to 4dvar_state.npz')
    parser.add_argument('--output-dir',  default=None,
                        help='Directory with 4dvar_state_iter_NNN.npz and '
                             'forcing_iter_NNN.nc4 snapshots')
    parser.add_argument('--forcing-dir', default=None,
                        help='Directory containing CO2_adjoint_forcing_*.nc4 files '
                             '(for obs time series plot)')
    parser.add_argument('--plot-dir',    default=None,
                        help='Where to save PNGs (default: same dir as state file)')
    parser.add_argument('--nlat', type=int, default=None,
                        help='Override grid rows (read from state file by default)')
    parser.add_argument('--nlon', type=int, default=None,
                        help='Override grid cols (read from state file by default)')
    args = parser.parse_args()

    plot_dir = args.plot_dir or os.path.dirname(os.path.abspath(args.state_file))
    os.makedirs(plot_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Load state
    # ------------------------------------------------------------------
    s  = np.load(args.state_file)
    it = int(s['iteration']) - 1        # completed iterations
    J  = float(s['J_prev'])

    nlat = args.nlat or int(s['nlat'])
    nlon = args.nlon or int(s['nlon'])

    sigma      = s['x_next'].reshape(nlat, nlon)
    sigma_prev = s['x_prev'].reshape(nlat, nlon)
    grad       = s['g_prev'].reshape(nlat, nlon)
    lats, lons = make_grid(nlat, nlon)

    print(f'State file      : {args.state_file}')
    print(f'Iterations done : {it}')
    print(f'J               : {J:.6e}')
    print(f'sigma           : min={sigma.min():.4f}  max={sigma.max():.4f}  '
          f'mean={sigma.mean():.4f}')
    print(f'|grad|_inf      : {np.abs(grad).max():.4e}')
    print()

    # ------------------------------------------------------------------
    # 1. Sigma map
    # ------------------------------------------------------------------
    plot_map(sigma, lats, lons,
             title=f'CO₂ flux scaling factor σ  (after iteration {it})',
             path=os.path.join(plot_dir, 'sigma_map.png'),
             cmap='RdBu_r', vmin=0.0, vmax=2.0,
             cbar_label='σ  [ ]')

    # ------------------------------------------------------------------
    # 2. Sigma departure from prior  (sigma - 1)
    # ------------------------------------------------------------------
    plot_map(sigma - 1.0, lats, lons,
             title=f'Departure from prior  σ − 1  (after iteration {it})',
             path=os.path.join(plot_dir, 'sigma_departure.png'),
             cmap='RdBu_r',
             cbar_label='σ − 1  [ ]')

    # ------------------------------------------------------------------
    # 3. Gradient map
    # ------------------------------------------------------------------
    plot_map(grad, lats, lons,
             title=f'Total gradient ∂J/∂σ  (iteration {it})',
             path=os.path.join(plot_dir, 'gradient_map.png'),
             cmap='RdBu_r',
             cbar_label='∂J/∂σ')

    # ------------------------------------------------------------------
    # 4. Last step  (sigma_next - sigma_prev)
    # ------------------------------------------------------------------
    if it >= 2:
        plot_map(sigma - sigma_prev, lats, lons,
                 title=f'Last optimizer step  Δσ = σ_next − σ_prev  (iteration {it})',
                 path=os.path.join(plot_dir, 'sigma_step.png'),
                 cmap='RdBu_r',
                 cbar_label='Δσ  [ ]')

    # ------------------------------------------------------------------
    # 5. Convergence history (requires per-iter snapshot files)
    # ------------------------------------------------------------------
    snap_dir = args.output_dir
    if snap_dir and os.path.isdir(snap_dir):
        snap_files = sorted(glob.glob(
            os.path.join(snap_dir, '4dvar_state_iter*.npz')))

        iters, Js, J_obss, g_infs = [], [], [], []
        for f in snap_files:
            sv  = np.load(f)
            i   = int(sv['iteration']) - 1
            if i < 1:
                continue
            iters.append(i)
            Js.append(float(sv['J_prev']))
            J_obss.append(float(sv['J_obs_prev']) if 'J_obs_prev' in sv else None)
            g_infs.append(float(np.abs(sv['g_prev']).max()))

        if len(iters) >= 2:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

            ax1.plot(iters, Js, 'k-o', markersize=5, linewidth=1.5, label='J (total)')
            if any(v is not None for v in J_obss):
                J_obss_clean = [v for v in J_obss if v is not None]
                iters_obs    = [iters[k] for k, v in enumerate(J_obss) if v is not None]
                ax1.plot(iters_obs, J_obss_clean,
                         'b--s', markersize=4, linewidth=1.2, label='J_obs')
            ax1.set_ylabel('Cost  J')
            ax1.set_title('4D-Var convergence history')
            ax1.legend(fontsize=9)
            ax1.grid(True, which='both', alpha=0.3)

            ax2.semilogy(iters, g_infs, 'g-o', markersize=5, linewidth=1.5)
            ax2.set_ylabel('||∇J||∞')
            ax2.set_xlabel('Iteration')
            ax2.grid(True, which='both', alpha=0.3)

            plt.tight_layout()
            conv_path = os.path.join(plot_dir, 'convergence.png')
            plt.savefig(conv_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f'  Saved: {conv_path}')
        elif snap_files:
            print('  Convergence plot needs at least 2 iterations — skipping.')
    elif snap_dir:
        print(f'  WARNING: --output-dir not found: {snap_dir}')

    # ------------------------------------------------------------------
    # 6. Animation of sigma evolution across all iterations
    # ------------------------------------------------------------------
    if snap_dir and os.path.isdir(snap_dir):
        snap_files = sorted(glob.glob(
            os.path.join(snap_dir, '4dvar_state_iter*.npz')))

        frames = []
        frame_iters = []
        for f in snap_files:
            sv = np.load(f)
            i  = int(sv['iteration']) - 1
            if i < 1:
                continue
            frames.append(sv['x_prev'].reshape(nlat, nlon))
            frame_iters.append(i)

        if len(frames) >= 2:
            lon2d, lat2d = np.meshgrid(lons, lats)

            if HAS_CARTOPY:
                proj = ccrs.PlateCarree()
                fig, ax = plt.subplots(figsize=(12, 5),
                                       subplot_kw={'projection': proj})
                _add_land(ax)
                im = ax.pcolormesh(lon2d, lat2d, frames[0],
                                   cmap='RdBu_r', vmin=0.0, vmax=2.0,
                                   shading='auto', transform=proj, zorder=1)
                ax.set_extent([-180, 180, -90, 90], crs=proj)
                _set_ticks(ax)
                ax.gridlines(alpha=0.3, linewidth=0.5, draw_labels=False)
            else:
                fig, ax = plt.subplots(figsize=(12, 5))
                im = ax.pcolormesh(lon2d, lat2d, frames[0],
                                   cmap='RdBu_r', vmin=0.0, vmax=2.0,
                                   shading='auto')
                ax.set_xlabel('Longitude')
                ax.set_ylabel('Latitude')
                ax.set_xticks(range(-180, 181, 60))
                ax.set_yticks(range(-90,   91, 30))
                ax.grid(True, alpha=0.3, linewidth=0.5)

            plt.colorbar(im, ax=ax, label='σ  [ ]',
                         fraction=0.046, pad=0.04)
            title      = ax.set_title('')
            stats_text = ax.text(0.01, 0.02, '', transform=ax.transAxes,
                                 fontsize=8, color='0.3', va='bottom')

            def update(frame_idx):
                data = frames[frame_idx]
                i    = frame_iters[frame_idx]
                im.set_array(data.ravel())
                title.set_text(
                    f'CO₂ flux scaling factor σ  (iteration {i})')
                stats_text.set_text(
                    f'min={data.min():.3f}  max={data.max():.3f}  '
                    f'mean={data.mean():.3f}  |inf|={np.abs(data).max():.3e}')
                return im, title, stats_text

            anim = animation.FuncAnimation(
                fig, update, frames=len(frames),
                interval=500, blit=False)

            gif_path = os.path.join(plot_dir, 'sigma_animation.gif')
            writer   = animation.PillowWriter(fps=2)
            anim.save(gif_path, writer=writer, dpi=100)
            plt.close()
            print(f'  Saved: {gif_path}')
        elif snap_files:
            print('  Animation needs at least 2 iterations — skipping.')

    # ------------------------------------------------------------------
    # 7–9. Forcing diagnostics  (requires forcing_iter_NNN.nc4 files)
    # ------------------------------------------------------------------
    if not HAS_XARRAY:
        pass   # warning already printed at import time
    elif snap_dir and os.path.isdir(snap_dir):
        forcing_files = sorted(glob.glob(
            os.path.join(snap_dir, 'forcing_iter*.nc4')))

        if not forcing_files:
            print('  No forcing_iter_NNN.nc4 files found — '
                  'skipping forcing diagnostic plots '
                  '(run with SAVE_DIAG=true to generate them).')
        else:
            data_first = _load_forcing_file(forcing_files[0])
            data_last  = _load_forcing_file(forcing_files[-1])
            iter_first = int(forcing_files[0].split('forcing_iter')[1].split('.')[0])
            iter_last  = int(forcing_files[-1].split('forcing_iter')[1].split('.')[0])
            print(f'  Forcing files found: {len(forcing_files)}  '
                  f'(iter {iter_first} … {iter_last})')

            # 7. Observation locations  (use first-iteration file; locations
            #    are the same every iteration since the obs dataset is fixed)
            if data_first is not None:
                plot_obs_locations(data_first[0], data_first[1],
                                   iter_first, plot_dir)

            # 8. Forcing profile: mean ± std vs model level, overlaid
            plot_forcing_profile(data_first, data_last,
                                 iter_first, iter_last, plot_dir)

            # 9. Forcing lat-lon map (vertically averaged), 3 panels
            plot_forcing_latlon(data_first, data_last,
                                iter_first, iter_last,
                                lats, lons, plot_dir)

    # ------------------------------------------------------------------
    # 10. Obs time series  (reads per-checkpoint forcing files directly;
    #     obs count is identical every iteration — no rerun needed)
    # ------------------------------------------------------------------
    forcing_dir = args.forcing_dir
    if forcing_dir and os.path.isdir(forcing_dir):
        plot_obs_time_series(forcing_dir, plot_dir)
    elif forcing_dir:
        print(f'  WARNING: --forcing-dir not found: {forcing_dir}')


if __name__ == '__main__':
    main()
