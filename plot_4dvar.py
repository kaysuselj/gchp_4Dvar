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
  sigma_animation.gif - sigma evolving across all iterations
                        (requires --output-dir with per-iter snapshots)
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
# Main
# -----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--state-file',  required=True,
                        help='Path to 4dvar_state.npz')
    parser.add_argument('--output-dir',  default=None,
                        help='Directory with 4dvar_state_iter_NNN.npz snapshots '
                             '(for convergence plot)')
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

        iters, Js, J_obss, J_bs, g_infs = [], [], [], [], []
        for f in snap_files:
            sv  = np.load(f)
            i   = int(sv['iteration']) - 1
            if i < 1:
                continue
            iters.append(i)
            Js.append(float(sv['J_prev']))
            J_obss.append(float(sv['J_obs_prev']) if 'J_obs_prev' in sv else None)
            J_bs.append(float(sv['J_b_prev'])   if 'J_b_prev'   in sv else None)
            g_infs.append(float(np.abs(sv['g_prev']).max()))

        if len(iters) >= 2:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

            ax1.semilogy(iters, Js, 'k-o', markersize=5, linewidth=1.5, label='J (total)')
            if any(v is not None for v in J_obss):
                J_obss_clean = [v for v in J_obss if v is not None]
                iters_obs    = [iters[k] for k, v in enumerate(J_obss) if v is not None]
                ax1.semilogy(iters_obs, J_obss_clean,
                             'b--s', markersize=4, linewidth=1.2, label='J_obs')
            if any(v is not None for v in J_bs):
                J_bs_clean = [v for v in J_bs if v is not None]
                iters_b    = [iters[k] for k, v in enumerate(J_bs) if v is not None]
                ax1.semilogy(iters_b, J_bs_clean,
                             'r--^', markersize=4, linewidth=1.2, label='J_b')
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


if __name__ == '__main__':
    main()
