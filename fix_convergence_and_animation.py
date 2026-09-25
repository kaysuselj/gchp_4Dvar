#!/usr/bin/env python3
"""
Regenerate convergence.png, convergence_log.png, and daily_forcing_animation.gif
for the osse-exp_constant_offset_v1 experiment with two fixes:

  1. Convergence plots: cap history at the max snap-file iteration (11), so the
     stale iteration-12 entry from the previous experiment in the log is dropped.

  2. Animation: derive the colour scale from iterations 2+ only, so that the
     ~10x-larger iteration-1 forcing does not saturate every subsequent frame.
"""

import os
import glob
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import xarray as xr

# ── Paths ──────────────────────────────────────────────────────────────────
HERE      = os.path.dirname(os.path.abspath(__file__))
SNAP_DIR  = os.path.join(HERE, '4dvar_output_osse')
PLOT_DIR  = os.path.join(HERE, 'plots_osse', '20260923_092236')

# ── 1. Convergence plots (cap at max snap-file iteration) ──────────────────

snap_files = sorted(glob.glob(os.path.join(SNAP_DIR, '4dvar_state_iter*.npz')))
if not snap_files:
    sys.exit('No 4dvar_state_iter*.npz found in ' + SNAP_DIR)

hist = {}
for f in snap_files:
    sv = np.load(f)
    i  = int(sv['iteration']) - 1
    if i < 1:
        continue
    hist[i] = {
        'J':     float(sv['J_prev']),
        'J_obs': float(sv['J_obs_prev']) if 'J_obs_prev' in sv.files else None,
        'g_inf': float(np.abs(sv['g_prev']).max()),
    }

# Cap: only keep iterations backed by a snap file (prevents log bleed-through)
max_snap = max(hist.keys())
hist = {k: v for k, v in hist.items() if k <= max_snap}

iters  = sorted(hist)
Js     = [hist[i]['J']     for i in iters]
J_obss = [hist[i]['J_obs'] for i in iters]
g_infs = [hist[i]['g_inf'] for i in iters]

print(f'Convergence history: iterations {iters[0]}–{iters[-1]}  ({len(iters)} points)')

# ── 1a. convergence.png ───────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

ax1.plot(iters, Js, 'k-o', markersize=5, linewidth=1.5, label='J (total)')
valid_obs = [(iters[k], v) for k, v in enumerate(J_obss) if v is not None]
if valid_obs:
    it_obs, jo = zip(*valid_obs)
    ax1.plot(it_obs, jo, 'b--s', markersize=4, linewidth=1.2, label='J_obs')
ax1.set_ylabel('Cost  J')
ax1.set_ylim(bottom=0)
ax1.set_title('4D-Var convergence history')
ax1.legend(fontsize=9)
ax1.grid(True, which='both', alpha=0.3)

iters_g  = [iters[k] for k, v in enumerate(g_infs) if v is not None]
g_infs_g = [v for v in g_infs if v is not None]
ax2.semilogy(iters_g, g_infs_g, 'g-o', markersize=5, linewidth=1.5)
ax2.set_ylabel('||∇J||∞')
ax2.set_xlabel('Iteration')
ax2.grid(True, which='both', alpha=0.3)

plt.tight_layout()
out = os.path.join(PLOT_DIR, 'convergence.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'  Saved: {out}')

# ── 1b. convergence_log.png ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
ax.semilogy(iters, Js, 'k-o', markersize=5, linewidth=1.5, label='J (total)')
if valid_obs:
    ax.semilogy(it_obs, jo, 'b--s', markersize=4, linewidth=1.2, label='J_obs')
ax.set_xlabel('Iteration')
ax.set_ylabel('Cost  J  (log scale)')
ax.set_title('4D-Var cost function  (log y-scale)')
ax.legend(fontsize=9)
ax.grid(True, which='both', alpha=0.3)
plt.tight_layout()
out = os.path.join(PLOT_DIR, 'convergence_log.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'  Saved: {out}')

# ── 2. daily_forcing_animation.gif (scale from iter 2+ only) ──────────────

df_files = sorted(glob.glob(os.path.join(SNAP_DIR, 'daily_forcing_iter*.nc4')))
frames, frame_iters, dates = [], [], None

for fpath in df_files:
    itr = int(os.path.basename(fpath)
              .split('daily_forcing_iter')[1].split('.')[0])
    try:
        ds = xr.open_dataset(fpath)
        mf = ds['mean_forcing'].values.astype(float)
        if dates is None:
            dates = ds['date'].values
        ds.close()
    except Exception as exc:
        print(f'  WARNING: could not read {fpath}: {exc}')
        continue
    frames.append(mf)
    frame_iters.append(itr)

if len(frames) < 2:
    sys.exit('Need at least 2 daily_forcing_iter*.nc4 files')

n_days, nlev = frames[0].shape
day_labels = [str(np.datetime_as_string(d, unit='D')) for d in dates]

# Colour scale: use iterations 2+ so iter-1's 10x-larger amplitude does not
# saturate every other frame.  Iteration 1 will be clipped, but all subsequent
# frames will show their true structure.
amax_all  = max(np.nanmax(np.abs(m)) for m in frames)
amax_rest = max(np.nanmax(np.abs(m)) for m in frames[1:])
amax = amax_rest        # iter 1 saturates; iters 2+ use the full scale
print(f'Animation colour scale: ±{amax:.3e}  '
      f'(iter-1 max was {amax_all:.3e}, now saturated)')

fig, ax = plt.subplots(figsize=(10, 6))
im = ax.imshow(frames[0].T, aspect='auto', origin='lower',
               cmap='RdBu_r', vmin=-amax, vmax=amax,
               extent=[-0.5, n_days - 0.5, 0.5, nlev + 0.5])
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='∂J/∂CO₂')
ax.set_xlabel('Day')
ax.set_ylabel('Model level  (1 = surface)')
ax.set_ylim(0.5, nlev + 0.5)
step = max(1, n_days // 10)
ax.set_xticks(range(0, n_days, step))
ax.set_xticklabels(day_labels[::step], rotation=30, ha='right', fontsize=8)

def _title(k):
    note = '  [scale saturated — iter-1 amplitude ~10×]' if frame_iters[k] == 1 else ''
    return (f'Daily mean adjoint forcing ∂J/∂CO₂ per model level  '
            f'[1/(kg CO₂/kg dry air)]  —  iteration {frame_iters[k]}{note}')

title = ax.set_title(_title(0), fontsize=10)
fig.tight_layout()

def update(k):
    im.set_array(frames[k].T)
    title.set_text(_title(k))
    return im, title

anim = animation.FuncAnimation(fig, update, frames=len(frames),
                                interval=700, blit=False)
gif_path = os.path.join(PLOT_DIR, 'daily_forcing_animation.gif')
anim.save(gif_path, writer=animation.PillowWriter(fps=1.5), dpi=100)
plt.close()
print(f'  Saved: {gif_path}')
