#!/usr/bin/env python3
"""
Plot 4D-Var optimizer behaviour, reconstructed from the per-iteration state
snapshots (4dvar_state_iterNNN.npz) an optimizer run writes into its output dir.

Four panels vs iteration: cost (J, J_obs, J_b), gradient norms
(|g|, |g_obs|, |g_b|), control-sigma range (min/mean/max with the 1±sigma_b
band), and step quality rho = dJ_actual/dJ_pred (reconstructed from consecutive
states).  g_obs/g_b are split using J_b = 0.5||x-1||^2/sigma_b^2.

Modes:
  default          one summary figure over all iterations.
  --per-iteration  one FRAME per iteration (cumulative through iter N), fixed
                   axes across frames so the sequence reads like a flip-book;
                   iteration 1 is skipped (rho undefined, J_b=0, sigma=1).

Usage:
  python plot_optimizer_history.py --output-dir <run>/4dvar_output_osse \
      --sigma-b 0.2 --plot-dir <run>/plots_osse/history [--per-iteration]
"""
import os
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_trajectory(output_dir, sigma_b):
    """Reconstruct per-iteration arrays from 4dvar_state_iterNNN.npz snapshots.

    Each snapshot stores the values AT that iteration: J_prev/J_obs_prev/
    J_b_prev, x_prev (=sigma evaluated), g_prev (=gradient).  rho[N] is the
    step quality of the move N-1 -> N.
    """
    files = sorted(glob.glob(os.path.join(output_dir, '4dvar_state_iter*.npz')))
    if not files:
        raise SystemExit(f'no 4dvar_state_iter*.npz in {output_dir}')
    N = len(files)
    d = dict(it=np.arange(1, N + 1),
             J=np.zeros(N), Jo=np.zeros(N), Jb=np.zeros(N),
             smin=np.zeros(N), smax=np.zeros(N), smean=np.zeros(N),
             ginf=np.zeros(N), goinf=np.zeros(N), gbinf=np.zeros(N),
             rho=np.full(N, np.nan))
    xprev = np.zeros((N, np.load(files[0])['x_prev'].size))
    gprev = np.zeros_like(xprev)
    for k, f in enumerate(files):
        s = np.load(f)
        d['J'][k]  = float(s['J_prev'])
        d['Jo'][k] = float(s['J_obs_prev'])
        d['Jb'][k] = float(s['J_b_prev'])
        x, g = s['x_prev'], s['g_prev']
        xprev[k], gprev[k] = x, g
        d['smin'][k], d['smax'][k], d['smean'][k] = x.min(), x.max(), x.mean()
        gb = (x - 1.0) / sigma_b ** 2
        go = g - gb
        d['ginf'][k]  = np.abs(g).max()
        d['goinf'][k] = np.abs(go).max()
        d['gbinf'][k] = np.abs(gb).max()
    for k in range(1, N):
        pred = float(np.dot(gprev[k - 1], xprev[k] - xprev[k - 1]))
        d['rho'][k] = (d['J'][k] - d['J'][k - 1]) / pred if pred != 0 else np.nan
    # running best-so-far cost and the iteration that holds it
    d['Jbest']  = np.minimum.accumulate(d['J'])
    d['itbest'] = np.array([int(d['it'][np.argmin(d['J'][:k + 1])])
                            for k in range(N)])
    return d, N


def _pos(a):
    a = np.asarray(a, float)
    return a[np.isfinite(a) & (a > 0)]


def draw(d, upto, sigma_b, title, limits=None):
    """Draw the behaviour figure (5 panels) using data through iteration `upto`.
    Returns (fig, it_best_now, J_best_now) for the current frame."""
    sl = slice(0, upto)
    it = d['it'][sl]
    it_best_now = int(d['itbest'][upto - 1])
    J_best_now  = float(d['Jbest'][upto - 1])
    fig, ax = plt.subplots(2, 3, figsize=(18, 9))

    a = ax[0, 0]
    a.semilogy(it, d['J'][sl], 'o-', label='J')
    a.semilogy(it, d['Jo'][sl], 's-', label='J_obs')
    a.semilogy(it, np.where(d['Jb'][sl] > 0, d['Jb'][sl], np.nan), '^-', label='J_b')
    a.set_title('Cost function'); a.set_xlabel('iteration'); a.set_ylabel('J')
    a.legend(); a.grid(True, which='both', alpha=.3)

    # NEW: running best-so-far cost, with the best iteration marked
    a = ax[0, 1]
    a.semilogy(it, d['J'][sl], 'o-', color='0.7', label='J (current)')
    a.semilogy(it, d['Jbest'][sl], 'D-', color='C2', label='J best-so-far')
    a.scatter([it_best_now], [J_best_now], s=140, facecolors='none',
              edgecolors='r', linewidths=2, zorder=5, label='current best')
    a.set_title(f'Best cost so far — iter {it_best_now}, J={J_best_now:.4e}')
    a.set_xlabel('iteration'); a.set_ylabel('J best'); a.legend(); a.grid(True, which='both', alpha=.3)

    a = ax[0, 2]
    a.semilogy(it, d['ginf'][sl], 'o-', label='|g|inf')
    a.semilogy(it, d['goinf'][sl], 's-', label='|g_obs|inf')
    a.semilogy(it, np.where(d['gbinf'][sl] > 0, d['gbinf'][sl], np.nan), '^-', label='|g_b|inf')
    a.set_title('Gradient norms'); a.set_xlabel('iteration'); a.set_ylabel('|grad|_inf')
    a.legend(); a.grid(True, which='both', alpha=.3)

    a = ax[1, 0]
    a.plot(it, d['smax'][sl], '^-', label='sigma max')
    a.plot(it, d['smean'][sl], 'o-', label='sigma mean')
    a.plot(it, d['smin'][sl], 'v-', label='sigma min')
    a.axhline(1, color='k', lw=.7, ls=':')
    a.axhline(1 + sigma_b, color='r', lw=.7, ls='--', label=f'1±sigma_b')
    a.axhline(max(1 - sigma_b, 0), color='r', lw=.7, ls='--')
    a.set_title('Control sigma range'); a.set_xlabel('iteration'); a.set_ylabel('sigma')
    a.legend(); a.grid(True, alpha=.3)

    a = ax[1, 1]
    a.axhline(1, color='g', lw=.7, ls=':', label='rho=1 (model exact)')
    a.axhline(0, color='r', lw=.7, ls='--', label='rho=0 (J unchanged)')
    a.plot(it, d['rho'][sl], 'o-')
    a.set_title('Step quality rho = dJ_actual / dJ_pred')
    a.set_xlabel('iteration'); a.set_ylabel('rho'); a.legend(); a.grid(True, alpha=.3)

    ax[1, 2].axis('off')   # spare slot

    if limits:
        for aa in (ax[0, 0], ax[0, 1], ax[0, 2], ax[1, 0], ax[1, 1]):
            aa.set_xlim(*limits['x'])
        ax[0, 0].set_ylim(*limits['J']); ax[0, 1].set_ylim(*limits['J'])
        ax[0, 2].set_ylim(*limits['g'])
        ax[1, 0].set_ylim(*limits['s']); ax[1, 1].set_ylim(*limits['rho'])

    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    return fig, it_best_now, J_best_now


def full_limits(d, N):
    Jall = np.concatenate([_pos(d['J']), _pos(d['Jo']), _pos(d['Jb'])])
    gall = np.concatenate([_pos(d['ginf']), _pos(d['goinf']), _pos(d['gbinf'])])
    sall = np.concatenate([d['smin'], d['smax']])
    rall = d['rho'][np.isfinite(d['rho'])]
    return dict(
        x=(0.5, N + 0.5),
        J=(Jall.min() * 0.6, Jall.max() * 1.8),
        g=(gall.min() * 0.6, gall.max() * 1.8),
        s=(min(sall.min(), 1 - 0.3), sall.max() * 1.08),
        rho=(min(-1.0, rall.min() * 1.1 if rall.size else -1.0),
             max(2.5, rall.max() * 1.1 if rall.size else 2.5)),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--output-dir', required=True,
                   help='dir with 4dvar_state_iterNNN.npz snapshots')
    p.add_argument('--sigma-b', type=float, default=0.2)
    p.add_argument('--plot-dir', required=True)
    p.add_argument('--per-iteration', action='store_true',
                   help='one frame per iteration (cumulative), fixed axes')
    p.add_argument('--start-iter', type=int, default=2,
                   help='first iteration to emit a frame for (default 2; '
                        'iteration 1 is trivial)')
    p.add_argument('--tag', default='',
                   help='label added to titles/filenames (e.g. experiment name)')
    args = p.parse_args()

    d, N = load_trajectory(args.output_dir, args.sigma_b)
    os.makedirs(args.plot_dir, exist_ok=True)
    lab = f' [{args.tag}]' if args.tag else ''

    if args.per_iteration:
        lim = full_limits(d, N)
        made = []
        for upto in range(max(args.start_iter, 2), N + 1):
            title = (f'4D-Var optimizer behaviour{lab} — step {upto}/{N} '
                     f'(sigma_b={args.sigma_b})')
            fig, itb, Jb = draw(d, upto, args.sigma_b, title, limits=lim)
            fig.suptitle(title + f'   |   best so far: iter {itb}, J={Jb:.4e}',
                         fontsize=13)
            fp = os.path.join(args.plot_dir, f'optimizer_behaviour_iter{upto:03d}.png')
            fig.savefig(fp, dpi=110); plt.close(fig); made.append(fp)
        print(f'Wrote {len(made)} per-iteration frames to {args.plot_dir}/')
        print(f'  first: {os.path.basename(made[0])}   last: {os.path.basename(made[-1])}')
    else:
        fig, itb, Jb = draw(d, N, args.sigma_b, '')
        fig.suptitle(f'4D-Var optimizer behaviour{lab} — {N} iterations '
                     f'(sigma_b={args.sigma_b})   |   best: iter {itb}, J={Jb:.4e}',
                     fontsize=13)
        fp = os.path.join(args.plot_dir, 'optimizer_behaviour.png')
        fig.savefig(fp, dpi=130); plt.close(fig)
        print(f'Wrote {fp}  (best: iter {itb}, J={Jb:.4e})')


if __name__ == '__main__':
    main()
