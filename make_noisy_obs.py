#!/usr/bin/env python3
"""Add seeded Gaussian observation noise to an OSSE semi-observation file (O9, 2026-10-21).

    make_noisy_obs.py --in <obs dir>/semi_obs_xco2_2016Q1.nc --out-dir <new obs dir> --seed 1 [--scale 1.0]

y_i = H(x_truth)_i + eps_i,  eps_i ~ N(0, (scale * sigma_i)^2), sigma_i = the sounding's OCO-2 L2 XCO2 uncertainty
(xco2_uncertainty, stored in mol/mol), i.e. the same sigma the cost function divides by (co2_adjoint_forcing_osse.py
reads it from the L2 files). With scale = 1 the noise is consistent with R: at the true fluxes E[J_obs] = N/2.
Writes <out-dir>/<same name> (xco2_truth noisy, xco2_truth_clean and xco2_noise added) and links the truth emissions
file of the input directory, so the new directory works as an OSSE obs directory (setup_exp.sh --obs-dir).
"""
import argparse, os
import numpy as np
import xarray as xr

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--in', dest='inp', required=True); p.add_argument('--out-dir', required=True)
p.add_argument('--seed', type=int, required=True); p.add_argument('--scale', type=float, default=1.0)
a = p.parse_args()
ds = xr.open_dataset(a.inp).load()
sig_ppm = ds['xco2_uncertainty'].values.astype(np.float64) * 1e6          # stored in mol/mol (attribute says ppm)
if not (0.05 < np.median(sig_ppm) < 10):
    raise SystemExit(f'unexpected xco2_uncertainty scale (median {np.median(sig_ppm):.3g} ppm after x1e6)')
eps = np.random.default_rng(a.seed).standard_normal(sig_ppm.size) * a.scale * sig_ppm
clean = ds['xco2_truth'].values.astype(np.float64)
out = ds.copy()
out['xco2_truth_clean'] = ('time', clean.astype(np.float32), {'units': 'ppm', 'long_name': 'truth XCO2 without noise'})
out['xco2_noise'] = ('time', eps.astype(np.float32), {'units': 'ppm', 'long_name': 'added observation noise'})
out['xco2_truth'] = ('time', (clean + eps).astype(np.float32),
                     dict(ds['xco2_truth'].attrs, long_name='truth XCO2 + noise (the OSSE observation y)'))
for v in ('xco2_l2', 'xco2_apriori', 'xco2_uncertainty'):                   # fix the mislabelled units in the copy
    out[v].attrs['units'] = 'mol/mol'
out.attrs.update(noise=f'Gaussian, seed {a.seed}, sigma = {a.scale} x xco2_uncertainty (per sounding)',
                 noise_source=os.path.abspath(a.inp))
os.makedirs(a.out_dir, exist_ok=True)
dst = os.path.join(a.out_dir, os.path.basename(a.inp))
out.to_netcdf(dst)
src_dir = os.path.dirname(os.path.abspath(a.inp))
for f in os.listdir(src_dir):
    if f.startswith('truth_emissions') and not os.path.exists(os.path.join(a.out_dir, f)):
        os.symlink(os.path.join(src_dir, f), os.path.join(a.out_dir, f))
z = eps / (a.scale * sig_ppm)
print(f'{dst}: {eps.size} soundings, noise mean {eps.mean():+.4f} ppm, std {eps.std():.4f} ppm, '
      f'normalized mean {z.mean():+.4f} std {z.std():.4f}; J_obs of the noise alone = {0.5 * (eps ** 2 / sig_ppm ** 2).sum():.1f} '
      f'(expected N/2 x scale^2 = {0.5 * eps.size * a.scale ** 2:.1f})')
