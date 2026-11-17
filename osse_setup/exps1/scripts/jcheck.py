import sys, numpy as np, xarray as xr
X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
obs = xr.open_dataset(f"{X}/obs_c24_1510_1705/semi_obs_xco2_201510_201705.nc")
u = obs["xco2_uncertainty"].values
u = u * 1e6 if np.nanmedian(u) < 1e-3 else u
for e in ("OSSE1_C24", "OSSE2_C24"):
    J = float(open(f"{X}/{e}/forcing_files_osse/J_value.txt").read())
    N = int(open(f"{X}/{e}/forcing_files_osse/N_obs.txt").read())
    rms_norm = np.sqrt(2 * J / N)
    print(f"{e}: J_obs {J:.4g}, N {N}, 2J/N {2*J/N:.1f}, rms normalized misfit {rms_norm:.1f} sigma "
          f"= ~{rms_norm * np.median(u):.2f} ppm at the median uncertainty {np.median(u):.2f} ppm")
