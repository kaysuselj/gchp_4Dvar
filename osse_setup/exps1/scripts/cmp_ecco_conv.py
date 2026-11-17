import numpy as np, netCDF4 as nc, sys
for m in range(1, 13):
    a = nc.Dataset("/tmp/ecco_check/2016/%02d.nc" % m); b = nc.Dataset("/nobackupp27/ksuselj1/gchp_14.8.0_adj/inputs/surface_fluxes_osse/Ocean/ECCO-Darwin-MON-v05/2016/%02d.nc" % m)
    for k in ("CO2_Flux", "Uncertainty", "lat", "lon", "time"):
        assert np.array_equal(np.asarray(a[k][:]), np.asarray(b[k][:]), equal_nan=True), (m, k)
print("converter check: ECCO-Darwin 2016 identical to the existing files (CO2_Flux, Uncertainty, lat, lon, time)")
