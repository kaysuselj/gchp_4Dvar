#!/usr/bin/env python3
"""exps1 (2026-11-09): copy J. Liu V4 monthly 4x5 flux files into gchp_14.8.0_adj/inputs/surface_fluxes_osse in the
layout of the existing converted files (e.g. TER/CLASS-CTEM/2016/01.nc): dims (time=1 unlimited, lat 46, lon 72), vars
lon/lat (float32, the source's values), time (int32, hours since 1900-01-01, proleptic_gregorian, first of the month
00:00), CO2_Flux and Uncertainty (float32, time x lat x lon). Values are copied unchanged; the script checks this.
    convert_jliu_monthly.py <src root (V4/...)> <dst root> <year> [<year> ...]
"""
import os, sys
from datetime import datetime
import numpy as np, netCDF4 as nc

src, dst, years = sys.argv[1], sys.argv[2], [int(y) for y in sys.argv[3:]]
for y in years:
    os.makedirs(f"{dst}/{y}", exist_ok=True)
    for m in range(1, 13):
        fi, fo = f"{src}/{y}/{m:02d}.nc", f"{dst}/{y}/{m:02d}.nc"
        a = nc.Dataset(fi)
        latn = "latitude" if "latitude" in a.variables else "lat"
        lonn = "longitude" if "longitude" in a.variables else "lon"
        lat, lon = a[latn][:].astype(np.float32), a[lonn][:].astype(np.float32)
        f = np.squeeze(a["CO2_Flux"][:]).astype(np.float32)
        u = np.squeeze(a["Uncertainty"][:]).astype(np.float32)
        assert f.shape == (46, 72) and u.shape == (46, 72), (fi, f.shape, u.shape)
        tmp = fo + ".tmp"
        with nc.Dataset(tmp, "w", format="NETCDF4") as o:
            o.createDimension("time", None); o.createDimension("lat", 46); o.createDimension("lon", 72)
            t = o.createVariable("time", "i4", ("time",))
            t.units = "hours since 1900-01-01 00:00:00"; t.calendar = "proleptic_gregorian"; t.long_name = "time"
            t[:] = int((datetime(y, m, 1) - datetime(1900, 1, 1)).total_seconds() // 3600)
            v = o.createVariable("lon", "f4", ("lon",)); v.units = "degrees east"; v[:] = lon
            v = o.createVariable("lat", "f4", ("lat",)); v.units = "degrees north"; v[:] = lat
            v = o.createVariable("CO2_Flux", "f4", ("time", "lat", "lon")); v.units = "Kg C/Km^2/sec"; v[0] = f
            v = o.createVariable("Uncertainty", "f4", ("time", "lat", "lon")); v.units = "factor"; v[0] = u
            for k in a.ncattrs():
                o.setncattr(k, a.getncattr(k))
            o.setncattr("history", f"{datetime.now():%F %T}: exps1 convert_jliu_monthly.py from {fi} (values unchanged)")
        b = nc.Dataset(tmp)
        assert np.array_equal(b["CO2_Flux"][0], f, equal_nan=True) and np.array_equal(b["Uncertainty"][0], u, equal_nan=True), fi
        # 4x5 GEOS-Chem grid: interior rows -86..86 every 4 deg; the polar rows are +-89 (land files) or +-90 (ocean
        # files, also in our ECCO-Darwin copy), kept as in the source
        la = b["lat"][:]
        assert np.allclose(la[1:-1], np.arange(-86, 87, 4)) and abs(abs(la[0]) - 89.5) <= 0.5 and np.allclose(b["lon"][:], np.arange(-180, 180, 5)), fi
        b.close(); os.replace(tmp, fo)
    print(f"{dst}/{y}: 12 months written, values identical to {src}/{y}")
