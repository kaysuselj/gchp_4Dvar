#!/usr/bin/env python3
"""exps1 (2026-11-10): J. Liu V4 daily 4x5 flux files -> inputs/surface_fluxes_osse/<flux>/YYYY/MM/DD.nc in the layout of the
existing converted files:
  fossil   (Fossilfuel/FF_regrid2): 24 hourly records, values unchanged; time units "hours since <day> 00:00:00"
  balbio   (Balbio/CARDAMOM-ECCO) : J. Liu's 8 three-hourly records repeated x3 to 24 hourly records (as the existing
                                    files: each hour keeps the value of its 3-hour block); time 0.5 .. 23.5 h
dims time / lat / lon, vars time (f8), lat / lon (f4), CO2_Flux (f4). Checks the result against the source.
    convert_jliu_daily.py fossil|balbio <src root> <dst root> <first YYYY-MM> <last YYYY-MM>
"""
import calendar, os, sys
import numpy as np, netCDF4 as nc

kind, src, dst, first, last = sys.argv[1:6]
assert kind in ("fossil", "balbio")
y0, m0 = map(int, first.split("-")); y1, m1 = map(int, last.split("-"))
y, m = y0, m0
while (y, m) <= (y1, m1):
    os.makedirs(f"{dst}/{y}/{m:02d}", exist_ok=True)
    for d in range(1, calendar.monthrange(y, m)[1] + 1):
        fi, fo = f"{src}/{y}/{m:02d}/{d:02d}.nc", f"{dst}/{y}/{m:02d}/{d:02d}.nc"
        a = nc.Dataset(fi)
        latn = "latitude" if "latitude" in a.variables else "lat"
        lonn = "longitude" if "longitude" in a.variables else "lon"
        lat, lon = np.asarray(a[latn][:], "f4"), np.asarray(a[lonn][:], "f4")
        f = np.asarray(a["CO2_Flux"][:], "f4")
        if kind == "balbio":
            assert f.shape == (8, 46, 72), (fi, f.shape)
            f = np.repeat(f, 3, axis=0)
        assert f.shape == (24, 46, 72), (fi, f.shape)
        tmp = fo + ".tmp"
        with nc.Dataset(tmp, "w", format="NETCDF4") as o:
            o.createDimension("time", 24); o.createDimension("lat", 46); o.createDimension("lon", 72)
            t = o.createVariable("time", "f8", ("time",))
            t.units = f"hours since {y}-{m:02d}-{d:02d} 00:00:00"; t.long_name = "time"; t[:] = np.arange(24) + 0.5
            v = o.createVariable("lat", "f4", ("lat",)); v.units = "degrees_north"; v[:] = lat
            v = o.createVariable("lon", "f4", ("lon",)); v.units = "degrees_east"; v[:] = lon
            v = o.createVariable("CO2_Flux", "f4", ("time", "lat", "lon")); v.units = "Kg C/Km^2/sec"; v[:] = f
        b = nc.Dataset(tmp)
        src_f = np.asarray(a["CO2_Flux"][:], "f4")
        got = np.asarray(b["CO2_Flux"][:])
        assert np.array_equal(got[::3] if kind == "balbio" else got, src_f, equal_nan=True), fi
        b.close(); os.replace(tmp, fo)
    print(f"{dst}/{y}/{m:02d}: {calendar.monthrange(y, m)[1]} days written ({kind})")
    y, m = (y + 1, 1) if m == 12 else (y, m + 1)
