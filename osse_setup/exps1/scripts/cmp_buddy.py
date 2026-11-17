import os, numpy as np, netCDF4 as nc
L = "/nobackupp17/jliu7/OCO2/L2#/OCO2-B10"
for s in ("LN", "LG", "OG"):
    for day in ("2016/07/15", "2017/03/15"):
        p, r = f"{L}/PSEUDO/TRENDY/ORCHIDEE-ECCO2/OCO2b10_29Aug2020_{s}/{day}.nc", f"{L}/OCO2b10_29Aug2020_{s}/{day}.nc"
        if not os.path.exists(r):
            print(s, day, "no real 50 km buddy file"); continue
        a, b = nc.Dataset(p), nc.Dataset(r)
        same = {k: bool(np.array_equal(np.asarray(a[k][:]), np.asarray(b[k][:])))
                for k in ("latitude", "longitude", "time", "pressure", "xCO2-averagingKernel", "xCO2-uncertainty", "CO2-apriori", "xCO2-apriori")}
        dx = (np.asarray(a["xCO2"][:]) - np.asarray(b["xCO2"][:])) * 1e6
        print(f"{s} {day}: n pseudo {a['latitude'].size}, real buddy {b['latitude'].size}; identical: {same}; xCO2 pseudo-real mean {dx.mean():+.2f} std {dx.std():.2f} ppm")
for s in ("LN", "LG", "OG"):
    d = f"{L}/OCO2b10_29Aug2020_{s}"
    print(s, "real buddy years:", sorted(os.listdir(d)) if os.path.isdir(d) else "missing",
          "| pseudo 2017 months:", sorted(os.listdir(f"{L}/PSEUDO/TRENDY/ORCHIDEE-ECCO2/OCO2b10_29Aug2020_{s}/2017")))
