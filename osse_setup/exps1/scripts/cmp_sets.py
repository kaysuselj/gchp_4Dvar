import sys, netCDF4 as nc, numpy as np
L = "/nobackupp17/jliu7/OCO2/L2#"
sets = {"pseudo": "OCO2-B10/PSEUDO/TRENDY/ORCHIDEE-ECCO2/OCO2b10_29Aug2020_%s",
        "B10": {"LN": "LandNadir", "LG": "LandGlint", "OG": "OceanGlint"}, "B11.2": "OCO2-B11.2/%s"}
day = "2016/07/15"
D = {}
for s in ("LN", "OG"):
    for name, pat in sets.items():
        p = pat[s] if isinstance(pat, dict) else pat % s
        d = nc.Dataset(f"{L}/{p}/{day}.nc")
        lat, lon, t = d["latitude"][:], d["longitude"][:], d["time"][:]
        u = d["xCO2-uncertainty"][:] * 1e6; x = d["xCO2"][:] * 1e6
        D[(name, s)] = (np.round(lat, 2), np.round(lon, 2), np.round(t, 2))
        title = getattr(d, "TITLE", "")[:45]
        print(f"{s} {name:6s} n={lat.size:4d}  uncertainty median {np.median(u):.2f} ppm  xCO2 median {np.median(x):.1f}  {title}")
    a = set(zip(*D[("pseudo", s)])); b = set(zip(*D[("B10", s)])); c = set(zip(*D[("B11.2", s)]))
    print(f"   {s}: pseudo soundings at the same lat/lon/time in real B10: {len(a & b)} of {len(a)}; in B11.2: {len(a & c)}")
