#!/usr/bin/env python3
"""exps1 (2026-11-09) OSSE2 background error, C24 (the prior error E files are written by setup_exp.sh / make_prior_err.py:
OSSE1 E = 1.3, OSSE2 E = 1 because its prior is a different model, set in the control set):
  sigma_b        : <exp>/sigma_b.npz (key sigma_b, the optimizer's x layout: field (TER, ocean) x control months
                   x C24 cell, stacked 6*24 x 24) = the files' Uncertainty, TER2/CLASS-CTEM for TER and MOM6-MON for the
                   ocean, 4x5 -> C24 by area-weighted sub-cell sampling (as 4dvar_realcase/scripts/make_realcase_inputs.py).
                   CLASS-CTEM Uncertainty is NaN in 31 cells where the TER flux is 0: filled with 0.30 (the file's floor);
                   these cells have no gradient, so sigma stays 1 there.
    make_exps1_inputs.py <exp dir> osse2 [<first YYYY-MM> <last YYYY-MM> [<cube res, 24> [<sigma_b scale, 1>]]]
"""
import os, sys
import numpy as np, netCDF4 as nc

sys.path.insert(0, "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_realcase/scripts")
from make_inputs import COORD, subcell_points   # noqa: E402
COORDS = {24: COORD, 48: "/nobackupp27/ksuselj1/gchp_14.8.0_adj/inputs/osse_base/coords/coord_c48.nc4"}

I = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/inputs/surface_fluxes_osse"
LAT_EDGES = np.concatenate([[-90.0], np.arange(-88.0, 89.0, 4.0), [90.0]])
exp, kind = sys.argv[1], sys.argv[2]
assert kind == "osse2", "only OSSE2 has a sigma_b file (OSSE1: constant SIGMA_B 0.2)"
# window (2026-11-10): first and last control month, default 2016-01 .. 2016-12
FIRST = tuple(int(v) for v in (sys.argv[3] if len(sys.argv) > 3 else "2016-01").split("-"))
LAST = tuple(int(v) for v in (sys.argv[4] if len(sys.argv) > 4 else "2016-12").split("-"))
# resolution and sigma_b scale (2026-11-11): C48 with sigma_b x 48/24 (J_b ~ constant across resolutions)
RES = int(sys.argv[5]) if len(sys.argv) > 5 else 24
SCALE = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0


def months(first, last):
    y, m = first
    while (y, m) <= last:
        yield y, m
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)


if kind == "osse2":
    c = nc.Dataset(COORDS[RES])
    lon, lat = subcell_points(c["corner_lons"][:], c["corner_lats"][:])
    j = np.clip(np.searchsorted(LAT_EDGES, lat, side="right") - 1, 0, 45)
    i = np.floor(((lon + 182.5) % 360.0) / 5.0).astype(int) % 72
    w = np.cos(np.radians(lat))
    to_c24 = lambda f: ((f[j, i] * w).sum(axis=(-1, -2)) / w.sum(axis=(-1, -2))).reshape(6 * RES, RES) * SCALE
    out = []
    for name, d in (("TER", "TER/CLASS-CTEM"), ("Ocean", "Ocean/MOM6-MON")):
        nfill = 0
        for y, m in months(FIRST, LAST):
            g = nc.Dataset(f"{I}/{d}/{y}/{m:02d}.nc")
            u = np.squeeze(g["Uncertainty"][:]).astype("f8"); f = np.squeeze(g["CO2_Flux"][:]).astype("f8")
            bad = ~np.isfinite(u) | (u <= 0)
            if bad.any():
                if np.abs(f[bad]).max() > 0:
                    raise SystemExit(f"{d} {y}-{m:02d}: Uncertainty missing where the flux is not 0")
                u[bad] = 0.30; nfill += int(bad.sum())
            out.append(to_c24(u))
        s = np.array(out[-len(list(months(FIRST, LAST))):])
        print(f"sigma_b {name} ({d}): min {s.min():.3f} median {np.median(s):.3f} max {s.max():.3f}; {nfill} 4x5 cell-months filled with 0.30")
    sb = np.array(out).ravel()
    np.savez(os.path.join(exp, "sigma_b.npz"), sigma_b=sb, fields=np.array(["TER", "Ocean"]),
             first=np.str_("%d-%02d" % FIRST), last=np.str_("%d-%02d" % LAST),
             source=np.str_(f"Uncertainty of TER2/CLASS-CTEM (TER) and MOM6-MON (ocean), 4x5 -> C{RES} sub-cell mean, x {SCALE}"))
    print(f"sigma_b: {sb.size} elements -> {exp}/sigma_b.npz")
