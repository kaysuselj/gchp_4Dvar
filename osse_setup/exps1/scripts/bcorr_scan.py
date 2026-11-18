"""(1) J_b of the TRUE correction (SF_true - 1) under diagonal and correlated B, for a range of length scales and tau:
how much cheaper does correlated B make the truth, field by field (OSSE 1 and OSSE 2). (2) Empirical spatial and temporal
correlation of the OSSE 2 prior error (SF_true - 1) per field (sanity check only, not used to tune B)."""
import sys, numpy as np, netCDF4 as nc
sys.path.insert(0, "/nobackupp27/ksuselj1/gchp_14.8.0_adj/gchp_4Dvar")
from b_corr import unit_vectors, great_circle_km, corr_time
X = "/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse/exps1"
t = nc.Dataset(f"{X}/obs_c24_1510_1705/truth_emissions_monthly_201510_201705.nc")
lon, lat = np.asarray(t["lons"][:], float).ravel(), np.asarray(t["lats"][:], float).ravel()
u = unit_vectors(lon, lat)
FLD = (("TER", "EmisCO2_NetTerrExch", 0, 4), ("ocean", "EmisCO2_Ocean", -4, 4))


def true_delta(osse, k):
    name, v, lo, hi = FLD[k]
    T = np.asarray(t[v][:], float).reshape(19, -1)
    P = np.asarray(nc.Dataset(f"{X}/{osse}_C24_DL/4dvar_output_osse/emissions_monthly_iter001.nc4")[v][:], float).reshape(19, -1)
    m = np.abs(P).max(0) > 0
    x = np.ones_like(P); nz = np.abs(P) > 0
    x[nz] = T[nz] / P[nz]
    x = np.clip(x, lo, hi)
    if osse == "OSSE2":
        sb = np.load(f"{X}/OSSE2_C24_DL/sigma_b.npz")["sigma_b"].reshape(2, 19, -1)[k]
    else:
        sb = np.full_like(x, 0.2)
    return ((x - 1) / sb)[:, m], m                         # (19, n_mask), in sigma_b units


def jb(Dl, Lkm, tau, d):
    """0.5 * vec(D)^T (C_t (x) C_h)^-1 vec(D), D (months, cells), via eigendecompositions (relative cutoff 1e-12)."""
    if Lkm <= 0 and tau <= 0:
        return 0.5 * np.sum(Dl ** 2), 1.0
    lh, Uh = np.linalg.eigh(np.exp(-d / Lkm) if Lkm > 0 else np.eye(d.shape[0]))
    lt, Ut = np.linalg.eigh(corr_time(Dl.shape[0], tau))
    A = Ut.T @ Dl @ Uh                                      # coefficients in the eigenbasis
    lam = lt[:, None] * lh[None, :]
    ok = lam > 1e-12 * lam.max()
    return 0.5 * np.sum(A[ok] ** 2 / lam[ok]), lh.max() / max(lh.min(), 1e-300)


print("(1) J_b of the true correction [diagonal B = L 0, tau 0]; in brackets: ratio to diagonal")
for osse in ("OSSE1", "OSSE2"):
    for k, (name, *_ ) in enumerate(FLD):
        Dl, m = true_delta(osse, k)
        d = great_circle_km(u[m])
        base, _ = jb(Dl, 0, 0, d)
        Ls = (250, 500, 1000, 2000, 4000) if k == 0 else (500, 1000, 2000, 4000, 8000)
        out = []
        for tau in (0.0, 1.0, 3.0):
            row = []
            for L in Ls:
                v, cond = jb(Dl, L, tau, d)
                row.append(f"{L}:{v:.3g}({v / base:.2f})")
            out.append(f"      tau {tau:g}: " + "  ".join(row))
        print(f"  {osse} {name:5s} diagonal {base:.4g} (n = {Dl.size} cell-months)")
        print("\n".join(out))

print("\n(2) OSSE 2 prior error (SF_true - 1) correlation vs distance (pairs within the same month, mean over months)")
bins = np.array([0, 1, 300, 600, 900, 1200, 1800, 2500, 3500, 5000])
for k, (name, *_ ) in enumerate(FLD):
    Dl, m = true_delta("OSSE2", k)
    Dl = Dl - Dl.mean(1, keepdims=True)                     # remove the monthly global mean (only the pattern)
    d = great_circle_km(u[m])
    corr = []
    for a, b in zip(bins[:-1], bins[1:]):
        sel = (d >= a) & (d < b)
        num = np.mean([np.sum(np.outer(x, x)[sel]) for x in Dl]); den = np.mean([np.sum(x ** 2) * sel.sum() / x.size for x in Dl])
        corr.append(num / den)
    lagc = [np.mean([np.corrcoef(Dl[i], Dl[i + l])[0, 1] for i in range(19 - l)]) for l in (1, 2, 3, 6, 12)]
    print(f"  {name:5s} distance bins km " + " ".join(f"{a}-{b}" for a, b in zip(bins[:-1], bins[1:])))
    print(f"        correlation     " + " ".join(f"{c:6.2f}" for c in corr))
    print(f"        e-folding ~ first bin with corr < 0.37: "
          f"{next((f'{a}-{b} km' for (a, b), c in zip(zip(bins[:-1], bins[1:]), corr) if c < 0.37), '> 5000 km')}")
    print(f"        month-lag correlation (1, 2, 3, 6, 12): " + " ".join(f"{c:.2f}" for c in lagc))
