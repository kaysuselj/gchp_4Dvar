import sys, numpy as np, netCDF4 as nc
dst, ref = sys.argv[1], sys.argv[2]
r = nc.Dataset(ref)
with nc.Dataset(dst, "a") as d:
    for k in sorted(set(r.variables) - set(d.variables)):
        v = r[k]
        ref_vals = np.asarray(v[:])
        dims = tuple(x for x in v.dimensions)
        for x in dims:
            assert x in d.dimensions, (k, x)
        w = d.createVariable(k, v.dtype, dims)
        for a in v.ncattrs():
            if a != "_FillValue":
                w.setncattr(a, v.getncattr(a))
        w[:] = ref_vals
        print(f"added {k} {dims} from the 2016 restart (max |value| {np.abs(ref_vals).max():.3g})")
a = nc.Dataset(dst)
print("variables now missing:", sorted(set(r.variables) - set(a.variables)),
      "| extra:", sorted(set(a.variables) - set(r.variables)),
      "| DELP_DRY max", float(np.asarray(a["DELP_DRY"][:]).max()), "| time", a["time"].units, a["time"][:])
