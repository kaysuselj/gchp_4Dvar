#!/usr/bin/env python3
"""
Compare two satellite track files to diagnose structural differences.

Usage:
    python compare_sat_track_files.py old_working_file.nc new_file.nc
"""

import sys
import xarray as xr

def compare_files(file1, file2):
    print(f"Comparing:")
    print(f"  File 1 (OLD): {file1}")
    print(f"  File 2 (NEW): {file2}")
    print("=" * 80)

    ds1 = xr.open_dataset(file1)
    ds2 = xr.open_dataset(file2)

    # Compare global attributes
    print("\n### GLOBAL ATTRIBUTES ###")
    print(f"File 1: {dict(ds1.attrs)}")
    print(f"File 2: {dict(ds2.attrs)}")

    # Compare dimensions
    print("\n### DIMENSIONS ###")
    print(f"File 1: {dict(ds1.dims)}")
    print(f"File 2: {dict(ds2.dims)}")

    # Compare variables
    print("\n### VARIABLES ###")
    print(f"File 1: {list(ds1.data_vars)}")
    print(f"File 2: {list(ds2.data_vars)}")

    # Compare coordinates
    print("\n### COORDINATES ###")
    print(f"File 1: {list(ds1.coords)}")
    print(f"File 2: {list(ds2.coords)}")

    # Compare each variable in detail
    all_vars = set(list(ds1.data_vars) + list(ds1.coords) +
                   list(ds2.data_vars) + list(ds2.coords))

    for var in sorted(all_vars):
        print(f"\n### VARIABLE: {var} ###")

        in_file1 = var in ds1
        in_file2 = var in ds2

        print(f"  In File 1: {in_file1}")
        print(f"  In File 2: {in_file2}")

        if in_file1 and in_file2:
            v1 = ds1[var]
            v2 = ds2[var]

            print(f"  Dimensions:")
            print(f"    File 1: {v1.dims}")
            print(f"    File 2: {v2.dims}")

            print(f"  Shape:")
            print(f"    File 1: {v1.shape}")
            print(f"    File 2: {v2.shape}")

            print(f"  Dtype:")
            print(f"    File 1: {v1.dtype}")
            print(f"    File 2: {v2.dtype}")

            print(f"  Attributes:")
            print(f"    File 1: {dict(v1.attrs)}")
            print(f"    File 2: {dict(v2.attrs)}")

            # Check if values are similar
            if v1.shape == v2.shape:
                import numpy as np
                try:
                    if np.allclose(v1.values, v2.values, rtol=1e-5, atol=1e-8, equal_nan=True):
                        print(f"  Values: IDENTICAL (within tolerance)")
                    else:
                        diff = np.abs(v1.values - v2.values)
                        print(f"  Values: DIFFERENT (max diff: {np.nanmax(diff)})")
                except:
                    print(f"  Values: Could not compare")

    ds1.close()
    ds2.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_sat_track_files.py old_file.nc new_file.nc")
        sys.exit(1)

    compare_files(sys.argv[1], sys.argv[2])
