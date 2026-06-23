#!/bin/bash
# Check coordinate names in all emission files

echo "=== Checking OCEAN_EXCH_SCALED ==="
ncdump -h /nobackupp17/lwu5/GCHP/ExtData/HEMCO/CO2/v2022-11/OCEAN/Scaled_Ocean_CO2_monthly.nc | grep -E "dimensions:|lon|lat" | head -20

echo ""
echo "=== Checking GPP_ORCHIDEE (2016/01.nc) ==="
ncdump -h /nobackupp17/jliu7/INVENTORY/V4/GPP/ORCHIDEE/2016/01.nc | grep -E "dimensions:|lon|lat" | head -20

echo ""
echo "=== Checking TER_CLASS_CTEM (2016/01.nc) ==="
ncdump -h /nobackupp17/jliu7/INVENTORY/V4/TER2/CLASS-CTEM/2016/01.nc | grep -E "dimensions:|lon|lat" | head -20
