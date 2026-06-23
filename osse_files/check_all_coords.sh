#!/bin/bash

# Check coordinate names in all OSSE emission files that will be read
# for simulation starting 2016-01-01

DEST_BASE="/nobackup/ksuselj1/gchp_14.5.3_adjoint_surfaceF/surface_fluxes_osse"

echo "========================================"
echo "Checking coordinates in OSSE files"
echo "========================================"
echo ""

echo "=== Fossil Fuel (2016/01/01.nc) ==="
if [ -f "${DEST_BASE}/Fossilfuel/FF_regrid2/2016/01/01.nc" ]; then
    ncdump -h "${DEST_BASE}/Fossilfuel/FF_regrid2/2016/01/01.nc" | grep -E "^\s*(float|double) (lon|lat|longitude|latitude)\(" | head -5
else
    echo "   FILE NOT FOUND - need to run fix_co2_fluxes.sh"
fi

echo ""
echo "=== Ocean ECCO-Darwin (2016/01.nc) ==="
if [ -f "${DEST_BASE}/Ocean/ECCO-Darwin-MON-v05/2016/01.nc" ]; then
    ncdump -h "${DEST_BASE}/Ocean/ECCO-Darwin-MON-v05/2016/01.nc" | grep -E "^\s*(float|double) (lon|lat|longitude|latitude)\(" | head -5
else
    echo "   FILE NOT FOUND - need to run fix_co2_fluxes.sh"
fi

echo ""
echo "=== NBE CARDAMOM (2016/01/01.nc) ==="
if [ -f "${DEST_BASE}/Balbio/CARDAMOM-ECCO/2016/01/01.nc" ]; then
    ncdump -h "${DEST_BASE}/Balbio/CARDAMOM-ECCO/2016/01/01.nc" | grep -E "^\s*(float|double) (lon|lat|longitude|latitude)\(" | head -5
else
    echo "   FILE NOT FOUND - need to run fix_co2_fluxes.sh"
fi

echo ""
echo "=== GPP ORCHIDEE (2016/01.nc) ==="
if [ -f "${DEST_BASE}/GPP/ORCHIDEE/2016/01.nc" ]; then
    ncdump -h "${DEST_BASE}/GPP/ORCHIDEE/2016/01.nc" | grep -E "^\s*(float|double) (lon|lat|longitude|latitude)\(" | head -5
else
    echo "   FILE NOT FOUND - need to run fix_co2_fluxes.sh"
fi

echo ""
echo "=== TER CLASS-CTEM (2016/01.nc) ==="
if [ -f "${DEST_BASE}/TER/CLASS-CTEM/2016/01.nc" ]; then
    ncdump -h "${DEST_BASE}/TER/CLASS-CTEM/2016/01.nc" | grep -E "^\s*(float|double) (lon|lat|longitude|latitude)\(" | head -5
else
    echo "   FILE NOT FOUND - need to run fix_co2_fluxes.sh"
fi

echo ""
echo "========================================"
echo "Summary"
echo "========================================"
echo "All files should show:"
echo "  float lon(lon) ;"
echo "  float lat(lat) ;"
echo ""
echo "If any show 'longitude' or 'latitude', run:"
echo "  bash osse_files/fix_co2_fluxes.sh 2016"
