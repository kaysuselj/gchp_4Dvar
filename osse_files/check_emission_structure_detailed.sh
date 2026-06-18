#!/bin/bash

# Detailed check for remaining questions

BASE_PATH="/nobackupp17/jliu7/INVENTORY/V4"

echo "========================================"
echo "Detailed Emission Structure Check"
echo "========================================"
echo ""

# Check Ocean 2016 directory contents
echo "1. OCEAN - What's inside 2016/:"
echo "   Path: ${BASE_PATH}/Ocean/ECCO-Darwin-MON-v05/2016/"
if [ -d "${BASE_PATH}/Ocean/ECCO-Darwin-MON-v05/2016" ]; then
    ls -1 "${BASE_PATH}/Ocean/ECCO-Darwin-MON-v05/2016/" | head -15
else
    echo "   2016 directory not found!"
fi
echo ""

# Check if NBE files are daily or hourly
echo "2. NBE - Checking file pattern in 2016/01/:"
echo "   Are files DD.nc (daily) or DDHH.nc (hourly)?"
if [ -d "${BASE_PATH}/Balbio/CARDAMOM-ECCO/2016/01" ]; then
    ls -1 "${BASE_PATH}/Balbio/CARDAMOM-ECCO/2016/01/" | head -5
    echo "   Checking 2016/01/01/ for hourly files:"
    if [ -d "${BASE_PATH}/Balbio/CARDAMOM-ECCO/2016/01/01" ]; then
        ls -1 "${BASE_PATH}/Balbio/CARDAMOM-ECCO/2016/01/01/" | head -5
    else
        echo "   No 01/ subdirectory, files are daily"
    fi
fi
echo ""

# Search for GPP
echo "3. Searching for GPP directories:"
echo "   Checking alternative paths:"
find /nobackupp17/jliu7/INVENTORY/V4/ -maxdepth 2 -type d -iname "*gpp*" 2>/dev/null
echo ""

# Search for TER
echo "4. Searching for TER directories:"
echo "   Checking alternative paths:"
find /nobackupp17/jliu7/INVENTORY/V4/ -maxdepth 2 -type d -iname "*ter*" 2>/dev/null
echo ""

# Check what's actually in V4/
echo "5. All directories in V4/:"
ls -1 /nobackupp17/jliu7/INVENTORY/V4/
