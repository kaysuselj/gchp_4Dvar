#!/bin/bash

# Check GPP and TER file structure inside year directories

BASE_PATH="/nobackupp17/jliu7/INVENTORY/V4"

echo "========================================"
echo "GPP and TER File Structure"
echo "========================================"
echo ""

# Check GPP/ORCHIDEE/2016/
echo "1. GPP (ORCHIDEE) - What's inside 2016/:"
echo "   Path: ${BASE_PATH}/GPP/ORCHIDEE/2016/"
if [ -d "${BASE_PATH}/GPP/ORCHIDEE/2016" ]; then
    ls -1 "${BASE_PATH}/GPP/ORCHIDEE/2016/" | head -15
else
    echo "   2016 directory not found!"
fi
echo ""

# Check TER2/CLASS-CTEM/2016/
echo "2. TER (CLASS-CTEM) - What's inside 2016/:"
echo "   Path: ${BASE_PATH}/TER2/CLASS-CTEM/2016/"
if [ -d "${BASE_PATH}/TER2/CLASS-CTEM/2016" ]; then
    ls -1 "${BASE_PATH}/TER2/CLASS-CTEM/2016/" | head -15
else
    echo "   2016 directory not found!"
fi
echo ""

echo "========================================"
echo "Summary:"
echo "========================================"
echo "Fossil: YYYY/MM/DD.nc (daily files)"
echo "Ocean: YYYY/MM.nc (monthly files 01-12)"
echo "NBE: YYYY/MM/DD.nc (daily files 01-31)"
echo "GPP: YYYY/?.nc (checking...)"
echo "TER: YYYY/?.nc (checking...)"
