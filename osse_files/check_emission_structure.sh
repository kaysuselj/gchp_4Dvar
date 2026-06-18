#!/bin/bash

# Script to check the directory structure of emission data sources
# Run this on the supercomputer to verify file organization

BASE_PATH="/nobackupp17/jliu7/INVENTORY/V4"

echo "========================================"
echo "Checking Emission Data Directory Structure"
echo "========================================"
echo ""

# 1. Ocean (monthly) - ECCO-Darwin-MON-v05
echo "1. OCEAN (ECCO-Darwin-MON-v05) - Monthly files:"
echo "   Path: ${BASE_PATH}/Ocean/ECCO-Darwin-MON-v05/"
if [ -d "${BASE_PATH}/Ocean/ECCO-Darwin-MON-v05/" ]; then
    echo "   First 10 files/directories:"
    ls -1 "${BASE_PATH}/Ocean/ECCO-Darwin-MON-v05/" | head -10
    echo "   Sample file from 2016:"
    ls -1 "${BASE_PATH}/Ocean/ECCO-Darwin-MON-v05/" | grep 2016 | head -5
else
    echo "   Directory not found!"
fi
echo ""

# 2. Biosphere (hourly) - CARDAMOM-ECCO
echo "2. NBE (CARDAMOM-ECCO) - Hourly files:"
echo "   Path: ${BASE_PATH}/Balbio/CARDAMOM-ECCO/"
if [ -d "${BASE_PATH}/Balbio/CARDAMOM-ECCO/" ]; then
    echo "   First 10 files/directories:"
    ls -1 "${BASE_PATH}/Balbio/CARDAMOM-ECCO/" | head -10
    echo "   Checking for year directories (2016/):"
    if [ -d "${BASE_PATH}/Balbio/CARDAMOM-ECCO/2016" ]; then
        echo "   Found 2016/ directory, checking inside:"
        ls -1 "${BASE_PATH}/Balbio/CARDAMOM-ECCO/2016/" | head -10
        if [ -d "${BASE_PATH}/Balbio/CARDAMOM-ECCO/2016/01" ]; then
            echo "   Found 2016/01/ directory, checking inside:"
            ls -1 "${BASE_PATH}/Balbio/CARDAMOM-ECCO/2016/01/" | head -10
        fi
    else
        echo "   No year directory found, files appear to be flat"
    fi
else
    echo "   Directory not found!"
fi
echo ""

# 3. GPP (monthly) - ORCHIDEE
echo "3. GPP (ORCHIDEE) - Monthly files:"
echo "   Path: ${BASE_PATH}/GPP/ORCHIDEE/"
if [ -d "${BASE_PATH}/GPP/ORCHIDEE/" ]; then
    echo "   First 10 files/directories:"
    ls -1 "${BASE_PATH}/GPP/ORCHIDEE/" | head -10
    echo "   Sample file from 2016:"
    ls -1 "${BASE_PATH}/GPP/ORCHIDEE/" | grep 2016 | head -5
else
    echo "   Directory not found!"
fi
echo ""

# 4. TER (monthly) - CLASS-CTEM
echo "4. TER (CLASS-CTEM) - Monthly files:"
echo "   Path: ${BASE_PATH}/TER2/CLASS-CTEM/"
if [ -d "${BASE_PATH}/TER2/CLASS-CTEM/" ]; then
    echo "   First 10 files/directories:"
    ls -1 "${BASE_PATH}/TER2/CLASS-CTEM/" | head -10
    echo "   Sample file from 2016:"
    ls -1 "${BASE_PATH}/TER2/CLASS-CTEM/" | grep 2016 | head -5
else
    echo "   Directory not found!"
fi
echo ""

echo "========================================"
echo "Summary of expected patterns:"
echo "========================================"
echo "Ocean (monthly): ECCO_Darwin_YYYYMM.nc"
echo "NBE (hourly): CARDAMOM_YYYYMMDDHH.nc or YYYY/MM/DD/HH.nc or YYYY/MM/DDTHH.nc"
echo "GPP (monthly): ORCHIDEE_GPP_YYYYMM.nc"
echo "TER (monthly): CLASS_CTEM_TER_YYYYMM.nc"
echo ""
echo "Check if files follow flat naming or directory structure like Fossil Fuel (YYYY/MM/DD.nc)"
