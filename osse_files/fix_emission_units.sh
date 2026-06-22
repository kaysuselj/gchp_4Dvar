#!/bin/bash

# Script to fix 'unit' -> 'units' attribute in CO2 emission files
# Also fixes time coordinate units to be CF-compliant
# Usage: bash fix_emission_units.sh YEAR
# Example: bash fix_emission_units.sh 2016

if [ $# -ne 1 ]; then
    echo "Usage: $0 YEAR"
    echo "Example: $0 2016"
    exit 1
fi

YEAR=$1

# Source and destination directories
SRC_BASE="/nobackupp17/jliu7/INVENTORY/V4"
DEST_BASE="/nobackup/ksuselj1/gchp_14.5.3_adjoint_surfaceF/surface_fluxes_osse"

echo "========================================"
echo "Fixing emission units attribute for year ${YEAR}"
echo "========================================"
echo ""

# Check if ncatted is available
if ! command -v ncatted &> /dev/null; then
    echo "ERROR: ncatted (NCO tools) not found. Load the NCO module:"
    echo "  module load nco"
    exit 1
fi

# -----------------------------------------------------------------------------
# 1. Fix Fossil Fuel files (daily: YYYY/MM/DD.nc)
# -----------------------------------------------------------------------------
echo "1. Processing Fossil Fuel emissions (daily files)..."
SRC_FF="${SRC_BASE}/Fossilfuel/FF_regrid2/${YEAR}"
DEST_FF="${DEST_BASE}/Fossilfuel/FF_regrid2/${YEAR}"

if [ ! -d "${SRC_FF}" ]; then
    echo "   ERROR: Source directory not found: ${SRC_FF}"
else
    # Create destination directory structure
    mkdir -p "${DEST_FF}"

    # Process all months
    for MONTH in $(seq -f "%02g" 1 12); do
        if [ ! -d "${SRC_FF}/${MONTH}" ]; then
            echo "   Skipping month ${MONTH} (directory not found)"
            continue
        fi

        mkdir -p "${DEST_FF}/${MONTH}"

        # Process all days in this month
        for DAY_FILE in ${SRC_FF}/${MONTH}/*.nc; do
            if [ ! -f "${DAY_FILE}" ]; then
                continue
            fi

            DAY=$(basename "${DAY_FILE}")
            DEST_FILE="${DEST_FF}/${MONTH}/${DAY}"

            # Copy file first
            cp "${DAY_FILE}" "${DEST_FILE}"

            # Fix CO2_Flux attribute: delete 'unit', add 'units'
            ncatted -O -a unit,CO2_Flux,d,, -a units,CO2_Flux,c,c,"Kg C/Km^2/sec" "${DEST_FILE}"

            # Fix time attribute: delete 'unit', add proper CF-compliant 'units' with reference time
            ncatted -O -a unit,time,d,, -a units,time,c,c,"hours since 2000-01-01 00:00:00" "${DEST_FILE}"

            echo "   Fixed: ${YEAR}/${MONTH}/${DAY}"
        done
    done

    FF_COUNT=$(find "${DEST_FF}" -name "*.nc" 2>/dev/null | wc -l)
    echo "   Completed: ${FF_COUNT} Fossil Fuel files processed"
fi

echo ""

# -----------------------------------------------------------------------------
# 2. Fix Ocean files (monthly: YYYY/MM.nc)
# -----------------------------------------------------------------------------
echo "2. Processing Ocean emissions (monthly files)..."
SRC_OCEAN="${SRC_BASE}/Ocean/ECCO-Darwin-MON-v05/${YEAR}"
DEST_OCEAN="${DEST_BASE}/Ocean/ECCO-Darwin-MON-v05/${YEAR}"

if [ ! -d "${SRC_OCEAN}" ]; then
    echo "   ERROR: Source directory not found: ${SRC_OCEAN}"
else
    # Create destination directory
    mkdir -p "${DEST_OCEAN}"

    # Process all monthly files
    for MONTH_FILE in ${SRC_OCEAN}/*.nc; do
        if [ ! -f "${MONTH_FILE}" ]; then
            echo "   No monthly files found in ${SRC_OCEAN}"
            break
        fi

        MONTH=$(basename "${MONTH_FILE}")
        DEST_FILE="${DEST_OCEAN}/${MONTH}"

        # Copy file first
        cp "${MONTH_FILE}" "${DEST_FILE}"

        # Fix attribute: delete 'unit', add 'units'
        ncatted -O -a unit,CO2_Flux,d,, -a units,CO2_Flux,c,c,"Kg C/km2/sec" "${DEST_FILE}"

        echo "   Fixed: ${YEAR}/${MONTH}"
    done

    OCEAN_COUNT=$(find "${DEST_OCEAN}" -name "*.nc" 2>/dev/null | wc -l)
    echo "   Completed: ${OCEAN_COUNT} Ocean files processed"
fi

echo ""

# -----------------------------------------------------------------------------
# 3. Fix NBE files (daily: YYYY/MM/DD.nc)
# -----------------------------------------------------------------------------
echo "3. Processing NBE emissions (daily files)..."
SRC_NBE="${SRC_BASE}/Balbio/CARDAMOM-ECCO/${YEAR}"
DEST_NBE="${DEST_BASE}/Balbio/CARDAMOM-ECCO/${YEAR}"

if [ ! -d "${SRC_NBE}" ]; then
    echo "   ERROR: Source directory not found: ${SRC_NBE}"
else
    # Create destination directory structure
    mkdir -p "${DEST_NBE}"

    # Process all months
    for MONTH in $(seq -f "%02g" 1 12); do
        if [ ! -d "${SRC_NBE}/${MONTH}" ]; then
            echo "   Skipping month ${MONTH} (directory not found)"
            continue
        fi

        mkdir -p "${DEST_NBE}/${MONTH}"

        # Process all days in this month
        for DAY_FILE in ${SRC_NBE}/${MONTH}/*.nc; do
            if [ ! -f "${DAY_FILE}" ]; then
                continue
            fi

            DAY=$(basename "${DAY_FILE}")
            DEST_FILE="${DEST_NBE}/${MONTH}/${DAY}"

            # Copy file first
            cp "${DAY_FILE}" "${DEST_FILE}"

            # Fix time units: change "hour" to proper CF format
            # NBE files already have units (not unit), just need proper reference time
            REF_DATE="${YEAR}-${MONTH}-01 00:00:00"
            ncatted -O -a units,time,o,c,"hours since ${REF_DATE}" "${DEST_FILE}"

            echo "   Fixed: ${YEAR}/${MONTH}/${DAY}"
        done
    done

    NBE_COUNT=$(find "${DEST_NBE}" -name "*.nc" 2>/dev/null | wc -l)
    echo "   Completed: ${NBE_COUNT} NBE files processed"
fi

echo ""

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "========================================"
echo "Summary"
echo "========================================"
echo "Output directory: ${DEST_BASE}"
echo ""
echo "Fixed emissions:"
echo "  - Fossil Fuel: ${DEST_BASE}/Fossilfuel/FF_regrid2/${YEAR}/"
echo "     * CO2_Flux: unit -> units"
echo "     * time: unit -> units with reference time"
echo "  - Ocean:       ${DEST_BASE}/Ocean/ECCO-Darwin-MON-v05/${YEAR}/"
echo "     * CO2_Flux: unit -> units"
echo "  - NBE:         ${DEST_BASE}/Balbio/CARDAMOM-ECCO/${YEAR}/"
echo "     * time: added proper reference time"
echo ""
echo "GPP, TER already have correct attributes - use original paths"
echo ""
echo "Next steps:"
echo "1. Update HEMCO_Config.rc and ExtData.rc to point to ${DEST_BASE}"
echo "2. Verify one file from each type:"
echo "   ncdump -h ${DEST_BASE}/Fossilfuel/FF_regrid2/${YEAR}/01/01.nc | grep -E 'time:units|CO2_Flux:units'"
echo "   ncdump -h ${DEST_BASE}/Ocean/ECCO-Darwin-MON-v05/${YEAR}/01.nc | grep 'CO2_Flux:units'"
echo "   ncdump -h ${DEST_BASE}/Balbio/CARDAMOM-ECCO/${YEAR}/01/01.nc | grep 'time:units'"
