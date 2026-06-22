#!/bin/bash

# Script to check all CO2 emission files for required attributes
# Verifies that variables and coordinates have proper 'units' attributes
# Usage: bash check_all_emission_files.sh YEAR [MONTH] [DAY]
# Examples:
#   bash check_all_emission_files.sh 2016         # Check all first available files
#   bash check_all_emission_files.sh 2016 01 01   # Check specific date

YEAR=${1:-2016}
MONTH=${2:-01}
DAY=${3:-01}

echo "========================================"
echo "Checking CO2 Emission Files for Year ${YEAR}"
echo "========================================"
echo ""

check_file() {
    local label=$1
    local file=$2
    local check_time=$3  # "yes" or "no"

    echo "=== ${label} ==="
    echo "File: ${file}"

    if [ ! -f "${file}" ]; then
        echo "❌ FILE NOT FOUND"
        echo ""
        return 1
    fi

    # Check dimensions
    echo "Dimensions:"
    ncdump -h "${file}" | grep "dimensions:" -A10 | grep "=" | head -5

    # Check CO2_Flux variable and units
    echo ""
    echo "CO2_Flux variable:"
    if ncdump -h "${file}" | grep -q "CO2_Flux"; then
        ncdump -h "${file}" | grep -A5 "CO2_Flux(" | head -6

        # Check for correct 'units' attribute
        if ncdump -h "${file}" | grep "CO2_Flux:" | grep -q "units ="; then
            echo "✓ CO2_Flux has 'units' attribute"
        elif ncdump -h "${file}" | grep "CO2_Flux:" | grep -q "unit ="; then
            echo "❌ CO2_Flux has 'unit' (singular) - needs to be 'units'"
        else
            echo "❌ CO2_Flux missing 'units' attribute"
        fi
    else
        echo "❌ CO2_Flux variable not found!"
    fi

    # Check time coordinate if required
    if [ "${check_time}" == "yes" ]; then
        echo ""
        echo "Time coordinate:"
        if ncdump -h "${file}" | grep -q "time("; then
            ncdump -h "${file}" | grep -A5 "time(" | head -6

            # Check for correct 'units' attribute with reference time
            if ncdump -h "${file}" | grep "time:" | grep -q 'units = "hours since'; then
                echo "✓ time has proper 'units' with reference time"
            elif ncdump -h "${file}" | grep "time:" | grep -q 'units = "hour"'; then
                echo "❌ time has 'units = \"hour\"' - needs reference time like 'hours since YYYY-MM-DD'"
            elif ncdump -h "${file}" | grep "time:" | grep -q 'unit ='; then
                echo "❌ time has 'unit' (singular) - needs 'units' with reference time"
            else
                echo "❌ time missing proper 'units' attribute"
            fi
        else
            echo "⚠️  No time dimension (file is constant/climatology)"
        fi
    fi

    echo ""
    echo "---"
    echo ""
}

# Check corrected files in surface_fluxes_osse
CORRECTED_BASE="/nobackup/ksuselj1/gchp_14.5.3_adjoint_surfaceF/surface_fluxes_osse"

# Check original files
ORIGINAL_BASE="/nobackupp17/jliu7/INVENTORY/V4"

echo "CHECKING CORRECTED FILES"
echo "========================================"
echo ""

# 1. Fossil Fuel (daily, with time)
check_file "1. Fossil Fuel (corrected)" \
    "${CORRECTED_BASE}/Fossilfuel/FF_regrid2/${YEAR}/${MONTH}/${DAY}.nc" \
    "yes"

# 2. Ocean (monthly, no time)
check_file "2. Ocean (corrected)" \
    "${CORRECTED_BASE}/Ocean/ECCO-Darwin-MON-v05/${YEAR}/${MONTH}.nc" \
    "no"

# 3. NBE (daily, with time)
check_file "3. NBE (corrected)" \
    "${CORRECTED_BASE}/Balbio/CARDAMOM-ECCO/${YEAR}/${MONTH}/${DAY}.nc" \
    "yes"

echo ""
echo "CHECKING ORIGINAL FILES (should already be correct)"
echo "========================================"
echo ""

# 4. GPP (monthly, no time)
check_file "4. GPP (original - should be OK)" \
    "${ORIGINAL_BASE}/GPP/ORCHIDEE/${YEAR}/${MONTH}.nc" \
    "no"

# 5. TER (monthly, no time)
check_file "5. TER (original - should be OK)" \
    "${ORIGINAL_BASE}/TER2/CLASS-CTEM/${YEAR}/${MONTH}.nc" \
    "no"

echo ""
echo "========================================"
echo "SUMMARY"
echo "========================================"
echo ""
echo "Files checked for ${YEAR}-${MONTH}-${DAY}"
echo ""
echo "If any file shows ❌, run fix_emission_units.sh to correct it"
echo ""
echo "Expected results:"
echo "  ✓ All variables have 'units' (not 'unit')"
echo "  ✓ Time coordinates (if present) have 'units = \"hours since YYYY-MM-DD\"'"
echo "  ⚠️  Files without time dimension are OK (constant/climatology)"
