#!/bin/bash

# Comprehensive script to fix all CO2 emission file issues for GCHP OSSE
# Fixes:
#   1. Fossil Fuel: unit->units for CO2_Flux and time, proper time reference
#   2. Ocean: unit->units for CO2_Flux, ADD time coordinate
#   3. NBE: proper time reference for each day
#   4. GPP: rename longitude->lon, latitude->lat
#   5. TER: rename longitude->lon, latitude->lat
#
# Usage: bash fix_co2_fluxes.sh YEAR
# Example: bash fix_co2_fluxes.sh 2016

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 YEAR [COMPONENT]"
    echo "  COMPONENT (optional): ALL (default) | FF | OCEAN | NBE | GPP | TER"
    echo "Example: $0 2016        # all components"
    echo "Example: $0 2016 NBE    # only balanced biosphere"
    exit 1
fi

YEAR=$1

# Optional 2nd arg: which component(s) to (re)process. Default ALL.
#   ALL | FF | OCEAN | NBE | GPP | TER
# Example: bash fix_co2_fluxes.sh 2016 NBE   # only rebuild the balanced biosphere
COMPONENT="${2:-ALL}"
want() { [[ "${COMPONENT}" == "ALL" || "${COMPONENT}" == "$1" ]]; }

# Source and destination directories
SRC_BASE="/nobackupp17/jliu7/INVENTORY/V4"
DEST_BASE="/nobackup/ksuselj1/gchp_14.5.3_adjoint_surfaceF/surface_fluxes_osse"

echo "========================================"
echo "Fixing emission files for year ${YEAR}"
echo "========================================"
echo ""

# Check if NCO tools are available (NBE is now pure-python and needs no NCO)
if want FF || want OCEAN || want GPP || want TER; then
    for tool in ncatted ncap2 ncrename ncecat; do
        if ! command -v ${tool} &> /dev/null; then
            echo "ERROR: ${tool} (NCO tools) not found. Load the NCO module:"
            echo "  module load nco"
            exit 1
        fi
    done
fi

# -----------------------------------------------------------------------------
# 1. Fix Fossil Fuel files (daily: YYYY/MM/DD.nc)
# -----------------------------------------------------------------------------
if want FF; then
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

            DAY=$(basename "${DAY_FILE}" .nc)
            DEST_FILE="${DEST_FF}/${MONTH}/${DAY}.nc"

            # Copy file first
            cp "${DAY_FILE}" "${DEST_FILE}"

            # Fix CO2_Flux attribute: delete 'unit', add 'units'
            ncatted -O -a unit,CO2_Flux,d,, -a units,CO2_Flux,c,c,"Kg C/Km^2/sec" "${DEST_FILE}"

            # Fix time attribute: delete 'unit', add proper CF-compliant 'units' with reference time FOR THIS SPECIFIC DAY
            REF_DATE="${YEAR}-${MONTH}-${DAY} 00:00:00"
            ncatted -O -a unit,time,d,, -a units,time,c,c,"hours since ${REF_DATE}" "${DEST_FILE}"

            echo "   Fixed: ${YEAR}/${MONTH}/${DAY}"
        done
    done

    FF_COUNT=$(find "${DEST_FF}" -name "*.nc" 2>/dev/null | wc -l)
    echo "   Completed: ${FF_COUNT} Fossil Fuel files processed"
fi
fi  # end want FF

echo ""

# -----------------------------------------------------------------------------
# 2. Fix Ocean files (monthly: YYYY/MM.nc) - ADD TIME COORDINATE
# -----------------------------------------------------------------------------
if want OCEAN; then
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
        MONTH_NUM=$(basename "${MONTH_FILE}" .nc | sed 's/^0*//') # Strip leading zeros for Python
        DEST_FILE="${DEST_OCEAN}/${MONTH}"

        # Copy file first
        cp "${MONTH_FILE}" "${DEST_FILE}"

        # Fix CO2_Flux attribute: delete 'unit', add 'units'
        ncatted -O -a unit,CO2_Flux,d,, -a units,CO2_Flux,c,c,"Kg C/km2/sec" "${DEST_FILE}"

        # Fix coordinate attributes: delete 'unit', add 'units'
        ncatted -O -a unit,lon,d,, -a units,lon,c,c,"degrees_east" "${DEST_FILE}"
        ncatted -O -a unit,lat,d,, -a units,lat,c,c,"degrees_north" "${DEST_FILE}"

        # ADD time coordinate (day 1 of month, hour 0 - matches HEMCO time spec 1-12/1/0)
        HOURS_SINCE_1900=$(python3 -c "from datetime import datetime; d=datetime(${YEAR},${MONTH_NUM},1,0,0,0); ref=datetime(1900,1,1); print(int((d-ref).total_seconds()/3600))")

        # ncecat adds time as record dimension to data variables (CO2_Flux becomes
        # CO2_Flux(time,lat,lon)); coordinate variables lon(lon)/lat(lat) are preserved
        ncecat -O -u time "${DEST_FILE}" "${DEST_FILE}.tmp"
        ncap2 -O -s "time[time]=${HOURS_SINCE_1900}" "${DEST_FILE}.tmp" "${DEST_FILE}.tmp"
        ncatted -O \
            -a units,time,c,c,"hours since 1900-01-01 00:00:00" \
            -a calendar,time,c,c,"proleptic_gregorian" \
            -a long_name,time,c,c,"time" \
            "${DEST_FILE}.tmp" "${DEST_FILE}"
        rm -f "${DEST_FILE}.tmp"

        echo "   Fixed: ${YEAR}/${MONTH}"
    done

    OCEAN_COUNT=$(find "${DEST_OCEAN}" -name "*.nc" 2>/dev/null | wc -l)
    echo "   Completed: ${OCEAN_COUNT} Ocean files processed"
fi
fi  # end want OCEAN

echo ""

# -----------------------------------------------------------------------------
# 3. Fix NBE files (daily: YYYY/MM/DD.nc)
# -----------------------------------------------------------------------------
if want NBE; then
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

            DAY=$(basename "${DAY_FILE}" .nc)
            DEST_FILE="${DEST_NBE}/${MONTH}/${DAY}.nc"
            REF_DATE="${YEAR}-${MONTH}-${DAY} 00:00:00"

            # Resample 3-hourly (8 records, centered at 1.5,4.5,..,22.5 h) to HOURLY
            # (24 records at 0.5,1.5,..,23.5 h) by holding each 3-hourly value over its
            # 3 hours. HEMCO only cycles the diurnal cycle of a per-day file when there is
            # one record per requested hour (as with the hourly Fossilfuel field); with 8
            # records it holds a single slice all day, so the balanced biosphere must be
            # expanded to 24 hourly records. Also renames longitude/latitude -> lon/lat.
            python3 - "${DAY_FILE}" "${DEST_FILE}" "${REF_DATE}" <<'PYEOF'
import sys, numpy as np, netCDF4 as nc
src, dst, ref = sys.argv[1], sys.argv[2], sys.argv[3]
si = nc.Dataset(src)
lonv = 'longitude' if 'longitude' in si.variables else 'lon'
latv = 'latitude'  if 'latitude'  in si.variables else 'lat'
lon = np.array(si.variables[lonv][:]); lat = np.array(si.variables[latv][:])
flux = np.array(si.variables['CO2_Flux'][:], dtype='f4')   # (nt, lat, lon)
si.close()
nt = flux.shape[0]
if nt == 24:                       # already hourly
    flux24 = flux
elif 24 % nt == 0:                 # e.g. 8 (3-hourly) -> repeat each 24/8=3x
    flux24 = np.repeat(flux, 24 // nt, axis=0)
else:
    raise SystemExit(f"unexpected record count {nt} in {src}")
times = np.arange(24, dtype='f8') + 0.5
do = nc.Dataset(dst, 'w', format='NETCDF4')
do.createDimension('time', 24); do.createDimension('lat', lat.size); do.createDimension('lon', lon.size)
tv = do.createVariable('time', 'f8', ('time',)); tv.units = 'hours since ' + ref; tv.long_name = 'time'; tv[:] = times
yv = do.createVariable('lat', 'f4', ('lat',)); yv.units = 'degrees_north'; yv[:] = lat
xv = do.createVariable('lon', 'f4', ('lon',)); xv.units = 'degrees_east'; xv[:] = lon
fv = do.createVariable('CO2_Flux', 'f4', ('time', 'lat', 'lon',)); fv.units = 'Kg C/Km^2/sec'; fv[:] = flux24
do.close()
PYEOF

            echo "   Fixed (resampled 3-hourly -> hourly): ${YEAR}/${MONTH}/${DAY}"
        done
    done

    NBE_COUNT=$(find "${DEST_NBE}" -name "*.nc" 2>/dev/null | wc -l)
    echo "   Completed: ${NBE_COUNT} NBE files processed"
fi
fi  # end want NBE

echo ""

# -----------------------------------------------------------------------------
# 4. Fix GPP files (monthly: YYYY/MM.nc) - RENAME COORDINATES
# -----------------------------------------------------------------------------
if want GPP; then
echo "4. Processing GPP emissions (monthly files)..."
SRC_GPP="${SRC_BASE}/GPP/ORCHIDEE/${YEAR}"
DEST_GPP="${DEST_BASE}/GPP/ORCHIDEE/${YEAR}"

if [ ! -d "${SRC_GPP}" ]; then
    echo "   ERROR: Source directory not found: ${SRC_GPP}"
else
    # Create destination directory
    mkdir -p "${DEST_GPP}"

    # Process all monthly files
    for MONTH_FILE in ${SRC_GPP}/*.nc; do
        if [ ! -f "${MONTH_FILE}" ]; then
            echo "   No monthly files found in ${SRC_GPP}"
            break
        fi

        MONTH=$(basename "${MONTH_FILE}")
        MONTH_NUM=$(basename "${MONTH_FILE}" .nc | sed 's/^0*//')
        DEST_FILE="${DEST_GPP}/${MONTH}"

        # Copy file first
        cp "${MONTH_FILE}" "${DEST_FILE}"

        # Rename longitude -> lon and latitude -> lat
        ncrename -O -v longitude,lon -v latitude,lat "${DEST_FILE}"

        # ADD time coordinate (day 1 of month, hour 0 - matches HEMCO time spec 1-12/1/0)
        HOURS_SINCE_1900=$(python3 -c "from datetime import datetime; d=datetime(${YEAR},${MONTH_NUM},1,0,0,0); ref=datetime(1900,1,1); print(int((d-ref).total_seconds()/3600))")

        # ncecat adds time as record dimension to data variables (CO2_Flux becomes
        # CO2_Flux(time,lat,lon)); coordinate variables lon(lon)/lat(lat) are preserved
        ncecat -O -u time "${DEST_FILE}" "${DEST_FILE}.tmp"
        ncap2 -O -s "time[time]=${HOURS_SINCE_1900}" "${DEST_FILE}.tmp" "${DEST_FILE}.tmp"
        ncatted -O \
            -a units,time,c,c,"hours since 1900-01-01 00:00:00" \
            -a calendar,time,c,c,"proleptic_gregorian" \
            -a long_name,time,c,c,"time" \
            "${DEST_FILE}.tmp" "${DEST_FILE}"
        rm -f "${DEST_FILE}.tmp"

        echo "   Fixed: ${YEAR}/${MONTH}"
    done

    GPP_COUNT=$(find "${DEST_GPP}" -name "*.nc" 2>/dev/null | wc -l)
    echo "   Completed: ${GPP_COUNT} GPP files processed"
fi
fi  # end want GPP

echo ""

# -----------------------------------------------------------------------------
# 5. Fix TER files (monthly: YYYY/MM.nc) - RENAME COORDINATES
# -----------------------------------------------------------------------------
if want TER; then
echo "5. Processing TER emissions (monthly files)..."
SRC_TER="${SRC_BASE}/TER2/CLASS-CTEM/${YEAR}"
DEST_TER="${DEST_BASE}/TER/CLASS-CTEM/${YEAR}"

if [ ! -d "${SRC_TER}" ]; then
    echo "   ERROR: Source directory not found: ${SRC_TER}"
else
    # Create destination directory
    mkdir -p "${DEST_TER}"

    # Process all monthly files
    for MONTH_FILE in ${SRC_TER}/*.nc; do
        if [ ! -f "${MONTH_FILE}" ]; then
            echo "   No monthly files found in ${SRC_TER}"
            break
        fi

        MONTH=$(basename "${MONTH_FILE}")
        MONTH_NUM=$(basename "${MONTH_FILE}" .nc | sed 's/^0*//')
        DEST_FILE="${DEST_TER}/${MONTH}"

        # Copy file first
        cp "${MONTH_FILE}" "${DEST_FILE}"

        # Rename longitude -> lon and latitude -> lat
        ncrename -O -v longitude,lon -v latitude,lat "${DEST_FILE}"

        # ADD time coordinate (day 1 of month, hour 0 - matches HEMCO time spec 1-12/1/0)
        HOURS_SINCE_1900=$(python3 -c "from datetime import datetime; d=datetime(${YEAR},${MONTH_NUM},1,0,0,0); ref=datetime(1900,1,1); print(int((d-ref).total_seconds()/3600))")

        # ncecat adds time as record dimension to data variables (CO2_Flux becomes
        # CO2_Flux(time,lat,lon)); coordinate variables lon(lon)/lat(lat) are preserved
        ncecat -O -u time "${DEST_FILE}" "${DEST_FILE}.tmp"
        ncap2 -O -s "time[time]=${HOURS_SINCE_1900}" "${DEST_FILE}.tmp" "${DEST_FILE}.tmp"
        ncatted -O \
            -a units,time,c,c,"hours since 1900-01-01 00:00:00" \
            -a calendar,time,c,c,"proleptic_gregorian" \
            -a long_name,time,c,c,"time" \
            "${DEST_FILE}.tmp" "${DEST_FILE}"
        rm -f "${DEST_FILE}.tmp"

        echo "   Fixed: ${YEAR}/${MONTH}"
    done

    TER_COUNT=$(find "${DEST_TER}" -name "*.nc" 2>/dev/null | wc -l)
    echo "   Completed: ${TER_COUNT} TER files processed"
fi
fi  # end want TER

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
echo "     * time: unit -> units with day-specific reference time"
echo "  - Ocean:       ${DEST_BASE}/Ocean/ECCO-Darwin-MON-v05/${YEAR}/"
echo "     * CO2_Flux: unit -> units"
echo "     * lon/lat: unit -> units"
echo "     * Added time coordinate (middle of month)"
echo "  - NBE:         ${DEST_BASE}/Balbio/CARDAMOM-ECCO/${YEAR}/"
echo "     * time: added day-specific reference time"
echo "  - GPP:         ${DEST_BASE}/GPP/ORCHIDEE/${YEAR}/"
echo "     * Renamed longitude -> lon, latitude -> lat"
echo "  - TER:         ${DEST_BASE}/TER/CLASS-CTEM/${YEAR}/"
echo "     * Renamed longitude -> lon, latitude -> lat"
echo ""
echo "Next steps:"
echo "1. HEMCO_Config.rc and ExtData.rc already updated to point to corrected paths"
echo "2. Verify with: ncdump -h ${DEST_GPP}/01.nc | grep -E 'float (lon|lat)'"
