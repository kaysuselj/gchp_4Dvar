#!/bin/bash
# exps1 (2026-11-11): download the MERRA-2 2x2.5 met for the C48 OSSEs (window 2015-10-01 .. 2017-05-01 plus one day at each
# end for the ExtData brackets) from the GEOS-Chem S3 bucket into the 2x2.5 tree the C48 control sets read
# (/nobackup/ksuselj1/geos_chem_new/ExtData/GEOS_2x2.5/MERRA2/YYYY/MM/MERRA2.YYYYMMDD.<coll>.2x25.nc4).
# Files that already exist with the server's size are skipped; every file is checked against the server's Content-Length.
#     get_met_2x25.sh <first YYYY-MM-DD> <last YYYY-MM-DD> [parallel downloads, default 8]
set -uo pipefail
D0=${1:?first day}; D1=${2:?last day}; NPAR=${3:-8}
URL=https://geos-chem.s3.amazonaws.com/GEOS_2x2.5/MERRA2
OUT=/nobackup/ksuselj1/geos_chem_new/ExtData/GEOS_2x2.5/MERRA2
COLLS="A1 A3cld A3dyn A3mstC A3mstE I3"
list=$(mktemp)
d=${D0}
while [ "$(date -d "${d}" +%Y%m%d)" -le "$(date -d "${D1}" +%Y%m%d)" ]; do
    y=$(date -d "${d}" +%Y); m=$(date -d "${d}" +%m); ymd=$(date -d "${d}" +%Y%m%d)
    for c in ${COLLS}; do echo "${y}/${m}/MERRA2.${ymd}.${c}.2x25.nc4"; done >> "${list}"
    d=$(date -d "${d} + 1 day" +%Y-%m-%d)
done
echo "$(date '+%F %T') ${D0} .. ${D1}: $(wc -l < "${list}") files, ${NPAR} parallel"
get_one() {
    local rel=$1 url=$2 out=$3
    local want have
    want=$(curl -sfI "${url}/${rel}" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
    [ -n "${want}" ] || { echo "MISSING_ON_SERVER ${rel}"; return 0; }
    have=$(stat -c %s "${out}/${rel}" 2>/dev/null || echo 0)
    [ "${have}" = "${want}" ] && { echo "SKIP ${rel}"; return 0; }
    mkdir -p "$(dirname "${out}/${rel}")"
    for try in 1 2 3; do
        curl -sf --retry 3 -o "${out}/${rel}.part" "${url}/${rel}" && \
            [ "$(stat -c %s "${out}/${rel}.part")" = "${want}" ] && { mv "${out}/${rel}.part" "${out}/${rel}"; echo "OK ${rel}"; return 0; }
        sleep 5
    done
    rm -f "${out}/${rel}.part"; echo "FAILED ${rel}"
}
export -f get_one
xargs -P "${NPAR}" -I{} bash -c 'get_one "$@"' _ {} "${URL}" "${OUT}" < "${list}" > "${list}.log"
for s in OK SKIP FAILED MISSING_ON_SERVER; do echo "  ${s}: $(grep -c "^${s} " "${list}.log")"; done
grep -E "^(FAILED|MISSING_ON_SERVER) " "${list}.log" | head -20
echo "$(date '+%F %T') done"
rm -f "${list}" "${list}.log"
