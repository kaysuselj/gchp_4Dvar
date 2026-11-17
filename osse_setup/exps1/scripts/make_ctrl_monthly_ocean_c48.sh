#!/bin/bash
# exps1 (2026-11-11): C48 version of control_files_148_ckpt_min_monthly_ocean_adjm (monthly TER + ocean control), built
# from the C48 OSSE set control_files_148_c48_min_adjm (annual TER control; C48 grid, MERRA-2 2x2.5 met, 2 x 24 ranks)
# with exactly the changes that turned the C24 annual set control_files_148_ckpt_min_adjm into the C24 monthly-ocean set:
#   HEMCO  : TER 750/751/752 (as before), ocean 751 -> forward 751/753/754, adjoint 751/753; scale factors 750/752 read
#            ../sigma_monthly/ and ../prior_err_monthly/ (monthly, 2015-2017) instead of the annual files; new 753
#            OCN_PRIOR_ERR (../prior_err_ocean_monthly/) and 754 Sigma_OCN (../sigma_ocean_monthly/)
#   ExtData: the same four entries, F%y4-%m2-01 (monthly) instead of constant
#   adjoint: GCHP.rc ADJ_HEMCO_SFLUX2_DIAGN: EmisCO2_Ocean, HISTORY.rc writes SurfaceFluxAdj2_CO2
# Output: 4dvar_osse/control_files_148_c48_min_monthly_ocean_adjm. make_ctrl_exps1.sh then turns it into the exps1 sets.
set -euo pipefail
B=/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse
SRC=${B}/control_files_148_c48_min_adjm
REF=${B}/control_files_148_ckpt_min_monthly_ocean_adjm
OUT=${B}/control_files_148_c48_min_monthly_ocean_adjm
rm -rf "${OUT}"; mkdir -p "${OUT}"
for D in forward adjoint; do
    rsync -a --exclude 'OutputDir/*' --exclude 'gcadj_*' --exclude '*.log' --exclude cap_restart "${SRC}/${D}/" "${OUT}/${D}/"
    W=${OUT}/${D}
    # HEMCO: scale-factor definitions (replace the annual 750/752 lines by the four monthly lines of the C24 set)
    python3 - "${W}/HEMCO_Config.rc" "${REF}/${D}/HEMCO_Config.rc" "${D}" <<'EOF'
import re, sys
p, ref, d = sys.argv[1:4]
s, r = open(p).read(), open(ref).read()
new = [l for l in r.splitlines() if re.match(r"^75[0-4] (Sigma_CO2|TER_PRIOR_ERR|OCN_PRIOR_ERR|Sigma_OCN) ", l)]
assert len(new) == 4, new
lines = s.splitlines()
idx = [i for i, l in enumerate(lines) if re.match(r"^75[02] (Sigma_CO2|TER_PRIOR_ERR) ", l)]
assert len(idx) == 2 and idx[1] == idx[0] + 1, idx
lines[idx[0]:idx[1] + 1] = new
s = "\n".join(lines) + "\n"
ocn = "751/753/754" if d == "forward" else "751/753"
s, k = re.subn(r"^(0 OCEAN_ECCO_DARWIN .* CO2 )751( 2 500)$", r"\g<1>" + ocn + r"\2", s, flags=re.M)
assert k == 1, "ocean line"
open(p, "w").write(s)
EOF
    # ExtData: same four monthly entries
    python3 - "${W}/ExtData.rc" "${REF}/${D}/ExtData.rc" <<'EOF'
import re, sys
p, ref = sys.argv[1:3]
s, r = open(p).read(), open(ref).read()
new = [l for l in r.splitlines() if re.match(r"^(Sigma_CO2|TER_PRIOR_ERR|OCN_PRIOR_ERR|Sigma_OCN) ", l)]
assert len(new) == 4, new
lines = s.splitlines()
idx = [i for i, l in enumerate(lines) if re.match(r"^(Sigma_CO2|TER_PRIOR_ERR) ", l)]
assert len(idx) == 2 and idx[1] == idx[0] + 1, idx
lines[idx[0]:idx[1] + 1] = new
open(p, "w").write("\n".join(lines) + "\n")
EOF
    if [ "${D}" = adjoint ]; then
        line=$(grep '^ADJ_HEMCO_SFLUX2_DIAGN:' "${REF}/adjoint/GCHP.rc")
        grep -q '^ADJ_HEMCO_SFLUX_SELECTOR:' "${W}/GCHP.rc" || { echo "no ADJ_HEMCO_SFLUX_SELECTOR in C48 adjoint GCHP.rc"; exit 1; }
        sed -i "/^ADJ_HEMCO_SFLUX_SELECTOR:/a ${line}" "${W}/GCHP.rc"
        # HISTORY: add SurfaceFluxAdj2_CO2 after SurfaceFluxAdj_CO2 in the Adjoint collection
        python3 - "${W}/HISTORY.rc" "${REF}/adjoint/HISTORY.rc" <<'EOF'
import sys
p, ref = sys.argv[1:3]
s, r = open(p).read().splitlines(), open(ref).read().splitlines()
add = [l for l in r if "'SurfaceFluxAdj2_CO2" in l]
i = [k for k, l in enumerate(s) if "'SurfaceFluxAdj_CO2" in l]
assert len(add) == 1 and len(i) == 1 and not any("SurfaceFluxAdj2_CO2" in l for l in s), (add, i)
s.insert(i[0] + 1, add[0])                 # the same line as in the C24 set, right after SurfaceFluxAdj_CO2
open(p, "w").write("\n".join(s) + "\n")
EOF
    fi
done
echo "$(date '+%F %T') exps1: C48 monthly TER + ocean set, from ${SRC} with the changes of ${REF} (make_ctrl_monthly_ocean_c48.sh)" > "${OUT}/CHANGES.txt"
# check: apart from resolution, the result must equal the C24 monthly-ocean set
for D in forward adjoint; do
    for f in HEMCO_Config.rc ExtData.rc GCHP.rc HISTORY.rc; do
        n=$(diff <(grep -v '^ *#' "${REF}/${D}/${f}" | sed -e 's/4x5/RES/g' -e 's/2x25/RES/g' -e 's/c24/cRES/g') \
                 <(grep -v '^ *#' "${OUT}/${D}/${f}" | sed -e 's/4x5/RES/g' -e 's/2x25/RES/g' -e 's/c48/cRES/g') | grep -c '^[<>]' || true)
        echo "${D}/${f}: ${n} lines differ from the C24 set (after 4x5/2x25 and c24/c48)"
    done
done
echo "control set ${OUT}"
