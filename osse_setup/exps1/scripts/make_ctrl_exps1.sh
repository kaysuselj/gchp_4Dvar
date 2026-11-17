#!/bin/bash
# exps1 control sets (2026-11-09), C24 only for now. Source = the O8 set control_files_148_ckpt_min_monthly_ocean_adjm
# (monthly TER + ocean control, ADJ_METHOD key, minimal checkpoints), changed to
#   - ADJ_METHOD: 3 (reverse-wind + simplified convection adjoint, no per-step checkpoints)
#   - slim ExtData.rc (44 met fields the CO2 simulation does not use -> /dev/null; bit-identical, scaling/scripts)
#   - TER from TER/ORCHIDEE or TER/CLASS-CTEM, ocean from ECCO-Darwin or MOM6-MON (inputs/surface_fluxes_osse)
#   - forward samples the LN + LG + OG track (exps1/obs_c24/track_file_lnlgog_2016.nc)
# The HEMCO container names stay TER_CLASS_CTEM / OCEAN_ECCO_DARWIN (the scripts and checks refer to them); only the
# file each one reads changes, and the README/CHANGES.txt of each set says which flux it is.
#     make_ctrl_exps1.sh <name> <TER: ORCHIDEE|CLASS-CTEM> <OCEAN: ECCO-Darwin-MON-v05|MOM6-MON> [truth|4dvar] [24|48] [rw|dl]
#   dl (2026-11-12, C24 only): discrete advection adjoint, source control_files_148_ckpt_min_monthly_ocean_DL_adjm,
#       ADJ_METHOD 4 (discrete advection + the same simplified convection adjoint as method 3; checkpoints TRAJ_PREADV)
#   truth: forward only, no sigma / prior-error factors (751 only), for the truth run
#   48 (2026-11-11): source = control_files_148_c48_min_monthly_ocean_adjm (C48, MERRA-2 2x2.5, 2 x 24 ranks), restart
#       exps1/inputs/GEOSChem.Restart.20151001_0000z.c48.noDELP.nc4, track exps1/obs_c24_1510_1705 (same soundings)
set -euo pipefail
NAME=${1:?name}; TER=${2:?TER model}; OCN=${3:?ocean model}; KIND=${4:-4dvar}; RES=${5:-24}; ADV=${6:-rw}
METHOD=3; [ "${ADV}" = dl ] && METHOD=4
B=/nobackupp27/ksuselj1/gchp_14.8.0_adj
SRC=${B}/4dvar_osse/control_files_148_ckpt_min_monthly_ocean_adjm
[ "${RES}" = 48 ] && SRC=${B}/4dvar_osse/control_files_148_c48_min_monthly_ocean_adjm
if [ "${ADV}" = dl ]; then
    [ "${RES}" = 24 ] || { echo "dl: C24 only"; exit 1; }
    SRC=${B}/4dvar_osse/control_files_148_ckpt_min_monthly_ocean_DL_adjm
fi
X=${B}/4dvar_osse/exps1
OUT=${X}/ctrl/${NAME}
I=${B}/inputs/surface_fluxes_osse
SLIM=${B}/scaling/scripts/slim_extdata.py
PY=/nobackup/ksuselj1/envs/gchp_4dvar/bin/python3
TRACK=${X}/obs_c24/track_file_lnlgog_2016.nc
[ "${RES}" = 48 ] && TRACK=${X}/obs_c24_1510_1705/track_file_lnlgog_201510_201705.nc
RST=${X}/inputs/GEOSChem.Restart.20151001_0000z.c${RES}.noDELP.nc4
for d in "${I}/TER/${TER}/2016" "${I}/Ocean/${OCN}/2016"; do [ -d "$d" ] || { echo "missing $d"; exit 1; }; done
[ -f "${TRACK}" ] || { echo "missing ${TRACK}"; exit 1; }
rm -rf "${OUT}"; mkdir -p "${OUT}"
PHASES="forward adjoint"; [ "${KIND}" = truth ] && PHASES=forward
for D in ${PHASES}; do
    W=${OUT}/${D}
    rsync -a --exclude 'OutputDir/*' --exclude 'gcadj_*' --exclude '*.log' --exclude cap_restart "${SRC}/${D}/" "${W}/"
    cd "${W}"
    sed -i "s/^ADJ_METHOD: .*/ADJ_METHOD: ${METHOD}/" GCHP.rc
    grep -q "^ADJ_METHOD: ${METHOD}\$" GCHP.rc || { echo "${D}: ADJ_METHOD not set"; exit 1; }
    # fluxes: same container, other file (HEMCO_Config.rc and ExtData.rc both carry the path)
    for f in HEMCO_Config.rc ExtData.rc; do
        sed -i -e "s|${I}/TER/CLASS-CTEM/|${I}/TER/${TER}/|g" -e "s|${I}/Ocean/ECCO-Darwin-MON-v05/|${I}/Ocean/${OCN}/|g" "$f"
        grep -q "${I}/TER/${TER}/" "$f" && grep -q "${I}/Ocean/${OCN}/" "$f" || { echo "${D}/${f}: flux paths not set"; exit 1; }
    done
    if [ "${KIND}" = truth ]; then
        # truth: unit conversion only (no sigma 750/754, no prior error 752/753)
        sed -i -E -e 's|(TER_CLASS_CTEM .* CO2 )750/751/752( 5 500)|\1751\2|' \
                  -e 's|(OCEAN_ECCO_DARWIN .* CO2 )751/753/754( 2 500)|\1751\2|' HEMCO_Config.rc
        grep -E '^0 (TER_CLASS_CTEM|OCEAN_ECCO_DARWIN) ' HEMCO_Config.rc | grep -q -E ' 75[02-4]' && { echo "truth: sigma/E factors left"; exit 1; }
        # and their ExtData imports (the 4D-Var set reads ../sigma_monthly etc.; control_files_148/truth has none)
        sed -i -E '/^(Sigma_CO2|TER_PRIOR_ERR|OCN_PRIOR_ERR|Sigma_OCN) /d' ExtData.rc
        sed -i -E '/^75[02-4] (Sigma_CO2|TER_PRIOR_ERR|OCN_PRIOR_ERR|Sigma_OCN) /d' HEMCO_Config.rc
        grep -q -E 'sigma_monthly|prior_err_monthly|prior_err_ocean_monthly|sigma_ocean_monthly' ExtData.rc HEMCO_Config.rc \
            && { echo "truth: sigma/E files still read"; exit 1; }
    fi
    # slim ExtData.rc
    cp ExtData.rc ExtData.rc.full
    ${PY} -I "${SLIM}" ExtData.rc.full ExtData.rc
    # LN + LG + OG track for the offline sampler (forward)
    if [ "${D}" = forward ]; then
        ln -sfn "${TRACK}" track_file_osse_2016.nc
    fi
    # fast offline sampler (identical output; the original takes ~1 h per month of LN + LG + OG soundings)
    [ -f post_forward.sh ] && sed -i "s|/nobackupp27/ksuselj1/gchp_14.8.0_adj/gchp_4Dvar/make_sat_track_offline.py|${X}/scripts/make_sat_track_offline_fast.py|" post_forward.sh
    if [ "${RES}" = 48 ]; then
        # 2015-10-01 C48 restart (the C24 spin-up state on the C48 grid, make_restart_cn.py)
        [ -f "${RST}" ] || { echo "missing ${RST}"; exit 1; }
        ln -sfn "${RST}" "Restarts/GEOSChem.Restart.20151001_0000z.c48.nc4"
        sed -i "s|^INITIAL_RESTART=.*|INITIAL_RESTART=./Restarts/GEOSChem.Restart.20151001_0000z.c48.nc4|" "runConfig_${D}.sh"
        grep -q '^GCHP.IM_WORLD: *48$' GCHP.rc || { echo "${D}: not C48"; exit 1; }
    fi
    {
        echo "$(date '+%F %T') exps1 control set ${NAME} (${KIND}): copy of ${SRC}/${D}"
        echo "  ADJ_METHOD ${METHOD} (${ADV} advection adjoint); slim ExtData.rc (full one kept as ExtData.rc.full)"
        echo "  TER = ${I}/TER/${TER} (HEMCO container TER_CLASS_CTEM), ocean = ${I}/Ocean/${OCN} (container OCEAN_ECCO_DARWIN)"
        [ "${D}" = forward ] && echo "  track_file_osse_2016.nc -> ${TRACK} (LN + LG + OG)"
        [ "${KIND}" = truth ] && echo "  truth: TER and ocean x 751 only (no sigma, no prior error)"
    } > CHANGES_exps1.txt
    echo "${D}: $(grep -E '^0 TER_CLASS_CTEM ' HEMCO_Config.rc | awk '{print $3, $11}')"
    echo "${D}: $(grep -E '^0 OCEAN_ECCO_DARWIN ' HEMCO_Config.rc | awk '{print $3, $11}')"
done
echo "control set ${OUT}"
