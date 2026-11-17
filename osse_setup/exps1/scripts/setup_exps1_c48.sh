#!/bin/bash
# exps1 C48 (2026-11-11): set up OSSE1_C48 (and OSSE2_C48) with setup_exp.sh, as setup_exps1_c24.sh but
#   C48 control sets exps1/ctrl/osse{1,2}_c48 (MERRA-2 2x2.5, 2 Broadwell nodes), C48 truth observations
#   exps1/obs_c48_1510_1705 (same soundings and noise as C24, sampled from the C48 truth), restart 2015-10-01 C48,
#   sigma_b rescaled with the cube size (x 48/24): OSSE1 0.4; OSSE2 = the files' Uncertainty x 2 (sigma_b.npz on C48),
#   sigma_b-scaled first step, at most 1 sigma_b per cell (FIRST_STEP_SIGMA_B=1.0), bounds TER [0, 4], ocean [-4, 4],
#   at most 60 iterations, early stop as C24, 5-day jobs in the long queue.
#     setup_exps1_c48.sh OSSE1|OSSE2|both [SUBMIT=1]
set -euo pipefail
WHICH=${1:?OSSE1, OSSE2 or both}
B=/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse
X=${B}/exps1
PY=/nobackup/ksuselj1/envs/gchp_4dvar/bin/python3
WIN=2015-10:2017-05; FIRST=2015-10; LAST=2017-04
OBSROOT=${X}/obs_c48_1510_1705
OSSE_SUBDIRS=OCO2b10_29Aug2020_LG:OCO2b10_29Aug2020_LN:OCO2b10_29Aug2020_OG
BOUNDS=0:4:-4:4
WALL=120:00:00
SCALE=2.0                                          # sigma_b x 48/24
cd "${B}"
source /nobackup/ksuselj1/envs/gchp_4dvar/bin/activate
[ -f "${OBSROOT}/semi_obs_xco2_201510_201705.nc" ] || { echo "C48 observations not made yet (truth_c48_1510_1705)"; exit 1; }
common=(--mode monthly --window "${WIN}" --res 48 --adj rw --nodes 2 --plot-every 5 --exps-dir "${X}" --obs-root "${OBSROOT}")
todo=(); [ "${WHICH}" = OSSE1 ] || [ "${WHICH}" = both ] && todo+=(OSSE1); [ "${WHICH}" = OSSE2 ] || [ "${WHICH}" = both ] && todo+=(OSSE2)
for o in "${todo[@]}"; do
    e=${o}_C48
    if [ "${o}" = OSSE1 ]; then
        MAXITX=60 ./setup_exp.sh "${e}" --prior-err "--type constant --value 1.3" --ocean "--type constant --value 1.3" \
            --ctrl "${X}/ctrl/osse1_c48" "${common[@]}"
        sed -i 's|^SIGMA_B=${SIGMA_B:-0.2}|SIGMA_B=${SIGMA_B:-0.4}|' "${X}/${e}/4dvar_optimizer.osse.${e}.run"
        grep -q '^SIGMA_B=${SIGMA_B:-0.4}' "${X}/${e}/4dvar_optimizer.osse.${e}.run" || { echo "${e}: SIGMA_B not set"; exit 1; }
    else
        MAXITX=60 ./setup_exp.sh "${e}" --prior-err "--type constant --value 1.0" --ocean "--type constant --value 1.0" \
            --ctrl "${X}/ctrl/osse2_c48" "${common[@]}"
        ${PY} -I "${X}/scripts/make_exps1_inputs.py" "${X}/${e}" osse2 "${FIRST}" "${LAST}" 48 "${SCALE}"
    fi
    R=${X}/${e}/4dvar_optimizer.osse.${e}.run
    grep -q '^MAX_ITER=${MAX_ITER:-60}' "${R}" || { echo "${e}: MAX_ITER is not 60"; exit 1; }
    grep -q '^PLATEAU_N=${PLATEAU_N:-3}' "${R}" && grep -q '^FTOL=${FTOL:-0.005}' "${R}" || { echo "${e}: stop rule not set"; exit 1; }
    grep -q '^N_FIELDS=${N_FIELDS:-2}' "${R}" || { echo "${e}: N_FIELDS is not 2"; exit 1; }
    grep -q '^T_START=${T_START:-2015-10-01}' "${R}" && grep -q '^    T_END_DEFAULT="2017-05-01"' "${R}" || { echo "${e}: window not set"; exit 1; }
    grep -q '^GCHP.IM_WORLD: *48$' "${X}/ctrl/$(echo "${o}" | tr 'A-Z' 'a-z')_c48/forward/GCHP.rc" || { echo "${e}: control set is not C48"; exit 1; }
done
V="SIGMA_BOUNDS_FIELD=${BOUNDS},OSSE_SUBDIRS=${OSSE_SUBDIRS},NODES=2,WALLTIME_S=432000,FIRST_STEP_SIGMA_B=1.0"
echo
for o in "${todo[@]}"; do
    e=${o}_C48; v="${V}"; [ "${o}" = OSSE2 ] && v="${V},SIGMA_B_FILE=${X}/${e}/sigma_b.npz"
    echo "${e}: qsub -q long -l walltime=${WALL} -v ${v}"
    if [ "${2:-}" = SUBMIT=1 ]; then
        cd "${X}/${e}" && /PBS/bin/qsub -q long -l walltime=${WALL} -v "${v}" "4dvar_optimizer.osse.${e}.run"
    fi
done
