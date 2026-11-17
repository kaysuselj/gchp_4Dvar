#!/bin/bash
# exps1 C24: set up OSSE1_C24 and OSSE2_C24 with setup_exp.sh and the exps1 options, then (SUBMIT=1) submit.
#   OSSE1: prior TER (ORCHIDEE) and ocean (ECCO-Darwin) x 1.3, sigma_b 0.2 for both fields
#   OSSE2: prior TER = CLASS-CTEM, ocean = MOM6 (E = 1), sigma_b = the files' Uncertainty (sigma_b.npz)
#   both : window 2015-10-01 .. 2017-05-01 (19 months, user 2026-11-10; was 2016), monthly TER + ocean control,
#          bounds TER [0, 4], ocean [-4, 4], LN + LG + OG observations with noise (exps1/obs_c24_1510_1705), ADJ_METHOD 3,
#          at most 60 iterations, early stop |dJ/J| < 0.005 on 3 consecutive accepted steps, 2 Broadwell nodes,
#          5-day jobs in the long queue (user 2026-11-10)
#     setup_exps1_c24.sh [SUBMIT=1]
#     setup_exps1_c24.sh DL [SUBMIT=1]   (2026-11-12) OSSE2_C24_DL only: OSSE2 with the discrete advection adjoint
#                                        (exps1/ctrl/osse2_c24_dl, ADJ_METHOD 4), sigma_b-scaled first step (1 sigma_b)
#     setup_exps1_c24.sh DL1 [SUBMIT=1]  (2026-11-12) OSSE1_C24_DL: OSSE1 with the discrete advection adjoint
#                                        (exps1/ctrl/osse1_c24_dl, ADJ_METHOD 4), sigma_b 0.2, same first step
set -euo pipefail
B=/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse
X=${B}/exps1
PY=/nobackup/ksuselj1/envs/gchp_4dvar/bin/python3
WIN=2015-10:2017-05; FIRST=2015-10; LAST=2017-04
OBSROOT=${X}/obs_c24_1510_1705
RST=${X}/inputs/GEOSChem.Restart.20151001_0000z.c24.noDELP.nc4
OSSE_SUBDIRS=OCO2b10_29Aug2020_LG:OCO2b10_29Aug2020_LN:OCO2b10_29Aug2020_OG
BOUNDS=0:4:-4:4
WALL=120:00:00
cd "${B}"
source /nobackup/ksuselj1/envs/gchp_4dvar/bin/activate
if [ "${1:-}" = DL ] || [ "${1:-}" = DL1 ]; then
    if [ "${1}" = DL ]; then o=osse2; e=OSSE2_C24_DL; E0=1.0; else o=osse1; e=OSSE1_C24_DL; E0=1.3; fi
    ln -sfn "${OBSROOT}/track_file_lnlgog_201510_201705.nc" "${X}/ctrl/${o}_c24_dl/forward/track_file_osse_2016.nc"
    MAXITX=60 ./setup_exp.sh "${e}" --prior-err "--type constant --value ${E0}" --ocean "--type constant --value ${E0}" \
        --ctrl "${X}/ctrl/${o}_c24_dl" --mode monthly --window "${WIN}" --res 24 --adj dl --nodes 2 --plot-every 5 \
        --exps-dir "${X}" --obs-root "${OBSROOT}"
    [ "${o}" = osse2 ] && ${PY} -I "${X}/scripts/make_exps1_inputs.py" "${X}/${e}" osse2 "${FIRST}" "${LAST}"
    R=${X}/${e}/4dvar_optimizer.osse.${e}.run
    grep -q '^MAX_ITER=${MAX_ITER:-60}' "${R}" && grep -q '^N_FIELDS=${N_FIELDS:-2}' "${R}" \
        && grep -q '^T_START=${T_START:-2015-10-01}' "${R}" || { echo "${e}: settings not set"; exit 1; }
    grep -q '^ADJ_METHOD: 4$' "${X}/ctrl/${o}_c24_dl/adjoint/GCHP.rc" || { echo "${e}: not ADJ_METHOD 4"; exit 1; }
    V="SIGMA_BOUNDS_FIELD=${BOUNDS},OSSE_SUBDIRS=${OSSE_SUBDIRS},NODES=2,WALLTIME_S=432000,FIRST_STEP_SIGMA_B=1.0"
    [ "${o}" = osse2 ] && V="${V},SIGMA_B_FILE=${X}/${e}/sigma_b.npz"
    echo "${e}: qsub -q long -l walltime=${WALL} -v ${V}"
    [ "${2:-}" = SUBMIT=1 ] && { cd "${X}/${e}" && /PBS/bin/qsub -q long -l walltime=${WALL} -v "${V}" "4dvar_optimizer.osse.${e}.run"; }
    exit 0
fi
# the 19-month observation files and the LN + LG + OG track for the forward sampler
[ -f "${OBSROOT}/semi_obs_xco2_201510_201705.nc" ] || { echo "observations not made yet (truth_c24_1510_1705)"; exit 1; }
for c in osse1_c24 osse2_c24; do
    ln -sfn "${OBSROOT}/track_file_lnlgog_201510_201705.nc" "${X}/ctrl/${c}/forward/track_file_osse_2016.nc"
done
common=(--mode monthly --window "${WIN}" --res 24 --adj rw --nodes 2 --plot-every 5 --exps-dir "${X}" --obs-root "${OBSROOT}")
MAXITX=60 ./setup_exp.sh OSSE1_C24 --prior-err "--type constant --value 1.3" --ocean "--type constant --value 1.3" \
    --ctrl "${X}/ctrl/osse1_c24" "${common[@]}"
MAXITX=60 ./setup_exp.sh OSSE2_C24 --prior-err "--type constant --value 1.0" --ocean "--type constant --value 1.0" \
    --ctrl "${X}/ctrl/osse2_c24" "${common[@]}"
${PY} -I "${X}/scripts/make_exps1_inputs.py" "${X}/OSSE2_C24" osse2 "${FIRST}" "${LAST}"
for e in OSSE1_C24 OSSE2_C24; do
    R=${X}/${e}/4dvar_optimizer.osse.${e}.run
    grep -q '^MAX_ITER=${MAX_ITER:-60}' "${R}" || { echo "${e}: MAX_ITER is not 60"; exit 1; }
    grep -q '^PLATEAU_N=${PLATEAU_N:-3}' "${R}" && grep -q '^FTOL=${FTOL:-0.005}' "${R}" || { echo "${e}: stop rule not set"; exit 1; }
    grep -q '^N_FIELDS=${N_FIELDS:-2}' "${R}" || { echo "${e}: N_FIELDS is not 2"; exit 1; }
    grep -q '^T_START=${T_START:-2015-10-01}' "${R}" && grep -q '^    T_END_DEFAULT="2017-05-01"' "${R}" || { echo "${e}: window not set"; exit 1; }
done
# forward initial CO2 = the 2015-10-01 spin-up restart (the driver sets the dates and, for the adjoint, its own restart)
for c in osse1_c24 osse2_c24; do
    for d in forward; do
        W=${X}/ctrl/${c}/${d}
        ln -sfn "${RST}" "${W}/Restarts/GEOSChem.Restart.20151001_0000z.c24.nc4"
        sed -i "s|^INITIAL_RESTART=.*|INITIAL_RESTART=./Restarts/GEOSChem.Restart.20151001_0000z.c24.nc4|" "${W}/runConfig_${d}.sh"
        grep -q '^INITIAL_RESTART=./Restarts/GEOSChem.Restart.20151001_0000z.c24.nc4' "${W}/runConfig_${d}.sh" \
            || { echo "${c}/${d}: INITIAL_RESTART not set"; exit 1; }
    done
done
V="SIGMA_BOUNDS_FIELD=${BOUNDS},OSSE_SUBDIRS=${OSSE_SUBDIRS},NODES=2,WALLTIME_S=432000"
echo
echo "OSSE1_C24: qsub -q long -l walltime=${WALL} -v ${V}"
echo "OSSE2_C24: same + SIGMA_B_FILE=${X}/OSSE2_C24/sigma_b.npz"
if [ "${1:-}" = SUBMIT=1 ]; then
    cd "${X}/OSSE1_C24"
    qsub -q long -l walltime=${WALL} -v "${V}" 4dvar_optimizer.osse.OSSE1_C24.run
    cd "${X}/OSSE2_C24"
    qsub -q long -l walltime=${WALL} -v "${V},SIGMA_B_FILE=${X}/OSSE2_C24/sigma_b.npz" 4dvar_optimizer.osse.OSSE2_C24.run
fi
