#!/bin/bash
# GCHP 14.8.0 shim with the interface of a 14.5.3 runConfig_*.sh: the 4D-Var driver
# (GCHP_VERSION=14.8) edits the lines below with sed and then runs this script, which applies
# them through setCommonRunSettings.sh, cap_restart and the restart link. Built by make_ctrl148.sh.
Start_Time="20160101 000000"
End_Time="20160401 000000"
Duration="00000300 000000"
TOTAL_CORES=24
NUM_NODES=1
NUM_CORES_PER_NODE=24
INITIAL_RESTART=./Restarts/GEOSChem.Restart.20160101_0000z.c24.nc4
PHASE=FORWARD
set -e
cd "$(dirname "$0")"
for f in GCHP.rc input.nml ExtData.rc HISTORY.rc HEMCO_Config.rc setCommonRunSettings.sh; do
    [ -s "$f" ] || { echo "ERROR: empty or missing $f"; exit 1; }
done
RL="Restarts/GEOSChem.Restart.${Start_Time:0:8}_${Start_Time:9:4}z.c24.nc4"
TGT=$(readlink -f "${INITIAL_RESTART}")
[ -f "${TGT}" ] || { echo "ERROR: restart ${INITIAL_RESTART} not found"; exit 1; }
[ "$(readlink -f "${RL}" 2>/dev/null)" = "${TGT}" ] || ln -sfn "${TGT}" "${RL}"
sed -i -e "s/^TOTAL_CORES=.*/TOTAL_CORES=${TOTAL_CORES}/" -e "s/^NUM_NODES=.*/NUM_NODES=${NUM_NODES}/" \
       -e "s/^NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=${NUM_CORES_PER_NODE}/" \
       -e "s/^Run_Duration=.*/Run_Duration=\"${Duration}\"/" -e "s/^Model_Phase=.*/Model_Phase=${PHASE}/" setCommonRunSettings.sh
echo "${Start_Time}" > cap_restart
set +eu; source setCommonRunSettings.sh; rc=$?; source setRestartLink.sh > /dev/null; set -eu
[ $rc -eq 0 ] || { echo "ERROR: setCommonRunSettings.sh failed"; exit 1; }
grep -q "^MODEL_PHASE: ${PHASE}" GCHP.rc || { echo "ERROR: MODEL_PHASE not ${PHASE}"; exit 1; }
echo "shim: ${PHASE} ${Start_Time} + ${Duration}, cores ${TOTAL_CORES}, restart ${TGT}"
