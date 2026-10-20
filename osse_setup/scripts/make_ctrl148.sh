#!/bin/bash
# Build the GCHP 14.8.0 control files for the 4dvar_final OSSE (same design as the 14.5.3
# control_files_co2): control_files_148/{truth,forward,adjoint}.
#
# Base: gchp_14.8.0_adj/tests_final/rundir_co2_template (CO2 only, OSSE flux set, 4x5 MERRA-2)
# via tests_final/common/make_rundir.sh; then the OSSE specifics of 14.5.3:
#   HEMCO   TER_CLASS_CTEM scale factors: truth 751, forward 750/751/752, adjoint 751/752;
#           750 Sigma_CO2 (../sigma_co2_osse.nc4), 752 TER_PRIOR_ERR (../prior_err_co2.nc4)
#   ExtData Sigma_CO2 and TER_PRIOR_ERR entries
#   History truth/forward: Emissions (monthly mean) + Inst30 (every-step instantaneous
#           SpeciesConcVV_CO2, Met_PMIDDRY, Met_PSC2WET, Met_PSC2DRY) for the offline
#           satellite sampler; adjoint: Adjoint (SpeciesAdj_CO2, SurfaceFluxAdj_CO2, native CS)
#   GCHP.rc adjoint keys as 14.5.3 (reverse_wind_pure, airmass, massflux forward, forcing from
#           file, net-terrestrial = category-5 flux diagnostic)
#   physics transport + convection + PBL mixing, no chemistry/deposition, 1800 s steps
#   restart 14.5.3 OSSE initial restart with DELP_DRY = 0 (no 14.8.0 first-step pressure
#           scaling; its adjoint is not ported, port notes 2.5)
#   runConfig_forward.sh / runConfig_adjoint.sh: shims that take the settings the 4D-Var
#           driver (GCHP_VERSION=14.8) writes and apply them through setCommonRunSettings.sh
#   post_forward.sh (truth, forward): offline satellite sampling -> GEOSChem.sat_track.*
set -euo pipefail
die() { echo "ERROR: $*" >&2; exit 1; }
R=/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse   # moved out of the 14.5.3 tree 2026-10-19
B=/nobackupp27/ksuselj1/gchp_14.8.0_adj/inputs/osse_base   # C24 CMS-Flux restart + OCO-2 track file (copies of the 14.5.3 control_files_co2 inputs)
T148=/home5/ksuselj1/nobackup/gchp_14.8.0_adj/tests_final
BIN=${BIN:-/home5/ksuselj1/nobackup/gchp_14.8.0_adj/code/build_adj/bin/gchp.adj_20261015}
G4D=/nobackupp27/ksuselj1/gchp_14.8.0_adj/gchp_4Dvar
OUT=${R}/control_files_148
RST=${OUT}/common/GEOSChem.Restart.20160101_0000z.c24.noDELP.nc4
PY=/nobackup/ksuselj1/envs/gchp_4dvar/bin/python3
mkdir -p ${OUT}/common
${PY} -I - ${B}/GEOSChem.Restart.20160101_0000z.c24.nc4 ${RST} <<'EOF'
import sys, shutil, netCDF4 as nc
shutil.copy2(sys.argv[1], sys.argv[2])
with nc.Dataset(sys.argv[2], 'r+') as d:
    d['DELP_DRY'][:] = 0.0
print('restart with DELP_DRY = 0:', sys.argv[2])
EOF
cp -p ${B}/track_file_osse_2016.nc ${OUT}/common/

setkey() { sed -i -E "/^#?[[:space:]]*$1:/d" GCHP.rc; sed -i "/^MODEL_PHASE:/a $1: $2" GCHP.rc; }
shim() {   # $1 = file name, $2 = Model_Phase
cat > "$1" <<EOF
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
PHASE=$2
set -e
cd "\$(dirname "\$0")"
for f in GCHP.rc input.nml ExtData.rc HISTORY.rc HEMCO_Config.rc setCommonRunSettings.sh; do
    [ -s "\$f" ] || { echo "ERROR: empty or missing \$f"; exit 1; }
done
RL="Restarts/GEOSChem.Restart.\${Start_Time:0:8}_\${Start_Time:9:4}z.c24.nc4"
TGT=\$(readlink -f "\${INITIAL_RESTART}")
[ -f "\${TGT}" ] || { echo "ERROR: restart \${INITIAL_RESTART} not found"; exit 1; }
[ "\$(readlink -f "\${RL}" 2>/dev/null)" = "\${TGT}" ] || ln -sfn "\${TGT}" "\${RL}"
sed -i -e "s/^TOTAL_CORES=.*/TOTAL_CORES=\${TOTAL_CORES}/" -e "s/^NUM_NODES=.*/NUM_NODES=\${NUM_NODES}/" \\
       -e "s/^NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=\${NUM_CORES_PER_NODE}/" \\
       -e "s/^Run_Duration=.*/Run_Duration=\"\${Duration}\"/" -e "s/^Model_Phase=.*/Model_Phase=\${PHASE}/" setCommonRunSettings.sh
echo "\${Start_Time}" > cap_restart
set +eu; source setCommonRunSettings.sh; rc=\$?; source setRestartLink.sh > /dev/null; set -eu
[ \$rc -eq 0 ] || { echo "ERROR: setCommonRunSettings.sh failed"; exit 1; }
grep -q "^MODEL_PHASE: \${PHASE}" GCHP.rc || { echo "ERROR: MODEL_PHASE not \${PHASE}"; exit 1; }
echo "shim: \${PHASE} \${Start_Time} + \${Duration}, cores \${TOTAL_CORES}, restart \${TGT}"
EOF
chmod +x "$1"
}
post_forward() {
cat > post_forward.sh <<EOF
#!/bin/bash
# GCHP 14.8.0: offline satellite-track sampling of the every-step Inst30 output (MAPL 2.69 has
# no track_file sampler); writes GEOSChem.sat_track.* exactly as the 14.5.3 MAPL sampler did.
set -e
cd "\$(dirname "\$0")"
${PY} -I ${G4D}/make_sat_track_offline.py --track track_file_osse_2016.nc \\
    --inst-glob 'OutputDir/GEOSChem.Inst30.*.nc4' --out-dir OutputDir
rm -f OutputDir/GEOSChem.Inst30.*.nc4
EOF
chmod +x post_forward.sh
}

for D in truth forward adjoint; do
    W=${OUT}/${D}
    rm -rf ${W}
    bash ${T148}/common/make_rundir.sh ${W} ${BIN} 90 ${RST} > /dev/null
    cd ${W}
    # physics: transport + convection + PBL mixing (make_rundir defaults), 1800 s
    PH=FORWARD; [ ${D} = adjoint ] && PH=ADJOINT
    sed -i "s/^Model_Phase=.*/Model_Phase=${PH}/" setCommonRunSettings.sh
    # ---- HEMCO: scale factors on TER and the sigma / prior-error fields
    case ${D} in truth) SF=751 ;; forward) SF=750/751/752 ;; adjoint) SF=751/752 ;; esac
    sed -i -E "s|^(0 TER_CLASS_CTEM .* CO2 +)751 5 500$|\1${SF} 5 500|" HEMCO_Config.rc
    grep -q "^0 TER_CLASS_CTEM .* ${SF} 5 500$" HEMCO_Config.rc || die "${D}: TER scale factor not set"
    if [ ${D} != truth ]; then
        sed -i '/^751 CO2_UNIT_CONV/a 750 Sigma_CO2 ../sigma_co2_osse.nc4 Sigma_CO2 2016/1/1/0 C xy 1 1\n752 TER_PRIOR_ERR ../prior_err_co2.nc4 Prior_Err 2016/1/1/0 C xy 1 1' HEMCO_Config.rc
        sed -i '/^TER_CLASS_CTEM /a Sigma_CO2  1  N  N  -  none  none  Sigma_CO2  ../sigma_co2_osse.nc4\nTER_PRIOR_ERR  1  N  N  -  none  none  Prior_Err  ../prior_err_co2.nc4' ExtData.rc
        grep -q "^752 TER_PRIOR_ERR" HEMCO_Config.rc && grep -q "^TER_PRIOR_ERR " ExtData.rc || die "${D}: sigma/prior-err entries"
    fi
    # ---- History
    sed -i -E "s/^([[:space:]]*)'(SpeciesConc|Emissions|ProdLoss|StateMet|StateMetLevEdge)',/\1#'\2',/" HISTORY.rc
    sed -i "s/^COLLECTIONS: *'SpeciesConc',/COLLECTIONS: 'Placeholder_none',/" HISTORY.rc
    grep -q "^COLLECTIONS: 'Placeholder_none'," HISTORY.rc || die "${D}: COLLECTIONS line not as expected"
    if [ ${D} = adjoint ]; then
        sed -i "s/^COLLECTIONS: 'Placeholder_none',/COLLECTIONS: 'Adjoint',/" HISTORY.rc
        cat >> HISTORY.rc <<'EOF'
#==============================================================================
# Adjoint (4D-Var): sensitivities on the native cubed sphere, every step, backwards in time
  Adjoint.template:       '%y4%m2%d2.nc4',
  Adjoint.format:         'CFIO',
  Adjoint.timestampStart: .true.
  Adjoint.frequency:      003000
  Adjoint.duration:       240000
  Adjoint.mode:           'instantaneous'
  Adjoint.backwards:      1
  Adjoint.fields:         'SpeciesAdj_CO2        ', 'GCHPchem',
                          'SurfaceFluxAdj_CO2    ', 'GCHPchem',
::
EOF
    else
        sed -i "s/^COLLECTIONS: 'Placeholder_none',/COLLECTIONS: 'Emissions',\n             'Inst30',/" HISTORY.rc
        sed -i -e "s/^  Emissions.monthly:.*/  Emissions.monthly:          1/" -e "s/^  Emissions.mode:.*/  Emissions.mode:             'time-averaged'/" HISTORY.rc
        grep -q "^  Emissions.monthly: *1" HISTORY.rc || sed -i "/^  Emissions.template:/a \  Emissions.monthly:          1" HISTORY.rc
        cat >> HISTORY.rc <<'EOF'
#==============================================================================
# Inst30: every-step instantaneous fields for the offline satellite sampler (post_forward.sh)
  Inst30.template:       '%y4%m2%d2_%h2%n2z.nc4',
  Inst30.format:         'CFIO',
  Inst30.timestampStart: .false.
  Inst30.frequency:      003000
  Inst30.duration:       240000
  Inst30.mode:           'instantaneous'
  Inst30.fields:         'SpeciesConcVV_CO2    ', 'GCHPchem',
                         'Met_PMIDDRY    ', 'GCHPchem',
                         'Met_PSC2WET     ', 'GCHPchem',
                         'Met_PSC2DRY     ', 'GCHPchem',
::
EOF
        cp -p ${OUT}/common/track_file_osse_2016.nc .
        # AutoUpdate_Diagnostics=ON overrides the Emissions frequency: monthly means as in 14.5.3
        sed -i 's/^Diag_Monthly=.*/Diag_Monthly="1"/' setCommonRunSettings.sh
        post_forward
    fi
    # ---- GCHP.rc keys
    if [ ${D} = adjoint ]; then
        for kv in "ADV_ADJ_METHOD reverse_wind_pure" "ADV_ADJ_WEIGHT airmass" "ADV_ADJ_MASSFLUX forward" \
                  "ADJOINT_FROM_FILE 1" "OCO2_FORCING_DIR ${OUT}/forcing_placeholder/" "FD_TYPE LAYER" "LFD 1" "FD_SPEC CO2" \
                  "ADJ_HEMCO_SFLUX_SELECTOR DIAGN" "ADJ_HEMCO_SFLUX_DIAGN EmisCO2_NetTerrExch"; do setkey ${kv% *} ${kv#* }; done
    else
        setkey FD_STEP 0
    fi
    # ---- shim and one settings pass (writes GCHP.rc for C24 so the driver reads the right CS_RES)
    if [ ${D} = adjoint ]; then shim runConfig_adjoint.sh ADJOINT; ./runConfig_adjoint.sh > settings.log 2>&1
    else shim runConfig_forward.sh FORWARD; ./runConfig_forward.sh > settings.log 2>&1; fi
    grep -q "^GCHP.IM_WORLD: *24" GCHP.rc || die "${D}: GCHP.rc not C24 after settings"
    [ ${D} = adjoint ] && sed -i "s/^\([\t ]*\)#'Adjoint',/\1'Adjoint',/" HISTORY.rc
    rm -f cap_restart
    printf '# GCHP 14.8.0 environment for the 4D-Var driver (sourced with set -eu active)\nset +eu\nsource /home5/ksuselj1/nobackup/gchp_14.8.0_adj/gchp_14.8.0_env.sh > /dev/null 2>&1\nset -eu\n' > gchp.env
    echo "${D}: $(readlink -f gchp) | TER ${SF} | $(sed -n '/^COLLECTIONS:/,/^ *::/p' HISTORY.rc | grep -v '#' | tr -d ' \n')"
done
${PY} -I ${G4D}/check_adjoint_flux_config.py ${OUT}/forward/HEMCO_Config.rc ${OUT}/adjoint/HEMCO_Config.rc || die "flux-consistency guard failed"
echo "control files in ${OUT}"
