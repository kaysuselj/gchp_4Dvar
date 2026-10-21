#!/bin/bash
# setup_exp.sh — create a GCHP 14.8.0 4D-Var OSSE experiment in gchp_14.8.0_adj/4dvar_osse/exps/<name>
# (moved out of the 14.5.3 tree 2026-10-19; 14.5.3 is stale, kept in gchp_14.5.3_adjoint_surfaceF/4dvar_final)
#
#   ./setup_exp.sh <name> --prior-err "<make_prior_err.py options>" [options]
#   ./setup_exp.sh <name> --prior-err-file <existing E file> [options]
#
# The truth (truth/, obs/) is shared by all experiments: unmodified TER.  Each
# experiment has its own prior error E (prior TER = truth TER * E), written to
# exps/<name>/prior_err_co2.nc4 and read by the forward AND adjoint runs (HEMCO
# factor 752).  The optimizer should find sigma = 1/E.
#
# Options:
#   --prior-err "<opts>"     options for gchp_4Dvar/make_prior_err.py (without --out), e.g.
#                            "--type constant --value 1.3"
#                            "--type zonal --bands=-30:30:1.3,30:90:0.8"
#                            "--type random --amplitude 0.3 --length-km 1500 --seed 1"
#   --prior-err-file <path>  use an existing E file instead
#   --obs satellite|6hourly|24hourly   observations (default satellite = OCO-2 track,
#                            per-sounding AK; 6hourly/24hourly = C24 grid, mean AK)
#   --nodes N                PBS nodes (default 2)
#   --gchp 14.8|14.8-fullckpt|14.5.3   model version. DEFAULT 14.8 (standard since 2026-10-16, see
#                            gchp_14.8.0_adj/STANDARD_SETUP.md): optimized 14.8.0 = control_files_148_ckpt_reduced/
#                            (binary gchp.ckpt_20261016, ADJ_CKPT_REDUCED 1; gradient bit-identical to full
#                            checkpoints), 14.8.0 truth observations obs148/, GCHP_VERSION=14.8 in the driver.
#                            14.8-fullckpt = the same with full checkpoints (control_files_148/, binary
#                            gchp.adj_20261015; used by E1.3_const_sat_148). 14.5.3 = legacy.
#   --ftol X                 stop when |dJ/J| between iterations < X (default 0.02 for 14.8, E1.3 reaches it
#                            near iteration 10; driver default 1e-8 = never for 14.5.3)
#   --plot-every N           per-iteration plots every N iterations (+ iteration 1, convergence, end)
#                            (default 5 for 14.8, 1 for 14.5.3)
#   --res 24|48              cubed-sphere resolution (default 24). 48 = control_files_148_c48/ (reduced
#                            checkpoints only, MERRA-2 2x2.5 met), C48 truth obs148_c48/, E on the C48 grid
#   --adj rw|dl              advection adjoint: rw = reverse_wind_pure (default), dl = discrete_linear (TAF adjoint of the
#                            operational scheme + linearized mass fixer): control_files_148_ckpt_reduced_DL (C24 only so far)
#   --mode annual|monthly    control vector: one sigma field for the window (default) or one per month
#                            (monthly: C24 RW or DL, native CS, control_files_148_ckpt_reduced_monthly[_ocean][_DL],
#                            3-month window Jan-Mar 2016, E per month in exps/<name>/prior_err_monthly/)
#   --month-factors f1,f2,f3 monthly: E_m = (E from --prior-err) x f_m, e.g. --prior-err "--type constant --value 1"
#                            --month-factors 1.3,1.0,0.8 (O7)
#   --window 2016Q1|2016     monthly: 3-month window Jan-Mar 2016 (default) or the full year 2016 (12 months,
#                            T_END 2017-01-01, truth obs148_2016/ from truth148_2016/run_truth148_2016.pbs)
#   --ocean "<opts>"         DEV3 (monthly only): ocean flux in the control as field 2 (prior ocean = ocean x E_ocn,
#                            E_ocn from make_prior_err.py <opts>, e.g. "--type constant --value 0.8"; same sigma_b):
#                            control_files_148_ckpt_reduced_monthly_ocean, exps/<name>/prior_err_ocean_monthly/
#   --ocean-month-factors f1,..  monthly factors for E_ocn (default 1 every month)
#   --obs-dir <dir>          observation directory under 4dvar_osse/ instead of the default (e.g. obs148_noise_s1 for O9:
#                            gchp_4Dvar/make_noisy_obs.py adds seeded noise to the truth observations)
#   --lin                    O10: linear advection (hord_tr 2 / kord_tr 17) in truth AND 4D-Var, DL adjoint (exact):
#                            control_files_148_LIN, observations obs148_lin (truth148_lin/run_truth148_lin.pbs); C24 only
#   --dry-run                checks only
#
# Then:  cd exps/<name>; qsub 4dvar_optimizer.osse.<name>.run   (qsub -v MAX_ITER=..,SIGMA_B=.. to override)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
G4D="$(cd "${ROOT}/../gchp_4Dvar" && pwd)"
TEMPLATE="${G4D}/4dvar_optimizer.osse.unified.run"
# one binary for all new runs (user 2026-10-20): RW and DL (ADV_ADJ_METHOD) and the ocean field (ADJ_HEMCO_SFLUX2_DIAGN)
GCHP_BIN="${GCHP_BIN:-/nobackupp27/ksuselj1/gchp_14.8.0_adj/code/build_disc/bin/gchp.ocn_20261020}"
PYTHON_ENV="${PYTHON_ENV:-/nobackup/ksuselj1/envs/gchp_4dvar}"
CTRL_FWD="${ROOT}/control_files/forward"
CTRL_ADJ="${ROOT}/control_files/adjoint"

NAME=""; PERR=""; PERR_FILE=""; OBS="satellite"; NODES=2; DRY=false; GV=14.8; FTOLX=""; PLOTX=""; RES=24; ADJ=rw; CMODE=annual; MFAC=""; LIN=false; WIN=2016Q1; OCN=""; OMFAC=""; OBSDIR_X=""
while [ $# -gt 0 ]; do
    case "$1" in
        --prior-err)      PERR="$2"; shift 2 ;;
        --prior-err-file) PERR_FILE="$(readlink -f "$2")"; shift 2 ;;
        --obs)            OBS="$2"; shift 2 ;;
        --nodes)          NODES="$2"; shift 2 ;;
        --gchp)           GV="$2"; shift 2 ;;
        --ftol)           FTOLX="$2"; shift 2 ;;
        --plot-every)     PLOTX="$2"; shift 2 ;;
        --res)            RES="$2"; shift 2 ;;
        --adj)            ADJ="$2"; shift 2 ;;
        --mode)           CMODE="$2"; shift 2 ;;
        --lin)            LIN=true; shift ;;
        --month-factors)  MFAC="$2"; shift 2 ;;
        --window)         WIN="$2"; shift 2 ;;
        --obs-dir)        OBSDIR_X="$2"; shift 2 ;;
        --ocean)          OCN="$2"; shift 2 ;;
        --ocean-month-factors) OMFAC="$2"; shift 2 ;;
        --dry-run)        DRY=true; shift ;;
        --*)              echo "unknown option $1"; exit 1 ;;
        *)                NAME="$1"; shift ;;
    esac
done
[ -n "${NAME}" ] || { sed -n '2,24p' "$0"; exit 1; }
{ [ -n "${PERR}" ] || [ -n "${PERR_FILE}" ]; } || { echo "ERROR: give --prior-err or --prior-err-file"; exit 1; }

OBSD=obs
case "${GV}" in
    14.5.3) echo "ERROR: 14.5.3 is stale; use gchp_14.5.3_adjoint_surfaceF/4dvar_final/setup_exp.sh"; exit 1 ;;
    14.8|14.8-fullckpt)
            CF=control_files_148_ckpt_reduced; [ "${GV}" = 14.8-fullckpt ] && CF=control_files_148
            CTRL_FWD="${ROOT}/${CF}/forward"; CTRL_ADJ="${ROOT}/${CF}/adjoint"; OBSD=obs148
            [ -n "${FTOLX}" ] || FTOLX=0.005; [ -n "${PLOTX}" ] || PLOTX=5
            [ "${OBS}" = satellite ] || { echo "ERROR: --gchp 14.8 has satellite observations only"; exit 1; } ;;
    *) echo "ERROR: --gchp must be 14.8, 14.8-fullckpt or 14.5.3"; exit 1 ;;
esac
case "${RES}" in
    24) ;;
    48) [ "${GV}" = 14.8 ] || { echo "ERROR: --res 48 only with --gchp 14.8 (reduced checkpoints)"; exit 1; }
        CF=control_files_148_c48; CTRL_FWD="${ROOT}/${CF}/forward"; CTRL_ADJ="${ROOT}/${CF}/adjoint"; OBSD=obs148_c48 ;;
    *) echo "ERROR: --res must be 24 or 48"; exit 1 ;;
esac
case "${ADJ}" in
    rw) ;;
    dl) [ "${GV}" = 14.8 ] || { echo "ERROR: --adj dl only with --gchp 14.8 (reduced checkpoints)"; exit 1; }
        CF=control_files_148_ckpt_reduced_DL; [ "${RES}" = 48 ] && CF=control_files_148_c48_DL
        CTRL_FWD="${ROOT}/${CF}/forward"; CTRL_ADJ="${ROOT}/${CF}/adjoint" ;;
    *) echo "ERROR: --adj must be rw or dl"; exit 1 ;;
esac
if [ "${LIN}" = true ]; then
    [ "${GV}" = 14.8 ] && [ "${RES}" = 24 ] && [ "${ADJ}" = dl ] && [ "${CMODE}" = annual ] || { echo "ERROR: --lin needs --gchp 14.8 --res 24 --adj dl (annual)"; exit 1; }
    CF=control_files_148_LIN; CTRL_FWD="${ROOT}/${CF}/forward"; CTRL_ADJ="${ROOT}/${CF}/adjoint"; OBSD=obs148_lin
fi
case "${CMODE}" in
    annual) [ -z "${MFAC}" ] || { echo "ERROR: --month-factors needs --mode monthly"; exit 1; } ;;
    monthly) [ "${GV}" = 14.8 ] && [ "${RES}" = 24 ] && [ "${LIN}" = false ] || { echo "ERROR: --mode monthly only with --gchp 14.8 --res 24 (rw or dl, no --lin) so far"; exit 1; }
             [ -z "${PERR_FILE}" ] || { echo "ERROR: --mode monthly needs --prior-err (E is written per month)"; exit 1; }
             CF=control_files_148_ckpt_reduced_monthly; [ -n "${OCN}" ] && CF=${CF}_ocean; [ "${ADJ}" = dl ] && CF=${CF}_DL
             CTRL_FWD="${ROOT}/${CF}/forward"; CTRL_ADJ="${ROOT}/${CF}/adjoint"
             NMONTHS=3 ;;
    *) echo "ERROR: --mode must be annual or monthly"; exit 1 ;;
esac
[ "${CMODE}" = monthly ] || { [ -z "${OCN}${OMFAC}" ] && [ "${WIN}" = 2016Q1 ]; } || { echo "ERROR: --ocean/--ocean-month-factors/--window need --mode monthly"; exit 1; }
[ -z "${OMFAC}" ] || [ -n "${OCN}" ] || { echo "ERROR: --ocean-month-factors needs --ocean"; exit 1; }
case "${WIN}" in
    2016Q1) TAG=2016Q1 ;;
    2016)   TAG=2016; NMONTHS=12; OBSD=obs148_2016 ;;
    *) echo "ERROR: --window must be 2016Q1 or 2016"; exit 1 ;;
esac
COORD="${ROOT}/../inputs/osse_base/coords/coord_c${RES}.nc4"
[ -z "${OBSDIR_X}" ] || { [ -d "${ROOT}/${OBSDIR_X}" ] || { echo "ERROR: --obs-dir ${OBSDIR_X} not found in ${ROOT}"; exit 1; }; OBSD="${OBSDIR_X}"; }
case "${OBS}" in
    satellite) MODE=real;       TRACK=track_file_osse_2016.nc;        SEMI="${ROOT}/${OBSD}/semi_obs_xco2_${TAG:-2016Q1}.nc" ;;
    6hourly)   MODE=simplified; TRACK=track_file_c24_6h_3months.nc;   SEMI="${ROOT}/obs/semi_obs_xco2_6h_simplified_2016Q1.nc" ;;
    24hourly)  MODE=simplified; TRACK=track_file_c24_24h_3months.nc;  SEMI="${ROOT}/obs/semi_obs_xco2_24h_simplified_2016Q1.nc" ;;
    *) echo "ERROR: --obs must be satellite, 6hourly or 24hourly"; exit 1 ;;
esac
EXP="${ROOT}/exps/${NAME}"

NFAIL=0
ok()   { echo "  OK   $1"; }
fail() { echo "  FAIL $1"; NFAIL=$((NFAIL + 1)); }
echo "Pre-flight checks"
[ -x "${GCHP_BIN}" ] && ok "binary ${GCHP_BIN} (replaces the control sets' gchp links)" || fail "GCHP_BIN ${GCHP_BIN} not executable"
for d in "${CTRL_FWD}" "${CTRL_ADJ}"; do
    [ -f "${d}/gchp.env" ] && ok "$(basename "$d")/gchp.env" || fail "no gchp.env in ${d}"
done
[ -e "${CTRL_FWD}/${TRACK}" ] && ok "track file ${TRACK}" || fail "track file ${TRACK} not in ${CTRL_FWD}"
[ -f "${SEMI}" ] && ok "observations $(basename "${SEMI}")" || fail "observations missing: ${SEMI} (truth/run_truth.pbs done?)"
[ -f "${ROOT}/${OBSD}/truth_emissions_monthly_${TAG:-2016Q1}.nc" ] && ok "truth emissions (${OBSD})" || fail "${OBSD}/truth_emissions_monthly_${TAG:-2016Q1}.nc missing"
if [ "${CMODE}" = monthly ]; then
    grep -q '^752 TER_PRIOR_ERR ../prior_err_monthly/' "${CTRL_FWD}/HEMCO_Config.rc" && ok "forward reads monthly E (factor 752)" || fail "forward HEMCO does not read ../prior_err_monthly/"
    grep -q '^Sigma_CO2 .*sigma_monthly' "${CTRL_ADJ}/ExtData.rc" && ok "adjoint reads monthly sigma" || fail "adjoint ExtData does not read ../sigma_monthly/"
    if [ -n "${OCN}" ]; then
        grep -q '^753 OCN_PRIOR_ERR ../prior_err_ocean_monthly/' "${CTRL_FWD}/HEMCO_Config.rc" && ok "forward reads monthly ocean E (753)" || fail "forward HEMCO does not read ../prior_err_ocean_monthly/"
        grep -q '^Sigma_OCN .*sigma_ocean_monthly' "${CTRL_ADJ}/ExtData.rc" && ok "adjoint reads monthly ocean sigma" || fail "adjoint ExtData does not read ../sigma_ocean_monthly/"
        grep -q '^ADJ_HEMCO_SFLUX2_DIAGN: *EmisCO2_Ocean' "${CTRL_ADJ}/GCHP.rc" && ok "adjoint field 2 = EmisCO2_Ocean" || fail "adjoint GCHP.rc has no ADJ_HEMCO_SFLUX2_DIAGN: EmisCO2_Ocean"
        grep -q "SurfaceFluxAdj2_CO2" "${CTRL_ADJ}/HISTORY.rc" && ok "adjoint writes SurfaceFluxAdj2_CO2" || fail "adjoint HISTORY.rc does not write SurfaceFluxAdj2_CO2"
    fi
else
    grep -q '^752 TER_PRIOR_ERR ../prior_err_co2.nc4' "${CTRL_FWD}/HEMCO_Config.rc" && ok "forward reads E (factor 752)" || fail "forward HEMCO does not read ../prior_err_co2.nc4"
fi
python3 "${G4D}/check_adjoint_flux_config.py" "${CTRL_FWD}/HEMCO_Config.rc" "${CTRL_ADJ}/HEMCO_Config.rc" >/dev/null \
    && ok "adjoint flux = forward flux without sigma" || fail "check_adjoint_flux_config.py failed"
python3 "${G4D}/set_adjoint_met_offsets.py" "${CTRL_FWD}/ExtData.rc" "${CTRL_ADJ}/ExtData.rc" --check >/dev/null \
    && ok "adjoint met offsets = forward" || fail "adjoint met offsets differ from forward"
[ -f "${PERR_FILE:-/dev/null}" ] || [ -n "${PERR}" ] || fail "E file not found: ${PERR_FILE}"
for d in "${CTRL_FWD}" "${CTRL_ADJ}"; do
    grep -q "^GCHP.IM_WORLD: *${RES}\$" "${d}/GCHP.rc" && ok "$(basename "$d") GCHP.rc is C${RES}" || fail "$(basename "$d") GCHP.rc is not C${RES}"
done
AM=reverse_wind_pure; [ "${ADJ}" = dl ] && AM=discrete_linear
[ "$(awk -F'[[:space:]:]+' '/^ADV_ADJ_METHOD:/{print $2}' "${CTRL_ADJ}/GCHP.rc" | tail -1)" = "${AM}" ] \
    && ok "adjoint ADV_ADJ_METHOD = ${AM}" || fail "adjoint GCHP.rc ADV_ADJ_METHOD is not ${AM}"
[ -f "${COORD}" ] && ok "C${RES} coordinates $(basename "${COORD}")" || fail "coordinate file missing: ${COORD}"
[ "${NFAIL}" -eq 0 ] || { echo "${NFAIL} check(s) failed"; exit 1; }
[ "${DRY}" = true ] && { echo "dry run: nothing created"; exit 0; }

mkdir -p "${EXP}"/{forward_run,adjoint,forcing_files_osse,4dvar_output_osse,logs}
source "${PYTHON_ENV}/bin/activate"
if [ "${CMODE}" = monthly ]; then
    # shellcheck disable=SC2086
    python3 -I "${G4D}/make_prior_err.py" --out-dir "${EXP}/prior_err_monthly" --n-months "${NMONTHS}" \
        ${MFAC:+--month-factors "${MFAC}"} --cs-res "${RES}" --coord-file "${COORD}" ${PERR}
    if [ -n "${OCN}" ]; then
        # shellcheck disable=SC2086
        python3 -I "${G4D}/make_prior_err.py" --out-dir "${EXP}/prior_err_ocean_monthly" --n-months "${NMONTHS}" \
            ${OMFAC:+--month-factors "${OMFAC}"} --cs-res "${RES}" --coord-file "${COORD}" ${OCN}
    fi
elif [ -n "${PERR_FILE}" ]; then
    cp -p "${PERR_FILE}" "${EXP}/prior_err_co2.nc4"; echo "  E copied from ${PERR_FILE}"
    python3 -I -c "import xarray as xr,sys;d=xr.open_dataset(sys.argv[1]);n=d.Prior_Err.size;sys.exit(0 if n==6*${RES}*${RES} else 'E has %d cells, C${RES} needs %d'%(n,6*${RES}*${RES}))" "${EXP}/prior_err_co2.nc4"
else
    # shellcheck disable=SC2086
    python3 -I "${G4D}/make_prior_err.py" --out "${EXP}/prior_err_co2.nc4" --cs-res "${RES}" --coord-file "${COORD}" ${PERR}
fi

RUN="${EXP}/4dvar_optimizer.osse.${NAME}.run"
cp "${TEMPLATE}" "${RUN}"; chmod +x "${RUN}"
JOB="gchp_${NAME:0:14}"
sed -i "s/^#PBS -N .*/#PBS -N ${JOB}/" "${RUN}"
sed -i "s/^#PBS -l select=[0-9]*:ncpus=24:/#PBS -l select=${NODES}:ncpus=24:/" "${RUN}"
sed -i "s/^#PBS -m bea$/#PBS -m ae/" "${RUN}"
sed -i "s|SCRIPT_DIR=\"\${SCRIPT_DIR:-\${BASE_DIR}}\"|SCRIPT_DIR=\"\${SCRIPT_DIR:-${G4D}}\"|" "${RUN}"
sed -i "s|SELF=\"\${SELF:-\${BASE_DIR}/4dvar_optimizer.osse.unified.run}\"|SELF=\"\${SELF:-${RUN}}\"|" "${RUN}"
sed -i "s|_OSSE4DVAR_DIR=\"\${OSSE_ROOT:-\${PARENT_DIR}/osse_4dvar}\"|_OSSE4DVAR_DIR=\"\${OSSE_ROOT:-${ROOT}}\"|" "${RUN}"
sed -i "s|OBS_TYPE=\${OBS_TYPE:-satellite}|OBS_TYPE=\${OBS_TYPE:-${OBS}}|" "${RUN}"
sed -i "s|CTRL_FWD_DIR=\"\${CTRL_FWD_DIR:-\${_OBS_CTRL_FWD:-\${_CTRL_FWD_DEFAULT}}}\"|CTRL_FWD_DIR=\"\${CTRL_FWD_DIR:-${CTRL_FWD}}\"|" "${RUN}"
sed -i "s|CTRL_ADJ_DIR=\"\${CTRL_ADJ_DIR:-\${_CTRL_ADJ_DEFAULT}}\"|CTRL_ADJ_DIR=\"\${CTRL_ADJ_DIR:-${CTRL_ADJ}}\"|" "${RUN}"
sed -i "s|^GCHP_BIN=\"\${GCHP_BIN:-}\"|GCHP_BIN=\"\${GCHP_BIN:-${GCHP_BIN}}\"|" "${RUN}"
sed -i "s|_WORK_DIR_OVERRIDE=\"\${WORK_DIR:-}\"|_WORK_DIR_OVERRIDE=\"\${WORK_DIR:-${EXP}}\"|" "${RUN}"
sed -i "s|    LOG_NAME=\"4dvar_optimizer_osse\"$|    LOG_NAME=\"4dvar_optimizer_osse_${NAME}\"|" "${RUN}"
sed -i "s|SAT_TRACK_FILE=\"\${SAT_TRACK_FILE:-\${_OBS_SAT_TRACK}}\"|SAT_TRACK_FILE=\"\${SAT_TRACK_FILE:-${TRACK}}\"|" "${RUN}"
sed -i "s|OBS_MODE=\"\${OBS_MODE:-\${_OBS_MODE}}\"|OBS_MODE=\"\${OBS_MODE:-${MODE}}\"|" "${RUN}"
sed -i "s|SEMI_OBS_FILE=\"\${SEMI_OBS_FILE:-\${_OBS_SEMI_OBS}}\"|SEMI_OBS_FILE=\"\${SEMI_OBS_FILE:-${SEMI}}\"|" "${RUN}"
[ "${MODE}" = simplified ] && \
    sed -i "s|SIMPLIFIED_AK_FILE=\"\${SIMPLIFIED_AK_FILE:-\${_OBS_SIMPL_AK}}\"|SIMPLIFIED_AK_FILE=\"\${SIMPLIFIED_AK_FILE:-${SEMI}}\"|" "${RUN}"
if [ "${GV}" != 14.5.3 ]; then
    sed -i 's|^GCHP_VERSION=${GCHP_VERSION:-14.5.3}|GCHP_VERSION=${GCHP_VERSION:-14.8}|' "${RUN}"
    sed -i "s|\${_OSSE4DVAR_DIR}/obs/|\${_OSSE4DVAR_DIR}/${OBSD}/|g" "${RUN}"
    grep -q '^GCHP_VERSION=${GCHP_VERSION:-14.8}' "${RUN}" || { echo "ERROR: could not set GCHP_VERSION=14.8"; exit 1; }
fi
if [ "${CMODE}" = monthly ]; then
    sed -i -e 's|^MODE=${MODE:-annual}|MODE=${MODE:-monthly}|' -e 's|^    GRID=${GRID:-latlon}|    GRID=${GRID:-cs}|' \
           -e '0,/^    T_END_DEFAULT=""/s||    T_END_DEFAULT="2016-04-01"|' "${RUN}"
    grep -q '^MODE=${MODE:-monthly}' "${RUN}" && grep -q '^    GRID=${GRID:-cs}' "${RUN}" || { echo "ERROR: could not set monthly defaults"; exit 1; }
    if [ "${WIN}" = 2016 ]; then
        sed -i -e 's|^    T_END_DEFAULT="2016-04-01"$|    T_END_DEFAULT="2017-01-01"|' \
               -e 's|truth_emissions_monthly_2016Q1\.nc|truth_emissions_monthly_2016.nc|g' "${RUN}"
        grep -q '^    T_END_DEFAULT="2017-01-01"' "${RUN}" && ! grep -q 'truth_emissions_monthly_2016Q1' "${RUN}" \
            || { echo "ERROR: could not set the 2016 window"; exit 1; }
    fi
    if [ -n "${OCN}" ]; then
        sed -i 's|^N_FIELDS=${N_FIELDS:-1}$|N_FIELDS=${N_FIELDS:-2}|' "${RUN}"
        grep -q '^N_FIELDS=${N_FIELDS:-2}$' "${RUN}" || { echo "ERROR: could not set N_FIELDS=2"; exit 1; }
    fi
fi
[ -n "${FTOLX}" ] && sed -i "s|^FTOL=\${FTOL:-1e-8}|FTOL=\${FTOL:-${FTOLX}}|" "${RUN}"
# stop rule (user 2026-10-21): at most 20 iterations, early exit when |dJ/J| < ftol (0.005) on 3 consecutive accepted steps
sed -i -e "s|^PLATEAU_N=\${PLATEAU_N:-1}|PLATEAU_N=\${PLATEAU_N:-${PLATEAUX:-3}}|" -e "s|^MAX_ITER=\${MAX_ITER:-20}|MAX_ITER=\${MAX_ITER:-${MAXITX:-20}}|" "${RUN}"
grep -q "^PLATEAU_N=\${PLATEAU_N:-${PLATEAUX:-3}}" "${RUN}" || { echo "ERROR: could not set PLATEAU_N"; exit 1; }
[ -n "${PLOTX}" ] && sed -i "s|^PLOT_EVERY=\${PLOT_EVERY:-1}|PLOT_EVERY=\${PLOT_EVERY:-${PLOTX}}|" "${RUN}"
if [ "${RES}" = 48 ]; then
    # C48 truth (92 days, 2 nodes) took 45 min incl. sampling; budgets 1800 s fwd / 3600 s adj per month
    # so the driver chains before the walltime. Refine from the C48 truth / first iteration timings.
    # (raised 2026-10-21 to the budgets the C48 RW OSSE needed: 7200 s forward / 9000 s adjoint per month)
    sed -i -e "s|^EST_MONTH_FWD_S=\${EST_MONTH_FWD_S:-[0-9]*}|EST_MONTH_FWD_S=\${EST_MONTH_FWD_S:-7200}|" \
           -e "s|^EST_MONTH_ADJ_S=\${EST_MONTH_ADJ_S:-[0-9]*}|EST_MONTH_ADJ_S=\${EST_MONTH_ADJ_S:-9000}|" "${RUN}"
    grep -q '^EST_MONTH_ADJ_S=${EST_MONTH_ADJ_S:-9000}' "${RUN}" || { echo "ERROR: could not set C48 time estimates"; exit 1; }
fi
for k in "GCHP_BIN:-${GCHP_BIN}" "SCRIPT_DIR:-${G4D}" "OSSE_ROOT:-${ROOT}" "WORK_DIR:-${EXP}" "SEMI_OBS_FILE:-${SEMI}" "CTRL_ADJ_DIR:-${CTRL_ADJ}"; do
    grep -qF "${k}" "${RUN}" || { echo "ERROR: could not set ${k%%:-*} in ${RUN} (template changed?)"; exit 1; }
done

{
    echo "# $(date '+%F %T')  setup of ${NAME}"
    printf '%q ' "$0" "${NAME}"; [ -n "${PERR}" ] && printf -- '--prior-err %q ' "${PERR}"
    [ -n "${PERR_FILE}" ] && printf -- '--prior-err-file %q ' "${PERR_FILE}"
    printf -- '--obs %s --nodes %s --gchp %s --res %s --adj %s' "${OBS}" "${NODES}" "${GV}" "${RES}" "${ADJ}"
    [ -n "${FTOLX}" ] && printf -- ' --ftol %s' "${FTOLX}"; [ -n "${PLOTX}" ] && printf -- ' --plot-every %s' "${PLOTX}"
    [ "${CMODE}" = monthly ] && printf -- ' --mode monthly --month-factors %s' "${MFAC}"
    [ -n "${OBSDIR_X}" ] && printf -- ' --obs-dir %s' "${OBSDIR_X}"
    [ "${CMODE}" = monthly ] && printf -- ' --window %s' "${WIN}"
    [ -n "${OCN}" ] && printf -- ' --ocean %q' "${OCN}"; [ -n "${OMFAC}" ] && printf -- ' --ocean-month-factors %s' "${OMFAC}"
    [ "${LIN}" = true ] && printf -- ' --lin'; echo
} > "${EXP}/setup_cmd.sh"
echo
echo "Experiment ready: ${EXP}"
if [ "${CMODE}" = monthly ]; then echo "  E: monthly (${WIN}, ${NMONTHS} months), factors ${MFAC} -> ${EXP}/prior_err_monthly/"
    [ -n "${OCN}" ] && echo "  E_ocean: ${OCN} factors ${OMFAC:-1} -> ${EXP}/prior_err_ocean_monthly/ (N_FIELDS=2)"; else
echo "  E: $(python3 -I -c "import xarray as xr,warnings;warnings.simplefilter('ignore');d=xr.open_dataset('${EXP}/prior_err_co2.nc4');print(d.attrs.get('recipe','?'),'| min %.3f max %.3f'%(float(d.Prior_Err.min()),float(d.Prior_Err.max())))")"
fi
echo "  obs: ${OBS} ($(basename "${SEMI}")), nodes ${NODES}, C${RES}, adjoint ${ADJ}, binary $(basename "${GCHP_BIN}")"
echo "  submit: cd ${EXP} && qsub 4dvar_optimizer.osse.${NAME}.run"
