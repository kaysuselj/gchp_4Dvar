#!/bin/bash
# Standalone script to regenerate 4D-Var diagnostic plots on the front end.
#
# Uses the same path layout as 4dvar_optimizer.run so it can be run at any
# point during or after an optimization without submitting a PBS job.
#
# Run from gchp_4Dvar/:
#   cd gchp_4Dvar && bash run_plots.sh
#
# Override any variable on the command line:
#   PLOT_DIR=/tmp/plots bash run_plots.sh
#   PYTHON_ENV=/my/env  bash run_plots.sh

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ============================================================
#  Paths  (must match 4dvar_optimizer.run)
# ============================================================
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # gchp_4Dvar/
WORK_DIR="$(cd "${BASE_DIR}/../4dvar_design" && pwd)"       # 4dvar_design/
SCRIPT_DIR="${BASE_DIR}"

FORCING_DIR="${WORK_DIR}/forcing_files"
STATE_FILE="${WORK_DIR}/4dvar_state.npz"
OUTPUT_DIR="${WORK_DIR}/4dvar_output"
PLOT_DIR="${PLOT_DIR:-${WORK_DIR}/plots}"

# ============================================================
#  Python environment
# ============================================================
PYTHON_ENV="${PYTHON_ENV:-/nobackup/ksuselj1/envs/gchp_4dvar}"
if [ -f "${PYTHON_ENV}/bin/activate" ]; then
    set +e; source "${PYTHON_ENV}/bin/activate"; set -e
else
    die "Python env not found: ${PYTHON_ENV}  (override with PYTHON_ENV=...)"
fi
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"

if ! python3 -c "import numpy, xarray" 2>/dev/null; then
    die "Python env is missing required packages (numpy xarray)"
fi

# ============================================================
#  Checks
# ============================================================
[ -f "${STATE_FILE}" ]      || die "State file not found: ${STATE_FILE}"
[ -d "${OUTPUT_DIR}" ]      || die "Output dir not found: ${OUTPUT_DIR}"

mkdir -p "${PLOT_DIR}"

# ============================================================
#  Run plots
# ============================================================
log "Generating plots → ${PLOT_DIR}"
log "  state-file  : ${STATE_FILE}"
log "  output-dir  : ${OUTPUT_DIR}"
log "  forcing-dir : ${FORCING_DIR}"

python3 "${SCRIPT_DIR}/plot_4dvar.py" \
    --state-file  "${STATE_FILE}" \
    --output-dir  "${OUTPUT_DIR}" \
    --forcing-dir "${FORCING_DIR}" \
    --plot-dir    "${PLOT_DIR}"

log "Done — plots saved to ${PLOT_DIR}"
