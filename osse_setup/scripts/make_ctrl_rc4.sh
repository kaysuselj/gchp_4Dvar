#!/bin/bash
# v1.0-rc4 control sets (2026-10-28): copies of the minimal-checkpoint rc3 sets with the single checkpoint key ADJ_CKPT
# in place of ADJ_CKPT_REDUCED / ADJ_CKPT_TRAJ (rc4 stops on the old keys; rc3 ignores ADJ_CKPT, so the rc3 sets stay as
# they are for the runs that still use them). Key translation in every GCHP.rc of a set:
#     ADJ_CKPT_REDUCED 2 + ADJ_CKPT_TRAJ 0 -> ADJ_CKPT: 1     (minimal, RW)
#     ADJ_CKPT_REDUCED 2 + ADJ_CKPT_TRAJ 1 -> ADJ_CKPT: 2     (minimal + CO2 before advection, DL)
#     adjoint of a DL set (ADV_ADJ_METHOD discrete_linear)    -> ADJ_CKPT: 2
#     ADJ_CKPT_REDUCED 1 -> ADJ_CKPT: 3,  ADJ_CKPT_REDUCED 0 -> ADJ_CKPT: 4
# Everything else (symlinks, HEMCO, ExtData, HISTORY) is copied unchanged.
#     ./make_ctrl_rc4.sh            -> control_files_148_*_min*_rc4 (existing rc4 sets are replaced)
set -euo pipefail
cd "$(dirname "$0")"
SETS="control_files_148_ckpt_min control_files_148_ckpt_min_DL control_files_148_c48_min control_files_148_c48_min_DL
      control_files_148_ckpt_min_monthly control_files_148_ckpt_min_monthly_DL control_files_148_ckpt_min_monthly_ocean
      control_files_148_ckpt_min_monthly_ocean_DL control_files_148_LIN_min"
for s in ${SETS}; do
    d=${s}_rc4; rm -rf "${d}"; cp -a "${s}" "${d}"
    for f in "${d}"/*/GCHP.rc; do
        red=$(awk -F'[[:space:]:]+' '/^ADJ_CKPT_REDUCED:/{print $2}' "${f}" | tail -1)
        trj=$(awk -F'[[:space:]:]+' '/^ADJ_CKPT_TRAJ:/{print $2}' "${f}" | tail -1)
        met=$(awk -F'[[:space:]:]+' '/^ADV_ADJ_METHOD:/{print $2}' "${f}" | tail -1)
        case "${red:-0}" in
            2) v=1; { [ "${trj:-0}" = 1 ] || [ "${met}" = discrete_linear ]; } && v=2 ;;
            1) v=3 ;;
            *) v=4 ;;
        esac
        [ -z "${red}" ] && [ -z "${trj}" ] && v=""                 # no checkpoint key (truth runs): leave the default
        sed -i '/^ADJ_CKPT_REDUCED:/d; /^ADJ_CKPT_TRAJ:/d' "${f}"
        grep -q '^ADJ_CKPT:' "${f}" && sed -i '/^ADJ_CKPT:/d' "${f}"
        if [ -n "${v}" ]; then
            printf '# adjoint checkpoints (v1.0-rc4): 1 minimal (RW), 2 minimal + CO2 before advection (DL), 3 reduced, 4 full\nADJ_CKPT: %s\n' "${v}" >> "${f}"
        fi
        echo "${f}: REDUCED=${red:--} TRAJ=${trj:--} METHOD=${met:--} -> ADJ_CKPT ${v:-(default)}"
    done
done
