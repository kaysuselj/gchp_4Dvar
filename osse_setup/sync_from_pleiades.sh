#!/bin/bash
# Mirror the GCHP 14.8.0 OSSE setup (gchp_14.8.0_adj/4dvar_osse on Pleiades) into this folder.
# The live copies stay in 4dvar_osse (running experiments reference those paths); this folder is
# what is versioned and shared. Text files only: binaries, symlinks (listed in LINKS.txt per
# control-file set), restarts, track and observation files (NetCDF data, see README.md) and logs
# are left out.
#   ./sync_from_pleiades.sh           then review `git status` and commit in gchp_4Dvar
set -euo pipefail
SRC=/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse
DST="$(cd "$(dirname "$0")" && pwd)"
EXCL=(--no-links --exclude='*.nc4' --exclude='*.nc' --exclude='*.log' --exclude='*.o[0-9]*'
      --exclude='OutputDir/' --exclude='Restarts/' --exclude='*.stockflux_*' --exclude='*.bak*' --exclude='*.pre_*'
      --exclude='build/' --exclude='CreateRunDirLogs/' --exclude='utils/simulation-status')

# scripts
mkdir -p "${DST}/scripts" "${DST}/truth" "${DST}/fd_check"
for f in setup_exp.sh make_ctrl148.sh make_ctrl148_c48.sh SPEEDUP_REPORT.md; do cp -p "${SRC}/${f}" "${DST}/scripts/"; done
cp -p "${SRC}/truth148/run_truth148.pbs" "${SRC}/truth148_c48/run_truth148_c48.pbs" "${SRC}/truth148_lin/run_truth148_lin.pbs" "${SRC}/truth148_2016/run_truth148_2016.pbs" "${DST}/truth/"
cp -p "${SRC}/fd_check/run_fd148.pbs" "${SRC}/fd_check/compare_fd.py" "${SRC}/fd_check/cells_E1.3.txt" "${DST}/fd_check/"

# control-file sets
mkdir -p "${DST}/control_files"
for c in "${SRC}"/control_files_148*; do
    n=$(basename "${c}")
    rsync -a -q --delete "${EXCL[@]}" "${c}/" "${DST}/control_files/${n}/"
    ( cd "${c}" && find . -type l -printf '%p -> %l\n' | sort ) > "${DST}/control_files/${n}/LINKS.txt"
    ( cd "${c}" && find . -type f \( -name '*.nc4' -o -name '*.nc' \) -not -path '*/OutputDir/*' -printf '%p  %s bytes\n' | sort ) \
        > "${DST}/control_files/${n}/DATA_FILES.txt"
done
echo "synced from ${SRC}: $(find "${DST}" -type f | wc -l) files, $(du -sh "${DST}" | cut -f1)"
