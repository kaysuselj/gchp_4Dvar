#!/bin/bash
# GCHP 14.8.0: offline satellite-track sampling of the every-step Inst30 output (MAPL 2.69 has
# no track_file sampler); writes GEOSChem.sat_track.* exactly as the 14.5.3 MAPL sampler did.
set -e
cd "$(dirname "$0")"
/nobackup/ksuselj1/envs/gchp_4dvar/bin/python3 -I /nobackupp27/ksuselj1/gchp_14.8.0_adj/gchp_4Dvar/make_sat_track_offline.py --track track_file_osse_2016.nc \
    --inst-glob 'OutputDir/GEOSChem.Inst30.*.nc4' --out-dir OutputDir
rm -f OutputDir/GEOSChem.Inst30.*.nc4
