# Simplified OSSE 4D-Var Pipeline Manual

**Purpose**: Run a GCHP CO2 surface flux 4D-Var OSSE using a fixed mean OCO-2 averaging
kernel instead of per-sounding L2# data.  The mean AK is averaged once from
ORCHIDEE-ECCO2, embedded in the semi-obs netCDF, and applied to every GCHP profile —
giving dense full-grid observations without any L2# dependency at run time.

---

## Overview

```
Truth run (GCHP, forward only)
   └─ sat_track output (GCHP CO2 + pressure at every C24 grid cell × timestep)
         │
         ▼
average_l2_ak.py           ← read ORCHIDEE-ECCO2 L2# once, average AK/pressure/apriori
   └─ mean_ak_Jan_Mar_2016.nc
         │
         ▼
make_semi_obs_xco2_simplified.py  ← apply mean AK to sat_track profiles
   └─ semi_obs_xco2_6h_simplified_2016Q1.nc  (xco2_truth, xco2_uncertainty, embedded AK)
         │
         ▼
setup_osse_exp.sh          ← create experiment directory, generate run script
   └─ osse_4dvar/<exp_name>/4dvar_optimizer.osse.<exp_name>.run
         │
         ▼
qsub (OBS_MODE=simplified) ← 4D-Var outer loop via co2_adjoint_forcing_osse.py
```

---

## Step 1 — Truth forward run (sat_track)

The truth run is already done in `osse_4dvar/semi-obs-6hours_3months/`.
If you need to re-run it:

```bash
cd osse_4dvar/semi-obs-6hours_3months
qsub run_semi_obs_forward.pbs
```

Outputs land in `forward_run/OutputDir/GEOSChem.sat_track.YYYYMMDD_0030z.nc4`.
MAPL writes the file with the **previous month's** timestamp:

| Data month | File name                              |
|------------|----------------------------------------|
| Jan 2016   | `GEOSChem.sat_track.20151201_0030z.nc4` |
| Feb 2016   | `GEOSChem.sat_track.20160101_0030z.nc4` |
| Mar 2016   | `GEOSChem.sat_track.20160201_0030z.nc4` |

MAPL also writes short-form duplicates (`201601_0030z.nc4` etc.) — ignore them;
the scripts use the canonical `YYYYMMDD01_` form to avoid double-counting.

---

## Step 2 — Compute mean averaging kernel

Run once per time window.  Reads ORCHIDEE-ECCO2 L2# data from:
`/nobackupp17/jliu7/OCO2/L2#/OCO2-B10/PSEUDO/TRENDY/ORCHIDEE-ECCO2/`

```bash
cd gchp_4Dvar/
source /nobackup/ksuselj1/envs/gchp_4dvar/bin/activate

python3 average_l2_ak.py \
    --t-start 2016-01-01 \
    --t-end   2016-04-01 \
    --out     mean_ak_Jan_Mar_2016.nc
```

Output variables (all in `mean_ak_Jan_Mar_2016.nc`, dimension `lev_sat=20`):

| Variable       | Units    | Description                              |
|----------------|----------|------------------------------------------|
| `mean_ak`      | 1        | Mean OCO-2 averaging kernel              |
| `mean_prs`     | hPa      | Mean satellite pressure levels (0 = above model top) |
| `mean_co2_apr` | mol/mol  | Mean CO2 a priori profile                |
| `mean_xco2_apr`| mol/mol  | Mean XCO2 a priori (scalar)              |
| `mean_xco2_std`| ppm      | Mean XCO2 uncertainty (= obs error std)  |

Typical values: ~15 000 soundings, mean_xco2_std ≈ 0.93 ppm, 20 valid pressure levels.

---

## Step 3 — Create simplified semi-observations

```bash
python3 make_semi_obs_xco2_simplified.py \
    --sat-track-dir osse_4dvar/semi-obs-6hours_3months/forward_run/OutputDir \
    --mean-ak-file  gchp_4Dvar/mean_ak_Jan_Mar_2016.nc \
    --t-start       2016-01-01 \
    --t-end         2016-04-01 \
    --out           osse_4dvar/semi-obs-6hours_3months/obs/semi_obs_xco2_6h_simplified_2016Q1.nc
```

Output netCDF (`semi_obs_xco2_6h_simplified_2016Q1.nc`):

| Variable            | Dimension | Units   | Description                                   |
|---------------------|-----------|---------|-----------------------------------------------|
| `xco2_truth`        | time      | ppm     | XCO2 truth from GCHP + mean AK               |
| `xco2_uncertainty`  | time      | ppm     | Fixed obs error (= mean_xco2_std, same for all) |
| `latitude`          | time      | deg N   |                                               |
| `longitude`         | time      | deg E   |                                               |
| `mean_ak`           | lev_sat   | 1       | Embedded mean AK (for forcing script)        |
| `mean_prs`          | lev_sat   | hPa     | Embedded mean pressure levels                |
| `mean_co2_apr`      | lev_sat   | mol/mol | Embedded mean CO2 a priori                   |
| `mean_xco2_apr`     | scalar    | mol/mol | Embedded mean XCO2 a priori                  |

The AK variables are embedded so **one file carries everything** — no separate
`--mean-ak-file` is needed at 4D-Var run time.

Expected obs count for Jan–Mar 2016: **15 459** (one profile per C24 cell per
sat_track pass, 6-hourly).

---

## Step 4 — Create 4D-Var experiment

```bash
cd osse_4dvar/

bash setup_osse_exp.sh <exp_name> \
    --semi-obs-file     semi-obs-6hours_3months/obs/semi_obs_xco2_6h_simplified_2016Q1.nc \
    --simplified-ak-file semi-obs-6hours_3months/obs/semi_obs_xco2_6h_simplified_2016Q1.nc
```

`--semi-obs-file` and `--simplified-ak-file` point to **the same file** (AK is embedded).

The script will:
1. Run pre-flight checks (control_files, gchp binary, semi-obs file, etc.)
2. Create `osse_4dvar/<exp_name>/` with subdirs: `forward_run/`, `adjoint/`,
   `forcing_files_osse/`, `4dvar_output_osse/`, `logs/`, `sat_track/`
3. Copy `gchp_4Dvar/4dvar_optimizer.osse.unified.run` to
   `osse_4dvar/<exp_name>/4dvar_optimizer.osse.<exp_name>.run` and patch:
   - `SCRIPT_DIR` → `gchp_4Dvar/` (Python scripts location)
   - `SELF` → this experiment's run script (for self-resubmission)
   - `CTRL_FWD_DIR` / `CTRL_ADJ_DIR` → `osse_4dvar/control_files/{forward_6hourlyobs,adjoint}`
   - `WORK_DIR` → the experiment directory itself (no git clone)
   - `OBS_MODE` → `simplified`
   - `SEMI_OBS_FILE` / `SIMPLIFIED_AK_FILE` → absolute path to semi-obs file

Available options:

```
--ctrl-fwd <name>          forward control_files subdir (default: forward_6hourlyobs)
--ctrl-adj <name>          adjoint control_files subdir (default: adjoint)
--semi-obs-file <path>     pre-processed semi-obs .nc file (required for OBS_MODE=simplified)
--simplified-ak-file <p>   same file as --semi-obs-file for simplified mode
--track-file <path>        override sat_track.rcx track file
--nodes N                  PBS node count (default: 2)
--dry-run                  check only, create nothing
```

---

## Step 5 — Submit the 4D-Var job

```bash
cd osse_4dvar/<exp_name>/

# Standard 3-month run (Jan–Mar 2016), up to 20 iterations:
qsub -v "T_END=2016-04-01" 4dvar_optimizer.osse.<exp_name>.run

# Test run: 2 iterations only
qsub -v "T_END=2016-04-01,MAX_ITER=2" 4dvar_optimizer.osse.<exp_name>.run

# Restart after walltime expiry (job self-resubmits automatically when CHAIN=true):
qsub -v "T_END=2016-04-01,RESTART=true" 4dvar_optimizer.osse.<exp_name>.run
```

Common `-v` overrides (all have baked-in defaults from `setup_osse_exp.sh`):

| Variable        | Default      | Description                              |
|-----------------|--------------|------------------------------------------|
| `T_END`         | 2016-02-01   | End of assimilation window (set to 2016-04-01 for 3 months) |
| `MAX_ITER`      | 20           | L-BFGS-B iteration limit                |
| `SIGMA_B`       | 0.2          | Background error std dev                |
| `NODES`         | 2            | PBS nodes (must also set `-l select=N`)  |
| `RESTART`       | false        | Resume from last completed phase         |
| `OBS_MODE`      | simplified   | Already baked in by setup_osse_exp.sh   |

---

## How the forcing works in simplified mode

`co2_adjoint_forcing_osse.py --simplified-ak-file <f> --semi-obs-file <f>` is called
by the run script during phase 2 (adjoint forcing) of each iteration.

For each GCHP sat_track profile in the assimilation window:
1. Interpolate GCHP CO2 (mol/mol) to `mean_prs[valid]` via `np.interp`
2. Compute: `xco2_hat = mean_xco2_apr + Σ(mean_ak × (co2_interp − mean_co2_apr))`  (ppm)
3. Look up truth: `xco2_truth` from the semi-obs file at the nearest time/location
4. Compute cost contribution: `J_obs += 0.5 × ((xco2_hat − xco2_truth) / xco2_std)^2`
5. Back-propagate: `force_profile = W @ (mean_ak × (xco2_hat − xco2_truth) / xco2_std^2)`
   where `W` is the `get_intmap` pressure-weighted interpolation matrix

Forcing files are written as `forcing_files_osse/CO2_adjoint_forcing_YYYYMMDD_HHMMz.nc4`.

---

## Directory layout after setup

```
osse_4dvar/
├── semi-obs-6hours_3months/
│   ├── forward_run/OutputDir/     ← truth GCHP sat_track files
│   └── obs/
│       └── semi_obs_xco2_6h_simplified_2016Q1.nc  ← simplified semi-obs
│
├── control_files/
│   ├── forward_6hourlyobs/        ← forward run template
│   └── adjoint/                   ← adjoint run template
│
└── <exp_name>/
    ├── 4dvar_optimizer.osse.<exp_name>.run  ← patched run script
    ├── setup_cmd.sh                          ← reproducibility record
    ├── forward_run/               ← GCHP forward run dir (filled by run script)
    ├── adjoint/                   ← GCHP adjoint run dir
    ├── forcing_files_osse/        ← CO2 adjoint forcing files
    ├── 4dvar_output_osse/         ← sigma, gradient, J history
    └── logs/                      ← optimizer log

gchp_4Dvar/
├── 4dvar_optimizer.osse.unified.run    ← canonical run script (template source)
├── co2_adjoint_forcing_osse.py         ← forcing + cost function (OBS_MODE aware)
├── average_l2_ak.py                    ← compute mean AK from L2#
├── make_semi_obs_xco2_simplified.py    ← apply mean AK to sat_track → semi-obs
├── mean_ak_Jan_Mar_2016.nc             ← mean AK output (Jan–Mar 2016)
└── setup_osse_exp.sh                   ← (symlink from osse_4dvar/)
```

---

## Switching between OBS_MODE=real and OBS_MODE=simplified

`setup_osse_exp.sh` automatically sets `OBS_MODE`:

- Pass `--simplified-ak-file` → `OBS_MODE=simplified`
- Omit `--simplified-ak-file` → `OBS_MODE=real` (standard OCO-2 track OSSE with per-sounding L2# AK)

For `OBS_MODE=real`, `--semi-obs-file` should point to the standard
`semi_obs_xco2_2016Q1.nc` file (produced by `make_semi_obs_xco2.py` with ORCHIDEE-ECCO2).

You can override `OBS_MODE` at submit time without regenerating the run script:
```bash
# Switch a simplified experiment to real mode for debugging:
qsub -v "T_END=2016-04-01,OBS_MODE=real,SEMI_OBS_FILE=/path/to/semi_obs_xco2_2016Q1.nc" \
     4dvar_optimizer.osse.<exp_name>.run
```

---

## Troubleshooting

**Semi-obs file has wrong obs count**: Re-check that `forward_run/OutputDir/` has exactly
3 canonical sat_track files (one per month). The script skips the short-form MAPL
duplicates automatically.

**`mean_xco2_std` is tiny (< 0.01)**: The L2# uncertainty is in mol/mol; `average_l2_ak.py`
converts to ppm (×1e6). Expected: ~0.93 ppm.

**`co2_adjoint_forcing_osse.py` finds 0 obs in window**: Check that `T_START`/`T_END`
match the semi-obs file time range and that `SEMI_OBS_FILE` is set correctly.

**`runConfig_forward.sh` domain decomp crash**: The control-files rsync runs on every
job start but the run script immediately re-patches `runConfig` for `NODES` and timing,
so this should not happen. If it does, check `NODES` is consistent between `#PBS -l select`
and the `-v NODES=` override.
