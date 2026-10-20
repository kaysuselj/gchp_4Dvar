# GCHP 14.8.0 CO2 4D-Var OSSE setup

Everything needed to set up and run the OSSEs of the GCHP CO2 adjoint paper, apart from the model binary
and the input data. The scripts and control files are a mirror of the live setup on NASA Pleiades
(`/nobackupp27/ksuselj1/gchp_14.8.0_adj/4dvar_osse`), refreshed with `sync_from_pleiades.sh`. The 4D-Var
driver and the Python tools it calls live one level up in this repository (`gchp_4Dvar/`).

The paths inside the scripts and control files are the Pleiades paths. To run elsewhere, change them to
your locations (see "Inputs" below); a relocatable version is planned for the code release.

## Components

| Where | What |
|---|---|
| GCHP model | `kaysuselj/GCHP`, tag `v1.0-rc1` (`git clone --recurse-submodules -b v1.0-rc1 https://github.com/kaysuselj/GCHP.git`), build with `-DMECH=carbon -DADJOINT=yes -DREVERSE_OPERATORS=yes` |
| `../4dvar_optimizer.osse.unified.run` | 4D-Var driver (PBS job, self-resubmitting): forward → J_obs and adjoint forcing → adjoint → L-BFGS-B step, per iteration |
| `../*.py` | optimizer steps (`lbfgsb_step*.py`), control-vector writers (`write_sigma*.py`), prior error E (`make_prior_err.py`), observation operator (`co2_adjoint_forcing_osse.py`, `make_sat_track_offline.py`), plots, Taylor test (`taylor_test.*`) |
| `scripts/setup_exp.sh` | creates one experiment in `exps/<name>`: E file, a copy of the driver bound to a control-file set, observations and resolution |
| `scripts/make_ctrl148.sh`, `make_ctrl148_c48.sh` | build the base control files (`control_files_148`, `control_files_148_c48`: truth, forward, adjoint run directories) from a GCHP run directory |
| `control_files/` | the control-file sets (GCHP.rc, HEMCO_Config.rc, ExtData.rc, HISTORY.rc, input.nml, run scripts) |
| `truth/` | truth runs (C24, C48): unmodified fluxes, sampled along the OCO-2 track → synthetic observations |
| `fd_check/` | finite-difference check of the 4D-Var gradient for a set of control cells (P3b) |

### Control-file sets

| Set | Use |
|---|---|
| `control_files_148` | base set (truth / forward / adjoint), C24, full adjoint checkpoints |
| `control_files_148_ckpt_reduced` | **standard** C24 OSSE set: reduced checkpoints (gradient bit-identical), RW advection adjoint |
| `control_files_148_ckpt_reduced_DL` | as standard, DL advection adjoint (`ADV_ADJ_METHOD: discrete_linear`) |
| `control_files_148_ckpt_reduced_monthly` | as standard, monthly control (sigma and E read per month) |
| `control_files_148_c48`, `control_files_148_c48_DL` | C48 (reduced checkpoints), RW / DL |
| `control_files_148_ckpt_full`, `control_files_148_DL` | full-checkpoint variants (reference / P3b) |

Each set has `LINKS.txt` (the symlinks of the run directories: model binary, met, chemistry and HEMCO
data directories) and `DATA_FILES.txt` (the NetCDF files that were not copied: restarts, OCO-2 track).

## Running an OSSE

1. Build GCHP (tag above) and make the control files: `scripts/make_ctrl148.sh` (C24) /
   `make_ctrl148_c48.sh` (C48); the derived sets are copies with a few changed lines (described in `CHANGES.txt` where present, otherwise in the set table above).
2. Truth: `qsub truth/run_truth148.pbs` → `obs148/semi_obs_xco2_2016Q1.nc` (XCO2 along the OCO-2 track with
   per-sounding averaging kernels) and `obs148/truth_emissions_monthly_2016Q1.nc`.
3. Experiment: e.g. constant prior error 1.3 (O1):
   `scripts/setup_exp.sh E1.3_const --prior-err "--type constant --value 1.3" --obs satellite`
   options: `--res 24|48`, `--adj rw|dl`, `--mode annual|monthly` (+ `--month-factors 1.3,1.0,0.8`), `--ftol`
4. `cd exps/E1.3_const && qsub 4dvar_optimizer.osse.E1.3_const.run`  (the job resubmits itself; it stops at
   `|dJ/J| < ftol`, default 0.02, or at MAX_ITER). Results: `4dvar_output_osse/` (state per iteration,
   `4dvar_history.csv`), `plots_osse/`.
5. Gradient checks: `fd_check/run_fd148.pbs` (finite differences, P3b) and `../taylor_test.pbs` (Taylor test, G1).

Design (all experiments): truth = unmodified fluxes; prior terrestrial flux = truth × E; control sigma
multiplies the prior terrestrial flux, so the optimum is sigma = 1/E; OCO-2 track, no noise unless stated;
σ_b = 0.2, Jan–Mar 2016.

## Inputs (not in this repository)

| Input | Pleiades location | Notes |
|---|---|---|
| MERRA-2 met, 4x5 (C24) and 2x2.5 (C48) | `/nobackupp17/lwu5/GCHP/ExtData/GEOS_*/MERRA2` | standard GEOS-Chem ExtData |
| GEOS-Chem chemistry / HEMCO data | `/nobackup/ksuselj1/geos_chem_new/ExtData/{CHEM_INPUTS,HEMCO}` | standard GEOS-Chem ExtData |
| surface fluxes: TER (CLASS-CTEM), fossil (ODIAC), ocean (ECCO-Darwin), NBE (CARDAMOM), GPP (ORCHIDEE) | `gchp_14.8.0_adj/inputs/surface_fluxes_osse/` | paths in HEMCO_Config.rc |
| initial CO2 restarts (DELP_DRY = 0), cubed-sphere coordinates, OCO-2 track file | `gchp_14.8.0_adj/inputs/osse_base/`, `control_files_148*/common/` | see DATA_FILES.txt |

These will be archived with the code release (Zenodo).
