# GCHP 4D-Var CO2 flux optimization — OSSE runners

This document covers the OSSE (Observing System Simulation Experiment) driver of
the GCHP 4D-Var CO2 surface-flux optimizer.

## Unified runner (`4dvar_optimizer.osse.unified.run`)

**`4dvar_optimizer.osse.unified.run`** is a single driver that covers both
control-variable layouts via a `MODE` switch:

| `MODE` | Control variable | σ writer / L-BFGS step | Control files | Work tree |
|--------|------------------|------------------------|---------------|-----------|
| `annual` (default) | one σ field over the window | `write_sigma.py` / `lbfgsb_step.py` | `control_files/{forward,adjoint}_osse` | `../4dvar_design_osse` |
| `monthly` | one σ field **per month** | `write_sigma_monthly.py` / `lbfgsb_step_monthly.py` | `control_files/{forward,adjoint}_osse_monthly` | `../4dvar_design_osse_monthly` |

```bash
cd gchp_4Dvar
qsub 4dvar_optimizer.osse.unified.run                      # annual (one-month), fresh
qsub -v MODE=monthly 4dvar_optimizer.osse.unified.run      # monthly (full-year), fresh
qsub -v RESTART=true[,MODE=monthly] 4dvar_optimizer.osse.unified.run   # resume
qsub -v MODE=annual,MAX_ITER=10,SIGMA_B=0.3 4dvar_optimizer.osse.unified.run
# short annual run wants less walltime (pass the matching WALLTIME_S):
qsub -l walltime=8:00:00 -v WALLTIME_S=28800 4dvar_optimizer.osse.unified.run
```

Key points of the unified driver:

- **`NMON`** (number of monthly fields, and the per-phase walltime estimate) is
  **derived from the window** `T_START..T_END`. An annual one-month run is
  `NMON=1`; the full-year monthly run is `NMON=12`.
- **Phase checkpointing + job chaining is always active** but is a no-op for
  short windows: each iteration runs as phases (`start → forward_done →
  forcing_done → adjoint_done`); before each GCHP phase the remaining walltime
  budget is checked, and a phase that will not fit triggers a self-resubmit with
  `RESTART=true`. A ~3 h annual iteration never resubmits; a ~36 h monthly
  iteration chains ~2 jobs. Disable with `CHAIN=false`.
- **Logs** are written under `<work tree>/logs/` (e.g.
  `4dvar_design_osse/logs/4dvar_optimizer_osse.log`).
- The observation window's `--t-end` is passed as `T00:00:00` in **both** modes
  (month-boundary-exclusive; the window's final month boundary is not double
  counted).

The two original drivers `4dvar_optimizer.osse.run` (one-month, single σ) and
`4dvar_optimizer.osse.monthly.run` (full-year, 12 monthly σ) have been **removed**
— the unified runner reproduces both. They remain available in the git history if
ever needed. The two sections below describe the behavior of each `MODE`.

For the base (non-OSSE) system see `README.md`; for verifying the adjoint
gradient that drives these optimizers see `adjoint_validation/README.md`.

## What "OSSE" means here

Instead of real OCO-2 retrievals, the cost function is evaluated against
**synthetic pseudo-observations** sampled from a known "truth" run
(ORCHIDEE-ECCO2). Because the truth is known, the recovered σ field can be
scored directly, which is what makes these runners useful for validating the
whole forward → forcing → adjoint → L-BFGS-B loop.

## Control variable and cost function

The control variable is a dimensionless scaling factor **σ** applied to the
terrestrial CO2 flux via HEMCO scale factor 750 (`EmisCO2_NetTerrExch`,
category 5). The prior is σ = 1 with a diagonal background
**B = σ_b² I** (default `σ_b = 0.2`). The cost function is

```
J(σ) = J_obs + J_b
     = (1/2) Σ_obs ((model−obs)/σ_obs)²  +  (1/2) Σ ((σ−1)/σ_b)²
```

The adjoint returns `dJ_obs/dσ`; L-BFGS-B adds the background term and takes the
step. σ scales a land flux, so the observation gradient is exactly zero over
ocean.

---

## `annual` mode — one-month, single σ

`MODE=annual` (the default) optimizes a single σ field. Each iteration reuses the
24-core allocation via `mpiexec`; for a one-month window the whole run fits in one
job (the phase/chaining machinery stays dormant).

**Iteration (repeated up to `MAX_ITER`):**

1. `write_sigma.py` writes `sigma_co2_osse.nc4` from the current L-BFGS-B state.
2. Forward GCHP simulation.
3. `co2_adjoint_forcing_osse.py` computes `J_obs` and the per-timestep adjoint
   forcing files from the OCO-2 pseudo-observation co-locations.
4. Move GCHP checkpoints from the forward to the adjoint run directory.
5. Adjoint GCHP simulation → `SurfaceFluxAdj_CO2` gradient.
6. `lbfgsb_step.py` takes one L-BFGS-B step and checks convergence.

**Layout** (paths relative to `gchp_4Dvar/`, the submit directory):

```
gchp_4Dvar/                         ← submit qsub from here
├── 4dvar_optimizer.osse.unified.run
├── write_sigma.py
├── co2_adjoint_forcing_osse.py
├── lbfgsb_step.py
├── plot_4dvar.py
└── control_files/
    ├── forward_osse/               ← forward run config (rsynced to run dir each job)
    └── adjoint_osse/               ← adjoint run config

../4dvar_design_osse/               ← execution area (WORK_DIR, sibling dir)
├── forward_run/  adjoint/          ← GCHP run directories (disposable)
├── sigma_co2_osse.nc4              ← written each iteration
├── 4dvar_state_osse.npz            ← L-BFGS-B state (persisted between iters)
├── forcing_files_osse/            ← J_value.txt + per-timestep forcing NetCDF
├── 4dvar_output_osse/             ← sigma_iter_NNN / J_iter_NNN snapshots
└── plots_osse/<YYYYMMDD_HHMMSS>/  ← per-run diagnostic figures
```

**Key parameters** (env-overridable via `qsub -v`):

| Variable | Default | Meaning |
|----------|---------|---------|
| `MAX_ITER` | 20 | Outer-loop iteration cap |
| `SIGMA_B` | 0.2 | Background std σ_b. **Do not change across a RESTART.** |
| `GTOL` / `FTOL` | 1e-5 / 1e-8 | L-BFGS-B convergence tolerances |
| `RESTART` | false | Resume from existing `4dvar_state_osse.npz` |
| `NLAT` / `NLON` | 46 / 72 | Control grid (lat-lon) dimensions |
| `SAVE_DIAG` | true | Write per-iteration snapshots + plots |
| `PYTHON_ENV` | `/nobackup/ksuselj1/envs/gchp_4dvar` | Conda env with numpy/pandas/xarray/netCDF4 |

Window is Jan 2016 (`T_START=2016-01-01`, `T_END=2016-02-01`); chemistry and
dynamic timesteps are 1800 s. Runtime ≈ 3 h/iteration (1.5 h forward + 1.5 h
adjoint); 20 iterations ≈ 60 h — set the `#PBS -l walltime` accordingly.

```bash
cd gchp_4Dvar
qsub 4dvar_optimizer.osse.unified.run                       # annual (default), fresh
qsub -v RESTART=true 4dvar_optimizer.osse.unified.run       # resume
qsub -v MAX_ITER=10,SIGMA_B=0.3 4dvar_optimizer.osse.unified.run
```

---

## `monthly` mode — full-year, 12 monthly σ

`MODE=monthly`: the control variable is **12 monthly σ fields** over a full-year
window (`2016-01-01 → 2017-01-01`). Each monthly field scales that month's TER
flux; the monthly gradient `g_m` is obtained by **differencing** the backward
`SurfaceFluxAdj_CO2` accumulator at the month boundaries
(`g_m = A(t_m) − A(t_{m+1})`, with `A(t_end) = 0` by construction).

**Why it is structured differently:** one iteration is ~36 h (12-month forward +
12-month adjoint) and does **not** fit in a single 24 h job. So the script
**checkpoints phases and chains jobs**:

```
phase sequence per iteration:
  start        → write σ, run forward GCHP, move checkpoints
  forward_done → compute J_obs + adjoint forcing (12 sat_track files joined)
  forcing_done → run adjoint GCHP, then free ~840 GB of gcadj checkpoints
  adjoint_done → L-BFGS-B step (lbfgsb_step_monthly.py), snapshot, plots, next iter
```

The last **completed** phase is recorded in `4dvar_progress_osse_monthly.txt`.
Before each GCHP phase the script checks the remaining walltime budget; if the
phase will not fit, it **resubmits itself** with `RESTART=true` and exits, and
the next job resumes at the recorded phase. A job killed mid-phase redoes only
that phase. Expect roughly two chained jobs per iteration at 24 h walltime.

Fully isolated from the one-month optimizer: separate work tree
(`../4dvar_design_osse_monthly`), control files (`control_files/*_osse_monthly`),
sigma dir, forcing dir, state file, progress file, log and output dirs.

**Additional / changed parameters:**

| Variable | Default | Meaning |
|----------|---------|---------|
| `T_START` / `T_END` | `2016-01-01` / `2017-01-01` | Window; `NMON` and `START_MONTH` are derived from it |
| `MAX_ITER` | 20 | **Global** iteration cap across all chained jobs |
| `CHAIN` | true | Resubmit self when walltime is spent (`CHAIN=false` to disable) |
| `MAX_CHAIN` | 60 | Safety cap on chained resubmissions |
| `WALLTIME_S` | 86400 | Must match the `#PBS -l walltime` |
| `SIGMA_B`, `GTOL`, `FTOL`, `RESTART` | as above | Same meaning as in `annual` mode |

Supporting scripts used in `monthly` mode: `write_sigma_monthly.py`,
`lbfgsb_step_monthly.py`, `check_monthly_sigma.py` (plumbing check that ExtData
applies the right monthly σ), and the OSSE-window `co2_adjoint_forcing_osse.py`
(joins all 12 monthly `sat_track` files).

```bash
cd gchp_4Dvar
qsub -v MODE=monthly 4dvar_optimizer.osse.unified.run                     # fresh start
qsub -v RESTART=true,MODE=monthly 4dvar_optimizer.osse.unified.run        # resume
qsub -v MODE=monthly,CHAIN=false 4dvar_optimizer.osse.unified.run         # single job, no self-chaining
# short shakedown (3-month window, one iteration, no chaining):
qsub -v MODE=monthly,T_END=2016-04-01,MAX_ITER=1,CHAIN=false 4dvar_optimizer.osse.unified.run
```

---

## Before trusting a production run: validate the adjoint

The L-BFGS-B step size is only as good as the adjoint gradient. Two layers of
validation live in this repo:

- **`adjoint_validation/`** — native cubed-sphere FD-vs-adjoint harness (the
  three-experiment suite). See `adjoint_validation/README.md`.
- **`fd_validate_gradient.run` / `fd_gradient_tools.py`** and
  **`fd_validate_multicell.run` / `fd_multicell_tools.py`** — lat-lon
  single-cell and multi-cell FD probes of the monthly gradient chain.

Do not start a production optimization until the adjoint gradient passes these
checks (or the known limiter/curvature caveats are explicitly accepted).

## Requirements

- GCHP with adjoint support (forward and adjoint executables staged into the
  `forward_run/` and `adjoint/` run directories).
- Python ≥ 3.9 with `numpy`, `pandas`, `xarray`, `netCDF4` (the `PYTHON_ENV`
  conda environment).
