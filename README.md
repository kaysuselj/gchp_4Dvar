# GCHP 4D-Var CO2 Surface Flux Optimizer

Batch-queue 4D-Var system for optimizing gridded CO2 surface flux scaling
factors using GCHP forward and adjoint simulations and OCO-2 observations.

## Overview

The optimizer runs an outer loop of up to `MAX_ITER` iterations entirely
within a single PBS allocation — no child `qsub` calls, no queue waits
between iterations. Each iteration consists of:

1. Write `sigma_co2.nc4` (scaling factor field) from the current L-BFGS-B state
2. Run the forward GCHP simulation
3. Compute the observation cost *J_obs* and adjoint forcing from OCO-2 co-locations
4. Move GCHP checkpoints from the forward to the adjoint run directory
5. Run the adjoint GCHP simulation
6. Take one L-BFGS-B step; check convergence

The **control variable** is σ(lat, lon) — a dimensionless scaling factor applied
to all active surface CO2 emission fields via HEMCO scale factor 750.
The prior is σ = 1 (unperturbed emissions) with background covariance
B = σ_b² I (diagonal, default σ_b = 0.5).

## Directory layout

```
gchp_4Dvar/                     ← submit qsub from here
├── 4dvar_optimizer.run          ← PBS job script (this system)
├── control_files/
│   ├── forward/                 ← canonical config files for the forward run
│   └── adjoint/                 ← canonical config files for the adjoint run
├── write_sigma.py               ← writes sigma_co2.nc4 from optimizer state
├── co2_adjoint_forcing.py       ← computes J_obs and adjoint forcing files
├── lbfgsb_step.py               ← one L-BFGS-B update step
└── co2_sat_compare_monthly.py   ← OCO-2 observation operator (imported by above)

4dvar_design/                    ← execution area (sibling directory)
├── forward_run/                 ← GCHP forward run directory
├── adjoint/                     ← GCHP adjoint run directory
├── sigma_co2.nc4                ← written by write_sigma.py each iteration
├── 4dvar_state.npz              ← L-BFGS-B state (persisted between iterations)
├── forcing_files/               ← J_value.txt + per-timestep forcing NetCDF
└── 4dvar_output/                ← 4dvar_state_iter_NNN.npz diagnostic snapshots
```

Config files are kept in `control_files/` and rsynced to the run directories
at the start of every optimizer invocation. Edit configs in `control_files/`;
the run directories are treated as disposable.

## Requirements

### Software
- GCHP with adjoint support (forward and adjoint executables in `4dvar_design/forward_run/` and `4dvar_design/adjoint/`)
- Python ≥ 3.9 environment with: `numpy`, `pandas`, `xarray`, `netCDF4`

### Python environment setup (one time, on a login node)
```bash
module load python3/3.11        # or whichever module is available
conda create -p /nobackup/<user>/envs/gchp_4dvar \
    python=3.11 numpy pandas xarray netCDF4 -c conda-forge
```

### HEMCO configuration
Both `control_files/forward/HEMCO_Config.rc` and
`control_files/adjoint/HEMCO_Config.rc` must include scale factor 750:
```
750 Sigma_CO2 ../sigma_co2.nc4 Sigma_CO2 2019/1/1/0 C xy 1 1
```
and the `/750` multiplier must be appended to every active surface CO2
emission line.

### ExtData.rc
Both `control_files/forward/ExtData.rc` and
`control_files/adjoint/ExtData.rc` must register the field so MAPL reads it:
```
Sigma_CO2  1  N  N  -  none  none  Sigma_CO2  ../sigma_co2.nc4
```

### HISTORY.rc (adjoint)
`control_files/adjoint/HISTORY.rc` must output `SurfaceFluxAdj_CO2` on a
lat-lon grid matching `NLAT × NLON`:
```
Adjoint.grid_label:   PC72x46-DC
Adjoint.conservative: 1
```

## Submitting

```bash
cd /path/to/gchp_4Dvar
qsub 4dvar_optimizer.run
```

### Configuration variables

All variables have defaults and can be overridden at submit time with
`qsub -v`:

| Variable | Default | Description |
|---|---|---|
| `MAX_ITER` | 20 | Maximum number of outer-loop iterations |
| `NLAT` | 46 | Latitude grid points of control variable (4° grid) |
| `NLON` | 72 | Longitude grid points of control variable (5° grid) |
| `SIGMA_B` | 0.5 | Background error standard deviation |
| `GTOL` | 1e-5 | Convergence: max absolute gradient threshold |
| `FTOL` | 1e-8 | Convergence: relative cost change threshold |
| `RESTART` | false | Set `true` to resume from an existing `4dvar_state.npz` |
| `PYTHON_ENV` | `/nobackup/ksuselj1/envs/gchp_4dvar` | Path to Python venv/conda env |

Examples:
```bash
# Fresh start, 10 iterations
qsub -v MAX_ITER=10 4dvar_optimizer.run

# Resume a previous run
qsub -v RESTART=true,MAX_ITER=10 4dvar_optimizer.run

# Coarser grid, stronger background
qsub -v NLAT=24,NLON=48,SIGMA_B=1.0 4dvar_optimizer.run
```

### Walltime guide

Each iteration takes approximately 3 hours (1.5 h forward + 1.5 h adjoint)
on 24 cores. The default walltime of 72 hours accommodates ~20 iterations.
Adjust `#PBS -l walltime` and `MAX_ITER` together.

## Cost function

$$J = J_\text{obs} + J_b = \frac{1}{2}\sum_i \frac{(H(\sigma)_i - y_i)^2}{\sigma_{o,i}^2} + \frac{1}{2}\frac{\|\sigma - 1\|^2}{\sigma_b^2}$$

- **H(σ)**: XCO2 column average from the forward simulation, co-located with OCO-2
- **y**: OCO-2 XCO2 observations
- **σ**: surface flux scaling factors (control variable, initialized to 1)
- **σ_b**: background error standard deviation (default 0.5)

The gradient ∂J/∂σ comes directly from `SurfaceFluxAdj_CO2` in the adjoint
output plus the analytic background gradient.

## Optimizer state

`4dvar_state.npz` stores the full L-BFGS-B state between iterations:

| Field | Description |
|---|---|
| `x_next` | Control vector for the next iteration |
| `x_prev` | Control vector from the previous iteration |
| `g_prev` | Gradient from the previous iteration |
| `J_prev` | Cost from the previous iteration |
| `s_hist` | L-BFGS-B displacement history (m × n) |
| `y_hist` | L-BFGS-B gradient-change history (m × n) |
| `m_used` | Number of valid history pairs |
| `iteration` | Current iteration counter |

Per-iteration snapshots are saved to `4dvar_output/4dvar_state_iter_NNN.npz`.

## Output files

| File | Description |
|---|---|
| `4dvar_design/4dvar_state.npz` | Final optimizer state |
| `4dvar_design/4dvar_output/4dvar_state_iter_NNN.npz` | Per-iteration state snapshots |
| `4dvar_design/forcing_files/J_value.txt` | Observation cost from the last iteration |
| `4dvar_design/forward_run/gchp.log` | Forward GCHP log (last iteration) |
| `4dvar_design/adjoint/gchp.log` | Adjoint GCHP log (last iteration) |

## Other files in this repository

| File | Description |
|---|---|
| `oco2_Baker.f` | Fortran satellite observation operator (original code from J. Liu) |
| `co2_sat_compare_monthly.py` | Python wrapper that applies the OCO-2 averaging kernel operator: x̂ = x_a + A_k (x_m − x_a) |
| `co2_biases.py` | CO2 bias diagnostics |
| `GCHP_LETKF_scale_factor.nc4` | Example emission scaling factor file (NOx LETKF, for reference) |
| `HEMCO_Config.rc` | Example HEMCO config showing how to apply an emission scale factor |
| `ExtData.rc` | Example ExtData.rc entry for registering a scale factor field with MAPL |
| `HISTORY.rc` | Example HISTORY.rc with `sat_track` collection (output at OCO-2 time/location) |
| `adjoint_forcing_description.pdf` | Technical description of the adjoint forcing computation |

### Data locations on Pleiades

| Data | Path |
|---|---|
| Monthly OCO-2 satellite track files | `/nobackup/lwu5/Carbon/sat_track/` |
| OCO-2 Level-2 observations (Baker B11.2) | `/nobackupp17/jliu7/OCO2/BAKER/B11.2_L2#/` |

## Troubleshooting

**Forward/adjoint GCHP fails** — check `4dvar_design/forward_run/gchp.log` or
`4dvar_design/adjoint/gchp.log`. The optimizer stops immediately on any
non-zero GCHP exit code.

**`sigma_co2.nc4` not found by HEMCO** — verify the `../sigma_co2.nc4` relative
path in both `HEMCO_Config.rc` files resolves correctly from the run directory,
and that the matching `Sigma_CO2` entry exists in `ExtData.rc`.

**Wrong gradient shape** — `SurfaceFluxAdj_CO2` shape must match `NLAT × NLON`.
Check `Adjoint.grid_label` in `control_files/adjoint/HISTORY.rc`.

**Resume fails to find state file** — ensure you pass `RESTART=true` *and* that
`4dvar_design/4dvar_state.npz` exists from a previous run.
