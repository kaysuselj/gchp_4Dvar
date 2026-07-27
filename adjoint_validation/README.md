# Adjoint validation: finite-difference vs. adjoint gradient

This directory holds the harness that verifies the GCHP CO2 adjoint gradient
against **finite differences (FD)** of the forward model, on the native
cubed-sphere grid. It is the tool used to decide whether the adjoint gradient
is trustworthy before running a production 4D-Var optimization
(see `../README_4dvar_osse.md`).

## What is being tested

The 4D-Var control variable is a per-cell scaling factor **σ** applied to the
terrestrial CO2 flux (HEMCO scale factor 750, `EmisCO2_NetTerrExch`). The cost
function is

```
J(σ) = (1/2) Σ_obs ( (model − obs) / σ_obs )²
```

The adjoint model returns the gradient **g = dJ/dσ** for every grid cell in a
single backward run. Because CO2 is a linear tracer, `J` is *exactly quadratic*
in σ, so a one-sided finite difference recovers the true gradient plus a known
curvature term:

```
[ J(σ+ε) − J(σ) ] / ε  =  g + (ε/2) H
```

Running both `+ε` and `−ε` (central difference) cancels the curvature and gives
the exact analytic gradient for a chosen cell. The test compares that FD value
to the adjoint's `g` for a set of fixed cells and reports slope, correlation,
and per-cell ratios.

A perfect adjoint gives **ratio ≈ 1.0** and **slope ≈ 1.0**, **corr ≈ 1.0**.
Deviations localize the error (e.g. cross-face seam cells, PPM limiter
non-differentiability, convection/PBL adjoint approximations).

## Files

| File | Role |
|------|------|
| `test_adjoint_vs_forward_exp1.run` | Experiment 1 PBS driver (both advection limiters OFF) |
| `test_adjoint_vs_forward_exp2.run` | Experiment 2 PBS driver (production limiters, transport only) |
| `test_adjoint_vs_forward_exp3.run` | Experiment 3 PBS driver (production limiters + convection/PBL) |
| `cells_fixed.txt` | The 13 validation cells (f, j, i, ε, tag), identical across experiments |
| `cs_fd_tools.py` | Native-CS helper: `write-sigma`, `check-p0`, `pick-cells`, `analyze` |
| `co2_adjoint_forcing_osse.py` | Computes `J_obs` and the adjoint forcing files (OSSE-window version) |
| `convert_satelite_tracks_osse.py` | Builds the OCO-2 pseudo-observation sat_track inputs |

`co2_adjoint_forcing_osse.py` and `convert_satelite_tracks_osse.py` here are the
**sub-month-window** variants (they accept `day_start`/`day_end` and are
import-safe) used by the tests; the optimizer at the repo root uses its own
copies. `cs_fd_tools.py` is test-only.

## Pipeline (phases)

Each run script is a single PBS job that walks five phases. Marker files in the
experiment's `results/` directory let a resubmission skip finished work
(`qsub -v RESTART=true ...`).

| Phase | What it does |
|-------|--------------|
| **SETUP** | Clone forward/adjoint run dirs + monthly control templates; apply this experiment's advection/physics settings |
| **P0** | Two 1-day forwards (base vs. unique-per-cell σ); asserts `pert == σ ⊙ base` exactly → proves ExtData delivers the CS σ **unregridded** and the face mapping is correct. Cheap gate that fails fast. |
| **P1** | Base full-window forward → checkpoints → `co2_adjoint_forcing_osse.py` → `J_base` |
| **P2** | Full-window adjoint → native-CS gradient `g(f,j,i)`; `pick-cells` writes `cells.txt` from `cells_fixed.txt` |
| **P3** | One perturbed forward per `cells.txt` entry: each cell at `±ε`, a global `±ε`, and a `+2ε` quadraticity probe at the first cell |
| **P4** | `analyze` → `results/report.txt`, `fd_scatter_cs.png`, `g_map_cs.png` |

**Exit codes:** `0` adjoint validated · `3` constant-factor error · `4`
spatially varying error · `5` incomplete (resubmit with `RESTART=true`).

## The three experiments (worked example)

The experiments turn on one additional source of adjoint error at a time, so the
differences between them attribute the error. All three share: Jan 2016 window
(`2016-01-01 → 2016-02-01`), C24, 24 cores, the same 13 cells, `ε = 0.05`,
`long` queue. Non-local mixing, dry deposition and wet deposition are OFF in all
three.

| | Purpose | `hord_tr` | `kord_tr` | Convection | PBL mixing | Walltime |
|--|---------|-----------|-----------|------------|------------|----------|
| **Exp 1** | Ideal linear baseline | **2** (linear PPM) | **17** (linear remap) | off | off | 18 h |
| **Exp 2** | Production limiters, transport only | **12** (monotone PPM) | **8** (monotone remap) | off | off | 18 h |
| **Exp 3** | Production limiters + physics | **12** | **8** | **on** | **on** | 24 h |

Interpretation:

- **Exp 1** — with both PPM limiters off the advection scheme is perfectly
  linear and differentiable, so the adjoint should be *exact*. This is the
  reference; any error here is a framework bug, not a limiter/physics artifact.
- **Exp 2 vs Exp 1** — isolates the **PPM flux/remap limiter** adjoint error.
  This is the production `hord_tr=12`/`kord_tr=8` transport path that the
  boundary fixes (Fix A/B/C in the fvdycore `adjoint_fixes_GCHP` branch) target.
- **Exp 3 vs Exp 2** — isolates the extra error contributed by the **convection
  and PBL-mixing adjoints**.

### Validation cells

`cells_fixed.txt` lists 13 land cells spread over the continents (3 NA, 2 SA,
2 AF, 2 EU, 3 AS, 1 AU). Cell **AS3** (`f=1 j=23 i=20`, N China) is intentionally
placed on the face-1 top-edge **seam** to exercise the cross-face boundary
adjoint (Fix B). σ scales the terrestrial flux, so the gradient is exactly zero
over ocean — cells are land-only by construction.

Format: `cell <f> <j> <i> <eps> <tag>` (one per line, `#` comments allowed).

## Running

```bash
# Submit each experiment FROM INSIDE its own run directory so that
# PBS_O_WORKDIR is correct.  (The scripts also hardcode an absolute TEST_DIR,
# so they are robust to the submission directory — see the note below.)
cd <run_dir_exp1> && qsub test_adjoint_vs_forward_exp1.run
cd <run_dir_exp2> && qsub test_adjoint_vs_forward_exp2.run
cd <run_dir_exp3> && qsub test_adjoint_vs_forward_exp3.run

# resume an incomplete run (P3 is the long phase):
qsub -v RESTART=true test_adjoint_vs_forward_exp1.run
```

Results land in each run directory's `results/`:

- `report.txt` — per-cell FD vs. adjoint, slope, correlation, verdict
- `fd_scatter_cs.png` — FD (x) vs. adjoint g (y), 1:1 line
- `g_map_cs.png` — native-CS adjoint sensitivity map

### Deployment note (important)

Each `*.run` script sets an **absolute** `TEST_DIR` near the top, e.g.

```bash
ROOT_DIR="/nobackupp27/ksuselj1/gchp_14.5.3_adjoint_surfaceF"
TEST_DIR="${ROOT_DIR}/tests_full_J/test_adjoint_vs_forward_1month_exp1"
```

This was hardened after jobs `24887282/83/84` died at SETUP: they were submitted
from the parent directory, so the old `pwd`-derived `ROOT_DIR` resolved one level
too high and the run-dir clone source did not exist. **If you deploy this
harness to a new location, edit `ROOT_DIR`/`TEST_DIR` in each script** (and place
`cs_fd_tools.py`, `cells_fixed.txt`, and the two helper scripts in `TEST_DIR`).
