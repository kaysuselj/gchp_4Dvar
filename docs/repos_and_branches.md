# Git repos & branches — adjoint / 4D-Var work

_Last verified: 2026-07-31 on Pleiades (`/nobackupp27/ksuselj1/gchp_14.5.3_adjoint_surfaceF`,
reached via `/home5/ksuselj1/nobackup` → `/nobackup/ksuselj1` → `/nobackupp27/ksuselj1`)._

## Repos that hold our work

| # | Repo / component | Path (under `gchp_14.5.3_adjoint_surfaceF/`) | Remote (`origin`) | Branch | HEAD | Backed up? |
|---|---|---|---|---|---|---|
| 1 | **GCHP** (superproject) | `GCHP/` | `github.com/kaysuselj/GCHP` (fork) | `adjoint_fixes_GCHP` | `7af2448` | ⚠️ **2 commits unpushed** |
| 2 | **FVdycoreCubed_GridComp** (submod) | `GCHP/src/GCHP_GridComp/FVdycoreCubed_GridComp` | `github.com/kaysuselj/FVdycoreCubed_GridComp.git` (fork) | `adjoint_fixes_GCHP` | `0992974` | ✅ clean, in sync |
| 3 | **fvdycore** — advection kernel (submod) | `…/FVdycoreCubed_GridComp/fvdycore` | `github.com/geoschem/GFDL_atmos_cubed_sphere.git` (**UPSTREAM, no fork**) | `adjoint_fixes_GCHP` | `25367d8` | ❌ **NOT on any remote** |
| 4 | **geos-chem** (adjoint, submod) | `GCHP/src/GCHP_GridComp/GEOSChem_GridComp/geos-chem` | `github.com/kaysuselj/geos-chem_adjoint.git` (fork) | `adjoint_fixes` | `cf1a58a3a` | ✅ pushed |
| 5 | **gchp_4Dvar** — optimizer/Python/control | `gchp_4Dvar/` | `github.com/kaysuselj/gchp_4Dvar.git` (fork) | `main` | `645f92b` | ⚠️ **14 files uncommitted** |
| 6 | **offline adjoint harness** | `tests_full_J/adjoint_dotprod_test/` | — | — | — | ❌ **not in git at all** |

**The advection adjoint lives in repo #3 (fvdycore):** `tp_core.F90` (horizontal PPM flux),
`fv_tracer2d.F90` (tracer driver / subcycling), `fv_mapz.F90` (vertical remap),
`boundary.F90` + `../tools/fv_mp_mod.F90` (cubed-sphere halo / edge / corner exchange).
Commit `25367d8` = "Fix A/B/C: adjoint boundary stencil + copy_corners_ad".

## Stock dependencies (unmodified — leave at their pinned submodule commits)

HEMCO (`geoschem/HEMCO`, `6ca374a`), MAPL, FMS, ESMA_cmake, and the other GCHP
submodules carry **no local changes**; they come straight from the GCHP submodule pins.

## Is everything committed?  No — 4 gaps, by risk

1. ❌ **fvdycore `25367d8` (the advection adjoint) is committed locally but pushed nowhere,
   and `origin` points at the upstream geoschem repo we can't push to.** Single point of
   failure — if this disk is lost the advection-adjoint work is gone. **Fix:** create a
   personal fork `kaysuselj/GFDL_atmos_cubed_sphere`, add it as a remote, and push
   `adjoint_fixes_GCHP`.
2. ⚠️ **gchp_4Dvar: 14 uncommitted files**, incl. today's production fixes
   (`4dvar_optimizer.osse.unified.run` JOB_SGMT + job-name, `lbfgsb_step_monthly.py`,
   `write_sigma_monthly.py`, control_files) plus untracked scratch run variants.
3. ⚠️ **GCHP superproject: 2 commits unpushed** (`7af2448`, `e5bdcc7`), the `geos-chem`
   submodule pointer is updated-but-uncommitted (` M …/geos-chem`), and there are untracked
   backup blobs (`src.BeforeAdvMod.tar`, `FVdycoreCubed_GridComp.tat`, `src/geos-chem`,
   `ADJOINT_ISCHEMTIME_FIX.md`).
4. ❌ **The offline harness (`tests_full_J/adjoint_dotprod_test/`) is not version-controlled.**
   It is the primary Mac dev vehicle — it should live in a repo so it can be cloned to the
   laptop and synced back.

## Packaging for laptop (Mac) development

For offline advection-adjoint dev on the Mac we only need the **source**, no ExtData:
- fvdycore sources: `tp_core.F90`, `fv_tracer2d.F90`, `fv_mapz.F90`, `boundary.F90`,
  `../tools/fv_mp_mod.F90`, plus the grid/array modules they use.
- the offline harness `tests_full_J/adjoint_dotprod_test/` (driver + stubs + build script).
- Build with `gfortran` (Homebrew `gcc`); no MPI/ESMF/FMS needed for the interior operator.

Cleanest path: get #3 and #6 into pushable repos (fork fvdycore; put the harness under git),
then `git clone` on the Mac. See the dev plan for the tier split
(interior arithmetic = Mac; full 6-face corner exchange = Pleiades).
