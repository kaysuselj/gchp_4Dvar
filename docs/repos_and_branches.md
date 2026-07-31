# Git repos & branches — adjoint / 4D-Var work

_Last verified: 2026-07-31 on Pleiades (`/nobackupp27/ksuselj1/gchp_14.5.3_adjoint_surfaceF`,
reached via `/home5/ksuselj1/nobackup` → `/nobackup/ksuselj1` → `/nobackupp27/ksuselj1`)._
_All work is now pushed to `kaysuselj` remotes (backup gaps closed 2026-07-31)._

## Repos that hold our work

| # | Repo / component | Path (under `gchp_14.5.3_adjoint_surfaceF/`) | Remote (`origin`, and extras) | Branch | HEAD | Backed up? |
|---|---|---|---|---|---|---|
| 1 | **GCHP** (superproject) | `GCHP/` | `github.com/kaysuselj/GCHP` (fork) | `adjoint_fixes_GCHP` | `95a5d02` | ✅ pushed |
| 2 | **FVdycoreCubed_GridComp** (submod) | `GCHP/src/GCHP_GridComp/FVdycoreCubed_GridComp` | `github.com/kaysuselj/FVdycoreCubed_GridComp.git` (fork) | `adjoint_fixes_GCHP` | `0992974` | ✅ in sync |
| 3 | **fvdycore** — advection kernel (submod) | `…/FVdycoreCubed_GridComp/fvdycore` | `origin`=upstream `geoschem/GFDL_atmos_cubed_sphere`; **`mine`=`github.com/kaysuselj/GFDL_atmos_cubed_sphere` (fork, created 2026-07-31)** | `adjoint_fixes_GCHP` | `25367d8` | ✅ pushed to `mine` |
| 4 | **geos-chem** (adjoint, submod) | `GCHP/src/GCHP_GridComp/GEOSChem_GridComp/geos-chem` | `github.com/kaysuselj/geos-chem_adjoint.git` (fork) | `adjoint_fixes` | `cf1a58a3a` | ✅ pushed |
| 5 | **gchp_4Dvar** — optimizer/Python/control | `gchp_4Dvar/` | `github.com/kaysuselj/gchp_4Dvar.git` (fork) | `main` | `8c582f4` | ✅ pushed |
| 6 | **offline adjoint harness** | `tests_full_J/adjoint_dotprod_test/` | **`github.com/kaysuselj/advection-adjoint-test` (new repo, private, 2026-07-31)** | `main` | `9ba7b15` | ✅ pushed |

**The advection adjoint lives in repo #3 (fvdycore):** `tp_core.F90` (horizontal PPM flux),
`fv_tracer2d.F90` (tracer driver / subcycling), `fv_mapz.F90` (vertical remap),
`boundary.F90` + `../tools/fv_mp_mod.F90` (cubed-sphere halo / edge / corner exchange).
Commit `25367d8` = "Fix A/B/C: adjoint boundary stencil + copy_corners_ad".
Push future fvdycore work with `git push mine adjoint_fixes_GCHP` (NOT `origin`, which is upstream).

## Stock dependencies (unmodified — leave at their pinned submodule commits)

HEMCO (`geoschem/HEMCO`, `6ca374a`), MAPL, FMS, ESMA_cmake, and the other GCHP
submodules carry **no local changes**; they come straight from the GCHP submodule pins.

## Backup status — all committed & pushed (2026-07-31)

Every repo above is pushed to a `kaysuselj` remote. Notes:
- fvdycore `25367d8` (the advection adjoint) was local-only until 2026-07-31; now on
  fork `kaysuselj/GFDL_atmos_cubed_sphere` via remote `mine`.
- The offline harness is now its own repo `kaysuselj/advection-adjoint-test` (private).
- **Deliberately NOT in git** (left untracked in the GCHP superproject): backup blobs
  `src.BeforeAdvMod.tar`, `FVdycoreCubed_GridComp.tat`, the stray `src/geos-chem` dir, and
  `ADJOINT_ISCHEMTIME_FIX.md`; in gchp_4Dvar: the 715 MB `gchp-env/` (now gitignored),
  and scratch run variants (`4dvar_optimizer.run.v1`, `.devel`, `_long`, `_old_long`,
  `adjoint_determinism_test.run`, `test_monthly_sigma.run`, `control_files/forward_osse_monthly_test/`).
  Review these if any are real work worth keeping.

## Packaging for laptop (Mac) development

For offline advection-adjoint dev on the Mac we only need **source**, no ExtData:
- Clone `kaysuselj/advection-adjoint-test` — it is self-contained (its own `tp_core.F90`
  working copy + stubs + driver). Build with `./build_mac.sh` (`brew install gcc` first).
- Source of truth for the kernel is fvdycore (`kaysuselj/GFDL_atmos_cubed_sphere`,
  `adjoint_fixes_GCHP`); port verified fixes from the harness copy back to fvdycore.
- No MPI/ESMF/FMS needed for the interior operator. Full 6-face MPI exchange = Pleiades.
