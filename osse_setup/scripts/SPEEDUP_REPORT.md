# Speeding up the GCHP 14.5.3 CO2 4D-Var — report (draft, 2026-10-15)

Status (2026-10-15 17:30): steps 1–3 done and verified. All test runs are 1 month (2016-01),
one 4D-Var iteration, C24, 48 cores (2 × bro_ele), OCO-2 track forcing, E = 1.3. Timings below
are from runs with a healthy `/nobackupp27` (the filesystem failed over between ≈ 10:00 and 16:00
on 2026-10-15; runs from that period are not used for timing).

## 1. Where the time goes (measured)

Baseline: `exps/E1.3_const_sat`, C24, 48 cores (2 × bro_ele), 3-month window, binary
`gchp.ref_20261014`. One iteration ≈ 3 h 20 min: forward 71 min, adjoint 109 min, rest
(forcing, optimizer) ≈ 20 min.

MAPL timer profile of one adjoint (iteration 4, 87 min of GCHP time):

| Component | Time | Share |
|---|---|---|
| Physics operators (advection, convection, PBL, emissions adjoints) | ≈ 25 s + DYNAMICS 60 s | < 2 % |
| ExtData (met, emissions input) | 18 min | 17 % |
| GCHPchem RUN, exclusive (checkpoint read, Push/Pop state, forcing load) | 71 min | ≈ 70 % |
| DO_CHEM exclusive | 12 min | ≈ 12 % |

The model is I/O-bound. Every 30-min step writes (forward) and reads (adjoint) one internal
checkpoint (45 MB, 41 variables) and one import checkpoint (4 MB): about 210 GB written and
read again per 3-month iteration, all through `/nobackupp27`.

Composition of the internal checkpoint: 15 `SPC_*` species fields (66 %), the full-chemistry
leftovers `KPPHvalue`, `AeroH2O_SNA`, `H2O2AfterChem`, `SO2AfterChem`, `ORVCSESQ`,
`STATE_PSC`, plus `DELP_DRY`, `BXHEIGHT` and 2-D fields. A CO2 adjoint step needs only the
2-D fields, `DELP_DRY`, `BXHEIGHT` and the advected species.

## 2. Step 1 — CO2-only control files

- `4dvar_final/control_files_co2/{forward,adjoint}`: made with GEOS-Chem's
  `run/shared/singleCarbonSpecies.sh CO2` on copies of `control_files` (transported species,
  HISTORY, HEMCO_Diagn and ExtData trimmed; flux-consistency and met-offset guards pass).
- Two CO2-only fixes were needed in the **adjoint** directory (found by the spot check):
  the adjoint HEMCO_Config.rc still listed `* SPC_CO/CH4/OCS ./Restarts/...` restart fields
  (HEMCO turns them into ExtData imports → "3 unfulfilled imports"), and the adjoint
  ExtData.rc still listed the same species from `gcchem_initial_restart.nc4`. Both removed
  (backups `*.bak_allspc_20261015`).
- Finding: CO2-only alone does **not** shrink the checkpoint (still 43 MB per step); the
  leftover full-chemistry 3-D fields dominate. Forward throughput is about the same.
- Spot check `exps/spot_allspc_1m` vs `exps/spot_co2_1m` (same binary `gchp.ref_20261014`):
  J 5401.253389 vs 5401.253506 (rel. 2e-8), XCO2 at the soundings max diff 3e-5 ppm, gradient max
  rel. diff 2e-7 (correlation 1.00000000), next σ max diff 4e-8. **Same result** — differences are
  float32 rounding of the different species/diagnostic sets, not physics.
- Wall time: all-species forward 20.0 min, adjoint 29.5 min; CO2-only forward 20.6 min, adjoint
  27.3 min. CO2-only alone saves almost nothing (the run is I/O-bound, and the internal state still
  carries all 40 carbon-simulation variables: 43.1 MB per checkpoint).

## 3. Step 2 — profiling timers and cheap fixes

- New MAPL timers `CKPT_READ` and `CKPT_WRITE` around the adjoint checkpoint read/write
  (GCHP_adj_dev geos-chem `Interfaces/GCHP/Chem_GridCompMod.F90`, uncommitted).
- `ADJ_CO2_SUM` debug allreduces: 6 scalar allreduces per step, negligible → left in.
- Lustre striping: run directories use the NAS default progressive layout; 45 MB files sit
  on one stripe. Not changed — reducing the volume (step 3) is the bigger lever.
- Timer results: see step 3. The timers show the checkpoint write takes 97 % of GCHPchem's run
  time in the forward (473 of 488 s per month) and the checkpoint read ≈ 90 % in the adjoint.

## 4. Step 3 — reduced checkpoints

- New GCHP.rc switch `ADJ_CKPT_REDUCED:` (default 0 = unchanged). With 1 the forward writes
  and the adjoint reads only the 2-D fields, `DELP_DRY`, `BXHEIGHT` and `SPC_<advected>`, via
  a subset ESMF state of references to the internal fields (`Adjoint_CkptSubset`). Forward
  and adjoint numerics are unchanged.
- MAPL needs three state attributes on the subset (`MAPL_GridTypeBits`,
  `MAPL_RestartRequired`, `POSITIVE`); the first two test builds failed without them.
  Note: a MAPL `FAIL` inside a component does not abort GCHP — the 4D-Var script catches it
  by grepping the log.
- Binaries: `gchp.ckpt_20261015` (timers), `gchp.ckpt3_20261015` (timers + working subset).
- Bit-identity: `spot_co2_1m` (old binary) = `ckpt_full_1m` (new binary, switch 0) =
  `ckpt_reduced_1m` (switch 1) — J, XCO2 at all soundings, gradient and next σ **bit-identical**.
- Checkpoint size: internal 43.1 MB → 6.3 MB per step (21 of 40 variables kept); with the 4 MB
  import checkpoint ≈ 47 → 10 MB per step, i.e. ≈ 210 → 45 GB per 3-month iteration.

| 1 month, CO2-only | full checkpoint | reduced checkpoint |
|---|---|---|
| forward wall time | 20.6 min | 10.9 min (−47 %) |
| adjoint wall time | 27.3 min | 18.3 min (−33 %) |
| forward + adjoint | 47.9 min | 29.2 min (−39 %) |
| forward CKPT_WRITE (max over ranks) | 474 s | 181 s |
| adjoint CKPT_READ | (not timed cleanly: outage) | 471–605 s |
| forward ExtData | 229 s | 203 s |

- Expected for the 3-month E1.3 iteration (scaling the 1-month ratios): forward 71 → ≈ 38 min,
  adjoint 109 → ≈ 73 min, iteration ≈ 3 h 20 → ≈ 2 h 10 min.
- Remaining costs after the change: checkpoint read in the adjoint (≈ 8 min/month), ExtData
  (≈ 3–5 min/month), the import checkpoint (AREA is constant and could be dropped; DryPLE/PLE are
  needed), job start-up and the Python forcing/step (≈ 4–6 min per iteration).
- Not applied to the running E1.3_const_sat (it stays on `gchp.ref_20261014`). It can be used for
  new 14.5.3 experiments with `control_files_co2_ckpt_reduced/` (binary `gchp.ckpt3_20261015`);
  per the 2026-10-15 decision the speed work continues on 14.8.0 (the change must be ported).

## 5. Aitken and higher resolution (C48, C90)

### 5.1 Why Aitken was slow on 2026-10-01
rom_ait / mil_ait runs were 8–10× slower than bro_ele. Attributed to I/O to `/nobackupp27`
from the Aitken nodes, but never measured directly. The binary (`-march=core-avx2`) runs on
AMD Rome/Milan, so a rebuild is not required.

### 5.2 How cost scales with resolution
| | C24 | C48 | C90 |
|---|---|---|---|
| Cells (relative) | 1 | 4 | 14 |
| Compute per iteration (relative; ×2 if the time step halves) | 1 | 4–8 | 15–30 |
| Reduced checkpoints per 3-month iteration | ≈ 50 GB | ≈ 200–400 GB | ≈ 1 TB+ |
| Met input with matching resolution (0.5×0.625 instead of 4×5) | – | ≈ 16× more | ≈ 16× more |

At C48/C90 compute grows enough that Aitken's faster cores, 128 cores per node and faster
network pay off — but only if checkpoint and met I/O stop going through a shared filesystem
every step. Otherwise the run is I/O-bound on any machine.

### 5.3 Proposed design
1. **Per-rank checkpoints in node memory.** Each MPI rank writes its own subdomain as a raw
   binary file to `/dev/shm` of its node in the forward and reads it back in the adjoint
   (replacing the MAPL checkpoint write/read in `Adjoint_StateRecord` / `Adjoint_StateRefresh`).
   No gather to writer ranks, no shared filesystem; cost per node stays nearly constant as
   resolution grows; memory spreads over nodes (C90 ≈ 1 TB over 8 Rome nodes = 125 GB/node of
   512 GB). Requires forward and adjoint in the same PBS job with the same node/rank layout
   (true for the 4D-Var script). Fallback if memory is short: checkpoint every k steps and
   recompute in between. Verification: bit-identical gradient against the MAPL checkpoints,
   as for step 3. Side benefit: the 4D-Var no longer depends on `/nobackupp27` health in
   its inner loop.
2. **Met staging.** Copy the window's met to node memory or to the filesystem with the best
   path to Aitken at job start; ask NAS support which filesystem that is once 5.4 gives numbers.
3. **History** limited to what the 4D-Var reads (sat_track, adjoint fields).
4. **Re-validate the adjoint at C48/C90** (FD gradient checks, advection tests): accuracy
   depends on resolution and Courant number; C24 results do not carry over automatically.

### 5.4 Aitken test plan (focus: C48 and C90)
1. C24 baseline on one rom_ait node (full vs reduced checkpoint, CKPT/ExtData timers) plus a
   plain file write/read test to `/nobackupp27` and `/dev/shm` → identifies the 1 October cause.
2. C48 and C90 forward + adjoint (1 month, one iteration) with reduced checkpoints on
   rom_ait, varying node count (C48: 1–2 nodes, C90: 2–8 nodes) → compute scaling and I/O share.
3. Implement per-rank memory checkpoints; repeat 2; check bit-identity.
4. Met staging if ExtData share is still large.
5. Pick node counts per resolution by time and SBU cost (Aitken costs more per node-hour,
   but replaces several Broadwell nodes).

### 5.5 14.5.3 vs 14.8.0
The checkpoint code (`Adjoint_StateRecord/Refresh`) and the 4D-Var structure are the same in
both, so the bottleneck pattern and every fix above apply to both — but the reduced-checkpoint
and per-rank changes must be ported to 14.8.0 separately. Absolute timings will differ:
14.8.0 uses MAPL 2.69 (History/ExtData/checkpoint I/O code changed since 14.5.3), ESMF 8.6.1,
Intel 2023.2.1 and its new advection, so C48/C90 numbers from one version do not transfer
exactly to the other. Decision (2026-10-15): finish the 14.8.0 adjoint validation, port and test the 4D-Var on
14.8.0 (the OSSE satellite sampling still has to be ported), and if it is fine adopt 14.8.0
and do all optimization (reduced and per-rank checkpoints, Aitken, C48/C90 scaling) on
14.8.0. On 14.5.3 only the running checkpoint tests are finished and reported here.

## 6. GCHP 14.8.0 results (2026-10-16)

Reduced-checkpoint option ported to 14.8.0 (same code; binary `gchp_14.8.0_adj/code/build_adj/bin/gchp.ckpt_20261016`).
1-month OSSE (E = 1.3, OCO-2 track), one iteration, C24, `exps/ckpt148_*`, `exps/scan148_n*`:

| run | forward | adjoint | total | gradient vs full checkpoints |
|---|---|---|---|---|
| 2 nodes, full checkpoints | 30.5 min | 31.6 min | 62.2 min | reference |
| 2 nodes, reduced | 12.7 min | 22.9 min | 35.6 min | bit-identical (J too) |
| 2 nodes, reduced (repeat) | 15.8 min | 27.1 min | 42.9 min | bit-identical |
| 1 node, reduced | 22.6 min | 22.6 min | 45.2 min | 7e-8 (summation order) |
| 4 nodes, reduced | 21.2 min | 31.0 min | 52.2 min | 1e-7 (summation order) |

- Reduced checkpoints: -40 to -45 % per iteration, bit-identical.
- C24 sweet spot is 2 nodes; 4 nodes is slower (MPI/I-O overhead); 1 node only ~15-25 % slower (best throughput per node when many experiments run).
- Run-to-run variability ~20 % (shared filesystem load).
- Stopping: E1.3 (14.5.3) has |dJ/J| < 3 % per iteration from iteration 10 and the flux error plateaus from iteration 3-5 -> `--ftol 0.02`; plots every 5 iterations (`--plot-every 5`).
- Recommended C24 OSSE setup on 14.8.0: reduced checkpoints, 2 nodes (1 node for many parallel experiments), --ftol 0.02, --plot-every 5:
  ~1.9 h per 3-month iteration and ~11 iterations = ~21 h per experiment, vs ~53 h now (20 x ~2.7 h).
- New in the driver (opt-in): `PLOT_EVERY`, `FTOL` via `setup_exp.sh --ftol/--plot-every`; finalize-hang watchdog (`run_gchp_148`) for 14.8.0 multi-node runs.
