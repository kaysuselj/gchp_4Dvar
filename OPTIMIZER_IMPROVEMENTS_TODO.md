# 4D-Var optimizer — improvement backlog

Context: analysis of the completed `3months_osse` monthly run (2026-08-05) showed the
bare unit-step L-BFGS **diverged** — J rose +20% over 20 iterations with 100×–1000×
excursions, sigma blew up to ~29, and the best iterate (iter 8, J≈1.24e4) was passed
through and discarded. See `plot_optimizer_history.py` output and the memory notes.

## Already DONE (2026-08-05) — under test in `3months_osse_improved_optimizer`
- **(a)** driver `EST_MONTH_FWD_S/ADJ_S` defaults 6000 → **1800/1200** (measured C24:
  ~18 min/mo fwd, ~12 min/mo adj). Old value over-budgeted long monthly iterations
  past 24 h and forced constant self-resubmit.
- **(2)** `lbfgsb_step_monthly.py`: deferred **backtracking line search** + **keep-best-
  iterate**, in a pure/unit-tested `lbfgs_backtrack_step()`. Accept a step only if J
  decreased (Armijo `--c1`, default 1e-4); on reject keep the anchor+direction, halve
  `alpha`, and DO NOT update the L-BFGS memory (each backtrack = the next outer
  fwd+adj); a non-descent direction (g·d≥0) resets memory to scaled steepest descent;
  stall after `--max-backtrack` (default 6) halvings → deliver best. Tracks
  `x_best/J_best/it_best`, writes `4dvar_state_best.npz`, and delivers best on
  convergence/stall. `--max-backtrack 0` restores the old pure unit step.
- **Diagnostics:** per-iteration `4dvar_output_osse/4dvar_history.csv` (rho, g·d, alpha,
  reject/reject_reason, n_reject_total, reset_memory, sigma eval vs proposal, bound
  counts, gamma, curvature, is_best, J_best/it_best) + console step ACCEPT/REJECT lines.
- **Plots:** `plot_optimizer_history.py` — 5-panel behaviour figure (cost, best-so-far,
  gradient norms, sigma range, rho) from the per-iteration state snapshots; `--per-
  iteration` emits one fixed-axes frame per step.

NOTE: annual `lbfgsb_step.py` has the diagnostics but NOT the backtracking safeguard yet.

---

## DEFERRED — do NOT start until the improved run is evaluated (user, 2026-08-05)

### 1. Deliver best at MAX_ITER (correctness gap)  [highest priority]
On convergence/stall the step already delivers `x_best`. But if the run stops by hitting
`MAX_ITER` (exit 0 each iter, driver loop ends by count), the driver's final sigma is the
last `x_next` — an **unevaluated proposal**, not the best evaluated iterate. Fix: after the
`while` loop in `4dvar_optimizer.osse.unified.run`, if `4dvar_state_best.npz` exists,
promote it (write final sigma from `x_best`, log `final = best iter N, J=…`). ~5–10 lines.

### 2. Relative gradient tolerance  [high]
`gtol=1e-5` is unreachable (|g|_inf ~ 1e2–1e3), so the run never declares gradient
convergence — it always hits MAX_ITER or stall. Store `g0 = |g|_inf` at iter 1 and stop
when `|g|_inf < gtol_rel * g0` (add `--gtol-rel`, default ~1e-3). Loosen `ftol` (~1e-3–1e-4)
too. Keep the absolute `gtol` as a fallback so nothing regresses.

### 3. Quadratic-interpolation step length instead of pure halving  [high value / cost]
Each backtrack is a full ~30-min fwd+adj, so halving is expensive. Fit a quadratic to
`(J_anchor, g·d, J_trial)` and jump to the predicted minimizer — typically 1 trial instead
of 2–4. Big wall-clock saving given the cost model.

### 4. Global backtrack budget per run
Cap total rejects (`n_reject_total`), not just per-direction, so a pathological run can't
burn unlimited wall. Exit-to-best when the budget is exhausted.

### 5. Mirror the safeguard to annual `lbfgsb_step.py`
It has the diagnostics but not the backtracking — port the identical
`lbfgs_backtrack_step()` for consistency once the monthly version is validated.

### 6. Re-sweep sigma_b  [item (c)]
The old sweep was confounded by divergence. With the trajectory stabilized, re-examine
`sigma_b` (0.2 gave J_b < 1% of J early; to bite earlier it needs to be *smaller*, ~0.05–0.1,
not larger). Use the new rho/history diagnostics to judge.

### 7. Revisit the observation error R  [physical]
`chi2 = 2·J_obs/N_obs ≈ 1.8–2.4` in the real run → residuals ~1.3–1.5σ, i.e. the assumed
obs error R is likely too tight. This is a physical dial affecting the whole obs/background
balance, independent of the optimizer mechanics.
