#!/usr/bin/env python3
"""exps1 (2026-11-11): redo iteration 1's L-BFGS step with the sigma_b-scaled first step, without rerunning the model.
Iteration 1 is evaluated at the prior (sigma = 1, J and gradient do not depend on the first-step rule). Its saved state
(4dvar_output_osse/4dvar_state_iter001.npz: x_anchor = 1, g_anchor = the iteration-1 gradient, m_used = 0) is taken, the
first step is recomputed with lbfgsb_step_monthly.lbfgsb_direction(..., sig_b_scale) and the per-field bounds, and the
state file is rewritten as after iteration 1 (iteration = 2, x_next = the new trial point). The driver then resumes at
iteration 2 from phase 'start'.
    redo_first_step.py <exp dir> <frac> <bounds lo:hi,lo:hi> [<sigma_b.npz> | <constant sigma_b>]
"""
import os, sys, shutil
import numpy as np
sys.path.insert(0, "/nobackupp27/ksuselj1/gchp_14.8.0_adj/gchp_4Dvar")
import lbfgsb_step_monthly as L

exp, frac, bounds, sb_arg = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
src = os.path.join(exp, "4dvar_output_osse", "4dvar_state_iter001.npz")
st0 = dict(np.load(src))
assert int(st0["iteration"]) == 2 and int(st0["m_used"]) == 0 and np.all(st0["x_anchor"] == 1.0), "not an iteration-1 state"
g = st0["g_anchor"].astype(np.float64); n = g.size
sb = np.load(sb_arg)["sigma_b"].astype(np.float64).ravel() if sb_arg.endswith(".npz") else np.full(n, float(sb_arg))
assert sb.size == n
pairs = [tuple(float(v) for v in b.split(":")) for b in bounds.split(",")]
blk = n // len(pairs)
lo = np.repeat([p[0] for p in pairs], blk); hi = np.repeat([p[1] for p in pairs], blk)
z = np.zeros_like(st0["s_hist"])
d = L.lbfgsb_direction(g, z, z, 0, frac, sb)
x_raw = 1.0 + d
x_next = np.clip(x_raw, lo, hi)
old = st0["x_next"]
print(f"{exp}: g.d {float(g @ d):.4e} (old {float(st0['g_anchor'] @ st0['d_cur']):.4e}); max|dx|/sigma_b {np.abs(d / sb).max():.3f}")
for k, name in enumerate(["TER", "ocean"][:len(pairs)]):
    s = slice(k * blk, (k + 1) * blk)
    print(f"  {name:5s}: new sigma {x_next[s].min():.3f}..{x_next[s].max():.3f} mean {x_next[s].mean():.4f} | "
          f"old {old[s].min():.3f}..{old[s].max():.3f} mean {old[s].mean():.4f} | clipped {int((x_raw[s] != x_next[s]).sum())}")
st = dict(st0)
st.update(x_next=x_next, d_cur=d, alpha=np.float64(1.0), n_backtrack=np.int64(0))
dst = os.path.join(exp, "4dvar_state_osse_monthly.npz")
shutil.copy2(dst, dst + ".before_redo_first_step")
np.savez(dst, **st)
with open(os.path.join(exp, "4dvar_progress_monthly.txt"), "w") as f:
    f.write("iter=2 phase=start\n")
# history: keep only iteration 1 (the iteration-2 row belongs to the old first step; saved in save_defaultstep_*)
h = os.path.join(exp, "4dvar_output_osse", "4dvar_history.csv")
rows = open(h).read().splitlines()
keep = [rows[0]] + [r for r in rows[1:] if r.split(",")[0] == "1"]
open(h, "w").write("\n".join(keep) + "\n")
print(f"  history: kept {len(keep) - 1} iteration-1 row(s), dropped {len(rows) - len(keep)}")
print(f"  state rewritten: {dst} (iteration 2, phase start); previous state kept as .before_redo_first_step")
