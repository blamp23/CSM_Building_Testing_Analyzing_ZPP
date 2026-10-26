"""
B2: Reconcile the oxphos contradiction.

Test: on v7 base_120, close ALL exchanges (lb=0), open ONLY the 17 MEDIUM
rows + glucose at -1 + O2 at -1000. Max ATP.

If 31.5 — the 2.0 result in s261 came from a different list. Show which
rows differ.
If 2.0 — ox-phos needs something non-MEDIUM. Open each closed exchange
singly on top and find which single additions restore ATP to ~31.5.
"""
import sys
from pathlib import Path
import numpy as np, scipy.sparse as sp, pandas as pd
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from apply_context_bounds import MEDIUM
from s250_atp_sweep import load_model, V7_BASE_DIR, ATP_RXN, O2_EX

OUT = REPO / "grace_v4/results/atp_sweep_v7"
GLU = "MAR09034"

base = load_model(V7_BASE_DIR / "baked_ocr_120_v7.mat", "m")
S = base["S"]; n_met, n_rxn = S.shape
is_ex = (np.diff(S.indptr) == 1)
idx = {r:i for i,r in enumerate(base["rxns"])}

def atp(extra_open):
    """Close all exchanges; open MEDIUM + glucose at -1 + O2 at -1000 +
    the given extra set at -1000 for uptake. Max ATP."""
    lb = base["lb"].copy(); ub = base["ub"].copy()
    for j in np.where(is_ex)[0]:
        lb[j] = 0.0; ub[j] = 1000.0
    for rid in MEDIUM:
        if rid in idx: lb[idx[rid]] = -1000.0
    for rid in extra_open:
        if rid in idx: lb[idx[rid]] = -1000.0
    if GLU in idx:
        lb[idx[GLU]] = -1.0; ub[idx[GLU]] = -1.0
    if O2_EX in idx:
        lb[idx[O2_EX]] = -1000.0; ub[idx[O2_EX]] = 0.0
    if "MAR00021" in idx:
        j = idx["MAR00021"]; lb[j] = 0.0; ub[j] = 1000.0
    ub[idx[ATP_RXN]] = 1e9
    env = gp.Env(empty=True); env.setParam("OutputFlag",0); env.start()
    gm = gp.Model("", env=env); v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val

# Step 1: MEDIUM only (no extra)
a0 = atp(set())
print(f"[B2.1] MEDIUM + glucose -1 + O2 -1000, nothing else: ATP = {a0:.3f}", flush=True)

if a0 > 30:
    print(f"       ox-phos works with MEDIUM alone. Need to re-check s261's primer baseline — probably a bounds-setup bug there.")
else:
    print(f"       ox-phos won't start on MEDIUM alone. Finding which single exchange restores to ~31.5...", flush=True)
    all_ex = [base["rxns"][j] for j in np.where(is_ex)[0]]
    all_ex_set = set(all_ex)
    closed = [r for r in all_ex_set - MEDIUM]
    print(f"       scanning {len(closed)} closed exchanges singly...", flush=True)
    rows = []
    for i, rid in enumerate(closed):
        a = atp({rid})
        if a > a0 + 0.5:
            col = S.getcol(idx[rid]).tocoo()
            met_i = int(col.row[0])
            rows.append(dict(mar=rid, met=base["metNames"][met_i],
                             formula=base["metFormulas"][met_i], atp=a))
    df = pd.DataFrame(rows).sort_values("atp", ascending=False)
    print(f"\n[B2.2] exchanges that singly restore ATP above {a0:.2f}+0.5:")
    print(df.head(50).to_string(index=False))
    df.to_csv(OUT / "oxphos_single_restorers.csv", index=False)
    print(f"\n[write] oxphos_single_restorers.csv  n_candidates={len(df)}")
