"""
Diagnose the base_120 glucose + O2=-1000 → 2107 ATP result.

Protocol:
 A. Reproduce 2107 (glucose lb=ub=-1, O2 lb=-1000, ATP ub=1e9, 17 MEDIUM
    open at lb=-1000, class-1 vitamins at -1000, everything else closed).
 B. Close each of the 17 MEDIUM exchanges one at a time (lb=0, ub=0),
    report resulting ATP + the met name + its formula.
 C. ATP-from-nothing control: glucose closed, everything else as A — must
    be 0 (confirms no free-energy cycle on its own).
 D. If no single MEDIUM collapses to ~32, run loopless MILP via the
    minimum-flux-magnitude flux-direction indicator formulation.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from s250_atp_sweep import (
    load_model, class1_vitamin_exchanges, prep_bounds,
    V7_BASE_DIR, ATP_RXN, O2_EX, MEDIUM,
)

OUT = REPO / "grace_v4/results/atp_sweep_v7"
GLU_EX = "MAR09034"

def atp_max(model, lb, ub, idx, atp_ub=1e9):
    S = model["S"]; n_met, n_rxn = S.shape
    ub = ub.copy(); ub[idx[ATP_RXN]] = atp_ub
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    status = gm.Status
    flux = np.array([v[j].X for j in range(n_rxn)]) if gm.Status == GRB.OPTIMAL else None
    gm.dispose(); env.dispose()
    return val, status, flux


# -------- load model --------------------------------------------------------
base = load_model(V7_BASE_DIR / "baked_ocr_120_v7.mat", "m")
idx = {r: i for i, r in enumerate(base["rxns"])}
cls1 = class1_vitamin_exchanges()

# -------- Part A: reproduce 2107 -------------------------------------------
lb, ub = prep_bounds(base, cls1)
lb[idx[GLU_EX]] = -1.0; ub[idx[GLU_EX]] = -1.0
lb[idx[O2_EX]]  = -1000.0; ub[idx[O2_EX]] = 0.0
atp_full, st, flux_full = atp_max(base, lb, ub, idx)
print(f"[A] reproduce 2107:    ATP={atp_full:.3f}  status={st}", flush=True)

# -------- Part C: ATP-from-nothing control ---------------------------------
lbC, ubC = prep_bounds(base, cls1)
lbC[idx[GLU_EX]] = 0.0; ubC[idx[GLU_EX]] = 0.0
lbC[idx[O2_EX]]  = -1000.0; ubC[idx[O2_EX]] = 0.0
atp_nothing, st_n, _ = atp_max(base, lbC, ubC, idx)
print(f"[C] ATP-from-nothing:  ATP={atp_nothing:.3f}  status={st_n}", flush=True)

# -------- Part B: close each MEDIUM one at a time --------------------------
mets = base["mets"]; metNames = base["metNames"]; metFormulas = base["metFormulas"]
rows = []
print("\n[B] per-MEDIUM-exchange closure (glucose -1, O2 -1000, ATP ub=1e9):", flush=True)
for mex in sorted(MEDIUM):
    if mex not in idx:
        continue
    lb2 = lb.copy(); ub2 = ub.copy()
    j = idx[mex]
    orig_lb, orig_ub = lb2[j], ub2[j]
    lb2[j] = 0.0; ub2[j] = 0.0
    atp_c, st_c, _ = atp_max(base, lb2, ub2, idx)
    # identify metabolite
    col = base["S"].getcol(j).tocoo()
    met_i = int(col.row[0]) if len(col.row) else -1
    met_name = metNames[met_i] if met_i >= 0 else "?"
    formula  = metFormulas[met_i] if met_i >= 0 else "?"
    rows.append(dict(exchange=mex, met_name=met_name, formula=formula,
                     orig_lb=orig_lb, atp_when_closed=atp_c, status=int(st_c)))
    print(f"  close {mex:10s}  {met_name:35s}  formula={formula:12s}  ATP={atp_c:.3f}", flush=True)

df = pd.DataFrame(rows).sort_values("atp_when_closed")
df["drop_from_full"] = atp_full - df["atp_when_closed"]
df.to_csv(OUT / "excess_atp_diagnosis_per_medium.csv", index=False)
print(f"\n[write] excess_atp_diagnosis_per_medium.csv (rows={len(df)})", flush=True)
print("\n=== ranked by atp_when_closed (lowest = this metabolite was the electron donor) ===")
print(df.to_string())

# -------- Part D: loopless MILP (only if no single MEDIUM collapses) -------
# Threshold: 'collapse' = ATP drops below 100 (textbook ≈ 32)
collapses = df[df.atp_when_closed < 100]
if len(collapses) >= 1:
    print(f"\n[D] skipping loopless MILP: {len(collapses)} MEDIUM exchange(s) alone drop ATP < 100",
          flush=True)
else:
    print("\n[D] no single MEDIUM collapse — running loopless (CFF-style) check", flush=True)
    # CycleFreeFlux: given an FBA solution v*, re-solve for v with same objective
    # value, same signs on non-zero fluxes, min sum |v|. Loops get removed.
    S = base["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("cff", env=env)
    # split v = vp - vn
    vp = gm.addMVar(n_rxn, lb=0.0, ub=np.where(ub > 0, ub, 0.0))
    vn = gm.addMVar(n_rxn, lb=0.0, ub=np.where(lb < 0, -lb, 0.0))
    # Force sign agreement with FBA: if v*>0, vn=0; if v*<0, vp=0; if 0, both 0
    for j in range(n_rxn):
        if flux_full[j] > 1e-7:
            vn[j].ub = 0.0
        elif flux_full[j] < -1e-7:
            vp[j].ub = 0.0
        else:
            vp[j].ub = 0.0; vn[j].ub = 0.0
    # S (vp - vn) = 0
    SS = sp.hstack([S, -S], format="csr")
    gm.addMConstr(SS, gp.MVar.fromlist(list(vp.tolist()) + list(vn.tolist())),
                  "=", np.zeros(n_met))
    # Force ATP objective to match FBA ATP
    gm.addConstr(vp[idx[ATP_RXN]] - vn[idx[ATP_RXN]] == atp_full)
    # min sum |v|
    one = np.ones(2 * n_rxn)
    gm.setMObjective(None, one, 0.0, sense=GRB.MINIMIZE)
    gm.optimize()
    if gm.Status == GRB.OPTIMAL:
        atp_cff = vp[idx[ATP_RXN]].X - vn[idx[ATP_RXN]].X
        print(f"  CFF-style loopless ATP: {atp_cff:.3f}  (if ~32, 2107 was a loop)", flush=True)
    else:
        print(f"  CFF solve failed: status={gm.Status}", flush=True)
    gm.dispose(); env.dispose()
