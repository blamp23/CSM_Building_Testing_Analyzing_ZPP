"""
B3: Set-gate, the right version.

Under final apply_context_bounds for LD_72 (any cell — pick LD_72):
  - Open together all exchanges that are free-supplied: lb = -1000 AND
    not in the 105 measured set (so neither MEDIUM nor measured).
  - Keep MEDIUM at -1000.
  - Close all 105 measured nutrients (lb = 0, ub = +1000).
  - O2 open at -1000 (already in MEDIUM).
  - Biomass unconstrained.
  - Max MAR03964 with ATP ub = 1e9.
  - If > 0.5, greedy-close those free-supplied rows one at a time
    (largest ATP drop first) until ATP = 0. Record the closed set.

Then repeat with a glucose primer at lb = ub = -1. Report ATP - 31.5
after each close, so the "offending as a 2nd-carbon-source" set
emerges distinct from the "offending as a 1st-carbon-source" set.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from apply_context_bounds import apply_context_bounds, MEDIUM
from s250_atp_sweep import load_model, V7_EXTRACTED, V7_BASE_DIR, ATP_RXN, O2_EX

OUT = REPO / "grace_v4/results/atp_sweep_v7"
GLU = "MAR09034"

cell = "LD_72"; cond = "LD"; hpf = 72
m = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
S = m["S"]; n_met, n_rxn = S.shape
is_ex = (np.diff(S.indptr) == 1)
idx = {r: i for i, r in enumerate(m["rxns"])}

lb_final, ub_final, _, _ = apply_context_bounds(m["lb"].copy(), m["ub"].copy(),
                                               m["rxns"], is_ex, cond, hpf)

phase2 = pd.read_csv(REPO / "grace_v4/data/phase2_exchange_bounds_long.csv")
measured_105 = set(phase2.ex_rxn.unique())

# Free-supplied set per user spec: lb = -1000 AND not in 105 AND not MEDIUM.
free_supplied = []
for j in np.where(is_ex)[0]:
    rid = m["rxns"][j]
    if abs(lb_final[j] - (-1000.0)) < 1e-6 and rid not in measured_105 and rid not in MEDIUM:
        free_supplied.append(rid)
free_supplied = sorted(free_supplied)
print(f"[info] free-supplied set (lb=-1000, not MEDIUM, not in-105): {len(free_supplied)}", flush=True)
for rid in free_supplied:
    j = idx[rid]
    col = S.getcol(j).tocoo(); met_i = int(col.row[0])
    print(f"   {rid}  {m['metNames'][met_i]:30s}  formula={m['metFormulas'][met_i]}", flush=True)


def gate_lp(free_open, glucose_primer=False):
    """Build the gate LP:
      - Close all exchanges (lb = 0, ub = 1000).
      - Open MEDIUM at -1000.
      - Open each rid in free_open at -1000.
      - Close all 105 measured (already closed by default).
      - If glucose_primer: lb[glucose] = ub = -1.
      - Biomass free, ATP ub = 1e9.
    """
    lb = m["lb"].copy(); ub = m["ub"].copy()
    for j in np.where(is_ex)[0]:
        lb[j] = 0.0; ub[j] = 1000.0
    for rid in MEDIUM:
        if rid in idx: lb[idx[rid]] = -1000.0
    for rid in free_open:
        if rid in idx: lb[idx[rid]] = -1000.0
    # O2 cap enforced
    if O2_EX in idx:
        lb[idx[O2_EX]] = -1000.0; ub[idx[O2_EX]] = 0.0
    if glucose_primer and GLU in idx:
        lb[idx[GLU]] = -1.0; ub[idx[GLU]] = -1.0
    if "MAR00021" in idx:
        j = idx["MAR00021"]; lb[j] = 0.0; ub[j] = 1000.0
    ub[idx[ATP_RXN]] = 1e9
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("", env=env); v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


# ---- Variant 1: no primer, ATP alone from the free-supplied set ------------
print("\n=== B3a: no glucose primer ===", flush=True)
all_open = set(free_supplied)
baseline = gate_lp(all_open, glucose_primer=False)
print(f"all free-supplied open: ATP = {baseline:.3f}", flush=True)
rows_a = []
if baseline > 0.5:
    cur = baseline
    step = 0
    while cur > 0.5 and all_open:
        step += 1
        # Try closing each row, find the one that drops ATP the most
        best_rid, best_atp = None, cur
        for rid in all_open:
            trial = all_open - {rid}
            v = gate_lp(trial, glucose_primer=False)
            if v < best_atp:
                best_atp = v; best_rid = rid
        if best_rid is None:
            print(f"  step {step}: no row drops ATP; stopping (ATP={cur:.2f})"); break
        all_open.discard(best_rid)
        drop = cur - best_atp
        cur = best_atp
        j = idx[best_rid]
        col = S.getcol(j).tocoo(); met_i = int(col.row[0])
        met_name = m["metNames"][met_i]
        rows_a.append(dict(step=step, closed=best_rid, met_name=met_name,
                           atp_after=cur, drop=drop))
        print(f"  step {step}: close {best_rid} ({met_name:30s})  drop={drop:+8.2f}  ATP→{cur:.2f}",
              flush=True)
else:
    print(f"  set-gate passes with no primer (ATP = {baseline:.2f}).")

# ---- Variant 2: with glucose primer -------------------------------------
print("\n=== B3b: with glucose primer at lb=ub=-1 ===", flush=True)
all_open_b = set(free_supplied)
baseline_b = gate_lp(all_open_b, glucose_primer=True)
baseline_no_free = gate_lp(set(), glucose_primer=True)
excess_b = baseline_b - baseline_no_free
print(f"all free-supplied + glucose: ATP = {baseline_b:.3f}", flush=True)
print(f"glucose-only (MEDIUM + glu, no free-supplied): ATP = {baseline_no_free:.3f}", flush=True)
print(f"excess over glucose-only (textbook): {excess_b:+8.3f}", flush=True)
rows_b = []
if excess_b > 0.5:
    cur = baseline_b
    step = 0
    while cur > baseline_no_free + 0.5 and all_open_b:
        step += 1
        best_rid, best_atp = None, cur
        for rid in all_open_b:
            trial = all_open_b - {rid}
            v = gate_lp(trial, glucose_primer=True)
            if v < best_atp:
                best_atp = v; best_rid = rid
        if best_rid is None:
            print(f"  step {step}: no row drops ATP; stopping (ATP={cur:.2f})"); break
        all_open_b.discard(best_rid)
        drop = cur - best_atp
        cur = best_atp
        excess = cur - baseline_no_free
        j = idx[best_rid]
        col = S.getcol(j).tocoo(); met_i = int(col.row[0])
        met_name = m["metNames"][met_i]
        rows_b.append(dict(step=step, closed=best_rid, met_name=met_name,
                           atp_after=cur, drop=drop, excess_over_textbook=excess))
        print(f"  step {step}: close {best_rid} ({met_name:30s})  drop={drop:+8.2f}  ATP→{cur:.2f}  excess={excess:+8.2f}",
              flush=True)
else:
    print(f"  set-gate passes with primer (excess = {excess_b:.2f}).")

pd.DataFrame(rows_a).to_csv(OUT / "set_gate_LD72_no_primer.csv", index=False)
pd.DataFrame(rows_b).to_csv(OUT / "set_gate_LD72_glucose_primer.csv", index=False)
print(f"\n[write] set_gate_LD72_no_primer.csv  rows={len(rows_a)}")
print(f"[write] set_gate_LD72_glucose_primer.csv  rows={len(rows_b)}")
print(f"\n=== Offending set under no primer: ===")
print([r["closed"] for r in rows_a])
print(f"\n=== Offending set under glucose primer: ===")
print([r["closed"] for r in rows_b])
