"""
Follow-up to s255: no single MEDIUM collapsed the 2107. Test whether it is
a thermodynamic loop by:
 (1) pFBA on the same LP — min sum|v| s.t. v[ATP]=2107. Inspect the sparse
     route for reactions with large flux and no external uptake path.
 (2) Pairwise MEDIUM closure for the top candidates — maybe TWO exchanges
     together explain the drop.
 (3) Proper loopless FBA (MILP) — if the loopless max is ~32, 2107 was a
     loop and the stoichiometric O2 cap is a legitimate workaround.
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

base = load_model(V7_BASE_DIR / "baked_ocr_120_v7.mat", "m")
idx = {r: i for i, r in enumerate(base["rxns"])}
cls1 = class1_vitamin_exchanges()
lb, ub = prep_bounds(base, cls1)
lb[idx[GLU_EX]] = -1.0; ub[idx[GLU_EX]] = -1.0
lb[idx[O2_EX]]  = -1000.0; ub[idx[O2_EX]] = 0.0
ub[idx[ATP_RXN]] = 1e9

S = base["S"]; n_met, n_rxn = S.shape

# ---------- Part A: Phase 1 FBA + Phase 2 pFBA ------------------------------
env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
gm = gp.Model("fba", env=env)
v = gm.addMVar(n_rxn, lb=lb, ub=ub)
gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
gm.optimize()
atp_star = gm.ObjVal
print(f"[A1] FBA ATP max = {atp_star:.3f}", flush=True)
gm.dispose()

gm2 = gp.Model("pfba", env=env)
vp = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(ub > 0, ub, 0.0))
vn = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(lb < 0, -lb, 0.0))
SS = sp.hstack([S, -S], format="csr")
gm2.addMConstr(SS, gp.MVar.fromlist(list(vp.tolist()) + list(vn.tolist())),
               "=", np.zeros(n_met))
gm2.addConstr(vp[idx[ATP_RXN]] - vn[idx[ATP_RXN]] == atp_star)
one = np.ones(2 * n_rxn)
gm2.setMObjective(None, one, 0.0, sense=GRB.MINIMIZE)
gm2.optimize()
flux = np.array([vp[j].X - vn[j].X for j in range(n_rxn)])
print(f"[A2] pFBA sum|v| = {gm2.ObjVal:.1f}  (sparse route size: {(np.abs(flux)>1e-6).sum()})",
      flush=True)
# What exchanges are active in pFBA?
is_ex = (np.diff(S.indptr) == 1)
ex_active = [(base["rxns"][j], float(flux[j]), base["metNames"][int(S.getcol(j).tocoo().row[0])])
             for j in np.where(is_ex)[0] if abs(flux[j]) > 1e-6]
ex_active.sort(key=lambda x: abs(x[1]), reverse=True)
print(f"[A3] pFBA active exchanges (top 20 by |v|):", flush=True)
for r, f, mn in ex_active[:20]:
    print(f"    {r:12s}  v={f:+9.3f}  {mn}", flush=True)
gm2.dispose()

# ---------- Part B: pairwise closure of top ATP-sensitive exchanges --------
# From s255: HCO3-, H+, NH3 each drop ~35-40. Try all pairs of top candidates.
print("\n[B] pairwise MEDIUM closure (candidates: HCO3-, H+, NH3, CO2):", flush=True)
candidates = ["MAR09078", "MAR09079", "MAR09073", "MAR09058"]
rows = []
# single closures first
for m1 in candidates:
    lb2 = lb.copy(); ub2 = ub.copy()
    lb2[idx[m1]] = 0.0; ub2[idx[m1]] = 0.0
    gm = gp.Model("", env=env)
    vv = gm.addMVar(n_rxn, lb=lb2, ub=ub2)
    gm.addMConstr(S.tocsr(), vv, "=", np.zeros(n_met))
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    atp_c = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    rows.append(dict(closed=m1, atp=atp_c))
    print(f"  single close {m1}: ATP={atp_c:.3f}", flush=True)
    gm.dispose()
# pairs
for i, m1 in enumerate(candidates):
    for m2 in candidates[i+1:]:
        lb2 = lb.copy(); ub2 = ub.copy()
        lb2[idx[m1]] = 0.0; ub2[idx[m1]] = 0.0
        lb2[idx[m2]] = 0.0; ub2[idx[m2]] = 0.0
        gm = gp.Model("", env=env)
        vv = gm.addMVar(n_rxn, lb=lb2, ub=ub2)
        gm.addMConstr(S.tocsr(), vv, "=", np.zeros(n_met))
        gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
        gm.optimize()
        atp_c = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
        rows.append(dict(closed=f"{m1}+{m2}", atp=atp_c))
        print(f"  pair close  {m1}+{m2}: ATP={atp_c:.3f}", flush=True)
        gm.dispose()
# triples
for i, m1 in enumerate(candidates):
    for j, m2 in enumerate(candidates[i+1:], i+1):
        for m3 in candidates[j+1:]:
            lb2 = lb.copy(); ub2 = ub.copy()
            for mm in (m1, m2, m3):
                lb2[idx[mm]] = 0.0; ub2[idx[mm]] = 0.0
            gm = gp.Model("", env=env)
            vv = gm.addMVar(n_rxn, lb=lb2, ub=ub2)
            gm.addMConstr(S.tocsr(), vv, "=", np.zeros(n_met))
            gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
            gm.optimize()
            atp_c = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
            rows.append(dict(closed=f"{m1}+{m2}+{m3}", atp=atp_c))
            print(f"  triple      {m1}+{m2}+{m3}: ATP={atp_c:.3f}", flush=True)
            gm.dispose()
# full 4
lb2 = lb.copy(); ub2 = ub.copy()
for mm in candidates: lb2[idx[mm]] = 0.0; ub2[idx[mm]] = 0.0
gm = gp.Model("", env=env)
vv = gm.addMVar(n_rxn, lb=lb2, ub=ub2)
gm.addMConstr(S.tocsr(), vv, "=", np.zeros(n_met))
gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
gm.optimize()
print(f"  quad closed: ATP={gm.ObjVal:.3f}", flush=True)
gm.dispose()
pd.DataFrame(rows).to_csv(OUT / "excess_atp_diagnosis_pairs.csv", index=False)

# ---------- Part C: proper loopless FBA (MILP) ------------------------------
# Add binary z_j for each INTERNAL reaction. For each internal j:
#    v_j >= lb_j * z_j and v_j <= ub_j * z_j  (so z=0 => v=0)
#    Then forbid thermodynamic loops: assign chemical potential g_j to each
#    internal metabolite; if z_j=1 and v_j>0, require sum_i S_ij * g_i <= -eps.
# This is Schellenberger's loopless COBRA formulation.
print("\n[C] loopless MILP (ll-FBA):", flush=True)
# Internal reaction mask: NOT a boundary exchange
is_internal = ~is_ex
n_int = int(is_internal.sum())
# For performance: approach ll-FBA via two-level test. First test ATP max
# under indicator-only (z=1 iff |v|>0, no loop constraint) baseline; then
# add loop-cutting on an iteratively-growing set. For now, use the simpler
# "flux must be in the right-null-space of internal S^T" formulation.
# Shortest proper: enforce a chemical-potential vector g for metabolites,
# and for every internal reaction j with |v_j|>0:
#   sign(v_j) * (S^T g)_j < 0
# We linearise with M=1000, delta=1e-6:
#   S^T g >= -M*(1 - z^+) + delta  if  v > 0 (v is in vp)
#   S^T g <=  M*(1 - z^-) - delta  if  v < 0 (v is in vn)
# where z^+/z^- indicator binaries.
Mcap = 1000.0; delta = 1e-6
gm = gp.Model("llfba", env=env)
# positive/negative flux variables
vp = gm.addMVar(n_rxn, lb=0.0, ub=np.where(ub > 0, ub, 0.0))
vn = gm.addMVar(n_rxn, lb=0.0, ub=np.where(lb < 0, -lb, 0.0))
# indicator binaries for internal reactions only
zp = gm.addMVar(n_rxn, vtype=GRB.BINARY)
zn = gm.addMVar(n_rxn, vtype=GRB.BINARY)
# metabolite potentials
g = gm.addMVar(n_met, lb=-Mcap, ub=Mcap)
# Steady state
SS = sp.hstack([S, -S], format="csr")
gm.addMConstr(SS, gp.MVar.fromlist(list(vp.tolist()) + list(vn.tolist())),
              "=", np.zeros(n_met))
# vp <= ub * zp, vn <= (-lb) * zn
for j in range(n_rxn):
    gm.addConstr(vp[j] <= float(np.where(ub[j] > 0, ub[j], 0.0)) * zp[j])
    gm.addConstr(vn[j] <= float(np.where(lb[j] < 0, -lb[j], 0.0)) * zn[j])
    gm.addConstr(zp[j] + zn[j] <= 1)   # exclusivity

# Thermodynamic constraint for internal reactions
int_rxns = [j for j in range(n_rxn) if is_internal[j]]
S_csr = S.tocsr()
for j in int_rxns:
    col = S.getcol(j).tocoo()
    stg = gp.quicksum(float(col.data[k]) * g[int(col.row[k])] for k in range(len(col.data)))
    gm.addConstr(stg <= -delta + Mcap * (1 - zp[j]))
    gm.addConstr(stg >=  delta - Mcap * (1 - zn[j]))

# Objective: max ATP
cvec = np.zeros(n_rxn); cvec[idx[ATP_RXN]] = 1.0
gm.setMObjective(None, np.concatenate([cvec, -cvec]), 0.0, sense=GRB.MAXIMIZE)
gm.setParam("TimeLimit", 300)
gm.setParam("MIPGap", 0.01)
gm.optimize()
if gm.Status in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL):
    atp_ll = vp[idx[ATP_RXN]].X - vn[idx[ATP_RXN]].X
    print(f"  loopless ATP = {atp_ll:.3f}   (status={gm.Status}, gap={gm.MIPGap:.3f})", flush=True)
    print(f"  if ~32: 2107 was a loop; stoichiometric O2 cap is legitimate workaround.", flush=True)
    print(f"  if >1000: not a loop; medium feeds extra oxidisable material.", flush=True)
else:
    print(f"  loopless MILP failed: status={gm.Status}", flush=True)
gm.dispose(); env.dispose()
