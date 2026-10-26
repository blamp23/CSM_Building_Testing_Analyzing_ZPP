"""
C3: Tier 0 biomass_max + top-3 exchange binders + 15-pair 2×2
     under OLD vs CORRECTED bounds.

OLD = β-carotene / inositol at lb=-1000 (pre-§11.19 behaviour)
NEW = β-carotene / inositol at demand bounds from s266

Uses apply_context_bounds for NEW. For OLD, we post-apply override lb=-1000
on the 4 free-supplied rows after calling apply_context_bounds, then re-run.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
import importlib
import apply_context_bounds as acb; importlib.reload(acb)
from s250_atp_sweep import load_model, V7_EXTRACTED, V7_BASE_DIR, CELLS, HPFS

OUT = REPO / "grace_v4/results/atp_sweep_v7"
BIOMASS = "MAR00021"
CONDS = ("BL", "D", "LD")

FREE = ["MAR09152","MAR09154","MAR09276","MAR09361"]


def fba(model, lb, ub, idx, with_duals=False):
    S = model["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag",0); env.start()
    gm = gp.Model("", env=env); v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[BIOMASS]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm.optimize()
    if gm.Status != GRB.OPTIMAL:
        gm.dispose(); env.dispose()
        return (np.nan, np.zeros(n_rxn)) if with_duals else np.nan
    bm = gm.ObjVal
    rcs = np.array([v[j].RC for j in range(n_rxn)]) if with_duals else None
    gm.dispose(); env.dispose()
    return (bm, rcs) if with_duals else bm


def top3_duals(model, lb, ub, idx, is_ex, n_top=3):
    bm, rcs = fba(model, lb, ub, idx, with_duals=True)
    hits = []
    for j in range(len(model["rxns"])):
        if not is_ex[j]: continue
        if model["rxns"][j] == BIOMASS: continue
        rc = rcs[j]
        if abs(rc) < 1e-9: continue
        if rc > 0:
            if lb[j] == 0.0 or abs(lb[j]) < 1e-12: continue
            hits.append((abs(rc * lb[j]), model["rxns"][j], "lb", lb[j], rc))
        else:
            if abs(ub[j]) < 1e-12: continue
            hits.append((abs(rc * ub[j]), model["rxns"][j], "ub", ub[j], rc))
    hits.sort(reverse=True)
    return bm, hits[:n_top]


def apply_old_bounds(model, cond, hpf):
    """Apply current apply_context_bounds then override free-supplied to -1000."""
    is_ex = (np.diff(model["S"].indptr) == 1)
    lb, ub, idx, _ = acb.apply_context_bounds(model["lb"].copy(), model["ub"].copy(),
                                              model["rxns"], is_ex, cond, hpf,
                                              skip_energy_gate=True)
    for rid in FREE:
        if rid in idx: lb[idx[rid]] = -1000.0
    return lb, ub, idx, is_ex


def apply_new_bounds(model, cond, hpf):
    is_ex = (np.diff(model["S"].indptr) == 1)
    lb, ub, idx, _ = acb.apply_context_bounds(model["lb"].copy(), model["ub"].copy(),
                                              model["rxns"], is_ex, cond, hpf,
                                              S=model["S"])  # gate on
    return lb, ub, idx, is_ex


# -------- PART A: per-cell biomass_max + top-3 ---------------------------
print("=== C3-A: biomass_max + top-3 exchange binders, OLD vs NEW ===\n", flush=True)
rows_a = []
model_cache = {}
def get_cell(cell):
    if cell not in model_cache:
        model_cache[cell] = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    return model_cache[cell]

for cell in CELLS:
    cond, h = cell.split("_"); h = int(h)
    m = get_cell(cell)
    lb_o, ub_o, idx, is_ex = apply_old_bounds(m, cond, h)
    bm_o, hits_o = top3_duals(m, lb_o, ub_o, idx, is_ex)
    lb_n, ub_n, _, _ = apply_new_bounds(m, cond, h)
    bm_n, hits_n = top3_duals(m, lb_n, ub_n, idx, is_ex)
    rec = dict(cell=cell, cond=cond, hpf=h,
               bm_old=bm_o, bm_new=bm_n, bm_delta=bm_n-bm_o,
               bm_pct_change=100.0*(bm_n-bm_o)/bm_o if bm_o else np.nan)
    for k in range(3):
        for tag, hits in [("old", hits_o), ("new", hits_n)]:
            if k < len(hits):
                c0,r0,s0,b0,rc0 = hits[k]
                rec[f"{tag}_top{k+1}"] = f"{r0}({s0}={b0:.2f}, rc={rc0:.3f}, contrib={c0:.2f})"
            else:
                rec[f"{tag}_top{k+1}"] = ""
    rows_a.append(rec)
    print(f"  {cell:10s} bm old={bm_o:.3f}  new={bm_n:.3f}  Δ={bm_n-bm_o:+.3f}  "
          f"old_top1={rec['old_top1']}  new_top1={rec['new_top1']}", flush=True)

tier0_a = pd.DataFrame(rows_a)
tier0_a.to_csv(OUT / "tier0_C3_biomass_max_before_after.csv", index=False)
print(f"\n[write] tier0_C3_biomass_max_before_after.csv", flush=True)


# -------- PART B: 15 F1 pairs 2x2 ---------------------------------------
print("\n=== C3-B: 15 F1 pairs 2×2 topology/envelope decomposition ===\n", flush=True)
# For each F1 pair (A, B same hpf, different cond):
# corners:  (topology_A, env_A), (topology_A, env_B), (topology_B, env_A), (topology_B, env_B)
# topology_X = model X, env_Y = apply bounds with cond_Y.
# We compute each corner under OLD and NEW.
F1_PAIRS = []
for h in HPFS:
    for a, b in [("BL","D"), ("BL","LD"), ("D","LD")]:
        F1_PAIRS.append((f"{a}_{h}", f"{b}_{h}"))

rows_b = []
for cA, cB in F1_PAIRS:
    condA, h = cA.split("_"); h = int(h); condB = cB.split("_")[0]
    mA, mB = get_cell(cA), get_cell(cB)
    def corner(model, env_cond, old_or_new):
        if old_or_new == "old":
            lb, ub, idx, _ = apply_old_bounds(model, env_cond, h)
        else:
            lb, ub, idx, _ = apply_new_bounds(model, env_cond, h)
        return fba(model, lb, ub, idx)
    rec = dict(pair=f"{cA}_vs_{cB}", hpf=h, condA=condA, condB=condB)
    for tag in ["old", "new"]:
        AA = corner(mA, condA, tag)   # A model, A env
        AB = corner(mA, condB, tag)   # A model, B env
        BA = corner(mB, condA, tag)   # B model, A env
        BB = corner(mB, condB, tag)   # B model, B env
        topo = (AA - BA)              # topology effect at A env
        env  = (AA - AB)              # env effect at A model
        inter = AA - AB - BA + BB
        total = AA - BB
        rec[f"AA_{tag}"] = AA
        rec[f"AB_{tag}"] = AB
        rec[f"BA_{tag}"] = BA
        rec[f"BB_{tag}"] = BB
        rec[f"topology_{tag}"] = topo
        rec[f"envelope_{tag}"] = env
        rec[f"interaction_{tag}"] = inter
        rec[f"total_gap_{tag}"] = total
    rows_b.append(rec)
    print(f"  {rec['pair']:12s}  topology old={rec['topology_old']:+.3f} new={rec['topology_new']:+.3f}  "
          f"envelope old={rec['envelope_old']:+.3f} new={rec['envelope_new']:+.3f}", flush=True)

tier0_b = pd.DataFrame(rows_b)
tier0_b.to_csv(OUT / "tier0_C3_2x2_before_after.csv", index=False)
print(f"\n[write] tier0_C3_2x2_before_after.csv  rows={len(tier0_b)}", flush=True)
