"""B3 extended: run the right set-gate across all 15 extracted cells + 5 base
models. Report per-cell offending sets + the union."""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from apply_context_bounds import apply_context_bounds, MEDIUM
from s250_atp_sweep import load_model, V7_EXTRACTED, V7_BASE_DIR, ATP_RXN, O2_EX, CELLS, HPFS, CONDS

OUT = REPO / "grace_v4/results/atp_sweep_v7"
GLU = "MAR09034"

phase2 = pd.read_csv(REPO / "grace_v4/data/phase2_exchange_bounds_long.csv")
MEAS = set(phase2.ex_rxn.unique())

def gate_lp(model, lb_f, ub_f, free_open, glucose_primer, atp_ub=1e9):
    S = model["S"]; n_met, n_rxn = S.shape
    is_ex = (np.diff(S.indptr) == 1)
    idx = {r:i for i,r in enumerate(model["rxns"])}
    lb = model["lb"].copy(); ub = model["ub"].copy()
    for j in np.where(is_ex)[0]:
        lb[j] = 0.0; ub[j] = 1000.0
    for rid in MEDIUM:
        if rid in idx: lb[idx[rid]] = -1000.0
    for rid in free_open:
        if rid in idx: lb[idx[rid]] = -1000.0
    if O2_EX in idx:
        lb[idx[O2_EX]] = -1000.0; ub[idx[O2_EX]] = 0.0
    if glucose_primer and GLU in idx:
        lb[idx[GLU]] = -1.0; ub[idx[GLU]] = -1.0
    if "MAR00021" in idx:
        j = idx["MAR00021"]; lb[j] = 0.0; ub[j] = 1000.0
    ub[idx[ATP_RXN]] = atp_ub
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("", env=env); v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


def run_cell(name, model_path, root_key, cond, hpf):
    print(f"\n[cell] {name}", flush=True)
    m = load_model(model_path, root_key)
    S = m["S"]
    is_ex = (np.diff(S.indptr) == 1)
    idx = {r:i for i,r in enumerate(m["rxns"])}
    lb_f, ub_f, _, _ = apply_context_bounds(m["lb"].copy(), m["ub"].copy(),
                                            m["rxns"], is_ex, cond, hpf)
    # Free-supplied: lb=-1000 ∧ not in 105 ∧ not MEDIUM
    free = []
    for j in np.where(is_ex)[0]:
        rid = m["rxns"][j]
        if abs(lb_f[j] - (-1000.0)) < 1e-6 and rid not in MEAS and rid not in MEDIUM:
            free.append(rid)
    free_set = set(free)
    print(f"  free-supplied set: {len(free_set)} rows: {sorted(free_set)}", flush=True)

    out = dict(cell=name, cond=cond, hpf=hpf, n_free=len(free_set))

    for primer_tag, prim in [("no_primer", False), ("primer", True)]:
        all_open = set(free_set)
        base_atp = gate_lp(m, lb_f, ub_f, all_open, prim)
        base_no_free = gate_lp(m, lb_f, ub_f, set(), prim)
        excess = base_atp - base_no_free
        print(f"  [{primer_tag}] all-open ATP={base_atp:.2f}  textbook={base_no_free:.2f}  excess={excess:+.2f}", flush=True)
        out[f"{primer_tag}_all_open_atp"] = base_atp
        out[f"{primer_tag}_textbook"] = base_no_free
        out[f"{primer_tag}_excess"] = excess
        picks = []
        cur = base_atp
        local_open = set(all_open)
        thresh = 0.5 if not prim else base_no_free + 0.5
        step = 0
        while cur > thresh and local_open and step < 20:
            step += 1
            best_rid, best_atp = None, cur
            for rid in local_open:
                trial = local_open - {rid}
                v = gate_lp(m, lb_f, ub_f, trial, prim)
                if v < best_atp - 1e-6:
                    best_atp = v; best_rid = rid
            if best_rid is None: break
            picks.append((best_rid, cur - best_atp, best_atp))
            local_open.discard(best_rid); cur = best_atp
            j = idx[best_rid]
            col = S.getcol(j).tocoo(); met_i = int(col.row[0])
            print(f"    step {step}: close {best_rid} ({m['metNames'][met_i]:20s}) → ATP={cur:.2f}", flush=True)
        out[f"{primer_tag}_closed"] = ";".join(p[0] for p in picks)
    return out


rows = []
# 15 extracted
for c in CELLS:
    cond, h = c.split("_"); h = int(h)
    rows.append(run_cell(c, V7_EXTRACTED / f"trans_rfastcormics_{c}.mat", "parent", cond, h))
# 5 base
for h in HPFS:
    rows.append(run_cell(f"base_{h}",
                         V7_BASE_DIR / f"baked_ocr_{h}_v7.mat", "m",
                         "LD", h))  # pick LD env to apply bounds (base is the same across cond/hpf except OCR)

df = pd.DataFrame(rows)
df.to_csv(OUT / "set_gate_all_cells.csv", index=False)
print(f"\n[write] set_gate_all_cells.csv  rows={len(df)}", flush=True)

print("\n=== Per-cell offending sets (no primer) ===")
print(df[["cell","no_primer_all_open_atp","no_primer_excess","no_primer_closed"]].to_string(index=False))
print("\n=== Per-cell offending sets (with primer) ===")
print(df[["cell","primer_all_open_atp","primer_textbook","primer_excess","primer_closed"]].to_string(index=False))

union_np = set()
for s in df["no_primer_closed"].str.split(";"): union_np.update(x for x in (s or []) if x)
union_p  = set()
for s in df["primer_closed"].str.split(";"): union_p.update(x for x in (s or []) if x)
print(f"\n[union] offending set no-primer across all cells: {sorted(union_np)}")
print(f"[union] offending set primer across all cells:    {sorted(union_p)}")
