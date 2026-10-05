"""
s150 -- Fix A for LD_72 polytope conditioning.

For every cell:
  1. Load extracted model + apply_context_bounds (Frame A).
  2. Plain LP FVA (loops allowed) via Gurobi.
  3. Tighten every reaction's bounds to [min - 1e-6, max + 1e-6].
  4. For reactions with FVA range < 1e-4, fix at midpoint.
  5. Count reactions falling in [1e-6, 1e-4) for the amendment.

Writes:
  v7_models_tight/trans_rfastcormics_{cell}.mat  (tightened copy)
  s150_fix_A_counts.csv                           (per-cell count of [1e-6, 1e-4) rxns)
  s150_fix_A_LD_72_fixed_rxns.csv                 (list of LD_72 reactions fixed at midpoint + subsystem)
"""
import sys, h5py, numpy as np, pandas as pd, shutil
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB
from pathlib import Path

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds

V7_MODELS = REPO / "reviewer_packet_v4_slim/v7_models"
V7_TIGHT  = REPO / "reviewer_packet_v4_slim/v7_models_tight"
V7_TIGHT.mkdir(exist_ok=True)
BIOMASS_RXN = "MAR00021"
MARGIN = 1e-4          # ±1e-6 was too tight for HiGHS to re-solve; 1e-4 is still small vs the fc envelope scale
FIX_THRESHOLD = 1e-4
CELLS = [f"{c}_{h}" for c in ("BL","D","LD") for h in (24,48,72,96,120)]
TARGET_CELLS = sys.argv[1:] if len(sys.argv) > 1 else CELLS


def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except: return ""
    return str(x)


def load_model(fp):
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        try:
            subs = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception:
            subs = [""] * len(rxns)
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        S = sp.csc_matrix((data, ir, jc), shape=(int(ir.max())+1, len(jc)-1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, subs=subs, S=S, lb=lb, ub=ub)


def fva_and_tighten(cell):
    cond, hpf = cell.split("_"); hpf = int(hpf)
    src = V7_MODELS / f"trans_rfastcormics_{cell}.mat"
    dst = V7_TIGHT / f"trans_rfastcormics_{cell}.mat"
    m = load_model(src)
    Snnz = np.diff(m["S"].tocsc().indptr); is_ex = Snnz == 1
    lb, ub, idx, is_ex = apply_context_bounds(m["lb"], m["ub"], m["rxns"], is_ex, cond, hpf)

    # Gurobi FVA model
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    mo = gp.Model("fva", env=env)
    mo.setParam("Threads", 4); mo.setParam("OutputFlag", 0)
    v = mo.addMVar(len(m["rxns"]), lb=lb, ub=ub)
    mo.addMConstr(m["S"].tocsr(), v, "=", np.zeros(m["S"].shape[0]))
    mo.setObjective(v[idx[BIOMASS_RXN]], GRB.MAXIMIZE); mo.optimize()
    bmax = mo.ObjVal
    v[idx[BIOMASS_RXN]].LB = 0.9 * bmax

    mn = np.zeros(len(m["rxns"])); mx = np.zeros(len(m["rxns"]))
    for j in range(len(m["rxns"])):
        mo.setObjective(v[j], GRB.MINIMIZE); mo.optimize()
        mn[j] = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
        mo.setObjective(v[j], GRB.MAXIMIZE); mo.optimize()
        mx[j] = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
    rng = mx - mn

    # tighten with margin on EVERY reaction (do NOT fix at midpoint).
    # Fixing coupled reactions at independent midpoints violates mass balance:
    # mn[j] and mx[j] come from independent LPs and may not be simultaneously
    # achievable. The ±MARGIN window gives the LP enough wiggle room to satisfy
    # coupling while still collapsing the slit directions MVE can't resolve.
    new_lb = mn - MARGIN
    new_ub = mx + MARGIN
    thin = rng < FIX_THRESHOLD   # kept for reporting only; no fixing applied
    # Keep biomass floor consistent with the FVA regime. Every other reaction was
    # tightened to [mn, mx] computed under biomass >= 0.9*bmax, so dropping biomass
    # back to pre-Frame-A (lb=0) makes the LP infeasible: tightened flux ranges
    # require biomass flux the relaxed biomass bound can't support.
    # Preprocess runs with --skip_apply_bounds, so what we save here is what the
    # sampler sees.
    new_lb[idx[BIOMASS_RXN]] = 0.9 * bmax
    new_ub[idx[BIOMASS_RXN]] = ub[idx[BIOMASS_RXN]]

    # count ranges in [1e-6, 1e-4)
    n_near_fixed = int(np.sum((rng >= 1e-6) & (rng < FIX_THRESHOLD)))
    n_fixed_below_1e6 = int(np.sum(rng < 1e-6))

    # write tightened copy
    shutil.copyfile(src, dst)
    with h5py.File(dst, "r+") as f:
        del f["parent/lb"]; f["parent/lb"] = new_lb.reshape(-1, 1)
        del f["parent/ub"]; f["parent/ub"] = new_ub.reshape(-1, 1)

    # if LD_72, also save the fixed reactions list
    if cell == "LD_72":
        rows = []
        for j in np.where(thin)[0]:
            rows.append(dict(rxn=m["rxns"][j], subsystem=m["subs"][j],
                             min=mn[j], max=mx[j], range=rng[j]))
        dfx = pd.DataFrame(rows).sort_values("range")
        dfx.to_csv(REPO / "reviewer_packet_v4_slim/s150_fix_A_LD_72_fixed_rxns.csv", index=False)
        print(f"  [LD_72] wrote s150_fix_A_LD_72_fixed_rxns.csv ({len(dfx)} fixed)")

    return n_near_fixed, n_fixed_below_1e6, int(np.sum(thin)), len(m["rxns"]), bmax


results = []
for cell in TARGET_CELLS:
    print(f"[{cell}] fva + tighten …", flush=True)
    n_near, n_sub1e6, n_fixed, n_total, bmax = fva_and_tighten(cell)
    results.append(dict(cell=cell, n_total=n_total,
                        n_sub_1e6=n_sub1e6,
                        n_in_1e6_1e4=n_near,
                        n_fixed_total=n_fixed,
                        biomass_max=bmax))
    print(f"  {cell}: total={n_total}  <1e-6={n_sub1e6}  [1e-6,1e-4)={n_near}  fixed={n_fixed}  bmax={bmax:.4f}", flush=True)

df = pd.DataFrame(results)
df.to_csv(REPO / "reviewer_packet_v4_slim/s150_fix_A_counts.csv", index=False)
print(f"\nwrote s150_fix_A_counts.csv")
print(df.to_string(index=False))
