"""
Which reactions carry β-carotene carbon in the FBA biomass_max solution
under OLD bounds? Trace active carotenoid pathway in each of the 4 cells
+ a base reference. Also compare under NEW bounds (should drop to ~zero
flux through the carotenoid pool).

Approach: FBA biomass_max under OLD bounds; pFBA for sparse flux; walk
from β-carotene-containing metabolites outward through the active
subgraph; report reactions with flux magnitudes.
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
from s250_atp_sweep import load_model, V7_EXTRACTED, V7_BASE_DIR

OUT = REPO / "grace_v4/results/atp_sweep_v7"
BIOMASS = "MAR00021"
FREE = ["MAR09152","MAR09154","MAR09276","MAR09361"]
CELLS_TO_RUN = ["LD_72","D_72","BL_96","LD_48"]


def apply_old(model, cond, hpf):
    is_ex = (np.diff(model["S"].indptr) == 1)
    lb, ub, idx, _ = acb.apply_context_bounds(model["lb"].copy(), model["ub"].copy(),
                                              model["rxns"], is_ex, cond, hpf,
                                              skip_energy_gate=True)
    for rid in FREE:
        if rid in idx: lb[idx[rid]] = -1000.0
    return lb, ub, idx


def pfba_biomax(model, lb, ub, idx):
    S = model["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag",0); env.start()
    gm = gp.Model("fba", env=env); v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[BIOMASS]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm.optimize()
    bm = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
    gm.dispose()
    gm2 = gp.Model("pfba", env=env)
    vp = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(ub>0, ub, 0.0))
    vn = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(lb<0, -lb, 0.0))
    SS = sp.hstack([S, -S], format="csr")
    gm2.addMConstr(SS, gp.MVar.fromlist(list(vp.tolist())+list(vn.tolist())),
                   "=", np.zeros(n_met))
    gm2.addConstr(vp[idx[BIOMASS]] - vn[idx[BIOMASS]] == bm)
    gm2.setMObjective(None, np.ones(2*n_rxn), 0.0, sense=GRB.MINIMIZE); gm2.optimize()
    flux = np.array([vp[j].X - vn[j].X for j in range(n_rxn)])
    gm2.dispose(); env.dispose()
    return bm, flux


rows = []
for cell in CELLS_TO_RUN:
    cond, h = cell.split("_"); h = int(h)
    print(f"\n[cell] {cell}", flush=True)
    m = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    idx = {r:i for i,r in enumerate(m["rxns"])}
    lb_o, ub_o, _ = apply_old(m, cond, h)
    bm, flux = pfba_biomax(m, lb_o, ub_o, idx)
    print(f"  biomass_max = {bm:.3f}", flush=True)
    # Carotenoid metabolite indices: any met name containing 'caroten', 'retinol', 'retinoid', 'retinal', 'retinoic'
    carotenoid_mets = []
    for mi, name in enumerate(m["metNames"]):
        nm = str(name).lower()
        if any(x in nm for x in ["caroten", "retinal", "retinol", "retinoic", "retinoid", "retinyl"]):
            carotenoid_mets.append(mi)
    print(f"  carotenoid/retinoid metabolite indices: {len(carotenoid_mets)}", flush=True)
    # Reactions that touch any carotenoid-related metabolite with |v| > 1e-6
    S = m["S"]
    for j in range(len(m["rxns"])):
        if abs(flux[j]) < 1e-6: continue
        col = S.getcol(j).tocoo()
        touches = [int(mi) for mi in col.row if int(mi) in carotenoid_mets]
        if not touches: continue
        rxn_name = m.get("rxnNames", [""]*len(m["rxns"]))[j] if "rxnNames" in m else ""
        subs = m["subs"][j] if j < len(m["subs"]) else ""
        gpr = str(m["grRules"][j]).replace("\x00","").strip() if m["grRules"][j] is not None else ""
        touched_names = [m["metNames"][mi] for mi in touches]
        rows.append(dict(cell=cell, biomass_max=bm, rxn=m["rxns"][j],
                         flux=float(flux[j]), subsystem=subs,
                         touched_metabolites=";".join(touched_names),
                         gpr=gpr))
        print(f"   {m['rxns'][j]:12s} v={flux[j]:+8.3f} subs={subs[:35]:35s} touches={touched_names}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "carotenoid_route_under_old.csv", index=False)
print(f"\n[write] carotenoid_route_under_old.csv rows={len(df)}", flush=True)
# Per-cell summary
print("\n=== Per-cell summary: # reactions carrying carotenoid carbon ===")
for cell in CELLS_TO_RUN:
    sub = df[df.cell == cell]
    print(f"  {cell}: {len(sub)} rxns with carotenoid/retinoid flux")
