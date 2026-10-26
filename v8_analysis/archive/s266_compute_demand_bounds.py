"""
Compute per-cell demand bounds for the 4 free-supplied class-1 exchanges
(medium_uptake rows that are not in the 105 measured set): inositol
MAR09361, β-carotene MAR09276, α-tocotrienol MAR09152, γ-tocotrienol
MAR09154.

Protocol (per cell):
  1. Load model, apply current apply_context_bounds.
  2. Max biomass under the current bounds (all 4 rows at lb=-1000) -> bm*.
  3. For each row, fix biomass >= 0.999 * bm*, then min v_row.
     v_demand = min v_row (most negative = largest uptake required).
  4. lb_corrected = max(-1.1 * |v_demand|, -1e-3).

Writes grace_v4/results/qc/medium_uptake_demand_bounds.csv with columns
cond, hpf, cell, mar, metabolite, biomass_max, fva_min_at_biomax,
lb_corrected.

Also runs for 5 base models (one per hpf) so the base can use the same
mechanism — we parameterise by (cond, hpf) with cond=base.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from apply_context_bounds import apply_context_bounds
from s250_atp_sweep import load_model, V7_EXTRACTED, V7_BASE_DIR, CELLS, HPFS

OUT_QC = REPO / "grace_v4/results/qc"
OUT_QC.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "MAR09361": "inositol",
    "MAR09276": "beta-carotene",
    "MAR09152": "alpha-tocotrienol",
    "MAR09154": "gamma-tocotrienol",
}
BIOMASS = "MAR00021"
FLOOR = -1e-3
SLACK = 1.1


def fva_min(model, lb, ub, idx, target_mar, bm_min):
    """Fix v[biomass] >= bm_min, minimise v[target_mar]. Returns min flux."""
    S = model["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    gm.addConstr(v[idx[BIOMASS]] >= bm_min)
    c = np.zeros(n_rxn); c[idx[target_mar]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MINIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


def biomass_max(model, lb, ub, idx):
    S = model["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[BIOMASS]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


def biomass_essential(model, lb, ub, idx, target_mar, floor=-1e-3, tol=0.01):
    """Does the row carry biomass-essential flow? Close it to floor, see if
    biomass_max drops by more than `tol` (fractional)."""
    if target_mar not in idx: return False
    bm0 = biomass_max(model, lb, ub, idx)
    if bm0 <= 1e-6: return False
    lb2 = lb.copy(); lb2[idx[target_mar]] = floor
    bm1 = biomass_max(model, lb2, ub, idx)
    return (bm0 - bm1) / bm0 > tol


rows = []
# 15 extracted
for cell in CELLS:
    cond, h = cell.split("_"); h = int(h)
    m = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    S = m["S"]; is_ex = (np.diff(S.indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(m["lb"].copy(), m["ub"].copy(),
                                          m["rxns"], is_ex, cond, h)
    bm = biomass_max(m, lb, ub, idx)
    print(f"[cell] {cell}  biomass_max={bm:.4f}", flush=True)
    for mar, met in TARGETS.items():
        if mar not in idx:
            rows.append(dict(cond=cond, hpf=h, cell=cell, mar=mar,
                             metabolite=met, biomass_max=bm,
                             fva_min_at_biomax=np.nan,
                             lb_corrected=np.nan, note="not in model"))
            continue
        fv = fva_min(m, lb, ub, idx, mar, 0.999 * bm)
        essential = biomass_essential(m, lb, ub, idx, mar)
        if essential:
            corrected = min(-SLACK * abs(fv), FLOOR) if np.isfinite(fv) else np.nan
            policy = "1.1x_FVA_min"
        else:
            corrected = FLOOR
            policy = "floor_biomass_insensitive"
        rows.append(dict(cond=cond, hpf=h, cell=cell, mar=mar,
                         metabolite=met, biomass_max=bm,
                         fva_min_at_biomax=fv,
                         biomass_essential=essential,
                         policy=policy,
                         lb_corrected=corrected, note=""))
        print(f"    {mar} ({met:20s}): fva_min={fv:+9.4f}  "
              f"ess={essential}  lb_corrected={corrected:+9.4f}  policy={policy}", flush=True)

# 5 base models (one per hpf, cond='base')
for h in HPFS:
    m = load_model(V7_BASE_DIR / f"baked_ocr_{h}_v7.mat", "m")
    S = m["S"]; is_ex = (np.diff(S.indptr) == 1)
    # Apply LD-env bounds to base (OCR cap per hpf is the main effect)
    lb, ub, idx, _ = apply_context_bounds(m["lb"].copy(), m["ub"].copy(),
                                          m["rxns"], is_ex, "LD", h)
    bm = biomass_max(m, lb, ub, idx)
    print(f"[cell] base_{h}  biomass_max={bm:.4f}", flush=True)
    for mar, met in TARGETS.items():
        if mar not in idx:
            rows.append(dict(cond="base", hpf=h, cell=f"base_{h}", mar=mar,
                             metabolite=met, biomass_max=bm,
                             fva_min_at_biomax=np.nan,
                             lb_corrected=np.nan, note="not in model"))
            continue
        fv = fva_min(m, lb, ub, idx, mar, 0.999 * bm)
        essential = biomass_essential(m, lb, ub, idx, mar)
        if essential:
            corrected = min(-SLACK * abs(fv), FLOOR) if np.isfinite(fv) else np.nan
            policy = "1.1x_FVA_min"
        else:
            corrected = FLOOR
            policy = "floor_biomass_insensitive"
        rows.append(dict(cond="base", hpf=h, cell=f"base_{h}", mar=mar,
                         metabolite=met, biomass_max=bm,
                         fva_min_at_biomax=fv,
                         biomass_essential=essential,
                         policy=policy,
                         lb_corrected=corrected, note=""))
        print(f"    {mar} ({met:20s}): fva_min={fv:+9.4f}  "
              f"ess={essential}  lb_corrected={corrected:+9.4f}  policy={policy}", flush=True)

df = pd.DataFrame(rows)
out_fp = OUT_QC / "medium_uptake_demand_bounds.csv"
df.to_csv(out_fp, index=False)
print(f"\n[write] {out_fp.relative_to(REPO)}  rows={len(df)}", flush=True)
print("\n=== Pivot: cell × metabolite, lb_corrected ===")
pivot = df[df.cond != "base"].pivot(index="cell", columns="metabolite", values="lb_corrected")
print(pivot.to_string())
