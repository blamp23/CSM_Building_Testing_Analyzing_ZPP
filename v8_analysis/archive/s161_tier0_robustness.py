"""
s161 -- Tier 0 robustness: perturb amino acid + glucose uptakes by 10%.

For each of 15 extracted cells × 4 uptake exchanges (valine, histidine,
isoleucine, glucose): FBA biomass_max at the normal bounds (`apply_context_bounds`),
then with lb(exchange) relaxed by 10% (lb_new = 1.10 * lb_old; uptake is lb < 0,
so this allows 10% more uptake). Report Δbiomass_max as %.

Expected: amino acid uptake relaxation lifts biomass_max; glucose does not.

Output: v7_final_addenda/tier0_robustness.csv
"""
import sys
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_final_addenda"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
EX = {"valine":    "MAR09046",
      "histidine": "MAR09038",
      "isoleucine":"MAR09039",
      "glucose":   "MAR09034"}
CONDS = ("BL","D","LD"); HPFS = (24,48,72,96,120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]


def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)


def _load(fp, root_key):
    with h5py.File(fp, "r") as f:
        m = f[root_key]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        data = np.asarray(m["S"]["data"]); ir = np.asarray(m["S"]["ir"]); jc = np.asarray(m["S"]["jc"])
        try: n_met_hdf = int(np.asarray(m["mets"]).shape[0])
        except Exception: n_met_hdf = 0
        n_met = max(int(ir.max())+1, n_met_hdf)
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc)-1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, S=S, lb=lb, ub=ub)


def fba(mdl, lb, ub, idx):
    S = mdl["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c_obj = np.zeros(n_rxn); c_obj[idx[BIOMASS_RXN]] = 1.0
    gm.setMObjective(None, c_obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    bm = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return bm


rows = []
for cell in CELLS:
    cond, hpf = cell.split("_"); hpf = int(hpf)
    ch = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    is_ex = (np.diff(ch["S"].tocsc().indptr) == 1)
    lb0, ub0, idx, _ = apply_context_bounds(ch["lb"].copy(), ch["ub"].copy(),
                                             ch["rxns"], is_ex, cond, hpf)
    base = fba(ch, lb0, ub0, idx)
    rec = {"cell": cell, "hpf": hpf, "condition": cond, "biomax_baseline": base}
    for name, rxn in EX.items():
        if rxn not in idx:
            rec[f"pct_{name}"] = np.nan
            rec[f"biomax_{name}_relaxed"] = np.nan
            continue
        j = idx[rxn]
        lb_new = lb0.copy()
        # relax by 10% in the direction of more uptake (lb more negative)
        lb_new[j] = 1.10 * lb0[j]
        bm = fba(ch, lb_new, ub0, idx)
        rec[f"biomax_{name}_relaxed"] = bm
        rec[f"pct_{name}"] = (bm - base) / base * 100 if base else np.nan
    rows.append(rec)
    print(f"  {cell}  base={base:.3f}  "
          f"val={rec['pct_valine']:+.3f}%  "
          f"his={rec['pct_histidine']:+.3f}%  "
          f"ile={rec['pct_isoleucine']:+.3f}%  "
          f"glc={rec['pct_glucose']:+.3f}%", flush=True)

df = pd.DataFrame(rows)
cols = (["cell","hpf","condition","biomax_baseline"]
        + [c for name in EX for c in (f"biomax_{name}_relaxed", f"pct_{name}")])
df = df[cols]
out_fp = OUT/"tier0_robustness.csv"
df.to_csv(out_fp, index=False)
print(f"\n[write] {out_fp}")

# verdicts
print("\n--- verdict ---")
for name in EX:
    sig = (df[f"pct_{name}"].abs() > 0.05).sum()
    mx = df[f"pct_{name}"].abs().max()
    print(f"  {name:11s}  max |Δ%| = {mx:+.3f}%  (nonzero in {sig}/15 cells)")

# flag glucose exceptions
glc_moved = df[df["pct_glucose"].abs() > 0.5]
if len(glc_moved):
    print(f"\n[EXCEPTION] glucose shifts biomass in these cells:")
    print(glc_moved[["cell","biomax_baseline","biomax_glucose_relaxed","pct_glucose"]].to_string(index=False))
else:
    print("\n[OK] glucose is nowhere near biomass-binding in any of the 15 cells (|Δ%| ≤ 0.5)")
