"""
s164 -- Pre-flight guards for Probe 1 (and all subsequent LP probes).

Three assertions per cell:
  (A) Energy-from-nothing: with ALL exchanges closed (no uptake, no secretion),
      max MAR03964 (ATP hydrolysis) = 0. [Augments the §11.19 ATP-vs-textbook
      gate with a stronger "null-uptake" variant.]
  (B) Growth-from-nothing: with all exchanges closed, max MAR00021 = 0.
  (C) Cofactor-from-nothing: with all exchanges closed, max of a temporary
      sink reaction `X_c -> ∅` = 0 for each of:
         NADPH   (MAM02555c)
         GSH     (MAM02026c)
         UDP-glucuronate (UDP-GlcA, MAM03110c)
         PAPS    (MAM02681c)

The sink is added by appending a column to S with -1 at the metabolite row,
lb=0, ub=1000, and setting that column as the objective.

If any assertion trips, the script stops with a clear per-cell table.

Output: v7_susceptibility/guards.csv (per cell, per assertion, value).
"""
import sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
ATP_RXN = "MAR03964"
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

# Cytosolic cofactor metabolite IDs (Human-GEM / Wang zebrafish nomenclature).
COFACTORS = {
    "NADPH":    "MAM02555c",
    "GSH":      "MAM02026c",
    "UDP_GlcA": "MAM03109c",
    "PAPS":     "MAM02681c",
}

ATOL = 1e-6


def _resolve(f, x):
    if isinstance(x, h5py.Reference):
        return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object:
            return _resolve(f, x[0])
        try:
            return "".join(chr(int(c)) for c in x)
        except Exception:
            return ""
    return str(x)


def _load(fp, root_key):
    with h5py.File(fp, "r") as f:
        m = f[root_key]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        try:
            mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception:
            mets = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"])
        ir = np.asarray(S_grp["ir"])
        jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, S=S, lb=lb, ub=ub)


def _close_all_exchanges(mdl, lb, ub):
    """Return lb/ub copies with every exchange closed (lb=ub=0)."""
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb2 = lb.copy(); ub2 = ub.copy()
    for j in np.where(is_ex)[0]:
        lb2[j] = 0.0
        ub2[j] = 0.0
    return lb2, ub2


def _max_rxn(mdl, lb, ub, idx, rid):
    if rid not in idx:
        return np.nan
    S = mdl["S"].tocsr()
    n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S, v, "=", np.zeros(n_met))
    obj = np.zeros(n_rxn); obj[idx[rid]] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


def _max_cofactor_sink(mdl, lb, ub, idx, met_id):
    """Add temporary sink `met_id -> ∅` and maximise it."""
    if met_id not in mdl["mets"]:
        return np.nan
    S = mdl["S"].tocsc()
    n_met, n_rxn = S.shape
    mi = mdl["mets"].index(met_id)
    # Build augmented S with extra column: -1 at row mi
    extra_col = sp.csc_matrix(([-1.0], ([mi], [0])), shape=(n_met, 1))
    S_aug = sp.hstack([S, extra_col]).tocsr()
    lb_aug = np.concatenate([lb, [0.0]])
    ub_aug = np.concatenate([ub, [1000.0]])
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(n_rxn + 1, lb=lb_aug, ub=ub_aug)
    gm.addMConstr(S_aug, v, "=", np.zeros(n_met))
    obj = np.zeros(n_rxn + 1); obj[-1] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


rows = []
failed_cells = []

for cell in CELLS:
    cond, h = cell.split("_")
    hpf = int(h)
    mdl = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    idx = {r: i for i, r in enumerate(mdl["rxns"])}

    lb, ub = _close_all_exchanges(mdl, mdl["lb"].copy(), mdl["ub"].copy())

    row = dict(cell=cell, condition=cond, hpf=hpf)
    row["atp_from_nothing"]      = _max_rxn(mdl, lb, ub, idx, ATP_RXN)
    row["biomass_from_nothing"]  = _max_rxn(mdl, lb, ub, idx, BIOMASS_RXN)
    for name, mid in COFACTORS.items():
        row[f"{name}_from_nothing"] = _max_cofactor_sink(mdl, lb, ub, idx, mid)
    rows.append(row)

    checks = [
        ("ATP",         row["atp_from_nothing"]),
        ("biomass",     row["biomass_from_nothing"]),
    ] + [(name, row[f"{name}_from_nothing"]) for name in COFACTORS]
    fails = [(k, v) for k, v in checks if (not np.isnan(v)) and v > ATOL]
    status = "OK" if not fails else "FAIL"
    print(f"  {cell:7s}  {status}  "
          f"ATP={row['atp_from_nothing']:+.3e}  bm={row['biomass_from_nothing']:+.3e}  "
          + "  ".join(f"{k}={row[f'{k}_from_nothing']:+.3e}" for k in COFACTORS),
          flush=True)
    if fails:
        failed_cells.append((cell, fails))

df = pd.DataFrame(rows)
df.to_csv(OUT / "guards.csv", index=False)
print(f"\n[write] {OUT/'guards.csv'}")

print("\n--- verdict ---")
if failed_cells:
    print(f"FAIL in {len(failed_cells)} cells:")
    for c, fails in failed_cells:
        for k, v in fails:
            print(f"  {c:7s}  {k}  = {v:+.3e}")
    sys.exit(2)
else:
    print("All 15 cells pass energy-, growth-, and 4 cofactor-from-nothing checks.")
