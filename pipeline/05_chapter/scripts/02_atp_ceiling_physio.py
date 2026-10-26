"""
s268 -- Maximum ATP maintenance flux under physiological OCR and measured
envelope. One number per cell under the three biomass regimes used by
Probe 1 (free, per-cell floor, matched floor). Reports the ATP ceiling
and the top-3 binding duals (which bound set the ceiling).

OCR cap and envelope are applied via functions/apply_context_bounds.py:
24/48/72/96/120 hpf → 97/173/227/273/313 µmol O2/g DW/hr on MAR09048.
Internal ±1000 Wang defaults widened to ±1e6 as in Probe 1.

Output: grace_v4/results/v7_susceptibility/atp_ceiling_physio_ocr.csv
"""
import sys
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FUNCTIONS, V7_MODELS as V7, SUSCEPTIBILITY as OUT, EXCHANGE_BOUNDS, LEDGER
sys.path.insert(0, str(FUNCTIONS))
from apply_context_bounds import apply_context_bounds, compute_biomax
OUT.mkdir(parents=True, exist_ok=True)

ATP_RXN = "MAR03964"
BIOMASS_RXN = "MAR00021"
O2_EX = "MAR09048"
OCR = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

# Load MEASURED/MEDIUM sets for dual classification
MEASURED_EX = set(pd.read_csv(EXCHANGE_BOUNDS).ex_rxn.unique())
_LED = pd.read_csv(LEDGER)
MEDIUM_EX = set(_LED.loc[_LED.action == "medium_uptake", "exchange_MAR_id"].dropna())
LEDGER_EX = set(_LED.loc[_LED.action.isin(["keep_exchange", "ysl_lipid_delivery"]),
                          "exchange_MAR_id"].dropna())
MEDIUM_HARD = {"MAR09047","MAR09058","MAR09072","MAR09073","MAR09074","MAR09076",
               "MAR09077","MAR09078","MAR09079","MAR09080","MAR09081","MAR09082",
               "MAR09148","MAR09150","MAR13066","MAR13072","MAR13073"}


def _classify_bound(rid, is_ex_var, lb, ub, tol=1e-6):
    if not is_ex_var:
        if (abs(ub - 1000.0) < tol) or (abs(lb + 1000.0) < tol):
            return "internal_default_1000"
        return "internal_other"
    if rid == O2_EX:                  return "OCR_cap"
    if rid in MEDIUM_HARD:            return "medium_hard"
    if rid in MEDIUM_EX:              return "medium_uptake_demand"
    if rid in MEASURED_EX:            return "measured_ex"
    if rid in LEDGER_EX:              return "ledger_keep_ex"
    return "exchange_other"


def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)


def _load(fp):
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = int(ir.max()) + 1 if len(ir) else 0
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, S=S, lb=lb, ub=ub)


def _max_atp(S, lb, ub, rxns, idx):
    atp_j = idx[ATP_RXN]
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    n_met, n_rxn = S.shape
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    obj = np.zeros(n_rxn); obj[atp_j] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    if gm.Status != GRB.OPTIMAL:
        gm.dispose(); env.dispose()
        return np.nan, []
    atp = gm.ObjVal
    v_opt = np.asarray(v.X); rc = np.asarray(v.RC)
    is_ex = (np.diff(S.tocsc().indptr) == 1)
    duals = []
    for j in range(n_rxn):
        if j == atp_j: continue
        if abs(rc[j]) < 1e-9: continue
        if abs(v_opt[j] - lb[j]) < 1e-6:
            bound_val = lb[j]
        elif abs(v_opt[j] - ub[j]) < 1e-6:
            bound_val = ub[j]
        else:
            continue
        score = abs(rc[j]) * abs(bound_val)
        rid = rxns[j]
        btype = _classify_bound(rid, bool(is_ex[j]), lb[j], ub[j])
        duals.append((score, rid, float(bound_val), btype))
    duals.sort(reverse=True)
    gm.dispose(); env.dispose()
    return float(atp), duals[:3]


# -------- compute biomax per cell --------
print("[prep] biomax per cell", flush=True)
cells_data = {}
biomax = {}
for cell in CELLS:
    cond, h = cell.split("_"); hpf = int(h)
    mdl = _load(V7 / f"trans_rfastcormics_{cell}.mat")
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex, cond, hpf,
        S=mdl["S"])
    non_ex = ~is_ex
    lb[(lb <= -999.9) & non_ex] = -1e6
    ub[(ub >= 999.9) & non_ex] = +1e6
    bm, _ = compute_biomax(mdl, lb, ub, idx)
    biomax[cell] = bm
    cells_data[cell] = (mdl, lb, ub, idx, is_ex)
    print(f"  {cell:7s}  biomax = {bm:.3f}  OCR_cap = {OCR[hpf]}", flush=True)

# -------- ATP ceiling under 3 regimes --------
rows = []

# (a) free biomass
print("\n[run] free biomass regime", flush=True)
for cell in CELLS:
    mdl, lb, ub, idx, _ = cells_data[cell]
    lb_r = lb.copy(); lb_r[idx[BIOMASS_RXN]] = 0.0  # no biomass floor
    atp, duals = _max_atp(mdl["S"], lb_r, ub, mdl["rxns"], idx)
    top1 = duals[0] if duals else (0, "", 0, "")
    rows.append(dict(cell=cell, regime="free", floor=0.0, partner="",
                     atp_max=atp, OCR_cap=OCR[int(cell.split("_")[1])],
                     dom_rid=top1[1], dom_bound=top1[2], dom_type=top1[3],
                     top3=[(round(d[0],2), d[1], d[2], d[3]) for d in duals]))
    print(f"  {cell:7s}  ATP_max = {atp:.3f}  dom = {top1[3]:25s} ({top1[1]})", flush=True)

# (b) per-cell floor (biomass >= 0.9 * biomax[cell])
print("\n[run] per-cell floor regime (biomass >= 0.9 * biomax)", flush=True)
for cell in CELLS:
    mdl, lb, ub, idx, _ = cells_data[cell]
    lb_r = lb.copy(); lb_r[idx[BIOMASS_RXN]] = 0.9 * biomax[cell]
    atp, duals = _max_atp(mdl["S"], lb_r, ub, mdl["rxns"], idx)
    top1 = duals[0] if duals else (0, "", 0, "")
    rows.append(dict(cell=cell, regime="per_cell_floor", floor=0.9*biomax[cell],
                     partner="", atp_max=atp,
                     OCR_cap=OCR[int(cell.split("_")[1])],
                     dom_rid=top1[1], dom_bound=top1[2], dom_type=top1[3],
                     top3=[(round(d[0],2), d[1], d[2], d[3]) for d in duals]))
    print(f"  {cell:7s}  ATP_max = {atp:.3f}  dom = {top1[3]:25s} ({top1[1]})", flush=True)

# (c) matched floor (biomass >= 0.9 * min(biomax[A], biomax[B]) per same-hpf pair)
print("\n[run] matched_floor regime (same-hpf pairs)", flush=True)
PAIRS = []
for hpf in HPFS:
    for A, B in [("BL","D"), ("BL","LD"), ("D","LD")]:
        PAIRS.append((f"{A}_{hpf}", f"{B}_{hpf}"))

for A, B in PAIRS:
    floor_val = 0.9 * min(biomax[A], biomax[B])
    for cell, partner in [(A, B), (B, A)]:
        mdl, lb, ub, idx, _ = cells_data[cell]
        lb_r = lb.copy(); lb_r[idx[BIOMASS_RXN]] = floor_val
        atp, duals = _max_atp(mdl["S"], lb_r, ub, mdl["rxns"], idx)
        top1 = duals[0] if duals else (0, "", 0, "")
        rows.append(dict(cell=cell, regime="matched_floor", floor=floor_val,
                         partner=partner, atp_max=atp,
                         OCR_cap=OCR[int(cell.split("_")[1])],
                         dom_rid=top1[1], dom_bound=top1[2], dom_type=top1[3],
                         top3=[(round(d[0],2), d[1], d[2], d[3]) for d in duals]))
        print(f"  {cell:7s} vs {partner:7s}  floor={floor_val:.3f}  "
              f"ATP = {atp:.3f}  dom = {top1[3]:25s} ({top1[1]})", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "atp_ceiling_physio_ocr.csv", index=False)
print(f"\n[write] {OUT/'atp_ceiling_physio_ocr.csv'}  ({len(df)} rows)")

# Summary table: matched_floor ATP per cell (take min across pairs for conservative)
mf = df[df.regime == "matched_floor"]
pivot = mf.groupby("cell").agg(
    atp_min=("atp_max", "min"),
    atp_mean=("atp_max", "mean"),
    atp_max=("atp_max", "max"),
    n_pairs=("atp_max", "count"),
    dom_type_mode=("dom_type", lambda s: s.mode().iloc[0] if len(s.mode()) else ""),
).reset_index()
print("\nMatched-floor ATP ceiling per cell:")
print(pivot.to_string(index=False))

# Cross-condition spread per hpf
print("\nMatched-floor ATP by (hpf, condition), mean across pairs:")
mf["hpf"] = mf.cell.str.split("_").str[1].astype(int)
mf["cond"] = mf.cell.str.split("_").str[0]
mat = mf.pivot_table(index="hpf", columns="cond", values="atp_max", aggfunc="mean")
print(mat.round(2))
for hpf, row in mat.iterrows():
    rng = row.max() - row.min()
    spread_pct = 100 * rng / row.min() if row.min() > 0 else 0
    print(f"  hpf={hpf}: min={row.min():.1f}  max={row.max():.1f}  spread={spread_pct:.1f}%  OCR={OCR[hpf]}")
