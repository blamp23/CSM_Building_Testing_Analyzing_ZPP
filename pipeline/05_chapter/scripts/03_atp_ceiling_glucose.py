"""
s269 -- Maximum ATP maintenance flux under physiological OCR with GLUCOSE
AS THE ONLY CARBON SOURCE.

Variant of s268 that closes every other carbon exchange in the ATP-sweep
panel (43 of 44) and leaves glucose (MAR09034) open at lb = -1000. The
physiological OCR cap stays on MAR09048. Three biomass regimes same as
s268: free / per-cell floor / matched floor.

Expected ratio: ~5 ATP per O2 under oxphos (glucose complete oxidation
~31.5 ATP per 6 O2). At OCR = 97 → ATP ≈ 509; at OCR = 313 → ATP ≈ 1643.

Deviations from that ratio tell us how efficiently each cell's retained
glycolysis/TCA/oxphos chain converts glucose to ATP.

Output: grace_v4/results/v7_susceptibility/atp_ceiling_glucose_only.csv
"""
import sys
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FUNCTIONS, V7_MODELS as V7, SUSCEPTIBILITY as OUT, ATP_SWEEP
sys.path.insert(0, str(FUNCTIONS))
from apply_context_bounds import apply_context_bounds, compute_biomax
OUT.mkdir(parents=True, exist_ok=True)

ATP_RXN = "MAR03964"
BIOMASS_RXN = "MAR00021"
GLUCOSE_EX = "MAR09034"
O2_EX = "MAR09048"
OCR = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

# All 44 substrates in the ATP sweep panel. We close every carbon exchange
# EXCEPT glucose. Bulk inorganics and vitamins remain open per
# apply_context_bounds.
PANEL = pd.read_csv(ATP_SWEEP / "atp_sweep_substrates.csv")
NON_GLUCOSE_CARBON_MARS = set(PANEL.loc[PANEL.mar != GLUCOSE_EX, "mar"])
print(f"[prep] closing {len(NON_GLUCOSE_CARBON_MARS)} non-glucose carbon exchanges")


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


def _max_atp(S, lb, ub, idx):
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
        return np.nan, np.nan, np.nan
    atp = gm.ObjVal
    v_opt = np.asarray(v.X)
    o2_uptake = -v_opt[idx[O2_EX]] if O2_EX in idx else np.nan
    glc_uptake = -v_opt[idx[GLUCOSE_EX]] if GLUCOSE_EX in idx else np.nan
    gm.dispose(); env.dispose()
    return float(atp), float(o2_uptake), float(glc_uptake)


# ---- prep: load each cell, build glucose-only bounds, compute biomax ----
print("[prep] biomax per cell under glucose-only + physiological OCR")
cells_data = {}
biomax = {}
for cell in CELLS:
    cond, h = cell.split("_"); hpf = int(h)
    mdl = _load(V7 / f"trans_rfastcormics_{cell}.mat")
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex, cond, hpf,
        S=mdl["S"])
    # Widen internals (same as Probe 1 / s268)
    non_ex = ~is_ex
    lb[(lb <= -999.9) & non_ex] = -1e6
    ub[(ub >= 999.9) & non_ex] = +1e6
    # Close all non-glucose carbon exchanges in the ATP sweep panel
    for mar in NON_GLUCOSE_CARBON_MARS:
        if mar in idx:
            j = idx[mar]; lb[j] = 0.0
    # Make sure glucose is open for uptake
    if GLUCOSE_EX in idx:
        j = idx[GLUCOSE_EX]
        lb[j] = -1000.0; ub[j] = 1000.0
    else:
        print(f"  {cell}: GLUCOSE EXCHANGE NOT IN MODEL — glucose was pruned", flush=True)
    bm, _ = compute_biomax(mdl, lb, ub, idx)
    biomax[cell] = bm
    cells_data[cell] = (mdl, lb, ub, idx)
    print(f"  {cell:7s}  biomax(glucose-only) = {bm:.3f}  OCR_cap = {OCR[hpf]}", flush=True)

# ---- ATP ceiling under 3 regimes ----
rows = []

print("\n[run] free biomass regime (glucose-only)")
for cell in CELLS:
    mdl, lb, ub, idx = cells_data[cell]
    lb_r = lb.copy(); lb_r[idx[BIOMASS_RXN]] = 0.0
    atp, o2, glc = _max_atp(mdl["S"], lb_r, ub, idx)
    hpf = int(cell.split("_")[1])
    atp_per_o2 = atp / o2 if (o2 and o2 > 1e-9) else np.nan
    atp_per_glc = atp / glc if (glc and glc > 1e-9) else np.nan
    rows.append(dict(cell=cell, regime="free", floor=0.0, partner="",
                     atp_max=atp, o2_uptake=o2, glucose_uptake=glc,
                     atp_per_o2=atp_per_o2, atp_per_glucose=atp_per_glc,
                     OCR_cap=OCR[hpf]))
    print(f"  {cell:7s}  ATP={atp:7.2f}  O2={o2:6.2f}  glc={glc:6.2f}  "
          f"ATP/O2={atp_per_o2:4.2f}  ATP/glc={atp_per_glc:4.2f}", flush=True)

print("\n[run] per-cell floor regime (biomass ≥ 0.9 × biomax)")
for cell in CELLS:
    mdl, lb, ub, idx = cells_data[cell]
    lb_r = lb.copy(); lb_r[idx[BIOMASS_RXN]] = 0.9 * biomax[cell]
    atp, o2, glc = _max_atp(mdl["S"], lb_r, ub, idx)
    hpf = int(cell.split("_")[1])
    atp_per_o2 = atp / o2 if (o2 and o2 > 1e-9) else np.nan
    atp_per_glc = atp / glc if (glc and glc > 1e-9) else np.nan
    rows.append(dict(cell=cell, regime="per_cell_floor", floor=0.9*biomax[cell],
                     partner="", atp_max=atp, o2_uptake=o2, glucose_uptake=glc,
                     atp_per_o2=atp_per_o2, atp_per_glucose=atp_per_glc,
                     OCR_cap=OCR[hpf]))
    print(f"  {cell:7s}  ATP={atp:7.2f}  O2={o2:6.2f}  glc={glc:6.2f}  "
          f"ATP/O2={atp_per_o2:4.2f}", flush=True)

print("\n[run] matched_floor regime (same-hpf pairs)")
PAIRS = [(f"{A}_{hpf}", f"{B}_{hpf}")
         for hpf in HPFS for A, B in [("BL","D"), ("BL","LD"), ("D","LD")]]
for A, B in PAIRS:
    floor_val = 0.9 * min(biomax[A], biomax[B])
    for cell, partner in [(A, B), (B, A)]:
        mdl, lb, ub, idx = cells_data[cell]
        lb_r = lb.copy(); lb_r[idx[BIOMASS_RXN]] = floor_val
        atp, o2, glc = _max_atp(mdl["S"], lb_r, ub, idx)
        hpf = int(cell.split("_")[1])
        atp_per_o2 = atp / o2 if (o2 and o2 > 1e-9) else np.nan
        atp_per_glc = atp / glc if (glc and glc > 1e-9) else np.nan
        rows.append(dict(cell=cell, regime="matched_floor", floor=floor_val,
                         partner=partner, atp_max=atp, o2_uptake=o2,
                         glucose_uptake=glc, atp_per_o2=atp_per_o2,
                         atp_per_glucose=atp_per_glc, OCR_cap=OCR[hpf]))

df = pd.DataFrame(rows)
df.to_csv(OUT / "atp_ceiling_glucose_only.csv", index=False)
print(f"\n[write] {OUT/'atp_ceiling_glucose_only.csv'}  ({len(df)} rows)")

# Summary
print("\n=== FREE regime (pure glucose + OCR capacity) ===")
free = df[df.regime == "free"].copy()
free["hpf"] = free.cell.str.split("_").str[1].astype(int)
free["cond"] = free.cell.str.split("_").str[0]
piv = free.pivot_table(index="hpf", columns="cond", values="atp_max")
print(piv.round(2))
print()
print("ATP/O2 ratio (free regime):")
piv2 = free.pivot_table(index="hpf", columns="cond", values="atp_per_o2")
print(piv2.round(2))
print()
print("Cross-cond spread per hpf:")
for hpf, row in piv.iterrows():
    rng = row.max() - row.min()
    spread = 100 * rng / row.min() if row.min() > 0 else 0
    print(f"  hpf={hpf}: BL={row['BL']:.1f} D={row['D']:.1f} LD={row['LD']:.1f}  spread={spread:.1f}%  OCR={OCR[hpf]}  expected≈{OCR[hpf]*5.25:.0f}")
