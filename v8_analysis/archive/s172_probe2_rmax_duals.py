"""
s172 -- What binds Probe 2 r_max? Top-3 duals at the H2O2-tolerance LP optimum
for a representative cell (BL_96, LD_96, D_96 at the matched_floor regime).

Catalases MAR03980/82/84 are blocked (as in s166). Internal ±1000 widened to
±1e6. H2O2 cytosol source MAM02041c pinned at r = r_max (bisection result
from probe2_ros_tolerance.csv).

Report: top-3 constraints by |RC × bound| at the LP optimum, with bound
classification. If the top binder is a measured_ex or medium exchange, the
value reflects the envelope; if it's biomass_floor, the value reflects the
matched-floor coupling; if it's an internal ±1e6 widened variable, that's
my introduced artefact.
"""
import sys, json
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds, compute_biomax

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"

BIOMASS_RXN = "MAR00021"
H2O2_C = "MAM02041c"
CATALASE_RXNS = ("MAR03980", "MAR03982", "MAR03984")

MEASURED_EX = set(pd.read_csv(REPO / "grace_v4/data/phase2_exchange_bounds_long.csv").ex_rxn.unique())
_LED = pd.read_csv(REPO / "reviewer_packet_v4_slim/s130_v7_ledger.csv")
MEDIUM_UP = set(_LED.loc[_LED.action == "medium_uptake", "exchange_MAR_id"].dropna())
LEDGER_KEEP = set(_LED.loc[_LED.action.isin(["keep_exchange", "ysl_lipid_delivery"]), "exchange_MAR_id"].dropna())
MEDIUM_HARD = {"MAR09047","MAR09058","MAR09072","MAR09073","MAR09074","MAR09076",
               "MAR09077","MAR09078","MAR09079","MAR09080","MAR09081","MAR09082",
               "MAR09148","MAR09150","MAR13066","MAR13072","MAR13073"}

CELLS = ["BL_96", "D_96", "LD_96", "BL_24", "D_24", "LD_24"]


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
        try: mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception: mets = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, S=S, lb=lb, ub=ub)


def _classify(rid, is_ex_var, lb, ub, side):
    """Classify the ACTIVE bound (side='lb' or 'ub')."""
    if rid == BIOMASS_RXN:
        return "biomass_floor" if side == "lb" else "biomass_cap"
    if rid in MEASURED_EX:   return "measured_ex"
    if rid in MEDIUM_UP:     return "medium_uptake_demand"
    if rid in MEDIUM_HARD:   return "medium_hard"
    if rid in LEDGER_KEEP:   return "ledger_keep_ex"
    val = lb if side == "lb" else ub
    if abs(abs(val) - 1e6) < 1e-3:
        return "internal_widened_1e6"
    if abs(abs(val) - 1000.0) < 1e-6:
        return "internal_default_1000"
    return "internal_other"


# r_max table from probe2_ros_tolerance.csv (per-cell regime value)
rmax_df = pd.read_csv(OUT / "probe2_ros_tolerance.csv")

rows = []
for cell in CELLS:
    cond, h = cell.split("_"); hpf = int(h)
    mdl = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat")
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex, cond, hpf,
        S=mdl["S"])
    non_ex = ~is_ex
    lb[(lb <= -999.9) & non_ex] = -1e6
    ub[(ub >= 999.9) & non_ex] = +1e6
    for rid in CATALASE_RXNS:
        if rid in idx:
            lb[idx[rid]] = 0.0; ub[idx[rid]] = 0.0
    bm, _ = compute_biomax(mdl, lb, ub, idx)

    sub = rmax_df[(rmax_df.cell == cell) & (rmax_df.regime == "per_cell")]
    if len(sub) == 0:
        print(f"  {cell}: no r_max row"); continue
    r_max = float(sub.iloc[0].r_max)

    # Build H2O2 source LP: pin ∅ -> MAM02041c at r_max; max biomass
    if H2O2_C not in mdl["mets"]:
        print(f"  {cell}: no H2O2 met"); continue
    mi = mdl["mets"].index(H2O2_C)
    S = mdl["S"].tocsc()
    extra = sp.csc_matrix(([1.0], ([mi], [0])), shape=(S.shape[0], 1))
    S_aug = sp.hstack([S, extra]).tocsc()
    n_rxn = S_aug.shape[1]
    sink_j = n_rxn - 1
    lb_aug = np.concatenate([lb, [r_max]])
    ub_aug = np.concatenate([ub, [r_max]])
    lb_aug[idx[BIOMASS_RXN]] = 0.9 * bm

    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(n_rxn, lb=lb_aug, ub=ub_aug)
    gm.addMConstr(S_aug.tocsr(), v, "=", np.zeros(S_aug.shape[0]))
    obj = np.zeros(n_rxn); obj[idx[BIOMASS_RXN]] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    if gm.Status != GRB.OPTIMAL:
        print(f"  {cell}: LP not optimal"); continue

    rc = np.asarray(v.RC)
    v_opt = np.asarray(v.X)
    rxns_aug = mdl["rxns"] + ["__SOURCE__"]
    is_ex_aug = np.concatenate([is_ex, [False]])
    duals = []
    for j in range(n_rxn):
        if j == sink_j: continue
        if abs(rc[j]) < 1e-9: continue
        if abs(v_opt[j] - lb_aug[j]) < 1e-6:
            side = "lb"; bval = lb_aug[j]
        elif abs(v_opt[j] - ub_aug[j]) < 1e-6:
            side = "ub"; bval = ub_aug[j]
        else:
            continue
        rid = rxns_aug[j]
        btype = _classify(rid, bool(is_ex_aug[j]), lb_aug[j], ub_aug[j], side)
        duals.append((abs(rc[j]) * abs(bval), rid, side, float(bval), float(rc[j]), btype))
    duals.sort(reverse=True)
    top3 = duals[:3]
    gm.dispose(); env.dispose()
    print(f"\n=== {cell}  r_max = {r_max:.1f}  biomax = {bm:.3f} ===")
    for s, r, sd, b, rcv, t in top3:
        print(f"  {r:12s} side={sd:2s} bound={b:+10.2f} rc={rcv:+.4e} score={s:9.2f} type={t}")
    rows.append(dict(cell=cell, r_max=r_max, biomax=bm,
                     top3=json.dumps([dict(score=round(s,3), rid=r, side=sd,
                                            bound=b, rc=round(rcv,6), type=t)
                                      for s, r, sd, b, rcv, t in top3])))

df = pd.DataFrame(rows)
df.to_csv(OUT / "probe2_rmax_duals.csv", index=False)
print(f"\n[write] {OUT/'probe2_rmax_duals.csv'}")
