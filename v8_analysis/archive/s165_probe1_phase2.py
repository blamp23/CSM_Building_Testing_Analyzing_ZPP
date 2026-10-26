"""
s165 -- Probe 1: Phase-II / redox cofactor supply capacity.

For each of 4 cytosolic cofactors (NADPH, GSH, UDP-glucuronate, PAPS), per
cell, maximise a temporary sink `X_c -> ∅` under three biomass regimes:

  (1) free_biomass      — no biomass lb (capacity ignoring growth cost)
  (2) per_cell_floor    — biomass ≥ 0.9 × biomax[cell]
  (3) matched_floor     — biomass ≥ 0.9 × min(biomax[A], biomax[B]) for
                          every same-hpf condition pair (BL-vs-D, BL-vs-LD,
                          D-vs-LD); the sink is maximised for each cell at
                          each of the 3 pair-matched floors it participates in

Primary comparison for the figure = matched_floor numbers; per_cell_floor
and free_biomass are sidecars.

Dual degeneracy / box-edge check: for any optimum within 1% of a ±1000
exchange box, flag and re-solve loopless via CycleFreeFlux-style minimisation
of |v| subject to v_sink = sink_optimum.

Output:
  probe1_phase2_capacity.csv   cell, cofactor, regime, floor_hpf, floor_src,
                               sink_max, biomass_at_opt, box_edge, loopless_sink
"""
import sys, time, json
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds, compute_biomax

# For dual classification -- load ledger + measured-ex set once.
import pandas as _pd
_BOUNDS_CSV = REPO / "grace_v4/data/phase2_exchange_bounds_long.csv"
if not _BOUNDS_CSV.exists():
    _BOUNDS_CSV = REPO / "data/phase2_exchange_bounds_long.csv"
if not _BOUNDS_CSV.exists():
    _BOUNDS_CSV = REPO / "results/qc/phase2_exchange_bounds_long.csv"
_LEDGER_CSV = REPO / "reviewer_packet_v4_slim/s130_v7_ledger.csv"
MEASURED_EX = set(_pd.read_csv(_BOUNDS_CSV).ex_rxn.unique())
_LED = _pd.read_csv(_LEDGER_CSV)
MEDIUM_EX = set(_LED.loc[_LED.action == "medium_uptake", "exchange_MAR_id"].dropna())
LEDGER_EX = set(_LED.loc[_LED.action.isin(["keep_exchange", "ysl_lipid_delivery"]),
                          "exchange_MAR_id"].dropna())
MEDIUM_HARD = {"MAR09047","MAR09058","MAR09072","MAR09073","MAR09074","MAR09076",
               "MAR09077","MAR09078","MAR09079","MAR09080","MAR09081","MAR09082",
               "MAR09148","MAR09150","MAR13066","MAR13072","MAR13073"}


def _classify_bound(rid, is_ex_var, lb, ub, tol=1e-6):
    """Return bound_type string for the dual-diagnostic row."""
    if not is_ex_var:
        # Non-exchange: check if ub or lb is at ±1000 default
        if (abs(ub - 1000.0) < tol) or (abs(lb + 1000.0) < tol):
            return "internal_default_1000"
        return "internal_other"
    if rid in MEDIUM_HARD:
        return "medium_hard"
    if rid in MEDIUM_EX:
        return "medium_uptake_demand"
    if rid in MEASURED_EX:
        return "measured_ex"
    if rid in LEDGER_EX:
        return "ledger_keep_ex"
    return "exchange_other"

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

COFACTORS = {
    "NADPH":    "MAM02555c",
    "GSH":      "MAM02026c",
    "UDP_GlcA": "MAM03109c",
    "PAPS":     "MAM02681c",
}
BOX = 1000.0
BOX_EDGE_FRAC = 0.99


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
        try: mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception: mets = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, S=S, lb=lb, ub=ub)


def _build_gurobi(S_aug, lb_aug, ub_aug):
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    n_met, n_rxn = S_aug.shape
    v = gm.addMVar(n_rxn, lb=lb_aug, ub=ub_aug)
    gm.addMConstr(S_aug, v, "=", np.zeros(n_met))
    return env, gm, v


def _augment(mdl, met_id):
    """Append sink column consuming met_id. Returns (S_aug_csr, sink_j, lb_aug, ub_aug)."""
    S = mdl["S"].tocsc()
    n_met, n_rxn = S.shape
    if met_id not in mdl["mets"]:
        return None, None, None, None
    mi = mdl["mets"].index(met_id)
    extra_col = sp.csc_matrix(([-1.0], ([mi], [0])), shape=(n_met, 1))
    S_aug = sp.hstack([S, extra_col]).tocsr()
    return S_aug, n_rxn  # sink column is at index n_rxn


def _maximise_sink(S_aug, sink_j, lb_aug, ub_aug, rxns=None, is_ex=None):
    env, gm, v = _build_gurobi(S_aug, lb_aug, ub_aug)
    obj = np.zeros(S_aug.shape[1]); obj[sink_j] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    if gm.Status != GRB.OPTIMAL:
        gm.dispose(); env.dispose()
        return np.nan, None, False, []
    sink = gm.ObjVal
    v_opt = np.asarray(v.X)
    at_box = ((np.abs(v_opt[:sink_j]) > BOX_EDGE_FRAC * BOX)).any()
    # Top-3 binding duals: for every variable at its bound, |RC| * |bound value|
    duals = []
    rc = np.asarray(v.RC)
    for j in range(S_aug.shape[1]):
        if j == sink_j: continue
        if abs(rc[j]) < 1e-9: continue
        # Determine which bound is active
        if abs(v_opt[j] - lb_aug[j]) < 1e-6:
            bound_val = lb_aug[j]; side = "lb"
        elif abs(v_opt[j] - ub_aug[j]) < 1e-6:
            bound_val = ub_aug[j]; side = "ub"
        else:
            continue
        score = abs(rc[j]) * abs(bound_val)
        rid = rxns[j] if (rxns is not None and j < len(rxns)) else f"col{j}"
        is_ex_var = bool(is_ex[j]) if (is_ex is not None and j < len(is_ex)) else False
        btype = _classify_bound(rid, is_ex_var, lb_aug[j], ub_aug[j])
        duals.append((score, rid, side, float(bound_val), float(rc[j]), btype))
    duals.sort(reverse=True)
    top3 = duals[:3]
    gm.dispose(); env.dispose()
    return sink, v_opt, bool(at_box), top3


def _loopless_sink(mdl, S_aug, sink_j, lb_aug, ub_aug, target_sink):
    """CycleFreeFlux-style: pin sink to target_sink, minimise sum |v| for
    internal non-exchange, non-sink reactions. Report the resulting sink
    (should equal target) and a flag if the LP is feasible."""
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    S = mdl["S"].tocsc()
    n_met_orig, n_rxn = S.shape
    # Pin sink to target
    lb2 = lb_aug.copy(); ub2 = ub_aug.copy()
    lb2[sink_j] = target_sink
    ub2[sink_j] = target_sink
    # Build LP with |v_i| via +/- split for internal reactions
    internal = ~is_ex
    int_idx = np.where(internal)[0]
    k = len(int_idx)
    # Augmented variable vector: [v (n_rxn+1), v_abs (k)]
    N = S_aug.shape[1] + k
    S_block = sp.hstack([S_aug, sp.csc_matrix((S_aug.shape[0], k))]).tocsr()
    # Add rows: v_i - v_abs_i <= 0, -v_i - v_abs_i <= 0 for i in internal
    # Use additional inequality rows; here we use setMObjective with abs trick
    # Simpler: add constraints row by row
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    lb_full = np.concatenate([lb2, np.zeros(k)])
    ub_full = np.concatenate([ub2, np.full(k, 1e4)])
    V = gm.addMVar(N, lb=lb_full, ub=ub_full)
    gm.addMConstr(S_aug, V[:S_aug.shape[1]], "=", np.zeros(S_aug.shape[0]))
    # |v_i| <= v_abs_i
    for a, j in enumerate(int_idx):
        gm.addConstr(V[j] - V[S_aug.shape[1] + a] <= 0)
        gm.addConstr(-V[j] - V[S_aug.shape[1] + a] <= 0)
    obj = np.zeros(N); obj[S_aug.shape[1]:] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MINIMIZE)
    gm.setParam("TimeLimit", 60)
    gm.optimize()
    if gm.Status not in (GRB.OPTIMAL, GRB.SUBOPTIMAL):
        gm.dispose(); env.dispose()
        return np.nan, False
    sink = V[sink_j].X
    gm.dispose(); env.dispose()
    return sink, True


# ---- per-cell biomax under final v7 bounds ----
print("[prep] computing biomax per cell (internal ±1000 widened to ±1e6)", flush=True)
biomax = {}
cell_cache = {}
for cell in CELLS:
    cond, h = cell.split("_"); hpf = int(h)
    mdl = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex, cond, hpf,
        S=mdl["S"],
    )
    # Widen all internal ±1000 Wang defaults to ±1e6 so no artefactual
    # internal cap can bind. Exchanges are untouched (envelope stays intact).
    non_ex = ~is_ex
    widen_lb = (lb <= -999.9) & non_ex
    widen_ub = (ub >= 999.9) & non_ex
    lb[widen_lb] = -1e6
    ub[widen_ub] = +1e6
    bm, _ = compute_biomax(mdl, lb, ub, idx)
    biomax[cell] = bm
    cell_cache[cell] = (mdl, lb, ub, idx)
    print(f"  {cell:7s}  biomax = {bm:.3f}  (widened {int(widen_lb.sum())} lb, {int(widen_ub.sum())} ub)", flush=True)

# ---- compute matched floors per hpf × pair ----
# For each hpf, for each pair (A,B), matched floor = 0.9 * min(biomax[A],biomax[B])
# A given cell may appear in up to 2 pairs at a given hpf.
matched_floors = {}  # (cell, partner) -> floor
for h in HPFS:
    cells_at_h = [f"{c}_{h}" for c in CONDS]
    for i in range(3):
        for j in range(i+1, 3):
            a, b = cells_at_h[i], cells_at_h[j]
            f = 0.9 * min(biomax[a], biomax[b])
            matched_floors[(a, b)] = f
            matched_floors[(b, a)] = f

# ---- run probe 1 ----
rows = []
t0 = time.time()
def _run_sink(cell, cond, hpf, cof_name, mdl, lb, ub, idx, S_aug, sink_j, lb_aug, ub_aug,
              regime, partner, partner_biomax, floor_val):
    rxns_aug = mdl["rxns"] + ["__SINK__"]
    is_ex_full = np.concatenate([(np.diff(mdl["S"].tocsc().indptr) == 1), [False]])
    sink_max, v_opt, box_edge, top3 = _maximise_sink(S_aug, sink_j, lb_aug, ub_aug,
                                                     rxns=rxns_aug, is_ex=is_ex_full)
    loopless = np.nan
    if box_edge and np.isfinite(sink_max):
        loopless, _ = _loopless_sink(mdl, S_aug, sink_j, lb_aug, ub_aug, sink_max)
    # Classify dominant bound
    dom_type = top3[0][5] if top3 else ""
    dom_rid  = top3[0][1] if top3 else ""
    dom_score = top3[0][0] if top3 else np.nan
    row = dict(cell=cell, condition=cond, hpf=hpf, cofactor=cof_name,
               regime=regime, floor=floor_val, partner=partner,
               partner_biomax=partner_biomax, sink_max=sink_max,
               box_edge=box_edge, loopless_sink=loopless,
               dom_bound_type=dom_type, dom_bound_rid=dom_rid,
               dom_score=dom_score,
               top3=json.dumps([dict(score=round(s,4), rid=r, side=sd, bound=b,
                                     rc=round(rc,6), type=t)
                                for s, r, sd, b, rc, t in top3]))
    return row


for cell in CELLS:
    cond, h = cell.split("_"); hpf = int(h)
    mdl, lb, ub, idx = cell_cache[cell]
    bm = biomax[cell]
    for cof_name, met_id in COFACTORS.items():
        aug = _augment(mdl, met_id)
        if aug is None or aug[0] is None:
            rows.append(dict(cell=cell, condition=cond, hpf=hpf, cofactor=cof_name,
                             regime="free", floor=0.0, partner="", partner_biomax=np.nan,
                             sink_max=np.nan, box_edge=False, loopless_sink=np.nan,
                             dom_bound_type="NA", dom_bound_rid="", dom_score=np.nan,
                             top3="[]"))
            continue
        S_aug, sink_j = aug
        lb_base = np.concatenate([lb, [0.0]])
        ub_base = np.concatenate([ub, [1000.0]])

        # Regime 1: free biomass
        lb_free = lb_base.copy(); lb_free[idx[BIOMASS_RXN]] = 0.0
        rows.append(_run_sink(cell, cond, hpf, cof_name, mdl, lb, ub, idx,
                              S_aug, sink_j, lb_free, ub_base,
                              "free", "", np.nan, 0.0))

        # Regime 2: per-cell floor
        lb_pc = lb_base.copy(); lb_pc[idx[BIOMASS_RXN]] = 0.9 * bm
        rows.append(_run_sink(cell, cond, hpf, cof_name, mdl, lb, ub, idx,
                              S_aug, sink_j, lb_pc, ub_base,
                              "per_cell_floor", "", np.nan, 0.9 * bm))

        # Regime 3: matched floor per partner
        for other in CONDS:
            partner = f"{other}_{hpf}"
            if partner == cell: continue
            mf = matched_floors[(cell, partner)]
            lb_mf = lb_base.copy(); lb_mf[idx[BIOMASS_RXN]] = mf
            rows.append(_run_sink(cell, cond, hpf, cof_name, mdl, lb, ub, idx,
                                  S_aug, sink_j, lb_mf, ub_base,
                                  "matched_floor", partner, biomax[partner], mf))
    dt = time.time() - t0
    print(f"  {cell:7s}  done ({dt:.1f}s)", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "probe1_phase2_capacity.csv", index=False)
print(f"\n[write] {OUT/'probe1_phase2_capacity.csv'}  ({len(df)} rows)")

# Verdict: pairwise matched_floor by hpf × pair direction
print("\n--- pairwise matched_floor sink_max (cell vs partner) ---")
mf = df[df.regime == "matched_floor"].copy()
for cof in COFACTORS:
    sub = mf[mf.cofactor == cof]
    piv = sub.pivot_table(index=["hpf", "partner"], columns="condition",
                          values="sink_max")
    print(f"\n[{cof}]")
    print(piv.round(3).to_string())

# Dominant bound classification summary
print("\n--- dominant bound type by regime ---")
counts = df.groupby(["regime", "dom_bound_type"]).size().unstack(fill_value=0)
print(counts.to_string())
# Fraction of matched_floor sinks whose dominant bound is an internal_default_1000
mf_bad = mf[mf.dom_bound_type == "internal_default_1000"]
print(f"\nmatched_floor sinks capped by internal_default_1000: {len(mf_bad)}/{len(mf)}")
if len(mf_bad):
    print(mf_bad[["cell","cofactor","partner","sink_max","dom_bound_rid"]].head(20).to_string(index=False))
