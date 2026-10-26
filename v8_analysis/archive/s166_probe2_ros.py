"""
s166 -- Probe 2: oxidative-stress tolerance.

For each cell, add a forced cytosolic H2O2 production `∅ -> MAM02041c` with
lb=ub=r. Bisection on r in [0, 1000]: find the largest r at which the
model remains feasible with biomass ≥ floor. Three regimes:

  (1) free        -- no biomass floor; r_max = pure stoichiometric limit
  (2) per_cell    -- biomass ≥ 0.9 × biomax[cell]
  (3) matched     -- biomass ≥ 0.9 × min(biomax[A], biomax[B]) per pair

At the tolerance limit (regime 'matched' with each partner), run FBA with
a secondary objective (min sum |v| across internal reactions, CFF'd), then
record:
  * H2O2 consumption by subsystem (ROS detox / Glutathione / Purine / other)
  * NADPH turnover (abs flux through G6PDH + 6PGD + NNT + ME + IDH, when present)
  * Biomass at r_max

Output:
  probe2_ros_tolerance.csv       cell, regime, partner, floor, r_max,
                                 h2o2_by_subsystem (JSON-ish cols), nadph_sources
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

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
H2O2_C = "MAM02041c"
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

# NADPH-source reactions (reaction IDs; abs flux summed as proxy)
NADPH_SOURCES = {
    "G6PDH": "MAR04306",   # typical glucose-6-P dehydrogenase in Human-GEM
    "PGD":   "MAR04304",
    "lactonase_6PG": "MAR04303",
    "PhosphoGlu_dehydro": "MAR08973",
    "NNT":   "MAR06975",
    "ME1":   "MAR03795",
    "ME2":   "MAR03794",
}

BISECT_TOL = 1e-2
R_HI = 1.0e5

# Catalase variants (disproportionation 2 H2O2 -> O2 + 2 H2O; no cofactor cost).
# Block these for Probe 2 so H2O2 clearance must route through NADPH-/GSH-coupled
# GPx or peroxiredoxin. Otherwise r_max is set by this "free vent" and tells us
# nothing about redox buffering.
CATALASE_RXNS = ("MAR03980", "MAR03982", "MAR03984")


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
        try: subs = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception: subs = [""] * len(rxns)
        try: mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception: mets = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, subs=subs, mets=mets, S=S, lb=lb, ub=ub)


def _augment_h2o2_source(mdl, r):
    """Append a production reaction ∅ -> MAM02041c with lb=ub=r.
    Returns (S_aug, src_j, lb_aug, ub_aug) or None if H2O2_C not in model."""
    S = mdl["S"].tocsc()
    n_met, n_rxn = S.shape
    if H2O2_C not in mdl["mets"]:
        return None
    mi = mdl["mets"].index(H2O2_C)
    # column: +1 at H2O2_C (produces)
    extra = sp.csc_matrix(([1.0], ([mi], [0])), shape=(n_met, 1))
    S_aug = sp.hstack([S, extra]).tocsc()
    return S_aug, n_rxn


def _feasible_at_r(mdl, lb, ub, idx, r, biomass_lb):
    """Build LP: v satisfies S*v=0, lb<=v<=ub, v_biomass >= biomass_lb,
    augmented with H2O2 source pinned at r. Return (feasible, v_opt or None)."""
    aug = _augment_h2o2_source(mdl, r)
    if aug is None:
        return False, None
    S_aug, src_j = aug
    n_rxn = S_aug.shape[1]
    lb_aug = np.concatenate([lb, [r]])
    ub_aug = np.concatenate([ub, [r]])
    lb_aug[idx[BIOMASS_RXN]] = biomass_lb
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(n_rxn, lb=lb_aug, ub=ub_aug)
    gm.addMConstr(S_aug.tocsr(), v, "=", np.zeros(S_aug.shape[0]))
    # Dummy objective (feasibility only)
    gm.setObjective(0.0, GRB.MINIMIZE)
    gm.optimize()
    ok = (gm.Status == GRB.OPTIMAL)
    v_opt = np.asarray(v.X) if ok else None
    gm.dispose(); env.dispose()
    return ok, v_opt


def _bisect_r_max(mdl, lb, ub, idx, biomass_lb):
    """Return largest r in [0, R_HI] with feasibility."""
    lo, hi = 0.0, R_HI
    # First check that lo is feasible
    ok0, _ = _feasible_at_r(mdl, lb, ub, idx, lo, biomass_lb)
    if not ok0:
        return 0.0
    ok_hi, _ = _feasible_at_r(mdl, lb, ub, idx, hi, biomass_lb)
    if ok_hi:
        return hi
    while hi - lo > BISECT_TOL:
        mid = 0.5 * (lo + hi)
        ok, _ = _feasible_at_r(mdl, lb, ub, idx, mid, biomass_lb)
        if ok:
            lo = mid
        else:
            hi = mid
    return lo


def _route_shares(mdl, lb, ub, idx, r_max, biomass_lb):
    """At r_max, solve a CFF-style LP: min sum|v| for internal reactions,
    subject to feasibility. Report H2O2 consumption by subsystem and NADPH
    source abs flux."""
    aug = _augment_h2o2_source(mdl, r_max)
    if aug is None:
        return {}, {}, np.nan
    S_aug, src_j = aug
    n_rxn = S_aug.shape[1]
    lb_aug = np.concatenate([lb, [r_max]])
    ub_aug = np.concatenate([ub, [r_max]])
    lb_aug[idx[BIOMASS_RXN]] = biomass_lb

    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    is_ex_aug = np.concatenate([is_ex, [False]])  # H2O2 source is internal
    internal = ~is_ex_aug
    int_idx = np.where(internal)[0]
    k = len(int_idx)
    N = n_rxn + k

    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    lb_full = np.concatenate([lb_aug, np.zeros(k)])
    ub_full = np.concatenate([ub_aug, np.full(k, 1e4)])
    V = gm.addMVar(N, lb=lb_full, ub=ub_full)
    gm.addMConstr(S_aug.tocsr(), V[:n_rxn], "=", np.zeros(S_aug.shape[0]))
    # Add v - v_abs <= 0 and -v - v_abs <= 0 for internal reactions
    # Build row by row via Python (slow but one-time per cell)
    for a, j in enumerate(int_idx):
        gm.addConstr(V[j] - V[n_rxn + a] <= 0)
        gm.addConstr(-V[j] - V[n_rxn + a] <= 0)
    obj = np.zeros(N); obj[n_rxn:] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MINIMIZE)
    gm.setParam("TimeLimit", 60)
    gm.optimize()
    if gm.Status not in (GRB.OPTIMAL, GRB.SUBOPTIMAL):
        gm.dispose(); env.dispose()
        return {}, {}, np.nan
    v_cff = np.asarray(V.X[:n_rxn])
    bm_at_opt = v_cff[idx[BIOMASS_RXN]]
    gm.dispose(); env.dispose()

    # H2O2 consumption by subsystem (sum of |flux| for reactions consuming H2O2_C)
    h2o2_rows = [i for i, mid in enumerate(mdl["mets"]) if mid == H2O2_C]
    Sc = mdl["S"].tocsc()
    h2o2_cons = {}
    for j in range(Sc.shape[1]):
        col = Sc[:, j]
        hits = [(int(r), float(vv)) for r, vv in zip(col.indices, col.data) if r in h2o2_rows]
        if not hits:
            continue
        # Reaction consumes H2O2 if any coefficient is negative and v > 0,
        # or any coefficient is positive and v < 0.
        v_j = v_cff[j]
        cons = sum(abs(vv * v_j) for r, vv in hits if (vv < 0 and v_j > 0) or (vv > 0 and v_j < 0))
        if cons < 1e-9:
            continue
        sub = mdl["subs"][j] if j < len(mdl["subs"]) else "other"
        h2o2_cons[sub] = h2o2_cons.get(sub, 0.0) + cons

    # NADPH source abs fluxes
    nadph_src = {}
    for name, rid in NADPH_SOURCES.items():
        if rid in idx:
            nadph_src[name] = float(abs(v_cff[idx[rid]]))
        else:
            nadph_src[name] = np.nan
    return h2o2_cons, nadph_src, bm_at_opt


# ---------- biomax prep ----------
print("[prep] computing biomax per cell (catalases blocked; internal ±1000 → ±1e6)", flush=True)
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
    # Widen internal ±1000 defaults to ±1e6
    non_ex = ~is_ex
    lb[(lb <= -999.9) & non_ex] = -1e6
    ub[(ub >= 999.9) & non_ex] = +1e6
    # Block catalase variants
    for rid in CATALASE_RXNS:
        if rid in idx:
            lb[idx[rid]] = 0.0
            ub[idx[rid]] = 0.0
    bm, _ = compute_biomax(mdl, lb, ub, idx)
    biomax[cell] = bm
    cell_cache[cell] = (mdl, lb, ub, idx)
    print(f"  {cell:7s} biomax={bm:.3f}", flush=True)

matched_floors = {}
for h in HPFS:
    cs = [f"{c}_{h}" for c in CONDS]
    for i in range(3):
        for j in range(i+1, 3):
            a, b = cs[i], cs[j]
            f = 0.9 * min(biomax[a], biomax[b])
            matched_floors[(a, b)] = f
            matched_floors[(b, a)] = f

# ---------- probe 2 ----------
rows = []
t0 = time.time()
for cell in CELLS:
    cond, h = cell.split("_"); hpf = int(h)
    mdl, lb, ub, idx = cell_cache[cell]

    # Regime 1: free biomass
    r1 = _bisect_r_max(mdl, lb, ub, idx, 0.0)
    rows.append(dict(cell=cell, condition=cond, hpf=hpf,
                     regime="free", partner="", floor=0.0, r_max=r1,
                     h2o2_shares="", nadph_sources="", biomass_at_opt=np.nan))

    # Regime 2: per-cell floor
    r2 = _bisect_r_max(mdl, lb, ub, idx, 0.9 * biomax[cell])
    rows.append(dict(cell=cell, condition=cond, hpf=hpf,
                     regime="per_cell", partner="", floor=0.9*biomax[cell], r_max=r2,
                     h2o2_shares="", nadph_sources="", biomass_at_opt=np.nan))

    # Regime 3: matched floor per partner (route shares recorded here)
    for other in CONDS:
        partner = f"{other}_{hpf}"
        if partner == cell: continue
        mf = matched_floors[(cell, partner)]
        r3 = _bisect_r_max(mdl, lb, ub, idx, mf)
        h2o2_cons, nadph_src, bm_at = _route_shares(mdl, lb, ub, idx, r3, mf) if r3 > 0 else ({}, {}, np.nan)
        rows.append(dict(cell=cell, condition=cond, hpf=hpf,
                         regime="matched", partner=partner, floor=mf, r_max=r3,
                         h2o2_shares=json.dumps({k: round(v, 3) for k, v in h2o2_cons.items()}),
                         nadph_sources=json.dumps({k: (round(v, 3) if np.isfinite(v) else None) for k, v in nadph_src.items()}),
                         biomass_at_opt=bm_at))
    dt = time.time() - t0
    print(f"  {cell:7s}  r_free={r1:8.3f}  r_pc={r2:8.3f}  "
          f"({dt:.1f}s elapsed)", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "probe2_ros_tolerance.csv", index=False)
print(f"\n[write] {OUT/'probe2_ros_tolerance.csv'}  ({len(df)} rows)")

print("\n--- verdict (matched; mean r_max over partners) ---")
mf = df[df.regime == "matched"].copy()
piv = mf.pivot_table(index="hpf", columns="condition", values="r_max", aggfunc="mean")
print(piv.round(3).to_string())
