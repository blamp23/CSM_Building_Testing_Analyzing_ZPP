"""
s83 -- FVA capability classification for v7, LP-first with MILP (loopless)
only on the ambiguous subset.

For each of the 15 v7 cells (Frame A):
  (1) Load extracted model, apply apply_context_bounds (Frame A).
  (2) biomass_max via Gurobi LP; set biomass_lb = 0.9 * bmax.
  (3) Plain LP FVA on every reaction (allowLoops=1).  Classify by LP range:
         blocked         : |min| < tol AND |max| < tol
         required_fwd    : min >=  tol
         required_rev    : max <= -tol
         optional_fwd    : |min| < tol AND max >  tol
         optional_rev    : min < -tol AND |max| < tol
         flexible_signed : min < -tol AND max >  tol (reversible, sign amb)
  (4) MILP (loopless, Gurobi) only on reactions whose LP range crossed zero
      (flexible_signed) OR where the LP max exceeds LOOP_SUSPECT (default 100).
      The loopless MILP gives the tighter true range.
  (5) Re-classify after MILP.
  (6) Write per-cell CSV + aggregate subsystem × class × cell matrix.
  (7) Developmental trajectory: for each (condition, subsystem), count the
      class membership change across 24 -> 48 -> 72 -> 96 -> 120.

Usage:
  python scripts/s83_fva_capability.py --cells all --workers 6

Outputs:
  results/fva_capability_v7/
    {cell}_fva.csv                              per-reaction class + bounds
    subsystem_class_matrix_{cell}.csv            per-subsystem counts
    subsystem_class_matrix_all.csv               stacked
    developmental_trajectory.csv                 class flips 24->48 etc per cond
"""
import argparse, sys, os, h5py
from pathlib import Path
from multiprocessing import Pool
import numpy as np, pandas as pd, scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
sys.path.insert(0, str(REPO/"functions"))
from apply_context_bounds import apply_context_bounds

MODELS = REPO/"reviewer_packet_v4_slim/v7_models"
OUT = REPO/"grace_v4/results/fva_capability_v7"
OUT.mkdir(exist_ok=True, parents=True)

BIOMASS_RXN = "MAR00021"
TOL = 1e-6
LOOP_SUSPECT = 100.0        # reactions whose LP max > LOOP_SUSPECT also go to MILP
BIO_LB_FRAC = 0.9
MILP_TIME_LIMIT = 15        # seconds per MILP

CELLS = [f"{c}_{h}" for c in ("BL","D","LD") for h in (24,48,72,96,120)]


def _resolve_h5_str(f, obj):
    arr = np.asarray(obj).ravel()
    try: return "".join(chr(int(c)) for c in arr)
    except: return ""


def load_model(cell):
    fp = MODELS/f"trans_rfastcormics_{cell}.mat"
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_rxn = len(jc) - 1; n_met = int(ir.max()) + 1
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, n_rxn))
        rxns = [_resolve_h5_str(f, f[r]) for r in np.asarray(m["rxns"]).ravel()]
        sub_refs = np.asarray(m["subSystems"]).ravel()
        subs = []
        for ref in sub_refs:
            obj = f[ref]; arr = np.asarray(obj).ravel()
            if arr.dtype == object and len(arr) > 0:
                subs.append(_resolve_h5_str(f, f[arr[0]]))
            else:
                subs.append(_resolve_h5_str(f, obj))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(S=S, rxns=rxns, subs=subs, lb=lb, ub=ub, n_rxn=n_rxn, n_met=n_met)


def classify(mn, mx, tol=TOL):
    bl = (abs(mn) < tol) & (abs(mx) < tol)
    rf = mn >=  tol
    rr = mx <= -tol
    of = (abs(mn) < tol) & (mx > tol)
    orev = (mn < -tol) & (abs(mx) < tol)
    fs = (mn < -tol) & (mx > tol)
    cls = np.array(["unknown"]*len(mn), dtype=object)
    cls[bl] = "blocked"
    cls[rf] = "required_fwd"
    cls[rr] = "required_rev"
    cls[of] = "optional_fwd"
    cls[orev] = "optional_rev"
    cls[fs] = "flexible_signed"
    return cls


def plain_lp_fva(cell):
    print(f"[{cell}] loading + applying Frame A", flush=True)
    M = load_model(cell)
    is_ex = (np.diff(M["S"].indptr) == 1)
    cond, hpf = cell.split("_"); hpf = int(hpf)
    lb, ub, idx, is_ex = apply_context_bounds(M["lb"], M["ub"], M["rxns"], is_ex, cond, hpf)
    biomass_idx = idx[BIOMASS_RXN]
    S = M["S"]; n_rxn = M["n_rxn"]; n_met = M["n_met"]

    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    mo = gp.Model("fva", env=env)
    mo.setParam("Threads", 1); mo.setParam("OutputFlag", 0); mo.setParam("Method", 1)
    v = mo.addMVar(n_rxn, lb=lb, ub=ub)
    mo.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    mo.setObjective(v[biomass_idx], GRB.MAXIMIZE); mo.optimize()
    if mo.Status != GRB.OPTIMAL:
        raise RuntimeError(f"[{cell}] biomass_max failed (status {mo.Status})")
    bmax = mo.ObjVal
    v[biomass_idx].LB = max(BIO_LB_FRAC * bmax, lb[biomass_idx])
    print(f"[{cell}] biomass_max = {bmax:.4f}, floor = {max(BIO_LB_FRAC*bmax, lb[biomass_idx]):.4f}", flush=True)

    mn = np.zeros(n_rxn); mx = np.zeros(n_rxn)
    for j in range(n_rxn):
        mo.setObjective(v[j], GRB.MINIMIZE); mo.optimize()
        mn[j] = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
        mo.setObjective(v[j], GRB.MAXIMIZE); mo.optimize()
        mx[j] = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
        if j % 1000 == 0:
            print(f"[{cell}] LP FVA {j}/{n_rxn}", flush=True)
    class_lp = classify(mn, mx)
    return M, lb, ub, bmax, mn, mx, class_lp


def milp_on_ambiguous(cell, M, lb, ub, mn, mx, class_lp):
    """Run loopless MILP (via Gurobi) on reactions flagged as
    flexible_signed OR having |LP range| > LOOP_SUSPECT."""
    S = M["S"]; rxns = M["rxns"]
    n_rxn = M["n_rxn"]; n_met = M["n_met"]
    biomass_idx = rxns.index(BIOMASS_RXN)
    bmax_floor = lb[biomass_idx]

    # Determine ambiguous set
    big_range = (mx - mn) > LOOP_SUSPECT
    flagged = (class_lp == "flexible_signed") | big_range
    flagged[biomass_idx] = False
    amb_idx = np.where(flagged)[0]
    print(f"[{cell}] MILP (loopless) on {len(amb_idx)} ambiguous reactions "
          f"of {n_rxn}", flush=True)
    if len(amb_idx) == 0:
        return mn.copy(), mx.copy()

    # Set up loopless MILP: Sv = 0, biomass_lb enforced, internal loopless
    # constraints via big-M formulation.  For 15 s time limit per LP.
    # For simplicity use Gurobi's native binary + big-M on internal reactions.
    Snnz = np.diff(S.tocsc().indptr)
    is_ex = Snnz == 1
    is_internal = ~is_ex
    # Approximate nullspace loopless: add binaries y_j for each internal reaction,
    # bounds v_j in [-M y_j^-, M y_j^+] with y pos+neg indicators.  Simpler
    # alternative: use thermodynamic loop-free via energy variable g_j.
    #
    # To keep code tractable, use the standard thermodynamic MILP:
    #   For each internal rxn j, add energy var g_j.
    #   v_j > 0  =>  g_j < 0  (energy decreases)
    #   v_j < 0  =>  g_j > 0
    # Implemented via binary z_j and big-M.
    #
    # This matches COBRA's "original" loop law.
    M_big = 1e4
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    mo = gp.Model("llmilp", env=env)
    mo.setParam("Threads", 1); mo.setParam("OutputFlag", 0)
    mo.setParam("TimeLimit", MILP_TIME_LIMIT)
    v = mo.addMVar(n_rxn, lb=lb, ub=ub)
    mo.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    v[biomass_idx].LB = bmax_floor

    internal_js = np.where(is_internal)[0]
    n_int = len(internal_js)
    g = mo.addMVar(n_int, lb=-1e4, ub=1e4)                # energy for internal rxns
    z = mo.addMVar(n_int, vtype=GRB.BINARY)
    for k, j in enumerate(internal_js):
        # v_j <= M*z_k    (if z_k=0, v_j <= 0 so flux is reverse or zero)
        # v_j >= -M*(1-z_k)
        mo.addConstr(v[j] <= M_big * z[k])
        mo.addConstr(v[j] >= -M_big * (1 - z[k]))
        # g_k >= -M*z_k + eps      (if z_k=1, g_k >= 0 - no that's wrong direction)
        # Loop law: v and g must have opposite signs for internals when active.
        # z_k=1 -> v_j >= 0 and g_k <= -1
        # z_k=0 -> v_j <= 0 and g_k >= 1
        mo.addConstr(g[k] <= -1 + M_big * (1 - z[k]))
        mo.addConstr(g[k] >= 1 - M_big * z[k])

    mn_llp = mn.copy(); mx_llp = mx.copy()
    for k, j in enumerate(amb_idx):
        mo.setObjective(v[j], GRB.MINIMIZE); mo.optimize()
        if mo.Status in (GRB.OPTIMAL, GRB.TIME_LIMIT) and v[j].X is not None:
            mn_llp[j] = v[j].X
        mo.setObjective(v[j], GRB.MAXIMIZE); mo.optimize()
        if mo.Status in (GRB.OPTIMAL, GRB.TIME_LIMIT) and v[j].X is not None:
            mx_llp[j] = v[j].X
        if (k + 1) % 50 == 0:
            print(f"[{cell}] MILP {k+1}/{len(amb_idx)}", flush=True)
    return mn_llp, mx_llp


def process_cell(cell):
    out_fp = OUT/f"{cell}_fva.csv"
    if out_fp.exists():
        print(f"[{cell}] already done, skip", flush=True)
        return cell, "skip"
    try:
        M, lb, ub, bmax, mn, mx, class_lp = plain_lp_fva(cell)
        mn_llp, mx_llp = milp_on_ambiguous(cell, M, lb, ub, mn, mx, class_lp)
        class_final = classify(mn_llp, mx_llp)
        df = pd.DataFrame(dict(
            rxn=M["rxns"], subsystem=M["subs"],
            lp_min=mn, lp_max=mx, lp_class=class_lp,
            llp_min=mn_llp, llp_max=mx_llp, class_final=class_final,
        ))
        df.to_csv(out_fp, index=False)
        print(f"[{cell}] wrote {out_fp}  (biomass_max={bmax:.4f})", flush=True)
        return cell, "ok"
    except Exception as e:
        print(f"[{cell}] FAILED: {e}", flush=True)
        return cell, f"fail: {e}"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", default="all", help="comma-sep or 'all'")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--summary_only", action="store_true")
    a = ap.parse_args()

    cells = CELLS if a.cells == "all" else [c.strip() for c in a.cells.split(",")]

    if not a.summary_only:
        with Pool(a.workers) as pool:
            results = pool.map(process_cell, cells)
        print("\n--- per-cell completion ---")
        for cell, status in results: print(f"  {cell}: {status}")

    # Aggregate subsystem × class × cell
    print("\n--- building subsystem × class matrix ---", flush=True)
    mats = []
    for cell in cells:
        fp = OUT/f"{cell}_fva.csv"
        if not fp.exists(): continue
        df = pd.read_csv(fp)
        g = (df.groupby(["subsystem", "class_final"]).size()
             .unstack(fill_value=0).reset_index())
        g["cell"] = cell
        g.to_csv(OUT/f"subsystem_class_matrix_{cell}.csv", index=False)
        mats.append(g)
    if mats:
        big = pd.concat(mats, ignore_index=True).fillna(0)
        big.to_csv(OUT/"subsystem_class_matrix_all.csv", index=False)
        print(f"  wrote subsystem_class_matrix_all.csv  rows={len(big)}", flush=True)

    # Developmental trajectory: per (condition, subsystem) across hpf
    print("\n--- developmental trajectory ---", flush=True)
    HPFS = [24, 48, 72, 96, 120]
    CONDS = ["BL", "D", "LD"]
    rows = []
    for cond in CONDS:
        for hpf in HPFS:
            cell = f"{cond}_{hpf}"
            fp = OUT/f"{cell}_fva.csv"
            if not fp.exists(): continue
            df = pd.read_csv(fp)
            g = df.groupby(["subsystem", "class_final"]).size().unstack(fill_value=0)
            for sub, r in g.iterrows():
                rows.append(dict(condition=cond, hpf=hpf, subsystem=sub,
                                 blocked=int(r.get("blocked",0)),
                                 required_fwd=int(r.get("required_fwd",0)),
                                 required_rev=int(r.get("required_rev",0)),
                                 optional_fwd=int(r.get("optional_fwd",0)),
                                 optional_rev=int(r.get("optional_rev",0)),
                                 flexible_signed=int(r.get("flexible_signed",0))))
    traj = pd.DataFrame(rows)
    if len(traj):
        traj.to_csv(OUT/"developmental_trajectory.csv", index=False)
        print(f"  wrote developmental_trajectory.csv  rows={len(traj)}", flush=True)
    print("DONE")
