"""
s192 -- CycleFreeFlux (Desouki 2015) post-processor for the lipid/F1/F2
swap-corner chain CSVs.

Mirrors s51b_cyclefreeflux_gurobi.py exactly: same LP formulation,
same SLACK_TOL=1e-6 on Sv constraints, same identity split
(vp - vn = v[internal]), same L1 objective, Gurobi dual simplex
(Method=1), warm start across samples via variable-bound update.

Per sample:
  - internal reactions: lb/ub set to the smaller / larger of [0, v0] by
    sign of v0 (sign-preserving, L1-minimising)
  - exchanges + biomass (MAR00021): pinned at v0 (lb = ub = v0)
  - Sv0 is right-hand side (± SLACK_TOL)

Per chain:
  - Overwrite the chain CSV in-place with corrected values
  - Log fraction of internal L1 removed (should land near 89.4 % per §2.9)
  - Hard abort on any LP failure (matches s51b policy)

Usage:
  python s192_cff_swap_corners.py --chain path/to/one_chain.csv
  python s192_cff_swap_corners.py --all              # all CSVs in chains/
  python s192_cff_swap_corners.py --all --jobs 3     # 3-wide parallel

Chain filename convention: {model_cell}_mdl__{env_cell}_env_ch{N}.csv
"""
import sys, time, argparse, re
from pathlib import Path
import numpy as np, pandas as pd
import h5py, scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB
from multiprocessing import Pool

REPO = Path(__file__).resolve().parents[2]
V7 = REPO / "reviewer_packet_v4_slim/v7_models"
CHAIN_DIR = REPO / "grace_v4/results/v7_final_addenda/lipid_relay_samples/chains"
STATS_FP = REPO / "grace_v4/results/v7_final_addenda/s192_cff_stats.csv"

BIOMASS_RXN = "MAR00021"
SLACK_TOL = 1e-6


def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)


def _load_model(fp):
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = int(ir.max()) + 1 if len(ir) else 0
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
    return rxns, S


_FN_RE = re.compile(r"^(?P<mdl>[A-Z]+_\d+)_mdl__(?P<env>[A-Z]+_\d+)_env_ch(?P<ch>\d+)\.csv$")


def cff_one_chain(chain_fp):
    """Apply CFF to one chain CSV in place. Returns stats dict."""
    chain_fp = Path(chain_fp)
    m = _FN_RE.match(chain_fp.name)
    if not m:
        raise ValueError(f"bad filename: {chain_fp.name}")
    model_cell = m.group("mdl")
    env_cell = m.group("env")
    chain_id = int(m.group("ch"))

    t_start = time.time()
    mat_fp = V7 / f"trans_rfastcormics_{model_cell}.mat"
    rxns, S = _load_model(mat_fp)
    n_rxn, n_met = len(rxns), S.shape[0]
    idx = {r: i for i, r in enumerate(rxns)}

    is_ex = (np.diff(S.tocsc().indptr) == 1)
    ex_or_bio = is_ex.copy(); ex_or_bio[idx[BIOMASS_RXN]] = True
    internal = ~ex_or_bio
    int_idx = np.where(internal)[0]
    n_int = int(internal.sum())

    df_in = pd.read_csv(chain_fp).set_index("rxn").reindex(rxns)
    V_in = df_in.values.astype(float)
    n_smp = V_in.shape[1]
    V_out = np.zeros_like(V_in)

    # Build CFF model once, reuse via bound updates
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    mo = gp.Model("cff", env=env)
    mo.setParam("Threads", 1)
    mo.setParam("OutputFlag", 0)
    mo.setParam("Method", 1)
    mo.setParam("FeasibilityTol", 1e-7)
    mo.setParam("OptimalityTol", 1e-7)

    v  = mo.addMVar(n_rxn, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="v")
    vp = mo.addMVar(n_int, lb=0.0, ub=GRB.INFINITY, name="vp")
    vn = mo.addMVar(n_int, lb=0.0, ub=GRB.INFINITY, name="vn")
    mo.setObjective(vp.sum() + vn.sum(), GRB.MINIMIZE)

    S_csr = S.tocsr()
    sv_le = mo.addMConstr(S_csr, v, '<', np.zeros(n_met))
    sv_ge = mo.addMConstr(S_csr, v, '>', np.zeros(n_met))

    I_int = sp.coo_matrix(
        (np.ones(n_int), (np.arange(n_int), int_idx)),
        shape=(n_int, n_rxn)).tocsr()
    mo.addConstr(I_int @ v - vp + vn == 0, name="l1_split")
    mo.update()

    S_dense = S.toarray()
    n_fail = 0
    t_last = time.time()
    for s in range(n_smp):
        v0 = V_in[:, s]
        lb_v = np.zeros(n_rxn); ub_v = np.zeros(n_rxn)
        pos_mask = v0 >= 0
        m_ipos = internal & pos_mask
        lb_v[m_ipos] = 0.0
        ub_v[m_ipos] = v0[m_ipos]
        m_ineg = internal & ~pos_mask
        lb_v[m_ineg] = v0[m_ineg]
        ub_v[m_ineg] = 0.0
        lb_v[ex_or_bio] = v0[ex_or_bio]
        ub_v[ex_or_bio] = v0[ex_or_bio]
        v.LB = lb_v
        v.UB = ub_v
        Sv0 = S_dense @ v0
        sv_le.RHS = Sv0 + SLACK_TOL
        sv_ge.RHS = Sv0 - SLACK_TOL
        mo.optimize()
        if mo.Status == GRB.OPTIMAL:
            V_out[:, s] = v.X
        else:
            n_fail += 1
            V_out[:, s] = np.nan
        if (s + 1) % 100 == 0 and __name__ == "__main__":
            dt = time.time() - t_last; t_last = time.time()
            rate = 100 / dt if dt > 0 else 0
            print(f"  {chain_fp.name}: {s+1}/{n_smp}  ({rate:.1f}/s, fails={n_fail})",
                  flush=True)

    mo.dispose(); env.dispose()

    if n_fail > 0:
        return dict(chain_fp=str(chain_fp), model_cell=model_cell,
                    env_cell=env_cell, chain=chain_id,
                    status=f"FAIL_{n_fail}_LPs", elapsed_min=(time.time()-t_start)/60)

    red_before = float(np.abs(V_in[internal, :]).sum())
    red_after  = float(np.abs(V_out[internal, :]).sum())
    removed_pct = 100 * (1 - red_after / red_before) if red_before > 0 else 0.0
    # Overwrite in-place
    df_out = pd.DataFrame(V_out, index=rxns, columns=df_in.columns)
    df_out.index.name = "rxn"
    df_out.to_csv(chain_fp, float_format="%.6g")
    return dict(
        chain_fp=str(chain_fp),
        model_cell=model_cell, env_cell=env_cell, chain=chain_id,
        status="OK",
        internal_L1_before=red_before,
        internal_L1_after=red_after,
        removed_pct=removed_pct,
        n_smp=n_smp,
        elapsed_min=(time.time() - t_start) / 60,
    )


def _worker(chain_fp):
    try:
        return cff_one_chain(chain_fp)
    except Exception as e:
        return dict(chain_fp=str(chain_fp), status=f"ERROR: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", help="apply CFF to one chain CSV")
    ap.add_argument("--all", action="store_true", help="apply CFF to every chain CSV in chains/")
    ap.add_argument("--jobs", type=int, default=1, help="parallel workers")
    args = ap.parse_args()
    if args.chain:
        print(cff_one_chain(args.chain))
        return
    if not args.all:
        ap.error("provide --chain or --all")
    chain_fps = sorted(CHAIN_DIR.glob("*_mdl__*_env_ch*.csv"))
    if not chain_fps:
        print(f"no chain CSVs found under {CHAIN_DIR}")
        return
    print(f"[s192] applying CFF to {len(chain_fps)} chain CSVs, {args.jobs} parallel")
    t0 = time.time()
    rows = []
    if args.jobs <= 1:
        for fp in chain_fps:
            r = _worker(fp); rows.append(r)
            print(f"  {Path(fp).name}: {r.get('status')}  "
                  f"removed={r.get('removed_pct',0):.2f}%  "
                  f"elapsed={r.get('elapsed_min',0):.1f}m", flush=True)
    else:
        with Pool(args.jobs) as pool:
            for r in pool.imap_unordered(_worker, chain_fps):
                rows.append(r)
                print(f"  {Path(r['chain_fp']).name}: {r.get('status')}  "
                      f"removed={r.get('removed_pct',0):.2f}%  "
                      f"elapsed={r.get('elapsed_min',0):.1f}m", flush=True)
    pd.DataFrame(rows).to_csv(STATS_FP, index=False)
    print(f"\n[write] {STATS_FP}  ({len(rows)} rows)  "
          f"total={(time.time()-t0)/60:.1f}m")


if __name__ == "__main__":
    main()
