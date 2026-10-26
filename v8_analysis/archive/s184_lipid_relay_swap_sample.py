"""
s184 -- lipid-relay swap-corner sampler, in two modes:

  --preprocess_only : build polytope locally (PolyRound simplify + round),
                      pickle to {out_dir}/polytopes/{tag}_polytope.pkl.
                      Needs Gurobi + PolyRound (local Mac).
  --chain_only N    : load cached polytope, run ONE hopsy chain, write
                      {out_dir}/chains/{tag}_ch{N}.csv.  No Gurobi needed
                      beyond what the preprocess already paid for.
                      Grace venv compatible (hopsy + numpy + pandas).

Four chains per corner can run in parallel across a Grace array:
  sbatch --array=1-56 run_lipid_relay.slurm
where 56 = 14 corners × 4 chains.

Default sampling: 500 samples × thinning 300 per chain. Final pool is
written locally by s186_lipid_relay_decompose.py which reads the four
chain CSVs per corner.

Internal ±1000 Wang defaults INTACT (match frozen pooled_v7). Matched
floor passed in via --floor.
"""
import sys, json, time, argparse, pickle
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import h5py

REPO = Path(__file__).resolve().parents[2]


def _import_poly_libs():
    """Preprocess-only imports. Needs Gurobi + PolyRound + apply_context_bounds.
    On Grace in chain-only mode these are never called."""
    # Try both local repo layout (parents[2]/functions) and Grace layout
    # (../functions relative to the batch directory).
    for cand in (REPO / "functions",
                 Path(__file__).resolve().parent.parent / "functions"):
        if (cand / "apply_context_bounds.py").exists():
            sys.path.insert(0, str(cand))
            break
    import gurobipy as gp
    from PolyRound.mutable_classes.polytope import Polytope
    from PolyRound.api import PolyRoundApi
    from apply_context_bounds import apply_context_bounds, compute_biomax
    return gp, Polytope, PolyRoundApi, apply_context_bounds, compute_biomax


def _import_hopsy():
    import hopsy
    return hopsy


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
        mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, S=S, lb=lb, ub=ub)


def preprocess(args):
    gp, Polytope, PolyRoundApi, apply_context_bounds, compute_biomax = _import_poly_libs()
    tag = f"{args.model_cell}_mdl__{args.env_cell}_env"
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    (out / "polytopes").mkdir(exist_ok=True)
    cache_fp = out / "polytopes" / f"{tag}_polytope.pkl"
    meta_fp  = out / "polytopes" / f"{tag}_meta.json"
    if cache_fp.exists() and not args.force:
        print(f"[cache hit] {cache_fp}", flush=True)
        return

    t0 = time.time()
    mdl_fp = REPO / f"reviewer_packet_v4_slim/v7_models/trans_rfastcormics_{args.model_cell}.mat"
    mdl = _load_model(mdl_fp)
    print(f"[prep] {tag} rxns={len(mdl['rxns'])}", flush=True)

    env_cond, env_hpf = args.env_cell.split("_"); env_hpf = int(env_hpf)
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex,
        env_cond, env_hpf, S=mdl["S"])
    # Internal ±1000 defaults INTACT. Pin matched floor.
    lb[idx["MAR00021"]] = args.floor
    bm, _ = compute_biomax(mdl, lb, ub, idx)
    if not np.isfinite(bm) or bm < args.floor - 1e-3:
        with open(meta_fp, "w") as fh:
            json.dump(dict(status="INFEASIBLE", biomax=bm, floor=args.floor,
                           model_cell=args.model_cell, env_cell=args.env_cell), fh)
        print(f"  INFEASIBLE at floor={args.floor} (biomax={bm}); skipping", flush=True)
        return

    n_met, n_rxn = mdl["S"].shape
    I = sp.eye(n_rxn, format="csr")
    A_poly = sp.vstack([I, -I]).toarray()
    b_poly = np.concatenate([ub, -lb])
    S_eq = mdl["S"].toarray()
    h_eq = np.zeros(n_met)
    polytope = Polytope(pd.DataFrame(A_poly), pd.Series(b_poly),
                        S=pd.DataFrame(S_eq), h=pd.Series(h_eq))

    try:
        polytope = PolyRoundApi.simplify_polytope(polytope)
        print(f"  simplify {polytope.A.shape}  ({(time.time()-t0)/60:.1f} min)", flush=True)
    except Exception as e:
        with open(meta_fp, "w") as fh:
            json.dump(dict(status=f"SIMPLIFY_FAIL: {e}", floor=args.floor,
                           model_cell=args.model_cell, env_cell=args.env_cell), fh)
        print(f"  SIMPLIFY_FAIL: {e}", flush=True)
        return
    try:
        polytope = PolyRoundApi.transform_polytope(polytope)
    except Exception as e:
        print(f"  TRANSFORM failed, using simplified polytope: {e}", flush=True)
    rounded = False
    try:
        polytope = PolyRoundApi.round_polytope(polytope)
        rounded = True
        print(f"  round    {polytope.A.shape}  ({(time.time()-t0)/60:.1f} min)", flush=True)
    except Exception as e:
        print(f"  ROUND failed, falling back to transformed polytope: {e}", flush=True)
        print(f"  (sampling will be less efficient but distribution is correct)", flush=True)

    with open(cache_fp, "wb") as fh:
        pickle.dump(dict(
            A=polytope.A.values, b=polytope.b.values,
            transformation=polytope.transformation.T.values
                if (hasattr(polytope, "transformation") and polytope.transformation is not None) else None,
            shift=polytope.shift.values
                if (hasattr(polytope, "shift") and polytope.shift is not None) else None,
            rxns=mdl["rxns"],
            n_rxn_full=n_rxn,
        ), fh)
    with open(meta_fp, "w") as fh:
        json.dump(dict(status="OK", biomax=bm, floor=args.floor,
                       model_cell=args.model_cell, env_cell=args.env_cell,
                       polytope_shape=list(polytope.A.shape),
                       n_rxn_full=n_rxn,
                       rounded=rounded,
                       elapsed_min=(time.time()-t0)/60), fh, indent=2)
    print(f"[prep] {tag} DONE in {(time.time()-t0)/60:.1f} min", flush=True)


def chain_only(args):
    hopsy = _import_hopsy()
    tag = f"{args.model_cell}_mdl__{args.env_cell}_env"
    out = Path(args.out_dir)
    (out / "chains").mkdir(parents=True, exist_ok=True)
    cache_fp = out / "polytopes" / f"{tag}_polytope.pkl"
    chain_fp = out / "chains" / f"{tag}_ch{args.chain_only}.csv"
    if chain_fp.exists() and not args.force:
        print(f"[chain hit] {chain_fp}", flush=True); return
    if not cache_fp.exists():
        raise FileNotFoundError(f"no polytope cache: {cache_fp}")
    with open(cache_fp, "rb") as fh:
        cache = pickle.load(fh)

    t0 = time.time()
    A_hp, b_hp = cache["A"], cache["b"]
    T = cache["transformation"]; shift = cache["shift"]
    poly_prob = hopsy.Problem(A=A_hp, b=b_hp)
    try:
        cc = np.asarray(hopsy.compute_chebyshev_center(poly_prob)).ravel()
    except Exception:
        cc = np.zeros(A_hp.shape[1])
    slack = b_hp - A_hp @ cc
    if (slack <= 0).any():
        for shrink in (0.99, 0.95, 0.9, 0.75, 0.5, 0.25, 0.1, 0.0):
            cc_try = cc * shrink
            if ((b_hp - A_hp @ cc_try) > 1e-9).all():
                cc = cc_try; break

    # Match the frozen s41d call style exactly: wrap mc and rng as single-
    # element lists and pass n_procs=1. hopsy.sample's C++ backend has a
    # faster path for the list form.
    seed = args.seed * 1000 + args.chain_only - 1
    rng = [hopsy.RandomNumberGenerator(seed=seed)]
    mc = [hopsy.MarkovChain(problem=poly_prob,
                            proposal=hopsy.BilliardWalkProposal,
                            starting_point=cc)]
    _, states = hopsy.sample(mc, rng, n_samples=args.n_samples,
                             thinning=args.thinning, n_procs=1)
    S_rounded = np.asarray(states[0])   # (n_samples, n_rxn_rounded)
    if T is not None:
        full = (S_rounded @ T) + (shift if shift is not None else 0)
    else:
        full = S_rounded
    # One CSV per chain: rows=reactions, columns = ch{N}_s{1..n_samples}
    out_df = pd.DataFrame(full.T, index=cache["rxns"],
                          columns=[f"ch{args.chain_only}_s{j+1}"
                                   for j in range(full.shape[0])])
    out_df.index.name = "rxn"
    out_df.to_csv(chain_fp)
    print(f"[chain] {tag} ch{args.chain_only} DONE in {(time.time()-t0)/60:.1f} min "
          f"({full.shape[0]} samples × {full.shape[1]} rxns)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model_cell", required=True)
    ap.add_argument("--env_cell", required=True)
    ap.add_argument("--floor", type=float, required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--preprocess_only", action="store_true")
    ap.add_argument("--chain_only", type=int, default=0, help="1..n_chains")
    ap.add_argument("--n_chains", type=int, default=4)
    ap.add_argument("--n_samples", type=int, default=500)
    ap.add_argument("--thinning", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if args.preprocess_only:
        preprocess(args)
    elif args.chain_only > 0:
        chain_only(args)
    else:
        preprocess(args)
        for ch in range(1, args.n_chains + 1):
            args.chain_only = ch
            chain_only(args)


if __name__ == "__main__":
    main()
