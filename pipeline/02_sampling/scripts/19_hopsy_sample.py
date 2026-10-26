"""
s41d -- DIRECT hopsy sampling of the FULL rFASTCORE polytope.

Grace-parallel version: supports three modes to enable per-chain SLURM tasks
- --preprocess_only: build + cache the rounded polytope, exit before sampling
- --chain_only N:    load cached polytope, sample ONE chain (N), write one CSV
- (default):        do everything (preprocess + all chains), Mac-style single-process

For maximum-parallel on Grace, use --preprocess_only first (15 SLURM tasks),
then --chain_only N (60 SLURM tasks) using seeds args.seed*1000 + N - 1.
Finally s41d_aggregate.py computes assertions + diagnostics + Rhat.
"""
import argparse, pickle, time
import numpy as np, pandas as pd
import h5py, hopsy, scipy.sparse as sp
import scipy.linalg
from pathlib import Path
from scipy.optimize import linprog

# SVD fallback (numpy 2.x default 'gesdd' fails on some stoichiometries; scipy 'gesvd' is robust)
_orig_svd = np.linalg.svd
def _robust_svd(a, *args, **kwargs):
    try:
        return _orig_svd(a, *args, **kwargs)
    except np.linalg.LinAlgError as e:
        if "did not converge" in str(e):
            print(f"    [SVD fallback] numpy SVD failed, using scipy lapack_driver='gesvd'", flush=True)
            full_matrices = kwargs.get("full_matrices", True)
            compute_uv = kwargs.get("compute_uv", True)
            return scipy.linalg.svd(a, full_matrices=full_matrices,
                                    compute_uv=compute_uv, lapack_driver="gesvd")
        raise
np.linalg.svd = _robust_svd

from PolyRound.mutable_classes.polytope import Polytope
from PolyRound.api import PolyRoundApi

OCR_MAP = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}
BOUND_EX = 10
BIO_LB_FRAC = 0.9
O2_EX = "MAR09048"
BIOMASS_RXN = "MAR00021"

p = argparse.ArgumentParser()
p.add_argument("--cell", default="LD_72")
p.add_argument("--n_samples", type=int, default=500)
p.add_argument("--n_chains", type=int, default=4)
p.add_argument("--thinning", type=int, default=300)
p.add_argument("--sampler", default="billiard_walk",
               choices=["uniform_coord_hit_and_run", "billiard_walk"])
p.add_argument("--seed", type=int, default=1)
p.add_argument("--biomass_lb", type=float, default=None)
p.add_argument("--outdir_name", default="hopsy_v4")
p.add_argument("--models_dir", default="models",
               help="dir under repo root containing trans_rfastcormics_{cell}.mat")
p.add_argument("--preprocess_only", action="store_true")
p.add_argument("--chain_only", type=int, default=None)
p.add_argument("--skip_simplify", action="store_true",
               help="skip PolyRound simplify_polytope; use if constraint_removal hits ValueError")
args = p.parse_args()

cond, hpf = args.cell.split("_"); hpf = int(hpf)

repo = Path.cwd()
mat_fp = repo/args.models_dir/f"trans_rfastcormics_{args.cell}.mat"
bnd_fp = repo/"data/phase2_exchange_bounds_long.csv"
outdir = repo/"results"/args.outdir_name
outdir.mkdir(parents=True, exist_ok=True)
cache_fp = outdir/f"{args.cell}_polytope_cache.pkl"
meta_fp  = outdir/f"{args.cell}_polytope_meta.pkl"

def load_v73(mat_fp):
    with h5py.File(mat_fp, "r") as f:
        m = f["parent"]
        def deref_str(ref):
            out = []
            for r in np.asarray(ref).ravel():
                obj = f[r]; arr = np.asarray(obj).ravel()
                out.append("".join(chr(int(c)) for c in arr))
            return out
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        return dict(
            S=sp.csc_matrix((data, ir, jc), shape=(int(ir.max())+1, len(jc)-1)),
            lb=np.asarray(m["lb"]).ravel().astype(float),
            ub=np.asarray(m["ub"]).ravel().astype(float),
            rxns=deref_str(m["rxns"]),
            mets=deref_str(m["mets"]))

def build_bounded_model():
    # 2026-09-30 CORRECTION: use apply_context_bounds. The prior version
    # reset every exchange to ±BOUND_EX before applying fc, which
    # discarded the v7 boundary (Part 1 secretion uncap, Part 2 closes,
    # ledger closes, medium). apply_context_bounds preserves stored
    # bounds and only edits (1) fc on the 105 measured, (2) OCR cap.
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "functions"))
    from apply_context_bounds import apply_context_bounds

    print(f"[{args.cell}] loading model...", flush=True)
    M = load_v73(mat_fp)
    S = M["S"]; lb = M["lb"].copy(); ub = M["ub"].copy()
    rxn_ids = M["rxns"]; n_rxn = len(rxn_ids); n_met = len(M["mets"])
    Snnz = np.diff(S.tocsc().indptr); is_ex = Snnz == 1
    lb, ub, idx, is_ex = apply_context_bounds(lb, ub, rxn_ids, is_ex, cond, hpf)

    c_vec = np.zeros(n_rxn); c_vec[idx[BIOMASS_RXN]] = -1.0
    S_dense = S.toarray()
    lp = linprog(c_vec, A_eq=S_dense, b_eq=np.zeros(n_met),
                 bounds=list(zip(lb, ub)), method="highs")
    bio_max_full = -lp.fun
    assert bio_max_full > 1e-6
    if args.biomass_lb is None:
        bio_lb = BIO_LB_FRAC * bio_max_full * (1 - 1e-4); frame = "A"
    else:
        bio_lb = args.biomass_lb; frame = "B"
    lb[idx[BIOMASS_RXN]] = bio_lb
    return dict(S=S, S_dense=S_dense, lb=lb, ub=ub, rxn_ids=rxn_ids, idx=idx,
                n_rxn=n_rxn, n_met=n_met, bio_max_full=bio_max_full,
                bio_lb=bio_lb, frame=frame)

def preprocess():
    M = build_bounded_model()
    print(f"[{args.cell}] biomass_max={M['bio_max_full']:.4f} lb={M['bio_lb']:.4f} frame={M['frame']}", flush=True)
    A_ineq = np.vstack([np.eye(M["n_rxn"]), -np.eye(M["n_rxn"])])
    b_ineq = np.concatenate([M["ub"], -M["lb"]])
    polyA = pd.DataFrame(A_ineq, columns=M["rxn_ids"])
    polyb = pd.Series(b_ineq)
    polyS = pd.DataFrame(M["S_dense"], columns=M["rxn_ids"])
    polyh = pd.Series(np.zeros(M["n_met"]))
    polytope = Polytope(polyA, polyb, S=polyS, h=polyh)
    if args.skip_simplify:
        print(f"[{args.cell}] SKIPPING simplify (--skip_simplify); polytope shape {polytope.A.shape}", flush=True)
    else:
        print(f"[{args.cell}] simplifying...", flush=True); t0=time.time()
        try:
            polytope = PolyRoundApi.simplify_polytope(polytope)
            print(f"[{args.cell}]   simplify done in {(time.time()-t0)/60:.1f}m -> {polytope.A.shape}", flush=True)
        except ValueError as e:
            # PolyRound's constraint_removal LP hit a non-optimal status
            # (numerical dead-end in glpk on some polytopes). Fall through
            # to transform+round without simplify -- gives a valid rounded
            # polytope with more constraints; sampling is slower per step
            # but still correct.
            print(f"[{args.cell}]   simplify FAILED after {(time.time()-t0)/60:.1f}m with ValueError; "
                  f"proceeding without constraint removal.", flush=True)
    print(f"[{args.cell}] transforming...", flush=True); t0=time.time()
    polytope = PolyRoundApi.transform_polytope(polytope)
    print(f"[{args.cell}]   transform done in {(time.time()-t0)/60:.1f}m -> reduced dim {polytope.A.shape[1]}", flush=True)
    print(f"[{args.cell}] rounding...", flush=True); t0=time.time()
    polytope = PolyRoundApi.round_polytope(polytope)
    print(f"[{args.cell}]   round done in {(time.time()-t0)/60:.1f}m", flush=True)
    with open(cache_fp,"wb") as fh: pickle.dump(polytope, fh)
    meta = dict(rxn_ids=M["rxn_ids"], idx=M["idx"], n_rxn=M["n_rxn"], n_met=M["n_met"],
                bio_max_full=M["bio_max_full"], bio_lb=M["bio_lb"], frame=M["frame"],
                lb=M["lb"], ub=M["ub"])
    with open(meta_fp,"wb") as fh: pickle.dump(meta, fh)
    print(f"[{args.cell}] cached to {cache_fp}", flush=True)

def sample_chain(chain_id):
    fp_out = outdir/f"flux_samples_{args.cell}_ch{chain_id}.csv"
    if fp_out.exists():
        print(f"[{args.cell}][ch{chain_id}] SKIP: {fp_out} already exists", flush=True)
        return
    with open(cache_fp,"rb") as fh: polytope = pickle.load(fh)
    with open(meta_fp,"rb") as fh: meta = pickle.load(fh)
    poly_prob = hopsy.Problem(A=polytope.A.values, b=polytope.b.values)
    proposal_cls = hopsy.BilliardWalkProposal if args.sampler == "billiard_walk" \
                   else hopsy.UniformCoordinateHitAndRunProposal
    mc = [hopsy.MarkovChain(problem=poly_prob, proposal=proposal_cls)]
    seed = args.seed * 1000 + chain_id - 1
    rng = [hopsy.RandomNumberGenerator(seed=seed)]
    print(f"[{args.cell}][ch{chain_id}] sampling {args.n_samples} samples thin={args.thinning} seed={seed}", flush=True)
    t0 = time.time()
    _, red = hopsy.sample(mc, rng, n_samples=args.n_samples, thinning=args.thinning, n_procs=1)
    print(f"[{args.cell}][ch{chain_id}] sampling done in {(time.time()-t0)/60:.1f}m", flush=True)
    T = polytope.transformation.values
    shift = polytope.shift.values
    full = red[0] @ T.T + shift
    df = pd.DataFrame(full.T, index=meta["rxn_ids"],
                      columns=[f"s{i+1}" for i in range(args.n_samples)])
    df.index.name = "rxn"
    df.to_csv(fp_out, float_format="%.6g")
    print(f"[{args.cell}][ch{chain_id}] wrote {fp_out}", flush=True)

def do_all():
    if not cache_fp.exists():
        preprocess()
    for ch in range(1, args.n_chains + 1):
        sample_chain(ch)
    print(f"[{args.cell}] all chains done. run s41d_aggregate.py --cell {args.cell} for assertions + Rhat.", flush=True)

if args.preprocess_only:
    preprocess()
elif args.chain_only is not None:
    sample_chain(args.chain_only)
else:
    do_all()
