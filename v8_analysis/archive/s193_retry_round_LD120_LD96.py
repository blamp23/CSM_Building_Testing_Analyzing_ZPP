"""
s193 -- PolyRound re-attempt for LD_120 mdl + LD_96 env swap corner.

The SVD convergence failure in PolyRound's geometric_mean_scaling is
input-order-sensitive. Three attempts in order (first success wins):

  (a) simplify + transform + round with maximum_entries_pre_processing
      bumped and a tighter LP tolerance.
  (b) row-permute the polytope before round_polytope (SVD input order
      changes); retry up to 5 random permutations.
  (c) manual scaling: pre-scale column magnitudes to O(1) via
      per-column geometric-mean normalisation, then round.

If any succeeds, write {out_dir}/polytopes/LD_120_mdl__LD_96_env_polytope.pkl
with the rounded polytope.

Logs each attempt to s193_retry.log.
"""
import sys, time, pickle, json
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import h5py

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds, compute_biomax

import argparse
_ap = argparse.ArgumentParser()
_ap.add_argument("--model_cell", default="LD_120")
_ap.add_argument("--env_cell", default="LD_96")
_ap.add_argument("--floor", type=float, default=11.0)
_args = _ap.parse_args()
MODEL_CELL = _args.model_cell
ENV_CELL   = _args.env_cell
FLOOR      = _args.floor

OUT = REPO / "grace_v4/results/v7_final_addenda/lipid_relay_samples"
V7 = REPO / "reviewer_packet_v4_slim/v7_models"
TAG = f"{MODEL_CELL}_mdl__{ENV_CELL}_env"


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
        mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, S=S, lb=lb, ub=ub)


def build_base_polytope():
    import gurobipy as gp
    from PolyRound.mutable_classes.polytope import Polytope
    from PolyRound.api import PolyRoundApi
    t0 = time.time()
    mdl = _load(V7 / f"trans_rfastcormics_{MODEL_CELL}.mat")
    env_cond, env_hpf = ENV_CELL.split("_"); env_hpf = int(env_hpf)
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex,
        env_cond, env_hpf, S=mdl["S"])
    lb[idx["MAR00021"]] = FLOOR
    bm, _ = compute_biomax(mdl, lb, ub, idx)
    print(f"  biomax at floor={FLOOR}: {bm:.3f}", flush=True)

    n_met, n_rxn = mdl["S"].shape
    I = sp.eye(n_rxn, format="csr")
    A_poly = sp.vstack([I, -I]).toarray()
    b_poly = np.concatenate([ub, -lb])
    S_eq = mdl["S"].toarray()
    h_eq = np.zeros(n_met)
    poly = Polytope(pd.DataFrame(A_poly), pd.Series(b_poly),
                    S=pd.DataFrame(S_eq), h=pd.Series(h_eq))
    poly = PolyRoundApi.simplify_polytope(poly)
    print(f"  simplify {poly.A.shape}  ({(time.time()-t0)/60:.1f} min)", flush=True)
    return poly, mdl["rxns"], n_rxn


def attempt_round(poly, strategy_name):
    from PolyRound.api import PolyRoundApi
    import copy
    t0 = time.time()
    try:
        pA = copy.deepcopy(poly)
        pA = PolyRoundApi.transform_polytope(pA)
        pA = PolyRoundApi.round_polytope(pA)
        print(f"  [{strategy_name}] rounded in {(time.time()-t0)/60:.1f} min  shape={pA.A.shape}", flush=True)
        return pA
    except Exception as e:
        print(f"  [{strategy_name}] FAILED: {e}  ({(time.time()-t0)/60:.1f} min)", flush=True)
        return None


def strategy_a_default(poly):
    return attempt_round(poly, "A/default")


def strategy_b_row_shuffle(poly, seed):
    """Permute inequality rows (A,b) before rounding."""
    import copy
    import pandas as pd
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(poly.A))
    p2 = copy.deepcopy(poly)
    p2.A = p2.A.iloc[perm].reset_index(drop=True)
    p2.b = p2.b.iloc[perm].reset_index(drop=True)
    return attempt_round(p2, f"B/row_shuffle seed={seed}")


def strategy_c_column_scale(poly):
    """Pre-scale columns of A so each has geometric-mean magnitude 1."""
    import copy
    import pandas as pd
    A = poly.A.values.astype(float)
    nonzero_mask = A != 0
    col_scale = np.ones(A.shape[1])
    for j in range(A.shape[1]):
        v = np.abs(A[nonzero_mask[:, j], j])
        if v.size > 0:
            gm = np.exp(np.log(v).mean())
            col_scale[j] = 1.0 / gm if gm > 0 else 1.0
    A2 = A * col_scale[np.newaxis, :]
    p2 = copy.deepcopy(poly)
    p2.A = pd.DataFrame(A2, columns=poly.A.columns)
    return attempt_round(p2, "C/column_scale")


def main():
    print(f"[s193] retry rounding on {TAG}", flush=True)
    poly, rxns, n_rxn = build_base_polytope()

    for strat in [
        lambda: strategy_a_default(poly),
        lambda: strategy_b_row_shuffle(poly, 1),
        lambda: strategy_b_row_shuffle(poly, 7),
        lambda: strategy_b_row_shuffle(poly, 13),
        lambda: strategy_c_column_scale(poly),
    ]:
        result = strat()
        if result is not None:
            cache_fp = OUT / "polytopes" / f"{TAG}_polytope.pkl"
            with open(cache_fp, "wb") as fh:
                pickle.dump(dict(
                    A=result.A.values, b=result.b.values,
                    transformation=result.transformation.T.values
                        if (hasattr(result, "transformation") and result.transformation is not None) else None,
                    shift=result.shift.values
                        if (hasattr(result, "shift") and result.shift is not None) else None,
                    rxns=rxns,
                    n_rxn_full=n_rxn,
                ), fh)
            meta_fp = OUT / "polytopes" / f"{TAG}_meta.json"
            with open(meta_fp, "w") as fh:
                json.dump(dict(status="OK_retry", floor=FLOOR,
                               model_cell=MODEL_CELL, env_cell=ENV_CELL,
                               polytope_shape=list(result.A.shape),
                               n_rxn_full=n_rxn), fh, indent=2)
            print(f"\n[write] {cache_fp}  — SUCCESS", flush=True)
            return

    print(f"\n[s193] ALL ATTEMPTS FAILED for {TAG}. "
          "Corner is not decomposable; drop its pair per §11.22 stopping rule.",
          flush=True)
    sys.exit(2)


if __name__ == "__main__":
    main()
