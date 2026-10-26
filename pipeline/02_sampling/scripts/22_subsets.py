"""
s52 -- Compute enzyme subsets AFTER sampling from S_full null space, purely
       for N_eff bookkeeping in the two-tier null test.

Per reviewer 2026-09-25: subsets never enter the sampling path. They only
label reactions for N_eff.

Usage:
    python s52_post_sample_subsets.py --cell LD_72
    python s52_post_sample_subsets.py --cell all
"""
import argparse, pickle
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from pathlib import Path
from scipy.linalg import null_space
from scipy.sparse.csgraph import connected_components

p = argparse.ArgumentParser()
p.add_argument("--cell", default="LD_72")
p.add_argument("--cos_thresh", type=float, default=0.9999)
p.add_argument("--models_dir", default="models")
p.add_argument("--outdir_name", default="subsets_v4")
args = p.parse_args()

repo = Path.cwd()
outdir = repo/"results"/args.outdir_name
outdir.mkdir(parents=True, exist_ok=True)

def process(cell):
    print(f"\n=== {cell} ===")
    mat_fp = repo/args.models_dir/f"trans_rfastcormics_{cell}.mat"
    with h5py.File(mat_fp,"r") as f:
        m = f["parent"]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        n_col = len(jc)-1
        S = sp.csc_matrix((data, ir, jc), shape=(int(ir.max())+1, n_col))
        rxns = []
        for r in np.asarray(m["rxns"]).ravel():
            arr = np.asarray(f[r]).ravel()
            rxns.append("".join(chr(int(c)) for c in arr))
    print(f"  {len(rxns)} rxns")

    # Null space of S
    N = null_space(S.toarray(), rcond=1e-10)
    print(f"  null-space dim: {N.shape[1]}")

    # Group reactions with proportional rows in N
    rn = np.linalg.norm(N, axis=1)
    act = rn > 1e-9
    Nu = np.zeros_like(N)
    Nu[act] = N[act] / rn[act, None]

    act_i = np.where(act)[0]
    Cmat = np.abs(Nu[act_i] @ Nu[act_i].T)
    adj = (Cmat > args.cos_thresh) & ~np.eye(len(act_i), dtype=bool)
    _, labels = connected_components(sp.csr_matrix(adj), directed=False)

    # Build subset membership + N_eff
    subs_map = pd.read_csv(repo/"data/rxn_subsystem_baked_base.csv")
    sub_of_rxn = dict(zip(subs_map.rxn, subs_map.subsystem))
    gpr_of_rxn = dict(zip(subs_map.rxn, subs_map.has_gpr))

    n_subsets_active = int(labels.max()+1) if len(labels) else 0
    subset_of = {}
    for pos_i, gl in enumerate(act_i):
        subset_of[gl] = int(labels[pos_i])
    next_sub = n_subsets_active
    for gl in np.where(~act)[0]:
        subset_of[gl] = next_sub; next_sub += 1

    members_of = {}
    for gl, sid in subset_of.items():
        members_of.setdefault(sid, []).append(gl)

    reps_of_subset = {sid: sorted(mm)[0] for sid, mm in members_of.items()}

    # expansion coefficients (v_m = c_m * v_rep) computed from N
    coeffs = np.zeros(len(rxns))
    for sid, mm in members_of.items():
        r = reps_of_subset[sid]
        Nr = N[r]; denom = Nr @ Nr
        if denom < 1e-16:
            for mg in mm:
                coeffs[mg] = 1.0 if mg == r else 0.0
            continue
        for mg in mm:
            coeffs[mg] = (N[mg] @ Nr) / denom

    # membership CSV
    rows = []
    for sid, mm in members_of.items():
        r = reps_of_subset[sid]
        for mg in mm:
            rows.append(dict(rxn=rxns[mg], representative=rxns[r], subset_id=sid,
                             subset_size=len(mm), expansion_coef=coeffs[mg],
                             subsystem=sub_of_rxn.get(rxns[mg], "UNMAPPED"),
                             has_gpr=int(gpr_of_rxn.get(rxns[mg], 0))))
    df = pd.DataFrame(rows)
    df.to_csv(outdir/f"{cell}_subset_membership.csv", index=False)

    # subsystem N vs N_eff table
    sub_neff = df.groupby("subsystem").agg(
        N=("rxn","count"),
        N_eff=("subset_id","nunique")).sort_values("N", ascending=False)
    sub_neff.to_csv(outdir/f"{cell}_subsystem_neff.csv")

    print(f"  wrote subset_membership and subsystem_neff for {cell} "
          f"({len(rxns)} -> {len(members_of)} subsets)")
    return True

if args.cell == "all":
    for c in ["BL_24","BL_48","BL_72","BL_96","BL_120",
              "D_24","D_48","D_72","D_96","D_120",
              "LD_24","LD_48","LD_72","LD_96","LD_120"]:
        try: process(c)
        except Exception as e: print(f"  {c}: ERROR {e}")
else:
    process(args.cell)
