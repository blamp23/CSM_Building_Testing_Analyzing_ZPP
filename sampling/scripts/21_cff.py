"""
s51 -- CycleFreeFlux (Desouki et al 2015) using scipy HiGHS (no Gurobi needed).

Per Desouki eq (2):
    minimize sum |v_i| over internal reactions
    subject to Sv = 0
               internal: 0 <= v_i <= v_i^(0) if v_i^(0) >= 0
                         v_i^(0) <= v_i <= 0 if v_i^(0) < 0
               exchange & biomass: pinned at v^(0)

Uses scipy.optimize.linprog with method='highs' (bundled). No solver license.

Usage:
    python s51_cyclefreeflux.py --cell LD_72
    python s51_cyclefreeflux.py --cell all
"""
import argparse, pickle, time
import numpy as np, pandas as pd
import h5py, scipy.sparse as sp
from pathlib import Path
from scipy.optimize import linprog

OCR_MAP = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}
BOUND_EX = 10
BIO_LB_FRAC = 0.9
O2_EX = "MAR09048"
BIOMASS_RXN = "MAR00021"

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
            rxns=deref_str(m["rxns"]))

def run_cell(cell, chains=(1,2,3,4)):
    repo = Path.cwd()
    print(f"\n=== CFF (HiGHS): {cell} ===", flush=True); t_all=time.time()
    cond, hpf = cell.split("_"); hpf = int(hpf)

    mat_fp = repo/args.models_dir/f"trans_rfastcormics_{cell}.mat"
    M = load_v73(mat_fp)
    S = M["S"]; lb = M["lb"].copy(); ub = M["ub"].copy()
    rxn_ids = M["rxns"]; n_rxn = len(rxn_ids); n_met = S.shape[0]
    idx = {r:i for i,r in enumerate(rxn_ids)}
    Snnz = np.diff(S.tocsc().indptr); is_ex = Snnz == 1

    # Apply s01/s17 bounds pipeline
    lb[is_ex] = -BOUND_EX; ub[is_ex] = BOUND_EX
    mb = pd.read_csv(repo/"data/phase2_exchange_bounds_long.csv")
    sel = (mb.condition == cond) & (mb.hpf == hpf)
    for rid, fc in zip(mb.loc[sel,"ex_rxn"], mb.loc[sel,"fc"]):
        if rid in idx:
            j = idx[rid]; fc = float(fc) if np.isfinite(fc) and fc>0 else 1.0
            lb[j] = -BOUND_EX*fc; ub[j] = BOUND_EX*fc
    lb[idx[O2_EX]] = -OCR_MAP[hpf]; ub[idx[O2_EX]] = 0

    # biomass_max + biomass_lb
    c = np.zeros(n_rxn); c[idx[BIOMASS_RXN]] = -1.0
    S_dense = S.toarray()
    lp = linprog(c, A_eq=S_dense, b_eq=np.zeros(n_met),
                 bounds=list(zip(lb,ub)), method="highs")
    bio_max = -lp.fun
    bio_lb = BIO_LB_FRAC * bio_max * (1 - 1e-4)
    lb[idx[BIOMASS_RXN]] = bio_lb
    print(f"  bio_max={bio_max:.3f} bio_lb={bio_lb:.3f}", flush=True)

    j_bio = idx[BIOMASS_RXN]
    ex_or_bio = is_ex.copy(); ex_or_bio[j_bio] = True
    internal = ~ex_or_bio
    n_int = int(internal.sum())
    int_idx = np.where(internal)[0]
    print(f"  {n_int} L1-min internal (biomass+exchanges pinned)", flush=True)

    # For HiGHS: reformulate |v_i| minimization as sum of vp + vn per internal,
    # with v_i = vp - vn, vp,vn >= 0. Variables: [v (n_rxn); vp (n_int); vn (n_int)]
    n_var = n_rxn + 2*n_int
    c_obj = np.zeros(n_var); c_obj[n_rxn:n_rxn+n_int] = 1.0; c_obj[n_rxn+n_int:] = 1.0

    # Constraint structure (FIX 2, 2026-09-28 per rerun failures):
    #   S*v - (S*v0) in [-tol, tol]      -> A_ub +/- slack   (2*n_met rows)
    #   v[int_j] - vp[k] + vn[k] = 0     -> A_eq              (n_int rows)
    # Slack replaces strict equality against v0's residual: HiGHS's default
    # primal feasibility tolerance (1e-7) was rejecting ~1% of samples where
    # S*v0 had entries close to numerical edge cases. Slack tol=1e-6 is
    # bigger than typical numerical noise but 1000x smaller than any flux we
    # care about, so removes false infeasibilities without loosening the
    # loopless problem.
    S_coo = S.tocoo()
    # A_ub: first n_met rows for  S*v - Sv0 <= tol;
    #       next n_met rows for -(S*v - Sv0) <= tol
    A_ub_data = np.concatenate([S_coo.data, -S_coo.data])
    A_ub_rows = np.concatenate([S_coo.row, S_coo.row + n_met])
    A_ub_cols = np.concatenate([S_coo.col, S_coo.col])
    A_ub = sp.csr_matrix((A_ub_data, (A_ub_rows, A_ub_cols)),
                          shape=(2*n_met, n_var))
    # A_eq: only the vp/vn identity
    A_eq_rows = []; A_eq_cols = []; A_eq_data = []
    for k, j in enumerate(int_idx):
        A_eq_rows.append(np.array([k, k, k]))
        A_eq_cols.append(np.array([j, n_rxn+k, n_rxn+n_int+k]))
        A_eq_data.append(np.array([1.0, -1.0, 1.0]))
    A_eq = sp.csr_matrix(
        (np.concatenate(A_eq_data),
         (np.concatenate(A_eq_rows), np.concatenate(A_eq_cols))),
        shape=(n_int, n_var))
    b_eq = np.zeros(n_int)
    SLACK_TOL = 1e-6

    # Loop over samples: build bounds from v^(0), solve, keep only v-part.
    # FIX (2026-09-28, per reviewer diagnostic): b_eq for the S*v block is
    # S*v0, not zero -- input samples satisfy S*v = 0 only to sampler LP
    # precision (residual ~1e-3), so an exact-zero constraint is infeasible.
    # The LP is now: minimize L1 of internal v s.t.
    #   S*v = S*v0  (i.e. preserve v0's residual exactly),
    #   exchanges + biomass pinned at v0,
    #   internal in sign-preserving [0, v0] or [v0, 0].
    # NO silent fallback: any LP failure aborts the cell with a status log.
    outdir = repo/"results"/args.outdir_out
    indir  = repo/"results"/args.outdir_in
    outdir.mkdir(parents=True, exist_ok=True)
    stats_rows = []

    # precompute S dense for O(n_met) matvec on each v0
    S_dense_full = S.toarray()

    for ch in chains:
        fp_in  = indir/f"flux_samples_{cell}_ch{ch}.csv"
        fp_out = outdir/f"flux_samples_{cell}_ch{ch}.csv"
        if not fp_in.exists():
            raise RuntimeError(f"chain {ch}: missing input {fp_in}")
        if fp_out.exists():
            print(f"  chain {ch}: OUTPUT EXISTS, skip", flush=True); continue
        df_in = pd.read_csv(fp_in, index_col=0).reindex(rxn_ids)
        V_in = df_in.values.astype(float)
        n_smp = V_in.shape[1]
        V_out = np.zeros_like(V_in)
        n_fail = 0
        fail_first_sample = -1
        fail_last_status = None
        fail_last_message = None
        t_ch = time.time()
        for s in range(n_smp):
            v0 = V_in[:, s]
            lb_v = np.zeros(n_rxn); ub_v = np.zeros(n_rxn)
            # non-biomass internals: sign-preserving [min(v0,0), max(v0,0)]
            for j in int_idx:
                if v0[j] >= 0:
                    lb_v[j] = 0.0; ub_v[j] = max(v0[j], 0.0)
                else:
                    lb_v[j] = min(v0[j], 0.0); ub_v[j] = 0.0
            # pin exchanges + biomass at v0
            for j in np.where(ex_or_bio)[0]:
                lb_v[j] = v0[j]; ub_v[j] = v0[j]
            bounds = list(zip(lb_v, ub_v)) + \
                     [(0.0, None)]*n_int + [(0.0, None)]*n_int
            # FIX 2: slack the S constraint by +/- SLACK_TOL around S*v0.
            Sv0 = S_dense_full @ v0
            b_ub = np.concatenate([Sv0 + SLACK_TOL, -Sv0 + SLACK_TOL])
            res = linprog(c_obj, A_ub=A_ub, b_ub=b_ub,
                          A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                          method="highs")
            if not res.success:
                n_fail += 1
                if fail_first_sample < 0:
                    fail_first_sample = s
                fail_last_status = res.status
                fail_last_message = res.message
                # Do NOT fall back. Record placeholder; will abort cell below.
                V_out[:, s] = np.nan
            else:
                V_out[:, s] = res.x[:n_rxn]
            if (s+1) % 100 == 0:
                el = time.time()-t_ch
                print(f"  chain {ch}: {s+1}/{n_smp} ({(s+1)/el:.1f}/s, ETA {(n_smp-s-1)/((s+1)/el):.0f}s, fails={n_fail})", flush=True)

        if n_fail > 0:
            # HARD ABORT: reviewer 2026-09-28: "count failures, log the status,
            # and abort the cell if any exceed zero."
            raise RuntimeError(
                f"[{cell} ch{ch}] CFF LP failed on {n_fail}/{n_smp} samples. "
                f"First failure at sample {fail_first_sample}. "
                f"Last status={fail_last_status} message={fail_last_message!r}. "
                f"NO fallback used. Cell aborted.")

        red_before = float(np.abs(V_in[internal,:]).sum())
        red_after  = float(np.abs(V_out[internal,:]).sum())
        print(f"  chain {ch} done in {(time.time()-t_ch)/60:.1f}m. "
              f"internal L1: {red_before:.4g} -> {red_after:.4g} "
              f"({100*(1-red_after/red_before):.2f}% removed) "
              f"[fails=0/{n_smp}]", flush=True)
        df_out = pd.DataFrame(V_out, index=rxn_ids, columns=df_in.columns)
        df_out.index.name = "rxn"
        df_out.to_csv(fp_out, float_format="%.6g")
        stats_rows.append(dict(cell=cell, chain=ch,
            internal_L1_before=red_before, internal_L1_after=red_after,
            removed_pct=100*(1-red_after/red_before),
            n_fail=0, n_smp=n_smp))
    if stats_rows:
        pd.DataFrame(stats_rows).to_csv(outdir/f"{cell}_loop_stats.csv", index=False)
    print(f"  {cell} total: {(time.time()-t_all)/60:.1f}m", flush=True)

p = argparse.ArgumentParser()
p.add_argument("--cell", default="LD_72")
p.add_argument("--outdir_in",  default="hopsy_v4")
p.add_argument("--outdir_out", default="hopsy_loopless")
p.add_argument("--models_dir", default="models")
p.add_argument("--chains", default="1,2,3,4",
               help="comma-separated chain ids; use to parallelize across SLURM tasks")
args = p.parse_args()
CHAINS = tuple(int(c) for c in args.chains.split(","))
if args.cell == "all":
    for c in ["BL_24","BL_48","BL_72","BL_96","BL_120",
              "D_24","D_48","D_72","D_96","D_120",
              "LD_24","LD_48","LD_72","LD_96","LD_120"]:
        try: run_cell(c, chains=CHAINS)
        except Exception as e: print(f"  {c}: ERROR {e}")
else:
    run_cell(args.cell, chains=CHAINS)
