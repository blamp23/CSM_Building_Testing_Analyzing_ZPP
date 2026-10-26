"""
s51b -- CycleFreeFlux (Desouki et al 2015) using Gurobi.

Faithful port of s51_cyclefreeflux.py: same LP formulation, same
slack-tol handling for Sv=Sv0, same hard-abort on any LP failure.
Only difference: Gurobi instead of scipy HiGHS, with warm-start reuse
across consecutive samples (Dual Simplex method 1).

Expected speedup: ~3-5x versus scipy HiGHS, driven mostly by warm
starts between samples whose RHS/bounds differ only slightly.

Usage:
    python s51b_cyclefreeflux_gurobi.py --cell BL_24
    python s51b_cyclefreeflux_gurobi.py --cell all
"""
import argparse, time
import numpy as np, pandas as pd
import h5py, scipy.sparse as sp
from pathlib import Path
import gurobipy as gp
from gurobipy import GRB

OCR_MAP = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}
BOUND_EX = 10
BIO_LB_FRAC = 0.9
O2_EX = "MAR09048"
BIOMASS_RXN = "MAR00021"
SLACK_TOL = 1e-6

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
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        return dict(
            S=sp.csc_matrix((data, ir, jc), shape=(int(ir.max())+1, len(jc)-1)),
            lb=np.asarray(m["lb"]).ravel().astype(float),
            ub=np.asarray(m["ub"]).ravel().astype(float),
            rxns=deref_str(m["rxns"]))


def run_cell(cell, chains=(1,2,3,4)):
    repo = Path.cwd()
    print(f"\n=== CFF (Gurobi): {cell} ===", flush=True); t_all=time.time()
    cond, hpf = cell.split("_"); hpf = int(hpf)

    mat_fp = repo/args.models_dir/f"trans_rfastcormics_{cell}.mat"
    M = load_v73(mat_fp)
    S = M["S"]; lb = M["lb"].copy(); ub = M["ub"].copy()
    rxn_ids = M["rxns"]; n_rxn = len(rxn_ids); n_met = S.shape[0]
    idx = {r:i for i,r in enumerate(rxn_ids)}
    Snnz = np.diff(S.tocsc().indptr); is_ex = Snnz == 1

    # Same bounds pipeline as s51 (only used for the biomass_max probe)
    lb[is_ex] = -BOUND_EX; ub[is_ex] = BOUND_EX
    mb = pd.read_csv(repo/"data/phase2_exchange_bounds_long.csv")
    sel = (mb.condition == cond) & (mb.hpf == hpf)
    for rid, fc in zip(mb.loc[sel,"ex_rxn"], mb.loc[sel,"fc"]):
        if rid in idx:
            j = idx[rid]; fc = float(fc) if np.isfinite(fc) and fc>0 else 1.0
            lb[j] = -BOUND_EX*fc; ub[j] = BOUND_EX*fc
    lb[idx[O2_EX]] = -OCR_MAP[hpf]; ub[idx[O2_EX]] = 0

    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()

    # ---------- biomass_max probe ----------
    mo_bm = gp.Model("biomax", env=env)
    mo_bm.setParam("Threads", 1); mo_bm.setParam("OutputFlag", 0)
    v_bm = mo_bm.addMVar(n_rxn, lb=lb, ub=ub)
    mo_bm.addMConstr(S.tocsr(), v_bm, '=', np.zeros(n_met))
    mo_bm.setObjective(v_bm[idx[BIOMASS_RXN]], GRB.MAXIMIZE)
    mo_bm.optimize()
    if mo_bm.Status != GRB.OPTIMAL:
        raise RuntimeError(f"[{cell}] biomass_max failed, status={mo_bm.Status}")
    bio_max = mo_bm.ObjVal
    bio_lb = BIO_LB_FRAC * bio_max * (1 - 1e-4)
    lb[idx[BIOMASS_RXN]] = bio_lb
    print(f"  bio_max={bio_max:.3f} bio_lb={bio_lb:.3f}", flush=True)

    j_bio = idx[BIOMASS_RXN]
    ex_or_bio = is_ex.copy(); ex_or_bio[j_bio] = True
    internal = ~ex_or_bio
    n_int = int(internal.sum())
    int_idx = np.where(internal)[0]
    print(f"  {n_int} L1-min internal (biomass + exchanges pinned)", flush=True)

    # ---------- Build CFF Gurobi model once (reused across chains + samples) ----------
    mo = gp.Model("cff", env=env)
    mo.setParam("Threads", 1)
    mo.setParam("OutputFlag", 0)
    mo.setParam("Method", 1)         # dual simplex -- best warm-start reuse under RHS changes
    mo.setParam("FeasibilityTol", 1e-7)
    mo.setParam("OptimalityTol", 1e-7)

    v  = mo.addMVar(n_rxn, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="v")
    vp = mo.addMVar(n_int, lb=0.0, ub=GRB.INFINITY, name="vp")
    vn = mo.addMVar(n_int, lb=0.0, ub=GRB.INFINITY, name="vn")

    # Objective: minimize sum(vp + vn)
    mo.setObjective(vp.sum() + vn.sum(), GRB.MINIMIZE)

    # Sv constraints -- RHS updated per sample
    S_csr = S.tocsr()
    sv_le = mo.addMConstr(S_csr, v, '<', np.zeros(n_met))
    sv_ge = mo.addMConstr(S_csr, v, '>', np.zeros(n_met))

    # Identity: v[int_idx[k]] - vp[k] + vn[k] = 0
    I_int = sp.coo_matrix((np.ones(n_int), (np.arange(n_int), int_idx)),
                          shape=(n_int, n_rxn)).tocsr()
    mo.addConstr(I_int @ v - vp + vn == 0, name="l1_split")

    mo.update()

    # Dense S for Sv0 matvec per sample
    S_dense_full = S.toarray()

    outdir = repo/"results"/args.outdir_out
    indir  = repo/"results"/args.outdir_in
    outdir.mkdir(parents=True, exist_ok=True)
    stats_rows = []

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
        t_ch = time.time()

        for s in range(n_smp):
            v0 = V_in[:, s]

            # Variable bounds on v: sign-preserving for internals, pinned for ex/bio
            lb_v = np.zeros(n_rxn); ub_v = np.zeros(n_rxn)
            pos_mask = v0 >= 0
            # internals with v0 >= 0: [0, v0]
            m_ipos = internal & pos_mask
            lb_v[m_ipos] = 0.0
            ub_v[m_ipos] = v0[m_ipos]      # already >= 0
            # internals with v0 < 0: [v0, 0]
            m_ineg = internal & ~pos_mask
            lb_v[m_ineg] = v0[m_ineg]      # already <= 0
            ub_v[m_ineg] = 0.0
            # exchanges + biomass: pinned at v0
            lb_v[ex_or_bio] = v0[ex_or_bio]
            ub_v[ex_or_bio] = v0[ex_or_bio]

            v.LB = lb_v
            v.UB = ub_v

            Sv0 = S_dense_full @ v0
            sv_le.RHS = Sv0 + SLACK_TOL
            sv_ge.RHS = Sv0 - SLACK_TOL

            mo.optimize()
            if mo.Status == GRB.OPTIMAL:
                V_out[:, s] = v.X
            else:
                n_fail += 1
                if fail_first_sample < 0:
                    fail_first_sample = s
                fail_last_status = mo.Status
                V_out[:, s] = np.nan

            if (s+1) % 100 == 0:
                el = time.time()-t_ch
                print(f"  chain {ch}: {s+1}/{n_smp} ({(s+1)/el:.1f}/s, ETA {(n_smp-s-1)/((s+1)/el):.0f}s, fails={n_fail})", flush=True)

        if n_fail > 0:
            raise RuntimeError(
                f"[{cell} ch{ch}] CFF LP failed on {n_fail}/{n_smp} samples. "
                f"First failure at sample {fail_first_sample}. "
                f"Last Gurobi status={fail_last_status}. NO fallback used. Cell aborted.")

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
p.add_argument("--cell", default="BL_24")
p.add_argument("--outdir_in",  default="hopsy_v7_local")
p.add_argument("--outdir_out", default="hopsy_v7_loopless_gurobi")
p.add_argument("--models_dir", default="../reviewer_packet_v4_slim/v7_models")
p.add_argument("--chains", default="1,2,3,4",
               help="comma-separated chain ids; use to parallelize across jobs")
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
