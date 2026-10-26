"""
s119 -- gamma-tail p + BH over pre-declared family, with the CORRECTED
Wasserstein statistic: drop exact zeros before log2 (they're solver tie-break
artifacts post-CFF, not biology). Also add a companion Fisher-on-zeros
column so zero-fraction differences are reported as a separate, weaker
Tier 1-style track.

Per reviewer 2026-09-28.

The Wasserstein statistic:
  For each subsystem in the family, restrict to carrying reps that also are
  non-loop-flagged (unchanged). For W1:
    for each rep, drop the sample-columns where AA[rep, s] == 0 OR AB[rep, s] == 0
    compute log2(|v| / mu) on the surviving values, quantile at 201 points,
    W1 = mean(|q_A - q_B|).
  Null is drawn from the same background as before (GPR-linked, non-exchange,
  non-loop, carry_both), same 2000-draw procedure, but the w1 within each
  draw also drops zeros before log2.

Fisher on zeros:
  For each subsystem, count total zeros across (reps in subsystem × samples)
  per condition. Test 2x2 (zero, nonzero) x (A, B) via scipy.stats.fisher_exact.
  Reported as p_fisher_zeros; carry as informational, not primary.
"""
import argparse, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.stats import gamma, fisher_exact
from statsmodels.stats.multitest import multipletests

MIN_CARRY_BOTH = 5
LOOP_FRAC_MAX  = 0.30
GAMMA_N_CUTOFF = 15
EXCLUDE_PATTERNS = ["Exchange", "Artificial", "Transport"]
NDRAW = 2000
Q_THRESH = 0.10

CONTRASTS = [
    ("D_24","LD_24","D24_vs_LD24","LD_24"),
    ("D_48","LD_48","D48_vs_LD48","LD_48"),
    ("D_72","LD_72","D72_vs_LD72","LD_72"),
    ("D_96","LD_96","D96_vs_LD96","LD_96"),
    ("D_120","LD_120","D120_vs_LD120","LD_120"),
    ("BL_24","LD_24","BL24_vs_LD24","LD_24"),
    ("BL_48","LD_48","BL48_vs_LD48","LD_48"),
    ("BL_72","LD_72","BL72_vs_LD72","LD_72"),
    ("BL_96","LD_96","BL96_vs_LD96","LD_96"),
    ("BL_120","LD_120","BL120_vs_LD120","LD_120"),
]


def w1_nonzero_pair(la_pool, lb_pool, q):
    """W1 between two 1D log-magnitude pools (already zero-dropped)."""
    if la_pool.size < 10 or lb_pool.size < 10:
        return np.nan
    return float(np.abs(np.quantile(la_pool, q) - np.quantile(lb_pool, q)).mean())


def build_log_pool(AA_slice, AB_slice, mu_slice):
    """Drop entries where AA==0 or AB==0 (sample-wise, per rep) then log2/mu.
    Both output pools are aligned index-wise post-drop; but for Wasserstein
    on marginals we treat them as independent samples of their own dist."""
    # For each rep row, take non-zero mask separately per side; pool ravel.
    mask_A = AA_slice > 0
    mask_B = AB_slice > 0
    la = np.log2(AA_slice[mask_A] / np.broadcast_to(mu_slice[:, None], AA_slice.shape)[mask_A])
    lb = np.log2(AB_slice[mask_B] / np.broadcast_to(mu_slice[:, None], AB_slice.shape)[mask_B])
    return la, lb


def run_contrast(A_csv, B_csv, M_csv, ndraw=NDRAW, carry=1e-4, loop=100.0, seed=0):
    XA = pd.read_csv(A_csv, index_col=0); XB = pd.read_csv(B_csv, index_col=0)
    mem = pd.read_csv(M_csv)
    common = XA.index.intersection(XB.index)
    AA = np.abs(XA.loc[common].values); AB = np.abs(XB.loc[common].values)
    medA = np.median(AA, 1); medB = np.median(AB, 1)
    carryA = medA > carry; carryB = medB > carry
    carry_both = carryA & carryB
    loop_mask = (AA.max(1) > loop) | (AB.max(1) > loop)

    rep_info = mem.drop_duplicates("representative").set_index("representative")
    gpr = rep_info.has_gpr.reindex(common).fillna(0).values == 1
    subsys = rep_info.subsystem.reindex(common).fillna("UNMAPPED")
    exch = subsys.str.contains("Exchange|Artificial|Transport", case=False).values

    mu = (AA.mean(1) + AB.mean(1)) / 2
    mu_safe = np.where(mu > 0, mu, 1.0)
    q = np.linspace(0, 1, 201)

    bg = np.where(gpr & ~exch & ~loop_mask & carry_both)[0]
    rng = np.random.default_rng(seed)
    rows = []

    for s in sorted(subsys.unique()):
        idx = np.where(subsys.values == s)[0]
        if len(idx) == 0 or s == "UNMAPPED":
            continue
        n_both = int(carry_both[idx].sum())
        n_loop = int((loop_mask & carry_both)[idx].sum())
        idx2 = idx[carry_both[idx] & ~loop_mask[idx]]
        row = dict(subsystem=s, N_rep=len(idx), N_tier2=len(idx2),
                   carry_both=n_both, loop_flagged=n_loop)
        # W1 on non-zero values only
        if len(idx2) >= 3 and len(bg) > len(idx2):
            la_obs, lb_obs = build_log_pool(AA[idx2], AB[idx2], mu_safe[idx2])
            obs = w1_nonzero_pair(la_obs, lb_obs, q)
            # zero fractions (informational)
            zA = int((AA[idx2] == 0).sum()); zB = int((AB[idx2] == 0).sum())
            totA = AA[idx2].size; totB = AB[idx2].size
            frac_A = zA / totA if totA else 0.0
            frac_B = zB / totB if totB else 0.0
            # Fisher on zero vs non-zero
            _, p_fisher_zeros = fisher_exact([[zA, totA - zA], [zB, totB - zB]])
            row.update(W1_obs=obs, zero_frac_A=frac_A, zero_frac_B=frac_B,
                       delta_zero_frac=frac_A - frac_B,
                       p_fisher_zeros=p_fisher_zeros)
            # Null
            if obs is not np.nan:
                null = np.empty(ndraw)
                for k in range(ndraw):
                    d = rng.choice(bg, len(idx2), replace=False)
                    la, lb = build_log_pool(AA[d], AB[d], mu_safe[d])
                    null[k] = w1_nonzero_pair(la, lb, q)
                null = null[np.isfinite(null)]
                if len(null) >= 30:
                    row.update(null_mean=float(np.mean(null)),
                               null_std=float(np.std(null)),
                               null_median=float(np.median(null)),
                               z=(obs - null.mean()) / null.std() if null.std() > 0 else np.nan,
                               p_emp=(1 + (null >= obs).sum()) / (len(null) + 1))
                    null_pos = null[null > 0]
                    if len(null_pos) >= 30:
                        try:
                            a_hat, loc_hat, sc_hat = gamma.fit(null_pos, floc=0)
                            row["gamma_shape"] = a_hat
                            row["gamma_scale"] = sc_hat
                            row["p_gamma"] = float(gamma.sf(obs, a_hat, loc=loc_hat, scale=sc_hat))
                        except Exception:
                            row["p_gamma"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def apply_family_and_bh(df, tag):
    df = df.copy()
    df["contrast"] = tag
    df["loop_frac"] = df["loop_flagged"] / df["N_tier2"].replace(0, np.nan)
    is_exch = df["subsystem"].astype(str).str.contains(
        "|".join(EXCLUDE_PATTERNS), case=False, na=False)
    df["in_family"] = (
        (df["carry_both"] >= MIN_CARRY_BOTH) &
        (df["loop_frac"] < LOOP_FRAC_MAX) &
        (~is_exch) &
        (df["N_tier2"] >= 3)
    )
    df["p_primary"] = np.where(
        df["N_tier2"].fillna(0) < GAMMA_N_CUTOFF,
        df.get("p_gamma", np.nan), df.get("p_emp", np.nan)
    )
    fam = df[df["in_family"] & df["p_primary"].notna()].copy()
    if len(fam):
        reject, qvals, _, _ = multipletests(fam["p_primary"].values, method="fdr_bh")
        fam["q_bh"] = qvals
        df = df.merge(fam[["subsystem", "q_bh"]], on="subsystem", how="left")
        df["family_N"] = len(fam)
    else:
        df["q_bh"] = np.nan; df["family_N"] = 0
    df["is_finding"] = df["in_family"] & (df["q_bh"] <= Q_THRESH)
    return df


def _process(work):
    tag, A, B, M, ndraw = work
    t0 = time.time()
    df = run_contrast(A, B, M, ndraw=ndraw)
    df = apply_family_and_bh(df, tag)
    n_fam = int(df["in_family"].sum())
    n_find = int(df["is_finding"].sum())
    dt = time.time() - t0
    print(f"[{tag}] done in {dt:.0f}s  family={n_fam}  findings(q<={Q_THRESH})={n_find}", flush=True)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps_dir", type=Path, required=True)
    ap.add_argument("--subsets_dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--contrasts", default="all")
    ap.add_argument("--ndraw", type=int, default=NDRAW)
    ap.add_argument("--workers", type=int, default=10)
    args = ap.parse_args()

    if args.contrasts != "all":
        selected = set(args.contrasts.split(","))
        chosen = [c for c in CONTRASTS if c[2] in selected]
    else:
        chosen = CONTRASTS

    work = []
    for A, B, tag, memb_cell in chosen:
        A_csv = args.reps_dir / f"reps_pooled_{A}.csv"
        B_csv = args.reps_dir / f"reps_pooled_{B}.csv"
        M_csv = args.subsets_dir / f"{memb_cell}_subset_membership.csv"
        if not (A_csv.exists() and B_csv.exists() and M_csv.exists()):
            print(f"  [{tag}] missing inputs; skip"); continue
        work.append((tag, A_csv, B_csv, M_csv, args.ndraw))

    print(f"submitting {len(work)} contrasts to {args.workers} workers ...", flush=True)
    t_all = time.time()
    all_dfs = []
    with Pool(processes=args.workers) as pool:
        for df in pool.imap_unordered(_process, work):
            all_dfs.append(df)
    print(f"all done in {time.time()-t_all:.0f}s", flush=True)

    out = pd.concat(all_dfs, ignore_index=True)
    cols_out = ["contrast", "subsystem", "N_rep", "N_tier2", "carry_both",
                "loop_flagged", "loop_frac", "in_family", "family_N",
                "zero_frac_A", "zero_frac_B", "delta_zero_frac",
                "p_fisher_zeros",
                "W1_obs", "null_median", "null_mean", "null_std", "z", "p_emp",
                "gamma_shape", "gamma_scale", "p_gamma",
                "p_primary", "q_bh", "is_finding"]
    for c in cols_out:
        if c not in out.columns: out[c] = np.nan
    out[cols_out].to_csv(args.out, index=False)
    print(f"\nwrote {args.out}")

    print(f"\n==== FINDINGS (in family, q_bh <= {Q_THRESH}, NON-ZERO W1 statistic) ====")
    f = out[out["is_finding"]].sort_values("q_bh")
    if len(f):
        show = f[["contrast","subsystem","N_tier2","loop_frac",
                  "zero_frac_A","zero_frac_B","W1_obs","z",
                  "p_primary","q_bh","p_fisher_zeros"]]
        print(show.to_string(index=False))
    else:
        print("(none)")

    # Also list subsystems with significant Fisher-on-zeros that DIDN'T make Tier 2
    print(f"\n==== ZERO-FRACTION HITS (in family, p_fisher_zeros <= 0.01, not in Tier 2 above) ====")
    zh = out[out["in_family"] & (out["p_fisher_zeros"] <= 0.01) & (~out["is_finding"].fillna(False))]
    if len(zh):
        show = zh[["contrast","subsystem","N_tier2","zero_frac_A","zero_frac_B",
                   "delta_zero_frac","p_fisher_zeros","z"]] \
              .sort_values("p_fisher_zeros")
        print(show.head(30).to_string(index=False))
    else:
        print("(none)")


if __name__ == "__main__":
    main()
