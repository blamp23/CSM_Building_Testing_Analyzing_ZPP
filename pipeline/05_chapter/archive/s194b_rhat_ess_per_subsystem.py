"""
s194b -- per (corner, subsystem) rank-R̂ and bulk-ESS on enzyme-subset
representatives, mirroring §2.8 convention.

For every (pair, subsystem) in s180_f1_relay + s181_f2_late_relay PAIRS
and (as the chapter's own anchor) the Tier A lipid F2 pairs from s186,
compute R̂ + ESS over the subsystem's reps in all four corners:
  - AA, BB: hopsy_v7_loopless/flux_samples_{cell}_ch{N}.csv (frozen, CFF'd)
  - AB, BA: lipid_relay_samples/chains/{model}_mdl__{env}_env_ch{N}.csv
            (CFF'd post-hoc by s192)

Report per (pair, subsystem, corner):
  rhat_median, rhat_max, rhat_frac_above_1p10
  ess_bulk_median, ess_bulk_min
  n_reps

Gate (per user, §2.8 wording): a pair's decomposition is permitted if
median R̂ < 1.05 AND median ESS > 100 in ALL FOUR corners.
Flag (not drop): any corner where > 10 % of reps exceed R̂ 1.10.

Writes s194b_rhat_ess_per_subsystem.csv.
"""
import sys, time, re
from pathlib import Path
from multiprocessing import Pool
import numpy as np, pandas as pd
from scipy.stats import norm, rankdata

REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
POOL_STD_CHAINS = R / "hopsy_v7_loopless"
SWAP_CHAIN_DIR = R / "v7_final_addenda/lipid_relay_samples/chains"
SUBSETS = R / "subsets_v7"
STATS_FP = R / "v7_final_addenda/s194b_rhat_ess_per_subsystem.csv"

# Thinning per sampling run. Frozen pooled_v7 and the default swap batch were
# both sampled at thinning 300. Overrides below record §11.22 re-samples.
DEFAULT_THINNING_STD = 300
DEFAULT_THINNING_SWAP = 300
THINNING_OVERRIDES_SWAP = {
    # (model_cell, env_cell): thinning
    ("D_24", "BL_24"): 1000,   # BA for BL_24↔D_24 — flagged by s194b first pass
}

PAIRS = [
    # (pair_label, source_set, A, B, subsystem, convention)
    # F1 pairs from s180
    ("F1_24_48", "F1", "BL_24", "D_24",  "Nucleotide metabolism"),
    ("F1_24_48", "F1", "BL_24", "D_24",  "Pyrimidine metabolism"),
    ("F1_24_48", "F1", "D_24",  "LD_24", "Nucleotide metabolism"),
    ("F1_24_48", "F1", "D_24",  "LD_24", "Pyrimidine metabolism"),
    ("F1_24_48", "F1", "BL_48", "LD_48", "Pyrimidine metabolism"),
    ("F1_24_48", "F1", "D_48",  "LD_48", "Nucleotide metabolism"),
    ("F1_24_48", "F1", "D_48",  "LD_48", "Pyrimidine metabolism"),
    # F2 late from s181
    ("F2_late",  "F2", "D_96",  "D_120",  "Nucleotide metabolism"),
    ("F2_late",  "F2", "D_96",  "D_120",  "Pyrimidine metabolism"),
    ("F2_late",  "F2", "D_96",  "D_120",  "Phenylalanine, tyrosine and tryptophan biosynthesis"),
    ("F2_late",  "F2", "LD_96", "LD_120", "Arginine and proline metabolism"),
    ("F2_late",  "F2", "LD_96", "LD_120", "Glycolysis / Gluconeogenesis"),
    ("F2_late",  "F2", "LD_96", "LD_120", "Phenylalanine, tyrosine and tryptophan biosynthesis"),
    ("F2_late",  "F2", "LD_96", "LD_120", "Valine, leucine, and isoleucine metabolism"),
]


def _rank_rhat_ess(x):
    """x shape (n_chains, n_draws) → (rank-R̂, bulk ESS)."""
    n_chains, n_draws = x.shape
    if n_draws < 4 or np.std(x) < 1e-12:
        return 1.0, float(n_chains * n_draws)
    pooled = x.flatten()
    ranks = rankdata(pooled)
    z = norm.ppf((ranks - 3.0/8.0) / (len(ranks) + 1.0/4.0))
    y = z.reshape(n_chains, n_draws)
    half = n_draws // 2
    y_split = y[:, :2*half].reshape(2 * n_chains, half)
    m = y_split.shape[0]; N = y_split.shape[1]
    chain_means = y_split.mean(axis=1)
    chain_vars  = y_split.var(axis=1, ddof=1)
    W = chain_vars.mean()
    B = N * chain_means.var(ddof=1) if m > 1 else 0.0
    var_hat = ((N - 1) / N) * W + B / N
    rhat = float(np.sqrt(var_hat / W)) if W > 0 else 1.0
    def _ac(ch):
        n = len(ch); mean = ch.mean(); var = ch.var()
        if var < 1e-20: return np.zeros(n)
        c = np.correlate(ch - mean, ch - mean, mode="full") / (var * n)
        return c[n-1:]
    rho = np.mean([_ac(y_split[i]) for i in range(m)], axis=0)
    tau = -1.0; k = 0
    while k + 1 < len(rho) and (rho[k] + rho[k+1]) > 0:
        tau += 2 * (rho[k] + rho[k+1]); k += 2
    tau = max(tau, 1.0)
    return rhat, float(m * N / tau)


def _reps(cell, subsystem):
    sub = pd.read_csv(SUBSETS / f"{cell}_subset_membership.csv")
    return sub.loc[sub.subsystem == subsystem, "representative"].unique().tolist()


def _std_chains(cell):
    """Return list of 4 pd.DataFrames (one per chain) for standard CFF'd cell."""
    frames = []
    for ch in (1, 2, 3, 4):
        fp = POOL_STD_CHAINS / f"flux_samples_{cell}_ch{ch}.csv"
        if not fp.exists():
            return None
        frames.append(pd.read_csv(fp, index_col=0))
    return frames


def _swap_chains(model_cell, env_cell):
    frames = []
    for ch in (1, 2, 3, 4):
        fp = SWAP_CHAIN_DIR / f"{model_cell}_mdl__{env_cell}_env_ch{ch}.csv"
        if not fp.exists():
            return None
        frames.append(pd.read_csv(fp, index_col=0))
    return frames


def _corner_stats(frames, reps, label):
    """frames = list of 4 chains (DFs indexed by rxn). Returns summary."""
    rep_in = [r for r in reps if r in frames[0].index]
    if not rep_in:
        return dict(corner=label, n_reps=0,
                    rhat_median=np.nan, rhat_max=np.nan,
                    rhat_frac_above_1p10=np.nan,
                    ess_bulk_median=np.nan, ess_bulk_min=np.nan)
    n_draws = min(f.shape[1] for f in frames)
    data = np.stack([f.loc[rep_in].values[:, :n_draws] for f in frames], axis=1)
    # data shape: (n_reps, 4 chains, n_draws)
    rhats = np.zeros(len(rep_in)); esss = np.zeros(len(rep_in))
    for i in range(len(rep_in)):
        rh, es = _rank_rhat_ess(data[i])
        rhats[i] = rh; esss[i] = es
    return dict(
        corner=label, n_reps=len(rep_in),
        rhat_median=float(np.median(rhats)),
        rhat_max=float(np.max(rhats)),
        rhat_frac_above_1p10=float((rhats > 1.10).mean()),
        ess_bulk_median=float(np.median(esss)),
        ess_bulk_min=float(np.min(esss)),
    )


def main():
    rows = []
    for pair_label, convention, A, B, subsystem in PAIRS:
        # Four corners × subsystem
        AA = _std_chains(A)
        BB = _std_chains(B)
        AB = _swap_chains(A, B)
        BA = _swap_chains(B, A)
        # (corner, frames, cell-for-reps, model_cell, env_cell, is_swap)
        corner_specs = [
            ("AA", AA, A, A, A, False),
            ("AB", AB, A, A, B, True),
            ("BA", BA, B, B, A, True),
            ("BB", BB, B, B, B, False),
        ]
        pair_rows = []
        for corner_label, frames, reps_cell, mdl_c, env_c, is_swap in corner_specs:
            if is_swap:
                thinning = THINNING_OVERRIDES_SWAP.get(
                    (mdl_c, env_c), DEFAULT_THINNING_SWAP)
            else:
                thinning = DEFAULT_THINNING_STD
            if frames is None:
                pair_rows.append(dict(corner=corner_label, thinning=thinning,
                                     model_cell=mdl_c, env_cell=env_c,
                                     n_reps=0,
                                     rhat_median=np.nan, rhat_max=np.nan,
                                     rhat_frac_above_1p10=np.nan,
                                     ess_bulk_median=np.nan, ess_bulk_min=np.nan,
                                     missing=True))
                continue
            reps = _reps(reps_cell, subsystem)
            s = _corner_stats(frames, reps, corner_label)
            s["thinning"] = thinning
            s["model_cell"] = mdl_c
            s["env_cell"] = env_c
            pair_rows.append(s)
        # One row per corner for this pair-subsystem
        for s in pair_rows:
            rows.append(dict(pair_label=pair_label, convention=convention,
                             A=A, B=B, subsystem=subsystem, **s))

    df = pd.DataFrame(rows)
    df.to_csv(STATS_FP, index=False)
    print(f"[write] {STATS_FP}  ({len(df)} rows)")
    print()

    # Gate: pair passes if median R̂ < 1.05 and median ESS > 100 in all four
    pair_gate = []
    for (A, B, subsystem), grp in df.groupby(["A", "B", "subsystem"]):
        all_pass = True
        notes = []
        for _, r in grp.iterrows():
            if pd.isna(r.rhat_median) or pd.isna(r.ess_bulk_median):
                all_pass = False
                notes.append(f"{r.corner}: missing")
                continue
            if not (r.rhat_median < 1.05 and r.ess_bulk_median > 100):
                all_pass = False
                notes.append(f"{r.corner}: R̂med={r.rhat_median:.3f} ESSmed={r.ess_bulk_median:.0f}")
            if r.rhat_frac_above_1p10 > 0.10:
                notes.append(f"{r.corner}: {r.rhat_frac_above_1p10*100:.0f}% reps above R̂ 1.10")
        pair_gate.append(dict(A=A, B=B, subsystem=subsystem,
                              gate="PASS" if all_pass else "FAIL",
                              notes="; ".join(notes)))
    pg = pd.DataFrame(pair_gate)
    pg.to_csv(R / "v7_final_addenda/s194b_pair_gate.csv", index=False)
    print(pg.to_string(index=False))
    print()
    print(f"PASS: {(pg.gate=='PASS').sum()} / {len(pg)}")
    print()
    print("Reference (§2.8 frozen run): median R̂ 1.013–1.023, median ESS 155–249.")


if __name__ == "__main__":
    main()
