"""
s194c -- resolve the BL_24↔D_24 nucleotide/pyrimidine R̂ flags.

Step 1  Parity check across the four corners
    From s194b_rhat_ess_per_subsystem.csv, pull AA, AB, BA, BB rows for
    BL_24↔D_24 Nucleotide and Pyrimidine. If AA/BB (frozen loopless)
    show a comparable rhat_frac_above_1p10 the flag is a property of
    24 hpf nucleotide metabolism; if only BA is high the swap corner
    is the problem.

Step 2  If BA is the outlier
    Load the BA chains (D_24 mdl + BL_24 env, CFF-corrected), compute
    per-representative R̂ and mean |v| within each subsystem, correlate,
    and report median R̂ of the top-|v|-quartile representatives.

Does not resample. Prints a verdict table and leaves the resample
decision to the caller.
"""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm, rankdata, spearmanr

REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
STATS = pd.read_csv(R / "v7_final_addenda/s194b_rhat_ess_per_subsystem.csv")
SWAP = R / "v7_final_addenda/lipid_relay_samples/chains"
POOL = R / "hopsy_v7_loopless"
SUBSETS = R / "subsets_v7"

PAIRS = [
    ("BL_24", "D_24", "Nucleotide metabolism"),
    ("BL_24", "D_24", "Pyrimidine metabolism"),
]


def _rank_rhat(x):
    n_chains, n_draws = x.shape
    if n_draws < 4 or np.std(x) < 1e-12:
        return 1.0
    pooled = x.flatten()
    ranks = rankdata(pooled)
    z = norm.ppf((ranks - 3.0/8.0) / (len(ranks) + 1.0/4.0))
    y = z.reshape(n_chains, n_draws)
    half = n_draws // 2
    y_split = y[:, :2*half].reshape(2 * n_chains, half)
    m, N = y_split.shape
    chain_means = y_split.mean(axis=1)
    chain_vars  = y_split.var(axis=1, ddof=1)
    W = chain_vars.mean()
    B = N * chain_means.var(ddof=1) if m > 1 else 0.0
    var_hat = ((N - 1) / N) * W + B / N
    return float(np.sqrt(var_hat / W)) if W > 0 else 1.0


def _load_chains_swap(model_cell, env_cell):
    frames = []
    for ch in (1, 2, 3, 4):
        fp = SWAP / f"{model_cell}_mdl__{env_cell}_env_ch{ch}.csv"
        frames.append(pd.read_csv(fp, index_col=0))
    return frames


def _load_chains_pool(cell):
    frames = []
    for ch in (1, 2, 3, 4):
        fp = POOL / f"flux_samples_{cell}_ch{ch}.csv"
        frames.append(pd.read_csv(fp, index_col=0))
    return frames


def _reps(cell, subsystem):
    sub = pd.read_csv(SUBSETS / f"{cell}_subset_membership.csv")
    return sub.loc[sub.subsystem == subsystem, "representative"].unique().tolist()


def main():
    print("Step 1 — four-corner parity check (rhat_frac_above_1p10, median R̂)")
    print("-" * 80)
    for A, B, subs in PAIRS:
        g = STATS[(STATS.A == A) & (STATS.B == B) & (STATS.subsystem == subs)]
        g = g.set_index("corner")
        print(f"\n{A}↔{B}  —  {subs}  (n_reps={int(g.n_reps.iloc[0])})")
        print(f"  corner    rhat_median   frac_above_1p10   ess_median   ess_min")
        for c in ["AA", "AB", "BA", "BB"]:
            if c in g.index:
                r = g.loc[c]
                print(f"  {c}        {r.rhat_median:.3f}         "
                      f"{r.rhat_frac_above_1p10*100:5.1f}%            "
                      f"{r.ess_bulk_median:6.0f}      {r.ess_bulk_min:5.0f}")

    print()
    print("Step 2 — Per-rep R̂ vs mean |v| in BA (D_24 mdl + BL_24 env)")
    print("-" * 80)
    ba_chains = _load_chains_swap("D_24", "BL_24")
    for A, B, subs in PAIRS:
        # Reps come from the model cell (D_24 for BA)
        reps = _reps("D_24", subs)
        rep_in = [r for r in reps if r in ba_chains[0].index]
        n_draws = min(f.shape[1] for f in ba_chains)
        data = np.stack([f.loc[rep_in].values[:, :n_draws] for f in ba_chains], axis=1)
        rhats = np.zeros(len(rep_in))
        for i in range(len(rep_in)):
            rhats[i] = _rank_rhat(data[i])
        mean_abs = np.abs(data).reshape(len(rep_in), -1).mean(axis=1)
        rho, p = spearmanr(mean_abs, rhats)
        q75 = np.quantile(mean_abs, 0.75)
        hi = mean_abs >= q75
        med_rhat_hi = float(np.median(rhats[hi]))
        med_rhat_lo = float(np.median(rhats[~hi]))
        n_hi_above = int((rhats[hi] > 1.10).sum())
        print(f"\n{subs}  —  BA  (n_reps={len(rep_in)})")
        print(f"  Spearman ρ(R̂, mean|v|) = {rho:+.3f}  (p={p:.2g})")
        print(f"  median R̂ of top-|v|-quartile ({hi.sum()} reps): {med_rhat_hi:.3f}  "
              f"({n_hi_above} above 1.10)")
        print(f"  median R̂ of bottom-three-quartiles ({(~hi).sum()} reps): {med_rhat_lo:.3f}")
        # Show individual high-flux offenders
        order = np.argsort(mean_abs)[::-1]
        print(f"  top-5 by mean |v|:")
        print(f"    {'rep':<12} {'mean|v|':>10}  {'R̂':>6}")
        for k in order[:5]:
            print(f"    {rep_in[k]:<12} {mean_abs[k]:10.3f}  {rhats[k]:6.3f}")

    print()
    print("Verdict logic:")
    print("  - If AA, BB, AB all show ~1.7% and BA shows 28%: swap corner is the problem.")
    print("  - If top-|v|-quartile median R̂ > 1.05 while bottom shows §2.8-like 1.01: "
          "high-flux reps drive the statistic and resample is warranted.")
    print("  - Otherwise the flag is informational.")


if __name__ == "__main__":
    main()
