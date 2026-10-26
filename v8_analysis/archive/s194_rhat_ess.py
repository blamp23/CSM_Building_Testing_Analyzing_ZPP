"""
s194 -- rank-R̂ + bulk-ESS per corner on CFF-corrected swap-corner chains.

Convention per chapter §2.8 (rank-normalised split-R̂, Vehtari et al.
2021). Each corner has 4 chains × 500 samples × ~16000 reactions. For
each reaction we have 2000 samples across 4 chains (split-R̂ splits each
chain into first/second half → 8 "chains" of 250).

Report per corner:
  - rhat_max_rxn, rhat_median_rxn      — across all reactions
  - ess_bulk_min_rxn, ess_bulk_median  — bulk ESS floor + typical

Gate at chapter standard (§2.8 convention): R̂ < 1.05 AND bulk-ESS > 100.
Corners failing the gate are flagged and excluded from decomposition.

Writes grace_v4/results/v7_final_addenda/s194_rhat_ess.csv.
"""
import sys, time, re
from pathlib import Path
from multiprocessing import Pool
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[2]
CHAIN_DIR = REPO / "grace_v4/results/v7_final_addenda/lipid_relay_samples/chains"
STATS_FP = REPO / "grace_v4/results/v7_final_addenda/s194_rhat_ess.csv"

_FN_RE = re.compile(r"^(?P<mdl>[A-Z]+_\d+)_mdl__(?P<env>[A-Z]+_\d+)_env_ch(?P<ch>\d+)\.csv$")


def _rank_rhat_ess(x):
    """
    x: array shape (n_chains, n_draws).
    Returns (rank-normalised split-R̂, bulk ESS).

    Minimal implementation (Vehtari+ 2021). We rank the pooled draws
    across all chains, apply the inverse normal CDF, then compute the
    standard split-R̂ and ESS on the rank-normalised series.
    """
    from scipy.stats import norm, rankdata
    n_chains, n_draws = x.shape
    if n_draws < 4 or np.std(x) < 1e-12:
        return 1.0, float(n_chains * n_draws)
    # rank-normalise
    pooled = x.flatten()
    ranks = rankdata(pooled)
    z = norm.ppf((ranks - 3.0/8.0) / (len(ranks) + 1.0/4.0))
    y = z.reshape(n_chains, n_draws)
    # split-R̂: split each chain in half
    half = n_draws // 2
    y_split = y[:, :2*half].reshape(2 * n_chains, half)
    m = y_split.shape[0]; N = y_split.shape[1]
    chain_means = y_split.mean(axis=1)
    chain_vars  = y_split.var(axis=1, ddof=1)
    W = chain_vars.mean()
    B = N * chain_means.var(ddof=1) if m > 1 else 0.0
    var_hat = ((N - 1) / N) * W + B / N
    rhat = float(np.sqrt(var_hat / W)) if W > 0 else 1.0

    # Bulk ESS (autocorrelation-based, standard formula)
    # Compute per-chain autocorrelations, average across chains, Geyer initial positive
    def _autocorr_per_chain(ch):
        n = len(ch); mean = ch.mean(); var = ch.var()
        if var < 1e-20: return np.zeros(n)
        c = np.correlate(ch - mean, ch - mean, mode="full") / (var * n)
        return c[n-1:]   # lags 0..n-1
    rho = np.mean([_autocorr_per_chain(y_split[i]) for i in range(m)], axis=0)
    # Geyer initial positive truncation: sum pairs
    tau = -1.0
    k = 0
    while k + 1 < len(rho) and (rho[k] + rho[k+1]) > 0:
        tau += 2 * (rho[k] + rho[k+1])
        k += 2
    tau = max(tau, 1.0)
    ess = float(m * N / tau)
    return rhat, ess


def _one_corner_group(group):
    """Compute rank-R̂ + ESS for one (model, env) corner.
    group: list of chain CSV paths for that corner."""
    group = sorted(group)
    frames = [pd.read_csv(fp).set_index("rxn") for fp in group]
    # align on union of reaction ids
    all_rxns = set(frames[0].index)
    for f in frames[1:]: all_rxns &= set(f.index)
    all_rxns = sorted(all_rxns)
    # Each frame is (n_rxn, n_samples_per_chain)
    n_chains = len(frames)
    n_draws = min(f.shape[1] for f in frames)
    # stack to (n_rxn, n_chains, n_draws)
    data = np.stack([f.loc[all_rxns].values[:, :n_draws] for f in frames], axis=1)
    # Compute R̂ + ESS per reaction
    rhats = np.zeros(len(all_rxns)); esss = np.zeros(len(all_rxns))
    for i in range(len(all_rxns)):
        rh, es = _rank_rhat_ess(data[i])
        rhats[i] = rh; esss[i] = es
    tag_match = _FN_RE.match(group[0].name)
    return dict(
        model_cell=tag_match.group("mdl"),
        env_cell=tag_match.group("env"),
        n_chains=n_chains,
        n_draws=n_draws,
        n_rxns=len(all_rxns),
        rhat_max=float(np.max(rhats)),
        rhat_median=float(np.median(rhats)),
        rhat_frac_above_1p05=float((rhats > 1.05).mean()),
        ess_bulk_min=float(np.min(esss)),
        ess_bulk_median=float(np.median(esss)),
        ess_bulk_frac_below_100=float((esss < 100).mean()),
    )


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=3)
    args = ap.parse_args()
    chain_fps = sorted(CHAIN_DIR.glob("*_mdl__*_env_ch*.csv"))
    # group by (model, env)
    groups = {}
    for fp in chain_fps:
        m = _FN_RE.match(fp.name)
        if not m: continue
        key = (m.group("mdl"), m.group("env"))
        groups.setdefault(key, []).append(fp)
    # filter to groups with ≥2 chains
    valid = [g for g in groups.values() if len(g) >= 2]
    print(f"[s194] {len(valid)} corners with ≥2 chains, {args.jobs} workers")
    t0 = time.time()
    rows = []
    with Pool(args.jobs) as pool:
        for r in pool.imap_unordered(_one_corner_group, valid):
            rows.append(r)
            gate = "PASS" if r["rhat_max"] < 1.05 and r["ess_bulk_min"] > 100 else "FAIL"
            print(f"  {r['model_cell']} mdl + {r['env_cell']} env  "
                  f"R̂max={r['rhat_max']:.3f}  ESSmin={r['ess_bulk_min']:.0f}  {gate}",
                  flush=True)
    df = pd.DataFrame(rows)
    df["gate"] = df.apply(lambda r: "PASS" if r.rhat_max < 1.05 and r.ess_bulk_min > 100 else "FAIL", axis=1)
    df.to_csv(STATS_FP, index=False)
    print(f"\n[write] {STATS_FP}  ({len(df)} corners)  total={(time.time()-t0)/60:.1f}m")
    print(f"PASS: {(df.gate=='PASS').sum()}   FAIL: {(df.gate=='FAIL').sum()}")


if __name__ == "__main__":
    main()
