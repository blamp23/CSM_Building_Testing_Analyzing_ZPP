"""
§11.14 method-of-record variant: MEAN per-rep normaliser + ZEROS DROPPED.

Differs from two_tier_subsystem_test_fast.py (§2 pre-reg as written) only in
how zeros are handled: §2 adds EPS=1e-30 to |v| so zero samples enter the
test at log2(ε) = -99.66. §11.14 corrects that to the v5 convention: samples
with |v| < 1e-10 are dropped and the zero fraction is handled by the
separate per-row Fisher-exact in s80.

Differs from two_tier_subsystem_test_fast_median.py (§11.12 sensitivity)
only in the normaliser: mu_r = mean of pooled non-zero values (not median).

Running all three isolates the two axes:
  v7 vs v7_mean_dropzero  → zero-handling effect (same normaliser)
  v7_median vs v7_mean_dropzero → normaliser effect (same zero handling)
"""
import argparse
import numpy as np, pandas as pd
from scipy.stats import fisher_exact

ap = argparse.ArgumentParser()
ap.add_argument("samples_a"); ap.add_argument("samples_b")
ap.add_argument("membership"); ap.add_argument("out")
ap.add_argument("--carry", type=float, default=1e-4)
ap.add_argument("--loop", type=float, default=100.0)
ap.add_argument("--ndraw", type=int, default=2000)
ap.add_argument("--min_carry", type=int, default=3)
ap.add_argument("--chunk", type=int, default=200)
a = ap.parse_args()

ZERO_TOL = 1e-10

XA = pd.read_csv(a.samples_a, index_col=0); XB = pd.read_csv(a.samples_b, index_col=0)
mem = pd.read_csv(a.membership)
common = XA.index.intersection(XB.index)
AA, AB = np.abs(XA.loc[common].values), np.abs(XB.loc[common].values)
medA, medB = np.median(AA, 1), np.median(AB, 1)
carryA, carryB = medA > a.carry, medB > a.carry
carry_both = carryA & carryB
loop = (AA.max(1) > a.loop) | (AB.max(1) > a.loop)

rep_info = mem.drop_duplicates("representative").set_index("representative")
gpr = rep_info.has_gpr.reindex(common).fillna(0).values == 1
subsys = rep_info.subsystem.reindex(common).fillna("UNMAPPED")
exch = subsys.str.contains("Exchange|Artificial|Transport", case=False).values

# §11.14 method of record: mu_r = mean of pooled NON-ZERO values; zeros dropped
# from the W1 input (NaN in log arrays, nanquantile below).
AB_concat = np.concatenate([AA, AB], axis=1)
mask_nz = AB_concat >= ZERO_TOL
with np.errstate(all="ignore"):
    masked = np.where(mask_nz, AB_concat, np.nan)
    mu = np.nanmean(masked, axis=1)
mu_safe = np.where(np.isfinite(mu) & (mu > 0), mu, 1.0)
mu_col = mu_safe[:, None]
log_AA = np.where(AA >= ZERO_TOL,
                   np.log2(np.maximum(AA, ZERO_TOL) / mu_col), np.nan)
log_AB = np.where(AB >= ZERO_TOL,
                   np.log2(np.maximum(AB, ZERO_TOL) / mu_col), np.nan)
q = np.linspace(0, 1, 201)

def w1_serial(idx):
    la = log_AA[idx].ravel(); lb = log_AB[idx].ravel()
    la = la[~np.isnan(la)]; lb = lb[~np.isnan(lb)]
    if len(la) < 10 or len(lb) < 10: return np.nan
    return np.abs(np.quantile(la, q) - np.quantile(lb, q)).mean()

def w1_batch(draws, chunk_size):
    n = draws.shape[0]
    out = np.empty(n)
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        d = draws[start:end]
        la = log_AA[d].reshape(end - start, -1)
        lb = log_AB[d].reshape(end - start, -1)
        with np.errstate(all="ignore"):
            qa = np.nanquantile(la, q, axis=1)
            qb = np.nanquantile(lb, q, axis=1)
        out[start:end] = np.abs(qa - qb).mean(axis=0)
    return out

bg = np.where(gpr & ~exch & ~loop & carry_both)[0]
rng = np.random.default_rng(0)
rows = []
for s in sorted(subsys.unique()):
    idx = np.where(subsys.values == s)[0]
    if len(idx) == 0 or s == "UNMAPPED":
        continue
    nA, nB = int(carryA[idx].sum()), int(carryB[idx].sum())
    n_both = int(carry_both[idx].sum())
    n_loop = int((loop & carry_both)[idx].sum())
    p_presence = fisher_exact([[nA, len(idx) - nA], [nB, len(idx) - nB]])[1]
    row = dict(subsystem=s, N_rep=len(idx), carry_A=nA, carry_B=nB,
               carry_both=n_both, loop_flagged=n_loop, p_presence=p_presence)
    idx2 = idx[carry_both[idx] & ~loop[idx]]
    if len(idx2) >= a.min_carry and len(bg) > len(idx2):
        obs = w1_serial(idx2)
        # SAME draws as reviewer's script: draw sequentially with the same rng
        draws = np.stack([rng.choice(bg, len(idx2), replace=False)
                          for _ in range(a.ndraw)])
        null = w1_batch(draws, a.chunk)
        null = null[np.isfinite(null)]
        if len(null) < 10 or not np.isfinite(obs):
            row.update(N_tier2=len(idx2), W1_obs=obs,
                       null_median=np.nan, z=np.nan, p_emp=np.nan)
        else:
            row.update(N_tier2=len(idx2), W1_obs=obs,
                       null_median=float(np.median(null)),
                       z=(obs - null.mean()) / null.std(),
                       p_emp=(1 + (null >= obs).sum()) / (a.ndraw + 1))
    rows.append(row)

out = pd.DataFrame(rows).sort_values("z", ascending=False, na_position="last")
out.to_csv(a.out, index=False)
print(f"reps in both cells: {len(common)}; carrying in A {carryA.sum()}, "
      f"in B {carryB.sum()}, in both {carry_both.sum()}; "
      f"loop-flagged {loop.sum()}; Tier-2 background size {len(bg)}")
print(out.head(15).to_string())
