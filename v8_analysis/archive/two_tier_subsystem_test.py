"""
Two-tier subsystem contrast on representative-level samples (A vs B, e.g. D_72 vs LD_72).

Tier 1 (presence): which representatives carry flux in each cell; Fisher exact test on
                   counts carrying-in-A vs carrying-in-B within the subsystem.
Tier 2 (redistribution): W1 on log2(|v| / pooled mean) restricted to representatives that
                   carry flux in BOTH cells, against a null of random sets of the same size
                   drawn from the flux-carrying, GPR-linked, non-exchange, non-loop background.
                   No epsilon anywhere: zero-flux representatives never enter Tier 2.

Usage:
    python two_tier_subsystem_test.py samples_A.csv samples_B.csv membership_A.csv out.csv \
        [--carry 1e-4] [--loop 100] [--ndraw 2000]

samples_*.csv   : rows = representatives (index 'rxn'), columns = samples
membership_A.csv: columns rxn, representative, subsystem, has_gpr (the subsystem map is taken
                  from this file; representatives are matched by id across cells)

Reviewer-provided script (unchanged).
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ap = argparse.ArgumentParser()
ap.add_argument("samples_a"); ap.add_argument("samples_b"); ap.add_argument("membership"); ap.add_argument("out")
ap.add_argument("--carry", type=float, default=1e-4, help="median |v| above this = carries flux")
ap.add_argument("--loop", type=float, default=100.0, help="max |v| above this = loop-flagged")
ap.add_argument("--ndraw", type=int, default=2000)
ap.add_argument("--min_carry", type=int, default=3, help="min reps carrying in both for Tier 2")
a = ap.parse_args()

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

# Tier 2 statistic: pooled log2 magnitude relative to per-rep mean over both cells, carrying reps only
# Bugfix 2026-09-28: post-CFF samples have exact zeros (L1 minimization
# pushes internals to 0 when possible). log2(0) = -inf propagates through
# quantile to NaN. EPS=1e-30 gives a finite floor (log2 -> ~-100) without
# affecting any non-zero flux (relative error < 1e-26 for medians ~1e-4).
EPS = 1e-30
mu = (AA.mean(1) + AB.mean(1)) / 2
mu_safe = np.where(mu > 0, mu, 1.0)
q = np.linspace(0, 1, 201)
def w1(idx):
    la = np.log2(AA[idx] / mu_safe[idx, None] + EPS).ravel()
    lb = np.log2(AB[idx] / mu_safe[idx, None] + EPS).ravel()
    return np.abs(np.quantile(la, q) - np.quantile(lb, q)).mean()

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
    # Tier 1: carrying vs not, A vs B
    p_presence = fisher_exact([[nA, len(idx) - nA], [nB, len(idx) - nB]])[1]
    row = dict(subsystem=s, N_rep=len(idx), carry_A=nA, carry_B=nB, carry_both=n_both,
               loop_flagged=n_loop, p_presence=p_presence)
    # Tier 2: redistribution among reps carrying in both, loop-flagged excluded
    idx2 = idx[carry_both[idx] & ~loop[idx]]
    if len(idx2) >= a.min_carry and len(bg) > len(idx2):
        obs = w1(idx2)
        null = np.array([w1(rng.choice(bg, len(idx2), replace=False)) for _ in range(a.ndraw)])
        row.update(N_tier2=len(idx2), W1_obs=obs, null_median=float(np.median(null)),
                   z=(obs - null.mean()) / null.std(),
                   p_emp=(1 + (null >= obs).sum()) / (a.ndraw + 1))  # never reports 0
    rows.append(row)

out = pd.DataFrame(rows).sort_values("z", ascending=False, na_position="last")
out.to_csv(a.out, index=False)
print(f"reps in both cells: {len(common)}; carrying in A {carryA.sum()}, in B {carryB.sum()}, in both {carry_both.sum()}; "
      f"loop-flagged {loop.sum()}; Tier-2 background size {len(bg)}")
print(out.head(15).to_string())
