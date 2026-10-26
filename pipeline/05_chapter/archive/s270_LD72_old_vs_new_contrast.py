"""
LD_72_old vs LD_72_new two-tier W1 contrast — method of record
(§11.14 variant: mean normaliser, zeros dropped).

Reads OLD pooled samples (existing) and NEW pooled samples (post-resample).
Runs the two_tier_subsystem_test_fast_mean_dropzero logic directly here
(importable and callable), BH-adjusts p_emp → q_bh, divides W1_obs by the
LD_72 split-chain 95th-pct noise floor per subsystem to get snr.

STOP CRITERION (per §11.19 spot-check):
  If any subsystem NOT in {Vitamin A metabolism, Retinol metabolism,
  Carotenoid biosynthesis} has q_bh ≤ 0.10 AND snr ≥ 2 → report + STOP.
Otherwise annotate retinoid rows with pre_correction_retinoid.
"""
import sys, argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "grace_v4/results/atp_sweep_v7"
NC = REPO / "grace_v4/results/findings_v7_mean_dropzero/negcontrol_per_subsystem.csv"

ap = argparse.ArgumentParser()
ap.add_argument("--old", default=str(REPO / "grace_v4/results/pooled_v7/LD_72_pooled.csv"))
ap.add_argument("--new", required=True, help="Path to new LD_72 pooled samples CSV")
ap.add_argument("--membership", default=str(REPO / "grace_v4/results/subsets_v7/LD_72_subset_membership.csv"))
ap.add_argument("--out", default=str(OUT / "LD72_old_vs_new_W1_subsystem.csv"))
ap.add_argument("--carry", type=float, default=1e-4)
ap.add_argument("--loop", type=float, default=100.0)
ap.add_argument("--ndraw", type=int, default=2000)
ap.add_argument("--min_carry", type=int, default=3)
ap.add_argument("--chunk", type=int, default=200)
a = ap.parse_args()

ZERO_TOL = 1e-10
RETINOID_SUBS = frozenset([
    "Vitamin A metabolism",
    "Retinol metabolism",
    "Carotenoid biosynthesis",
])

print(f"[load] OLD samples: {a.old}", flush=True)
XA = pd.read_csv(a.old, index_col=0)
print(f"[load] NEW samples: {a.new}", flush=True)
XB = pd.read_csv(a.new, index_col=0)
print(f"[load] membership: {a.membership}", flush=True)
mem = pd.read_csv(a.membership)

common = XA.index.intersection(XB.index)
print(f"OLD reps: {len(XA)}  NEW reps: {len(XB)}  common: {len(common)}", flush=True)
AA, AB = np.abs(XA.loc[common].values), np.abs(XB.loc[common].values)
medA, medB = np.median(AA, 1), np.median(AB, 1)
carryA, carryB = medA > a.carry, medB > a.carry
carry_both = carryA & carryB
loop = (AA.max(1) > a.loop) | (AB.max(1) > a.loop)

rep_info = mem.drop_duplicates("representative").set_index("representative")
gpr = rep_info.has_gpr.reindex(common).fillna(0).values == 1
subsys = rep_info.subsystem.reindex(common).fillna("UNMAPPED")
exch = subsys.str.contains("Exchange|Artificial|Transport", case=False).values

# §11.14 method of record: mean normaliser, zeros dropped
AB_concat = np.concatenate([AA, AB], axis=1)
mask_nz = AB_concat >= ZERO_TOL
with np.errstate(all="ignore"):
    masked = np.where(mask_nz, AB_concat, np.nan)
    mu = np.nanmean(masked, axis=1)
mu_safe = np.where(np.isfinite(mu) & (mu > 0), mu, 1.0)
mu_col = mu_safe[:, None]
log_AA = np.where(AA >= ZERO_TOL, np.log2(np.maximum(AA, ZERO_TOL) / mu_col), np.nan)
log_AB = np.where(AB >= ZERO_TOL, np.log2(np.maximum(AB, ZERO_TOL) / mu_col), np.nan)
q = np.linspace(0, 1, 201)

def w1_serial(idx):
    la = log_AA[idx].ravel(); lb = log_AB[idx].ravel()
    la = la[~np.isnan(la)]; lb = lb[~np.isnan(lb)]
    if len(la) < 10 or len(lb) < 10: return np.nan
    return np.abs(np.quantile(la, q) - np.quantile(lb, q)).mean()

def w1_batch(draws, chunk_size):
    n = draws.shape[0]; out = np.empty(n)
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
print(f"Tier-2 background size: {len(bg)}", flush=True)
rng = np.random.default_rng(0)

rows = []
for s in sorted(subsys.unique()):
    idx = np.where(subsys.values == s)[0]
    if len(idx) == 0 or s == "UNMAPPED": continue
    nA, nB = int(carryA[idx].sum()), int(carryB[idx].sum())
    n_both = int(carry_both[idx].sum())
    n_loop = int((loop & carry_both)[idx].sum())
    p_presence = fisher_exact([[nA, len(idx) - nA], [nB, len(idx) - nB]])[1]
    row = dict(subsystem=s, N_rep=len(idx), carry_A=nA, carry_B=nB,
               carry_both=n_both, loop_flagged=n_loop, p_presence=p_presence)
    idx2 = idx[carry_both[idx] & ~loop[idx]]
    if len(idx2) >= a.min_carry and len(bg) > len(idx2):
        obs = w1_serial(idx2)
        draws = np.stack([rng.choice(bg, len(idx2), replace=False)
                          for _ in range(a.ndraw)])
        null = w1_batch(draws, a.chunk)
        null = null[np.isfinite(null)]
        if len(null) < 10 or not np.isfinite(obs):
            row.update(N_tier2=len(idx2), W1_obs=obs, null_median=np.nan,
                       z=np.nan, p_emp=np.nan)
        else:
            row.update(N_tier2=len(idx2), W1_obs=obs,
                       null_median=float(np.median(null)),
                       z=(obs - null.mean()) / null.std(),
                       p_emp=(1 + (null >= obs).sum()) / (a.ndraw + 1))
    rows.append(row)

df = pd.DataFrame(rows)
# BH-FDR
pv = df.p_emp.values
good = np.isfinite(pv)
q_bh = np.full(len(pv), np.nan)
if good.sum():
    pg = pv[good]
    order = np.argsort(pg)
    ranks = np.empty_like(order); ranks[order] = np.arange(1, len(pg) + 1)
    q_sorted = pg[order] * len(pg) / np.arange(1, len(pg) + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_bh[np.where(good)[0]] = q_sorted[np.argsort(order)]
df["q_bh"] = q_bh

# snr: LD_72 split-chain W1 per subsystem as the floor (per user spec:
# "LD_72 split-chain null as the floor"). Both members of this contrast are
# LD_72, so using LD_72-specific noise is the most apt comparison.
nc = pd.read_csv(NC)
ld72_nc = nc[nc.cell == "LD_72"].set_index("subsystem").W1_obs
# Fallback for any subsystem not in LD_72 negcontrol: use across-cell p95
across_p95 = nc.groupby("subsystem").W1_obs.quantile(0.95)
df["noise_floor_LD72"] = df.subsystem.map(ld72_nc)
df["noise_floor_fallback_p95"] = df.subsystem.map(across_p95).fillna(np.median(across_p95))
df["noise_floor"] = df.noise_floor_LD72.fillna(df.noise_floor_fallback_p95)
df["snr"] = df.W1_obs / df.noise_floor

df["is_retinoid"] = df.subsystem.isin(RETINOID_SUBS)
df = df.sort_values("q_bh", ascending=True, na_position="last")
df.to_csv(a.out, index=False)
print(f"\n[write] {a.out}  rows={len(df)}", flush=True)

print("\n=== Subsystems with q_bh ≤ 0.10 ===")
sig = df[df.q_bh <= 0.10].copy()
if sig.empty:
    print("  (none)")
else:
    print(sig[["subsystem","N_tier2","W1_obs","noise_floor","snr","q_bh","z","is_retinoid"]].to_string(index=False))

print("\n=== STOP CRITERION ===")
outside = sig[(~sig.is_retinoid) & (sig.snr >= 2.0)]
if not outside.empty:
    print(f"  TRIGGERED. Non-retinoid subsystems with q ≤ 0.10 AND snr ≥ 2:")
    print(outside[["subsystem","snr","q_bh","z"]].to_string(index=False))
    print("  >>> GLOBAL RESAMPLE REQUIRED")
else:
    retinoid_sig = sig[sig.is_retinoid]
    print(f"  Not triggered. Non-retinoid snr≥2 count: 0.")
    if len(retinoid_sig):
        print(f"  Retinoid subsystems with q ≤ 0.10 ({len(retinoid_sig)}): mark as pre_correction_retinoid.")
        for s in retinoid_sig.subsystem:
            print(f"    - {s}")
    print("  >>> SAFE TO RETAIN v7 SAMPLES; annotate retinoid rows pre_correction_retinoid.")
