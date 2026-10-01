"""
s120 -- split-chain negative control on v5 samples using the NEW non-zero W1
statistic. Establishes the actual |z| noise floor under the corrected pipeline.

LD_72 pooled reps have 2000 columns = 4 chains x 500 samples in order
(ch1: s1..s500, ch2: s501..s1000, ch3: s1001..s1500, ch4: s1501..s2000).
Split into halfA (cols 1..1000, chains 1+2) and halfB (cols 1001..2000, chains 3+4),
run two-tier with drop-zeros W1 as if they were two conditions. Every z should be ~0.
"""
import numpy as np, pandas as pd
from pathlib import Path
from s119_gamma_bh_nonzero import (run_contrast, apply_family_and_bh, Q_THRESH)
import tempfile, os

REPS_DIR = Path("v5_pull/results/hopsy_loopless_v5_reps")
SUBS_DIR = Path("subsets")

# Read LD_72 pooled reps, split by columns
fp = REPS_DIR / "reps_pooled_LD_72.csv"
df = pd.read_csv(fp, index_col=0)
n_cols = df.shape[1]
assert n_cols == 2000, f"expected 2000 cols, got {n_cols}"
halfA = df.iloc[:, :1000]
halfB = df.iloc[:, 1000:]
print(f"halfA shape: {halfA.shape}, halfB shape: {halfB.shape}")

# Rename columns to sample IDs
halfA.columns = [f"s{i+1}" for i in range(1000)]
halfB.columns = [f"s{i+1}" for i in range(1000)]

with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as fa:
    halfA.to_csv(fa.name, float_format="%.5g")
    A_csv = fa.name
with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as fb:
    halfB.to_csv(fb.name, float_format="%.5g")
    B_csv = fb.name

M_csv = SUBS_DIR / "LD_72_subset_membership.csv"

print("running neg control two-tier with non-zero W1 statistic ...")
df_out = run_contrast(A_csv, B_csv, M_csv, ndraw=2000)
df_out = apply_family_and_bh(df_out, "LD72_halfA_vs_halfB_NEG")
df_out.to_csv("s120_neg_control_nonzero.csv", index=False)

os.unlink(A_csv); os.unlink(B_csv)

# Summary
fam = df_out[df_out["in_family"]]
z = fam["z"].dropna()
print("\n==== NEG CONTROL SUMMARY (v5 samples, non-zero W1 statistic) ====")
print(f"in-family subsystems tested: {len(fam)}")
print(f"|z| median              : {z.abs().median():.3f}")
print(f"|z| 90th pctile         : {z.abs().quantile(0.90):.3f}")
print(f"|z| 95th pctile         : {z.abs().quantile(0.95):.3f}")
print(f"|z| max                 : {z.abs().max():.3f}")
print(f"count |z|>2             : {int((z.abs()>2).sum())}")
print(f"count |z|>3             : {int((z.abs()>3).sum())}")
print(f"count |z|>5             : {int((z.abs()>5).sum())}")
print(f"n findings (q_bh<={Q_THRESH}): {int(df_out['is_finding'].fillna(False).sum())}")

top = fam.reindex(fam['z'].abs().sort_values(ascending=False).index).head(10)
print("\ntop 10 |z| in neg control:")
print(top[["subsystem","N_tier2","zero_frac_A","zero_frac_B",
          "delta_zero_frac","W1_obs","z","p_primary","q_bh"]].to_string(index=False))
