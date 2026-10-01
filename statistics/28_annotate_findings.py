"""
s124 -- annotate the 36 v5 findings with the flags the reviewer asked for:
  - asymmetric_zeros: True if |delta_zero_frac| >= 0.02 (arbitrary but
    matches the eicosanoid/prostaglandin/pentose-glucuronate cluster)
  - neg_control_z: the |z| for that same subsystem in the LD_72 halfA vs
    halfB neg control (blank if subsystem not in neg control family)
  - robustness_tier:
      "A" (build a paragraph on): mid-N (>=15), symmetric zeros, appears in
          both a D and BL contrast OR at two adjacent hpf
      "B" (report, less headline weight): mid-N + symmetric, single contrast
      "C" (report with caveat): small-N (<15) OR asymmetric zeros
"""
import pandas as pd
import numpy as np
from pathlib import Path

ASYMMETRIC_THRESH = 0.02
NEG_SUSPECT_Z = 2.0

df = pd.read_csv("two_tier_v5_nonzero_gamma_bh.csv")
neg = pd.read_csv("s120_neg_control_nonzero.csv")

# neg control z per subsystem
neg_z = neg.set_index("subsystem")["z"].to_dict()
df["neg_control_z"] = df["subsystem"].map(neg_z)
df["neg_control_absz"] = df["neg_control_z"].abs()

df["asymmetric_zeros"] = df["delta_zero_frac"].abs() >= ASYMMETRIC_THRESH

# Findings only
fin = df[df["is_finding"].fillna(False)].copy()

# Identify subsystems appearing in >=2 contrasts (any hpf, any A cell)
by_subs = fin.groupby("subsystem")["contrast"].nunique().to_dict()
fin["n_contrasts_at_this_subsystem"] = fin["subsystem"].map(by_subs)

# Also flag subsystems that appear at BOTH a D- and BL-vs-LD contrast
def has_both_D_and_BL(subs):
    contrasts = fin[fin.subsystem == subs].contrast.tolist()
    has_D = any(c.startswith("D") and "_vs_LD" in c for c in contrasts)
    has_BL = any(c.startswith("BL") and "_vs_LD" in c for c in contrasts)
    return has_D and has_BL
fin["has_D_and_BL"] = fin["subsystem"].apply(has_both_D_and_BL)

# Adjacent hpf flag - do two findings for same subsystem span hpfs 24/48, 48/72, 72/96, 96/120?
def hpf_of(contrast):
    A, _ = contrast.split("_vs_")
    for i, c in enumerate(A):
        if c.isdigit(): return int(A[i:])
    return None

def has_adjacent_hpf(subs):
    hpfs = sorted(set(hpf_of(c) for c in fin[fin.subsystem == subs].contrast if hpf_of(c) is not None))
    adjacent_pairs = [(24,48),(48,72),(72,96),(96,120)]
    for a,b in adjacent_pairs:
        if a in hpfs and b in hpfs: return True
    return False
fin["has_adjacent_hpf"] = fin["subsystem"].apply(has_adjacent_hpf)

# Robustness tier
def tier(row):
    if row.N_tier2 < 15 or row.asymmetric_zeros:
        return "C"
    if row.has_D_and_BL or row.has_adjacent_hpf:
        return "A"
    return "B"
fin["robustness_tier"] = fin.apply(tier, axis=1)

# Neg control suspect flag
fin["neg_control_suspect"] = fin["neg_control_absz"] >= NEG_SUSPECT_Z

# Sort: tier A first, then by q_bh within tier
fin = fin.sort_values(["robustness_tier", "q_bh"])

cols_show = ["robustness_tier", "contrast", "subsystem", "N_tier2",
             "zero_frac_A", "zero_frac_B", "asymmetric_zeros",
             "W1_obs", "z", "p_primary", "q_bh",
             "neg_control_absz", "neg_control_suspect",
             "n_contrasts_at_this_subsystem", "has_D_and_BL", "has_adjacent_hpf",
             "p_fisher_zeros"]
fin[cols_show].to_csv("s124_findings_annotated.csv", index=False)

print("==== FINDINGS BY ROBUSTNESS TIER ====\n")
for t in ["A", "B", "C"]:
    sub = fin[fin.robustness_tier == t]
    print(f"\n### Tier {t} ({len(sub)} findings)")
    if t == "A": print("(paragraph-quality: mid-N, symmetric zeros, appears in both D and BL vs LD or at adjacent hpf)")
    elif t == "B": print("(report-quality: mid-N, symmetric, single contrast, non-suspect)")
    else: print("(caveated: small-N or asymmetric-zero-driven)")
    print(sub[["contrast","subsystem","N_tier2","z","q_bh","asymmetric_zeros",
               "neg_control_suspect"]].to_string(index=False))

print(f"\n==== NEG-CONTROL-SUSPECT (subsystem's neg control |z| >= {NEG_SUSPECT_Z}) ====\n")
suspect = fin[fin.neg_control_suspect]
print(suspect[["contrast","subsystem","z","q_bh","neg_control_absz"]].to_string(index=False))

print(f"\ntotal findings: {len(fin)}, tiers A/B/C: {sum(fin.robustness_tier=='A')}/{sum(fin.robustness_tier=='B')}/{sum(fin.robustness_tier=='C')}")
