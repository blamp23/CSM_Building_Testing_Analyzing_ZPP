"""
Report step for split-chain check (task #105). Reads the three two-tier CSVs
produced by s60_split_chain.py and writes cross_half_headline.csv +
neg_control_summary.txt.
"""
from pathlib import Path
import numpy as np, pandas as pd

OUT = Path("results/split_chain_v4")

# --- headline cross-half comparison ---
a = pd.read_csv(OUT / "two_tier_D72_vs_LD72_halfA.csv")
b = pd.read_csv(OUT / "two_tier_D72_vs_LD72_halfB.csv")
m = (a[["subsystem", "z", "N_tier2", "loop_flagged"]]
     .rename(columns={"z": "z_halfA", "N_tier2": "N_halfA", "loop_flagged": "loops_halfA"})
     .merge(b[["subsystem", "z", "N_tier2", "loop_flagged"]]
              .rename(columns={"z": "z_halfB", "N_tier2": "N_halfB",
                               "loop_flagged": "loops_halfB"}),
            on="subsystem", how="outer"))

# splice in z_full from the packet's full run if we have it on disk
for candidate in ("results/two_tier_v4/two_tier_D72_vs_LD72.csv",
                  "reviewer_packet_v4/two_tier_D72_vs_LD72.csv"):
    fp = Path(candidate)
    if fp.exists():
        full = pd.read_csv(fp)[["subsystem", "z"]].rename(columns={"z": "z_full"})
        m = m.merge(full, on="subsystem", how="left")
        break
if "z_full" not in m.columns:
    m["z_full"] = np.nan

m["delta_halves"] = m["z_halfA"] - m["z_halfB"]
m = m.sort_values("z_halfA", ascending=False, na_position="last")
m.to_csv(OUT / "cross_half_headline.csv", index=False)

print("==== HEADLINE STABILITY (D_72 vs LD_72) ====")
top = m.dropna(subset=["z_halfA"]).head(20)
print(top[["subsystem", "z_halfA", "z_halfB", "z_full", "delta_halves"]]
      .to_string(index=False))
print()

# reviewer's decision rule: stable if |delta| small vs signal
if m["z_halfA"].notna().any():
    top_finite = m.dropna(subset=["z_halfA", "z_halfB"])
    r = np.corrcoef(top_finite["z_halfA"], top_finite["z_halfB"])[0, 1]
    med_abs_delta = top_finite["delta_halves"].abs().median()
    max_abs_delta = top_finite["delta_halves"].abs().max()
    print(f"corr(z_halfA, z_halfB)         : {r:.3f}")
    print(f"median |z_halfA - z_halfB|      : {med_abs_delta:.3f}")
    print(f"max    |z_halfA - z_halfB|      : {max_abs_delta:.3f}")
    print()

# --- negative control summary ---
neg = pd.read_csv(OUT / "two_tier_LD72_halfA_vs_halfB.csv")
neg_z = neg["z"].dropna()

lines = []
lines.append("==== NEG CONTROL (LD_72 halfA vs LD_72 halfB) ====")
lines.append(f"subsystems tested       : {len(neg_z)}")
lines.append(f"|z| median              : {neg_z.abs().median():.3f}")
lines.append(f"|z| 95th pctile         : {neg_z.abs().quantile(0.95):.3f}")
lines.append(f"|z| max                 : {neg_z.abs().max():.3f}")
lines.append(f"count |z|>2             : {(neg_z.abs()>2).sum()}")
lines.append(f"count |z|>3             : {(neg_z.abs()>3).sum()}")
lines.append("")
lines.append("Top 10 |z| in neg control:")
top_neg = neg.reindex(neg["z"].abs().sort_values(ascending=False).index)
lines.append(top_neg[["subsystem", "N_tier2", "z", "p_emp", "loop_flagged"]]
             .head(10).to_string(index=False))

text = "\n".join(lines)
print(text)
(OUT / "neg_control_summary.txt").write_text(text + "\n")
print(f"\nwrote {OUT/'cross_half_headline.csv'} and {OUT/'neg_control_summary.txt'}")
