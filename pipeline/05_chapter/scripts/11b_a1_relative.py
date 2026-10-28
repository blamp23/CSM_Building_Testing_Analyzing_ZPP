"""
11b_a1_relative.py — remake A1 as 100% stacked composition (relative abundance).

Reads SA_counts.csv (written by 11_master_figures.py). Produces:
    ppt_figures/A/A1_relative.{pdf,png}
"""
import re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FIGURES

OUT_A = FIGURES.parent / "ppt_figures" / "A"
sa_counts = pd.read_csv(OUT_A / "SA_counts.csv")

HPFS = (24, 48, 72, 96, 120)
CLASS_ORDER = ["Lipid", "Amino acid", "Sterol / bile", "Nucleotide",
               "Carbohydrate", "Cofactor / vitamin", "Transport / other"]
CLASS_COLORS = {
    "Lipid":              "#2CA02C",
    "Amino acid":         "#FF7F0E",
    "Sterol / bile":      "#D62728",
    "Nucleotide":         "#9467BD",
    "Carbohydrate":       "#8C564B",
    "Cofactor / vitamin": "#E377C2",
    "Transport / other":  "#7F7F7F",
}

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 12,
    "axes.titlesize": 14, "axes.labelsize": 13,
    "xtick.labelsize": 12, "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})


# Pool to F1 by stage and F2 by transition
def stage_from_f1(c):
    m = re.match(r"[A-Z]+_(\d+)_vs_[A-Z]+_\d+", c)
    return int(m.group(1)) if m else None
def trans_from_f2(c):
    m = re.match(r"[A-Z]+_(\d+)_vs_[A-Z]+_(\d+)", c)
    return (int(m.group(1)), int(m.group(2))) if m else None

f1 = sa_counts[sa_counts.family == "F1"].copy()
f1["stage"] = f1.contrast.map(stage_from_f1)
f1_pooled = f1.groupby("stage")[CLASS_ORDER].sum()

f2 = sa_counts[sa_counts.family == "F2"].copy()
f2["trans"] = f2.contrast.map(trans_from_f2)
f2_pooled = f2.groupby("trans")[CLASS_ORDER].sum()
f2_pooled.index = [f"{a}→{b}" for (a, b) in f2_pooled.index]

pool_df = pd.concat([f1_pooled, f2_pooled])
totals = pool_df.sum(axis=1)
pct = pool_df.div(totals, axis=0) * 100  # row sums to 100

labels = [f"{h} hpf" for h in f1_pooled.index] + list(f2_pooled.index)

# ----- plot -----
fig, ax = plt.subplots(figsize=(14, 7))
f1_x = np.arange(len(f1_pooled))
f2_x = np.arange(len(f2_pooled)) + len(f1_pooled) + 1
all_x = np.concatenate([f1_x, f2_x])

bottom = np.zeros(len(pct))
for cls in CLASS_ORDER:
    heights = pct[cls].values
    ax.bar(all_x, heights, bottom=bottom, color=CLASS_COLORS[cls],
            edgecolor="white", linewidth=0.5, label=cls, width=0.75)
    # In-cell percentage label if segment > 4%
    for xi, h, b in zip(all_x, heights, bottom):
        if h >= 4.0:
            ax.text(xi, b + h/2, f"{h:.0f}%", ha="center", va="center",
                     fontsize=9, color="white", weight="bold")
    bottom += heights

# total n above each bar
for xi, n in zip(all_x, totals.values):
    ax.text(xi, 102, f"n={int(n)}", ha="center", fontsize=10, color="#333333")

ax.set_xticks(all_x)
ax.set_xticklabels(labels, fontsize=11)
ax.axvline((f1_x[-1] + f2_x[0]) / 2, color="#aaaaaa", linestyle="--", alpha=0.6)
ax.text(f1_x.mean(), -8, "F1 — cross-regime at each stage\n(pooled over 3 contrasts)",
         ha="center", fontsize=11, color="#333333")
ax.text(f2_x.mean(), -8, "F2 — adjacent-stage within regime\n(pooled over 3 regimes)",
         ha="center", fontsize=11, color="#333333")

ax.set_ylim(0, 108)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
ax.set_ylabel("Share of notable representatives\nby metabolic class")
ax.set_title("A1  Where the differences live  (relative composition)",
              loc="left", fontsize=15, weight="bold", pad=12)
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=10,
           frameon=True, title="Class")
ax.grid(axis="y", linestyle=":", alpha=0.3)

fig.tight_layout()
fig.savefig(OUT_A / "A1_relative.pdf")
fig.savefig(OUT_A / "A1_relative.png", dpi=200)
plt.close(fig)
print(f"wrote {OUT_A}/A1_relative.{{pdf,png}}")
