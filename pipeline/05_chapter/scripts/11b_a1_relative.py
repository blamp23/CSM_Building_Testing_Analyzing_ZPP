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

# Denominator = total tested representatives (summed over contrasts in the group).
# Each bar height = (notable in that class) / (total tested across pooled contrasts) × 100
f1 = sa_counts[sa_counts.family == "F1"].copy()
f1["stage"] = f1.contrast.map(stage_from_f1)
f1_notable = f1.groupby("stage")[CLASS_ORDER].sum()
f1_tested = f1.groupby("stage").total_tested.sum()

f2 = sa_counts[sa_counts.family == "F2"].copy()
f2["trans"] = f2.contrast.map(trans_from_f2)
f2_notable = f2.groupby("trans")[CLASS_ORDER].sum()
f2_tested = f2.groupby("trans").total_tested.sum()
f2_notable.index = [f"{a}→{b}" for (a, b) in f2_notable.index]
f2_tested.index = [f"{a}→{b}" for (a, b) in f2_tested.index]

notable_df = pd.concat([f1_notable, f2_notable])
tested_s = pd.concat([f1_tested, f2_tested])

# pct of total tested
pct = notable_df.div(tested_s, axis=0) * 100
# "not notable" fills the remainder
pct["Not notable"] = 100 - pct.sum(axis=1)

labels = [f"{h} hpf" for h in f1_notable.index] + list(f2_notable.index)

# ----- plot: broken y-axis (chop out 40-85%) -----
fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(14, 7), sharex=True,
    gridspec_kw={"height_ratios": [1, 4], "hspace": 0.05})

f1_x = np.arange(len(f1_notable))
f2_x = np.arange(len(f2_notable)) + len(f1_notable) + 1
all_x = np.concatenate([f1_x, f2_x])

plot_order = CLASS_ORDER + ["Not notable"]
plot_colors = {**CLASS_COLORS, "Not notable": "#e8e8e8"}

def _stack_bars(ax, show_labels_in_grey=True):
    bottom = np.zeros(len(pct))
    for cls in plot_order:
        heights = pct[cls].values
        bars = ax.bar(all_x, heights, bottom=bottom, color=plot_colors[cls],
                      edgecolor="white", linewidth=0.5, label=cls, width=0.75)
        for xi, h, b in zip(all_x, heights, bottom):
            cx = xi; cy = b + h/2
            if cls == "Not notable":
                if show_labels_in_grey and h >= 10 and cy >= 85:
                    ax.text(cx, cy, f"{h:.0f}%", ha="center", va="center",
                             fontsize=9, color="#555555")
            else:
                if h >= 2.5:
                    ax.text(cx, cy, f"{h:.1f}%", ha="center", va="center",
                             fontsize=9, color="white", weight="bold")
        bottom += heights

_stack_bars(ax_top, show_labels_in_grey=True)
_stack_bars(ax_bot, show_labels_in_grey=False)

# totals annotation on the TOP axis only (above the 100% level)
for xi, lbl in zip(all_x, labels):
    n_notable = int(notable_df.iloc[list(labels).index(lbl)].sum())
    n_tested = int(tested_s.iloc[list(labels).index(lbl)])
    ax_top.text(xi, 102, f"{n_notable}/{n_tested}", ha="center",
                 fontsize=9, color="#333333")

# y-limits: top shows 85-105, bottom shows 0-35
ax_top.set_ylim(85, 105)
ax_bot.set_ylim(0, 35)
ax_top.set_yticks([90, 100])
ax_top.set_yticklabels(["90%", "100%"])
ax_bot.set_yticks([0, 10, 20, 30])
ax_bot.set_yticklabels(["0%", "10%", "20%", "30%"])

# hide shared spines + draw diagonal break marks
ax_top.spines["bottom"].set_visible(False)
ax_bot.spines["top"].set_visible(False)
ax_top.tick_params(labelbottom=False, bottom=False)

d = 0.012  # diagonal size
kwargs = dict(transform=ax_top.transAxes, color="k", clip_on=False,
               linewidth=1.2)
ax_top.plot((-d, +d), (-d*4, +d*4), **kwargs)
ax_top.plot((1 - d, 1 + d), (-d*4, +d*4), **kwargs)
kwargs.update(transform=ax_bot.transAxes)
ax_bot.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax_bot.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

# x-ticks + dividers on bottom
ax_bot.set_xticks(all_x)
ax_bot.set_xticklabels(labels, fontsize=11)
ax_bot.axvline((f1_x[-1] + f2_x[0]) / 2, color="#aaaaaa", linestyle="--", alpha=0.6)
ax_top.axvline((f1_x[-1] + f2_x[0]) / 2, color="#aaaaaa", linestyle="--", alpha=0.6)
ax_bot.text(f1_x.mean(), -5.5, "F1 — cross-regime at each stage\n(pooled over 3 contrasts)",
             ha="center", fontsize=11, color="#333333")
ax_bot.text(f2_x.mean(), -5.5, "F2 — adjacent-stage within regime\n(pooled over 3 regimes)",
             ha="center", fontsize=11, color="#333333")

# shared y-label spanning both axes
fig.text(0.015, 0.5, "Share of total tested representatives\n(notable split by metabolic class)",
          rotation=90, va="center", ha="center", fontsize=12)

ax_top.grid(axis="y", linestyle=":", alpha=0.3)
ax_bot.grid(axis="y", linestyle=":", alpha=0.3)

# title on top axis
ax_top.set_title("A1  Where the differences live  (fraction of all tested representatives)\n"
                  "Top number = notable / tested   ·   Grey = not notable (|centered log$_2$FC| < 1)   ·   Y-axis break at 35-85%",
                  loc="left", fontsize=13, weight="bold", pad=10)

# legend on the right side
ax_top.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=10,
               frameon=True, title="Class")

fig.tight_layout(rect=[0.03, 0, 1, 1])
fig.savefig(OUT_A / "A1_relative.pdf")
fig.savefig(OUT_A / "A1_relative.png", dpi=200)
plt.close(fig)
print(f"wrote {OUT_A}/A1_relative.{{pdf,png}}")
