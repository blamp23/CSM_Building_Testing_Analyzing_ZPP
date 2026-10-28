"""
12b_clock_dotplot.py — ordering-first clock-panel figure for Finding 6.

Heatmap colours encode level; the finding is ordering (which regime holds the
highest value per gene). A dot plot makes the ordering directly visible: 3
dots per gene, LD in colour, BL + D in a shared neutral so the eye reads them
as "constants." Rows sorted by LD margin descending; the 9 panel genes
(margin ≥ 0.5) have a maroon strip to the left.

Writes:
    ppt_figures/finding6/6c_clock_dotplot.{pdf,png}
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import TABLES, FIGURES

OUT = FIGURES.parent / "ppt_figures" / "finding6"
OUT.mkdir(parents=True, exist_ok=True)

LD_COLOR = "#F6B52C"       # orange — the cycling regime
CONST_COLOR = "#555555"    # neutral — BL and D grouped as "constants"
PANEL_COLOR = "#800020"    # maroon strip for the panel-9 rows
BL_MARKER = "o"
D_MARKER = "s"
LD_MARKER = "D"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 12,
    "axes.titlesize": 14, "axes.labelsize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})

# ---------- load clock table ----------
clk = pd.read_csv(TABLES / "table_6_clock_gene_zfpkm_120.csv")
clk = clk.dropna(subset=["BL_120", "D_120", "LD_120"]).copy()
clk["ld_margin"] = clk.LD_120 - clk[["BL_120", "D_120"]].max(axis=1)
# ordering label
def _order(row):
    bl, d, ld = row.BL_120, row.D_120, row.LD_120
    if bl <= d < ld: return "BL ≤ D < LD"
    if d < bl < ld: return "D < BL < LD"
    if bl > d and bl > ld: return "BL highest"
    if d > bl and d > ld: return "D highest"
    return "mixed"
clk["ordering"] = clk.apply(_order, axis=1)
clk = clk.sort_values("ld_margin", ascending=False).reset_index(drop=True)
print(f"[load] {len(clk)} clock genes with zFPKM data")

PANEL_MASK = clk.ld_margin >= 0.5
n_panel = PANEL_MASK.sum()
print(f"[load] {n_panel} genes in panel (LD margin ≥ 0.5)")

# ---------- plot ----------
fig, ax = plt.subplots(figsize=(10, max(7, 0.35 * len(clk) + 1.5)))

y = np.arange(len(clk))[::-1]  # top-to-bottom = highest LD margin first

# panel strip (maroon rectangle on the left margin for panel 9 rows)
for i, (_, row) in enumerate(clk.iterrows()):
    y_pos = y[i]
    if PANEL_MASK.iloc[i]:
        ax.axhspan(y_pos - 0.4, y_pos + 0.4, xmin=-0.02, xmax=-0.005,
                    color=PANEL_COLOR, clip_on=False, zorder=5)
    # faint connecting line BL -> LD span to show the range
    vals = [row.BL_120, row.D_120, row.LD_120]
    ax.plot([min(vals), max(vals)], [y_pos, y_pos], "-", color="#dddddd",
             linewidth=1.2, zorder=1)

# dots: BL (circle), D (square), LD (diamond)
ax.scatter(clk.BL_120, y, s=90, marker=BL_MARKER, color=CONST_COLOR,
            edgecolor="black", linewidth=0.8, label="BL (constant)", zorder=5)
ax.scatter(clk.D_120, y, s=90, marker=D_MARKER, color=CONST_COLOR,
            edgecolor="black", linewidth=0.8, label="D  (constant)", zorder=5)
ax.scatter(clk.LD_120, y, s=130, marker=LD_MARKER, color=LD_COLOR,
            edgecolor="black", linewidth=1.0, label="LD (cycling)", zorder=10)

# reference lines
ax.axvline(0, color="black", linewidth=0.8, linestyle="-", alpha=0.4)
ax.axvline(-3, color="#999999", linewidth=0.8, linestyle=":", alpha=0.7)
ax.text(0.08, len(clk) - 0.4, "typical\nexpressed", fontsize=8.5, color="#666666",
         ha="left", va="top")
ax.text(-2.9, len(clk) - 0.4, "noise\nfloor (−3)", fontsize=8.5, color="#999999",
         ha="left", va="top")

# horizontal divider between panel 9 and the rest
if n_panel > 0 and n_panel < len(clk):
    y_divider = y[n_panel - 1] - 0.5
    ax.axhline(y_divider, color="#999999", linewidth=1.0, linestyle="--",
                alpha=0.5)

# y-tick labels with ordering annotation
y_labels = []
for _, row in clk.iterrows():
    sym = row.symbol
    ord_short = {"BL ≤ D < LD": "★", "D < BL < LD": "◆",
                  "BL highest": "", "D highest": "", "mixed": ""}.get(row.ordering, "")
    y_labels.append(f"{sym}  {ord_short}")
ax.set_yticks(y)
ax.set_yticklabels(y_labels, fontsize=10)
ax.set_xlim(-3.5, 3.5)

ax.set_xlabel("zFPKM at 120 hpf     ←  below typical     above typical  →",
               fontsize=11)
ax.set_title("Figure 6c — Clock panel ordering at 120 hpf: LD sits to the right of both constants\n"
              "LD highest on 17/23 genes; by ≥ 0.5 zFPKM on 9 (★ 5 BL ≤ D < LD   ◆ 4 D < BL < LD).\n"
              "BL and D swap places gene to gene; LD does not.",
              fontsize=11, weight="bold", loc="left", pad=12)
ax.legend(loc="lower right", frameon=True, fontsize=10)
ax.grid(axis="x", linestyle=":", alpha=0.3)

fig.tight_layout()
fig.savefig(OUT / "6c_clock_dotplot.pdf")
fig.savefig(OUT / "6c_clock_dotplot.png", dpi=200)
plt.close(fig)
print(f"wrote {OUT}/6c_clock_dotplot.{{pdf,png}}")
