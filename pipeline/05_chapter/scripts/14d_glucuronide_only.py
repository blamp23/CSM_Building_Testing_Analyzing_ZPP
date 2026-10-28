"""Standalone glucuronide clearance trajectory figure.

Writes:
  ppt_figures/capacity/C_glucuronide_only.{pdf,png}
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import TABLES, FIGURES

OUT = FIGURES.parent / "ppt_figures" / "capacity"

REGIMES = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
REGIME_COLORS = {"BL": "#1976D2", "D": "#333333", "LD": "#F6B52C"}

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 12,
    "axes.titlesize": 14, "axes.labelsize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})

xeno = pd.read_csv(TABLES / "xenobiotic_capacity.csv")

fig, ax = plt.subplots(figsize=(10, 6))
for reg in REGIMES:
    sub = xeno[xeno.condition == reg].sort_values("hpf")
    ax.plot(sub.hpf, sub.max_X_glucuronide, "-o", color=REGIME_COLORS[reg],
             label=reg, linewidth=2.8, markersize=12,
             markeredgecolor="white", markeredgewidth=1.4)
    # value labels at each point
    for _, row in sub.iterrows():
        dy = {0: 8, 1: -16, 2: 8}[list(REGIMES).index(reg)]
        ax.annotate(f"{row.max_X_glucuronide:.0f}",
                     xy=(row.hpf, row.max_X_glucuronide),
                     xytext=(0, dy), textcoords="offset points",
                     ha="center", fontsize=9, color=REGIME_COLORS[reg],
                     weight="bold")

ax.set_xticks(list(HPFS))
ax.set_xlabel("hpf")
ax.set_ylabel("max X clearance via glucuronide route\n(mmol / gDW / hr)")
ax.set_title("Glucuronide route — rises with developmental OCR",
              loc="left", fontsize=14, weight="bold")
ax.legend(title="Regime", loc="upper left", fontsize=11, frameon=True)
ax.grid(True, linestyle=":", alpha=0.3)

fig.tight_layout()
fig.savefig(OUT / "C_glucuronide_only.pdf")
fig.savefig(OUT / "C_glucuronide_only.png", dpi=200)
plt.close(fig)
print(f"wrote {OUT}/C_glucuronide_only.{{pdf,png}}")
