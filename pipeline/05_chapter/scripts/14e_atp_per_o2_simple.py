"""Standalone ATP-per-O2 figure. Clean, no busy title.

Writes: ppt_figures/capacity/C_atp_per_o2.{pdf,png}
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import SUSCEPTIBILITY, FIGURES

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

phys = pd.read_csv(SUSCEPTIBILITY / "atp_ceiling_physio_ocr.csv")
phys["hpf"] = phys.cell.str.extract(r"_(\d+)$").astype(int)
phys["condition"] = phys.cell.str.extract(r"^([A-Z]+)_")
free = phys[phys.regime == "free"].copy()
free["atp_per_o2"] = free.atp_max / free.OCR_cap

fig, ax = plt.subplots(figsize=(10, 5.5))
for reg in REGIMES:
    sub = free[free.condition == reg].sort_values("hpf")
    ax.plot(sub.hpf, sub.atp_per_o2, "-o", color=REGIME_COLORS[reg],
             label=reg, linewidth=2.8, markersize=12,
             markeredgecolor="white", markeredgewidth=1.4)

ax.axhline(5.25, color="#800020", linewidth=1.5, linestyle="--")
ax.text(123, 5.25, " complete ox (5.25)", color="#800020",
         fontsize=10, va="center")

ax.set_xticks(list(HPFS))
ax.set_xlabel("hpf")
ax.set_ylabel("ATP per O$_2$")
ax.set_title("Oxidative share of ATP grows with development",
              loc="left", fontsize=14, weight="bold")
ax.legend(title="Regime", loc="upper right", fontsize=11, frameon=True)
ax.grid(True, linestyle=":", alpha=0.3)

fig.tight_layout()
fig.savefig(OUT / "C_atp_per_o2.pdf")
fig.savefig(OUT / "C_atp_per_o2.png", dpi=200)
plt.close(fig)
print(f"wrote {OUT}/C_atp_per_o2.{{pdf,png}}")
