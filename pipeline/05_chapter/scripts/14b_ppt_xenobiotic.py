"""
14b_ppt_xenobiotic.py — presentation-ready figures for the xenobiotic clearance probe.

Covers: generic X with Phase I hydroxylation (NADPH + O2) + Phase II conjugation
(UDP-glucuronate OR glutathione) + export; max X uptake at 0.9 × biomass_max.
A cofactor-drain test, not a toxicity model.

Writes to chapter_assets/ppt_figures/capacity/:
    C5_xenobiotic.{pdf,png}  — 2×2: glucuronide trajectory + GSH trajectory +
                                     BL_96 outlier bars + UDP-glucuronate check
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import TABLES, DATA, FIGURES

OUT = FIGURES.parent / "ppt_figures" / "capacity"
OUT.mkdir(parents=True, exist_ok=True)

REGIMES = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
REGIME_COLORS = {"BL": "#1976D2", "D": "#333333", "LD": "#F6B52C"}
UDP_GLCA_RXN = "MAR09864"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 12,
    "axes.titlesize": 14, "axes.labelsize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})


# ---------- load ----------
xeno = pd.read_csv(TABLES / "xenobiotic_capacity.csv")
bounds = pd.read_csv(DATA / "phase2_exchange_bounds_long.csv")

# UDP-glucuronate fc per cell (envelope value)
udp = bounds[bounds.ex_rxn == UDP_GLCA_RXN][["condition", "hpf", "fc"]].copy()
udp["cell"] = udp.condition + "_" + udp.hpf.astype(str)
xeno = xeno.merge(udp[["cell", "fc"]].rename(columns={"fc": "udp_glca_fc"}),
                   on="cell", how="left")

# correlation glucuronide ceiling vs UDP-glucuronate fc
r = xeno[["max_X_glucuronide", "udp_glca_fc"]].corr().iloc[0, 1]
print(f"[corr] glucuronide ceiling vs UDP-glucuronate envelope fc: r = {r:+.2f}")


# ---------- plot ----------
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.3,
                       left=0.08, right=0.96, top=0.90, bottom=0.08)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])


def _trajectory(ax, column, title, ylabel):
    for reg in REGIMES:
        sub = xeno[xeno.condition == reg].sort_values("hpf")
        ax.plot(sub.hpf, sub[column], "-o", color=REGIME_COLORS[reg],
                 label=reg, linewidth=2.5, markersize=10,
                 markeredgecolor="white", markeredgewidth=1.2)
    ax.set_xticks(list(HPFS))
    ax.set_xlabel("hpf")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=13, weight="bold", pad=6)
    ax.grid(True, linestyle=":", alpha=0.3)
    ax.legend(title="Regime", loc="upper left", fontsize=10, frameon=True)


# A — glucuronide trajectory
_trajectory(ax_a, "max_X_glucuronide",
             "A  Glucuronide route — rises with development; BL_96 outlier",
             "max X clearance (mmol / gDW / hr)")
bl96 = xeno[xeno.cell == "BL_96"]
if len(bl96):
    y_bl96 = float(bl96.max_X_glucuronide.values[0])
    ax_a.annotate(f"BL_96 = {y_bl96:.0f}\n(biomax low → less\ngrowth competes for C)",
                   xy=(96, y_bl96), xytext=(72, y_bl96 + 10),
                   fontsize=10, color="black",
                   arrowprops=dict(arrowstyle="->", color="black", lw=1))

# B — GSH trajectory
_trajectory(ax_b, "max_X_GSH",
             "B  Glutathione route — flat to falling",
             "max X clearance (mmol / gDW / hr)")

# C — BL_96 bar comparison (biomax + ceiling side-by-side at 96 hpf)
cells_96 = xeno[xeno.hpf == 96].sort_values("condition")
x_pos = np.arange(len(cells_96))
w = 0.4
# biomax (left y-axis), clearance (right y-axis)
ax_c2 = ax_c.twinx()
colors96 = [REGIME_COLORS[c] for c in cells_96.condition]
bars1 = ax_c.bar(x_pos - w/2, cells_96.biomax, w, color=colors96, edgecolor="black",
                   linewidth=0.8, label="biomass_max")
bars2 = ax_c2.bar(x_pos + w/2, cells_96.max_X_glucuronide, w, color=colors96,
                   edgecolor="black", linewidth=0.8, alpha=0.55,
                   label="glucuronide clearance")
ax_c.set_xticks(x_pos)
ax_c.set_xticklabels(cells_96.cell, fontsize=11)
ax_c.set_ylabel("biomass_max")
ax_c2.set_ylabel("glucuronide clearance\n(mmol / gDW / hr)")
ax_c.set_title("C  At 96 hpf, BL has lower biomax → higher clearance\n"
                 "(same isoleucine topology as §3.1 / §3.2)",
                 loc="left", fontsize=13, weight="bold", pad=6)
for b, val in zip(bars1, cells_96.biomax):
    ax_c.text(b.get_x() + b.get_width()/2, val + 0.3, f"{val:.1f}",
                ha="center", fontsize=9)
for b, val in zip(bars2, cells_96.max_X_glucuronide):
    ax_c2.text(b.get_x() + b.get_width()/2, val + 3, f"{val:.0f}",
                 ha="center", fontsize=9)
ax_c.grid(axis="y", linestyle=":", alpha=0.3)

# D — UDP-glucuronate correlation check
ax_d.scatter(xeno.udp_glca_fc, xeno.max_X_glucuronide,
              s=140, c=[REGIME_COLORS[c] for c in xeno.condition],
              edgecolor="black", linewidth=0.8, zorder=5)
for _, row in xeno.iterrows():
    ax_d.annotate(row.cell, (row.udp_glca_fc, row.max_X_glucuronide),
                   xytext=(5, 5), textcoords="offset points", fontsize=8,
                   alpha=0.7)
ax_d.set_xlabel("UDP-glucuronate envelope fold-change (measured)")
ax_d.set_ylabel("glucuronide clearance ceiling")
ax_d.set_title(f"D  UDP-glucuronate fc vs ceiling — r = {r:+.2f}\n"
                 "ceiling does not track the measured input → cofactor supplied internally",
                 loc="left", fontsize=13, weight="bold", pad=6)
ax_d.grid(True, linestyle=":", alpha=0.3)
# regime legend
ax_d.legend(handles=[Line2D([0],[0], marker="o", color="w",
                              markerfacecolor=REGIME_COLORS[r], markersize=11,
                              label=r, markeredgecolor="black") for r in REGIMES],
             title="Regime", loc="upper left", fontsize=10, frameon=True)

fig.suptitle("Figure C5 — Xenobiotic clearance: a cofactor-drain probe, not a toxicity model",
              fontsize=16, weight="bold", y=0.97)
fig.savefig(OUT / "C5_xenobiotic.pdf")
fig.savefig(OUT / "C5_xenobiotic.png", dpi=200)
plt.close(fig)
print(f"wrote {OUT}/C5_xenobiotic.{{pdf,png}}")
