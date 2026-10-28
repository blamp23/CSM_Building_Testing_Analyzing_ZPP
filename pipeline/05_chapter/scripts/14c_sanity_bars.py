"""
14c_sanity_bars.py — two simple capability-probe figures.

Writes:
  ppt_figures/capacity/C_sanity_glucose_palmitate.{pdf,png}
  ppt_figures/capacity/C_ceiling_rises_with_O2.{pdf,png}
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import SUSCEPTIBILITY, TABLES, FIGURES

OUT = FIGURES.parent / "ppt_figures" / "capacity"

REGIMES = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
REGIME_COLORS = {"BL": "#1976D2", "D": "#333333", "LD": "#F6B52C"}
BASE_COLOR = "#999999"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 12,
    "axes.titlesize": 14, "axes.labelsize": 13,
    "xtick.labelsize": 10, "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})


# ============================================================
# C_sanity_glucose_palmitate — bar graph
# ============================================================
swp = pd.read_csv(TABLES / "table_S7_substrate_yields.csv")

model_order = [f"base_{h}" for h in HPFS] + [f"{r}_{h}" for r in REGIMES for h in HPFS]

fig, axes = plt.subplots(1, 2, figsize=(16, 5.5))

for ax, substrate in [
    (axes[0], "glucose"),
    (axes[1], "FA 16:0 (palmitate)"),
]:
    # fuzzy-match palmitate in substrate names
    if substrate == "FA 16:0 (palmitate)":
        rows = swp[swp.substrate.str.lower().str.contains("palmit|16:0") &
                     ~swp.substrate.str.lower().str.contains("coa|carnitine")]
    else:
        rows = swp[swp.substrate == substrate]
    data = rows.set_index("model").atp_aerobic.reindex(model_order)
    colors = []
    for m in model_order:
        if m.startswith("base_"): colors.append(BASE_COLOR)
        else:
            reg = m.split("_")[0]
            colors.append(REGIME_COLORS[reg])
    x = np.arange(len(model_order))
    bars = ax.bar(x, data.values, color=colors, edgecolor="black", linewidth=0.6)
    # reference line = base model yield (what we're comparing extracted models against)
    ref_val = data["base_24"] if "base_24" in data.index else data.dropna().iloc[0]
    ax.axhline(ref_val, color="#800020", linewidth=1.4, linestyle="--")
    ax.text(len(x) - 0.5, ref_val, f"  base = {ref_val:.1f}",
             color="#800020", fontsize=10, va="center")

    # label bars with their value above
    for b, v in zip(bars, data.values):
        if pd.notna(v):
            ax.text(b.get_x() + b.get_width()/2, v + max(data.max(), 1)*0.015,
                     f"{v:.1f}", ha="center", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(model_order, rotation=45, ha="right", fontsize=9)
    # dividers between base and each regime
    for sep in [len(HPFS), 2*len(HPFS), 3*len(HPFS)]:
        ax.axvline(sep - 0.5, color="#aaaaaa", linewidth=0.6, linestyle=":")
    ax.set_ylabel("ATP yield per mol substrate")
    ax.set_title(substrate, loc="left", fontsize=14, weight="bold")
    ax.grid(axis="y", linestyle=":", alpha=0.3)

fig.tight_layout()
fig.savefig(OUT / "C_sanity_glucose_palmitate.pdf")
fig.savefig(OUT / "C_sanity_glucose_palmitate.png", dpi=200)
plt.close(fig)
print(f"wrote {OUT}/C_sanity_glucose_palmitate.{{pdf,png}}")


# ============================================================
# C_ceiling_rises_with_O2 — line graph
# ============================================================
phys = pd.read_csv(SUSCEPTIBILITY / "atp_ceiling_physio_ocr.csv")
phys["hpf"] = phys.cell.str.extract(r"_(\d+)$").astype(int)
phys["condition"] = phys.cell.str.extract(r"^([A-Z]+)_")
free = phys[phys.regime == "free"]

fig, ax = plt.subplots(figsize=(11, 6.5))
for reg in REGIMES:
    sub = free[free.condition == reg].sort_values("hpf")
    ax.plot(sub.hpf, sub.atp_max, "-o", color=REGIME_COLORS[reg],
             label=reg, linewidth=2.8, markersize=12,
             markeredgecolor="white", markeredgewidth=1.4)

# per-stage value labels + spread annotations
for h in HPFS:
    stage = free[free.hpf == h].set_index("condition").atp_max.reindex(list(REGIMES))
    hi, lo = stage.max(), stage.min()
    spread_pct = 100 * (hi - lo) / stage.mean()
    # value labels, offset per regime to avoid overlap
    for i, reg in enumerate(REGIMES):
        v = stage[reg]
        dy = {0: 32, 1: -40, 2: 65}[i]  # BL above, D below, LD above-above
        ax.annotate(f"{v:.0f}", xy=(h, v), xytext=(0, dy),
                     textcoords="offset points", ha="center", fontsize=9,
                     color=REGIME_COLORS[reg], weight="bold")
    # spread beneath x axis
    ax.text(h, lo - (hi - lo) * 1.3 - 100, f"spread\n{spread_pct:.1f}%",
             ha="center", fontsize=9, color="#555555")

ax.set_xticks(list(HPFS))
ax.set_xlabel("hpf", labelpad=35)
ax.set_ylabel("ATP maintenance ceiling\n(mmol / gDW / hr)")
ax.set_title("ATP ceiling rises with the O$_2$ cap; cross-regime spread 1-8%",
              loc="left", fontsize=14, weight="bold")
ax.legend(title="Regime", loc="upper left", fontsize=11, frameon=True)
ax.grid(True, linestyle=":", alpha=0.3)

# second axis: OCR cap
ax2 = ax.twinx()
ocr_vals = [97, 173, 227, 273, 313]
ax2.plot(HPFS, ocr_vals, "--", color="#800020", linewidth=1.8, alpha=0.7,
          label="O$_2$ cap")
ax2.set_ylabel("O$_2$ cap (µmol O$_2$ / gDW / hr)", color="#800020")
ax2.tick_params(axis="y", colors="#800020")
ax2.spines["right"].set_visible(True)
ax2.spines["right"].set_color("#800020")

fig.tight_layout()
fig.savefig(OUT / "C_ceiling_rises_with_O2.pdf")
fig.savefig(OUT / "C_ceiling_rises_with_O2.png", dpi=200)
plt.close(fig)
print(f"wrote {OUT}/C_ceiling_rises_with_O2.{{pdf,png}}")
