"""
14_ppt_capacity.py — presentation-ready figures for the ATP capacity-probe slide.

Covers: ATP yield sweep (sanity check), ATP ceiling under physiological OCR
(free + per-cell + matched floor variants), glucose-only control, and the
"flat baseline across regimes" conclusion.

Writes to chapter_assets/ppt_figures/capacity/:
    C1_substrate_yields.{pdf,png}    — Table S7 heatmap (44 substrates × 20 models)
    C2_atp_ceiling.{pdf,png}         — free + per-cell floor side-by-side; stage ramp
    C3_glucose_only.{pdf,png}        — full envelope vs glucose-only, per cell
    C4_atp_per_o2.{pdf,png}          — ATP per O2 falling with development
    C_combined.{pdf,png}             — one-slide summary
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import SUSCEPTIBILITY, TABLES, FIGURES

OUT = FIGURES.parent / "ppt_figures" / "capacity"
OUT.mkdir(parents=True, exist_ok=True)

REGIMES = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
REGIME_COLORS = {"BL": "#1976D2", "D": "#333333", "LD": "#F6B52C"}
OCR_CAP = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 12,
    "axes.titlesize": 14, "axes.labelsize": 12,
    "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})


def _save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)
    print(f"  wrote {name}.{{pdf,png}}")


# ---------- load data ----------
print("[load]")
phys = pd.read_csv(SUSCEPTIBILITY / "atp_ceiling_physio_ocr.csv")
glo = pd.read_csv(SUSCEPTIBILITY / "atp_ceiling_glucose_only.csv")
swp = pd.read_csv(TABLES / "table_S7_substrate_yields.csv")

phys["hpf"] = phys.cell.str.extract(r"_(\d+)$").astype(int)
phys["condition"] = phys.cell.str.extract(r"^([A-Z]+)_")
glo["hpf"] = glo.cell.str.extract(r"_(\d+)$").astype(int)
glo["condition"] = glo.cell.str.extract(r"^([A-Z]+)_")


# ============================================================
# C1 — substrate yield sanity-check heatmap
# ============================================================
print("\n[C1] substrate yields heatmap")

# pivot: substrates x models (atp_aerobic), mark unusable as NaN
mat = swp.pivot_table(index="substrate", columns="model",
                        values="atp_aerobic", aggfunc="first")
# column order: base per hpf | BL per hpf | D per hpf | LD per hpf
cols = [f"base_{h}" for h in HPFS] + \
       [f"{reg}_{h}" for reg in REGIMES for h in HPFS]
cols = [c for c in cols if c in mat.columns]
mat = mat[cols]

# sort substrates: usable in all 15 extracted first (central carbon),
# then variable, then always-unusable
extracted = [c for c in cols if not c.startswith("base_")]
usable_counts = (mat[extracted] > 0.5).sum(axis=1)
mat = mat.loc[usable_counts.sort_values(ascending=False).index]

# diverging centered at 31.5 (textbook glucose yield), capped at 150
vmax = 150
cmap = plt.cm.viridis.copy()
cmap.set_bad("#dddddd")

fig_c1, ax_c1 = plt.subplots(figsize=(14, 11))
display = mat.values.astype(float).copy()
display[display < 0.5] = np.nan  # unusable → grey
im = ax_c1.imshow(display, aspect="auto", cmap=cmap, vmin=0, vmax=vmax,
                   interpolation="nearest")

ax_c1.set_xticks(range(len(cols)))
ax_c1.set_xticklabels([c.replace("_", "\n") for c in cols], fontsize=8, rotation=0)
ax_c1.set_yticks(range(len(mat)))
ax_c1.set_yticklabels(mat.index, fontsize=9)
# vertical lines between model groups
for sep in [len(HPFS), 2*len(HPFS), 3*len(HPFS)]:
    ax_c1.axvline(sep - 0.5, color="black", linewidth=1.4)
# group labels on top
block_centers = {"base": 2, "BL": 7, "D": 12, "LD": 17}
for lbl, x in block_centers.items():
    if x < len(cols):
        ax_c1.text(x, -1.5, lbl, ha="center", fontsize=12, weight="bold")

cbar = fig_c1.colorbar(im, ax=ax_c1, pad=0.015, fraction=0.03,
                         shrink=0.8, extend="max")
cbar.set_label("Aerobic ATP yield  (per mol substrate)\n"
                "grey = not usable (yield < 0.5 ATP)",
                fontsize=10)
cbar.ax.axhline(31.5, color="white", linewidth=1.5)
cbar.ax.text(1.5, 31.5, " glucose\n complete\n oxidation (31.5)",
              fontsize=8, va="center", color="black")

ax_c1.set_title("Figure C1 — ATP yield sanity check: all 15 extracted models return 31.5 ATP per glucose\n"
                 "9 substrates not usable in extracted (acetate, glycerol, ketone bodies, DHA, EPA); "
                 "His, Tyr, Met, Trp, Lys fall short in some cells",
                 fontsize=12, weight="bold", loc="left", pad=18)
fig_c1.tight_layout()
_save(fig_c1, "C1_substrate_yields")


# ============================================================
# C2 — ATP ceiling: free biomass (primary) + per-cell floor
# ============================================================
print("\n[C2] ATP ceiling")

fig_c2, axes_c2 = plt.subplots(1, 2, figsize=(14, 5.5), sharex=True)

def _ceiling_lines(ax, df_subset, title):
    for reg in REGIMES:
        sub = df_subset[df_subset.condition == reg].sort_values("hpf")
        ax.plot(sub.hpf, sub.atp_max, "-o", color=REGIME_COLORS[reg],
                 label=reg, linewidth=2.5, markersize=10,
                 markeredgecolor="white", markeredgewidth=1.2)
    # overlay OCR cap scaled to show the ratio
    ax.set_xticks(list(HPFS))
    ax.set_xlabel("hpf")
    ax.set_ylabel("ATP maintenance ceiling\n(mmol / gDW / hr)")
    ax.set_title(title, fontsize=13, weight="bold", loc="left")
    ax.grid(True, linestyle=":", alpha=0.3)
    ax.legend(title="Regime", loc="upper left", fontsize=10, frameon=True)

_ceiling_lines(axes_c2[0], phys[phys.regime == "free"],
                 "A  Free biomass (primary)\nOCR-bound in all 15 cells; regimes within 8% at every stage")
_ceiling_lines(axes_c2[1], phys[phys.regime == "per_cell_floor"],
                 "B  Per-cell biomass floor (0.9 × biomass_max)\nBiomass pin becomes binding; BL_96 isoleucine topology visible")

# annotate BL_96 outlier on panel B
bl_96_percell = phys[(phys.regime == "per_cell_floor") & (phys.cell == "BL_96")]
if len(bl_96_percell):
    axes_c2[1].annotate(f"BL_96 pinned low\n({bl_96_percell.atp_max.values[0]:.0f})",
                          xy=(96, bl_96_percell.atp_max.values[0]),
                          xytext=(72, bl_96_percell.atp_max.values[0] + 200),
                          fontsize=9, color="black",
                          arrowprops=dict(arrowstyle="->", color="black", lw=1))

fig_c2.suptitle("Figure C2 — ATP ceiling scales with developmental OCR cap, not with light regime",
                 fontsize=14, weight="bold", y=1.00)
fig_c2.tight_layout(rect=[0, 0, 1, 0.95])
_save(fig_c2, "C2_atp_ceiling")


# ============================================================
# C3 — glucose-only control vs full envelope
# ============================================================
print("\n[C3] glucose-only vs full envelope")

fig_c3, ax_c3 = plt.subplots(figsize=(12, 6))

cells = [f"{r}_{h}" for h in HPFS for r in REGIMES]
x = np.arange(len(cells))
w = 0.38

# free biomass values from phys and from glo
phys_free = phys[phys.regime == "free"].set_index("cell").atp_max.to_dict()
glo_free = glo[glo.regime == "free"].set_index("cell").atp_max.to_dict()

phys_vals = [phys_free.get(c, 0) for c in cells]
glo_vals = [glo_free.get(c, 0) for c in cells]

ax_c3.bar(x - w/2, phys_vals, w, color="#1976D2", edgecolor="black",
           linewidth=0.8, label="Full envelope")
ax_c3.bar(x + w/2, glo_vals, w, color="#F6B52C", edgecolor="black",
           linewidth=0.8, label="Glucose only (envelope + OCR kept)")
ax_c3.set_xticks(x)
ax_c3.set_xticklabels(cells, rotation=45, ha="right", fontsize=10)
# stage block separators
for s in range(1, len(HPFS)):
    ax_c3.axvline(s*3 - 0.5, color="#999999", linewidth=0.6, linestyle=":")
ax_c3.set_ylabel("ATP maintenance ceiling\n(mmol / gDW / hr)")
ax_c3.set_title("Figure C3 — Glucose-only control: ceiling rises 1.4-1.6× but ATP per glucose < 3\n"
                 "Oxygen, not substrate, sets the ceiling",
                 fontsize=12, weight="bold", loc="left", pad=10)
ax_c3.legend(loc="upper left", fontsize=11, frameon=True)
ax_c3.grid(axis="y", linestyle=":", alpha=0.3)
fig_c3.tight_layout()
_save(fig_c3, "C3_glucose_only")


# ============================================================
# C4 — ATP per O2 falling with hpf
# ============================================================
print("\n[C4] ATP per O2 vs hpf")

phys_free_df = phys[phys.regime == "free"].copy()
# ATP/O2 from the free-biomass ceiling and the OCR cap
# At the ceiling, oxygen uptake equals the OCR cap (since dual is binding)
phys_free_df["atp_per_o2"] = phys_free_df.atp_max / phys_free_df.OCR_cap

fig_c4, ax_c4 = plt.subplots(figsize=(10, 6))
for reg in REGIMES:
    sub = phys_free_df[phys_free_df.condition == reg].sort_values("hpf")
    ax_c4.plot(sub.hpf, sub.atp_per_o2, "-o", color=REGIME_COLORS[reg],
                label=reg, linewidth=2.5, markersize=11,
                markeredgecolor="white", markeredgewidth=1.2)

# reference line at 5.25 = glucose complete-oxidation ratio
ax_c4.axhline(5.25, color="#800020", linewidth=1.5, linestyle="--")
ax_c4.text(123, 5.25, "glucose\ncomplete-\noxidation (5.25)",
            fontsize=10, va="center", color="#800020")

ax_c4.set_xticks(list(HPFS))
ax_c4.set_xlabel("hpf")
ax_c4.set_ylabel("ATP per O$_2$  (free-biomass ceiling / OCR cap)")
ax_c4.set_title("Figure C4 — ATP per O$_2$ falls from ~10.5 at 24 hpf to below 7 at 120 hpf\n"
                 "Excess above 5.25 = non-oxidative ATP (substrate-level phosphorylation + fermentation)",
                 fontsize=12, weight="bold", loc="left", pad=10)
ax_c4.legend(title="Regime", loc="upper right", fontsize=11, frameon=True)
ax_c4.grid(True, linestyle=":", alpha=0.3)
ax_c4.set_xlim(18, 135)
fig_c4.tight_layout()
_save(fig_c4, "C4_atp_per_o2")


# ============================================================
# C_combined — one-slide summary
# ============================================================
print("\n[C_combined]")

fig_cc = plt.figure(figsize=(18, 11))
gs = fig_cc.add_gridspec(2, 2, hspace=0.5, wspace=0.3,
                           left=0.08, right=0.95, top=0.91, bottom=0.08)

# A top-left: free biomass ceiling (the "flat baseline" headline)
ax_a = fig_cc.add_subplot(gs[0, 0])
_ceiling_lines(ax_a, phys[phys.regime == "free"],
                 "A  Free-biomass ATP ceiling (primary)")
ax_a.set_title("A  Free-biomass ATP ceiling\n2.1× rise 24→120 against 3.2× cap; regimes within 8% at every stage",
                 fontsize=12, weight="bold", loc="left")

# B top-right: per-cell floor
ax_b = fig_cc.add_subplot(gs[0, 1])
_ceiling_lines(ax_b, phys[phys.regime == "per_cell_floor"],
                 "B  Per-cell biomass floor")
ax_b.set_title("B  Per-cell biomass floor (0.9 × biomass_max)\nBL_96 pinned low — isoleucine topology (§3.1)",
                 fontsize=12, weight="bold", loc="left")

# C bottom-left: ATP per O2
ax_c = fig_cc.add_subplot(gs[1, 0])
for reg in REGIMES:
    sub = phys_free_df[phys_free_df.condition == reg].sort_values("hpf")
    ax_c.plot(sub.hpf, sub.atp_per_o2, "-o", color=REGIME_COLORS[reg],
               label=reg, linewidth=2.5, markersize=10,
               markeredgecolor="white", markeredgewidth=1.0)
ax_c.axhline(5.25, color="#800020", linewidth=1.5, linestyle="--")
ax_c.text(123, 5.25, "complete-ox (5.25)", fontsize=9, va="center", color="#800020")
ax_c.set_xticks(list(HPFS))
ax_c.set_xlabel("hpf")
ax_c.set_ylabel("ATP per O$_2$")
ax_c.set_title("C  ATP per O$_2$ falls — non-oxidative share shrinks with development",
                 fontsize=12, weight="bold", loc="left")
ax_c.legend(title="Regime", loc="upper right", fontsize=10, frameon=True)
ax_c.grid(True, linestyle=":", alpha=0.3)
ax_c.set_xlim(18, 135)

# D bottom-right: glucose-only control
ax_d = fig_cc.add_subplot(gs[1, 1])
ax_d.bar(x - w/2, phys_vals, w, color="#1976D2", edgecolor="black",
           linewidth=0.6, label="Full envelope")
ax_d.bar(x + w/2, glo_vals, w, color="#F6B52C", edgecolor="black",
           linewidth=0.6, label="Glucose only")
ax_d.set_xticks(x)
ax_d.set_xticklabels(cells, rotation=45, ha="right", fontsize=8)
for s in range(1, len(HPFS)):
    ax_d.axvline(s*3 - 0.5, color="#999999", linewidth=0.6, linestyle=":")
ax_d.set_ylabel("ATP ceiling\n(mmol / gDW / hr)")
ax_d.set_title("D  Glucose-only control — oxygen, not substrate, sets the ceiling\n(opening glucose to 1000 raises ceiling only 1.4-1.6×)",
                 fontsize=12, weight="bold", loc="left")
ax_d.legend(loc="upper left", fontsize=10, frameon=True)
ax_d.grid(axis="y", linestyle=":", alpha=0.3)

fig_cc.suptitle("Figure C — ATP capacity probe: bulk energetic capacity is set by developmental stage, not by light history",
                 fontsize=15, weight="bold", y=0.97)
_save(fig_cc, "C_combined")


print("\n[done]  outputs in", OUT)
