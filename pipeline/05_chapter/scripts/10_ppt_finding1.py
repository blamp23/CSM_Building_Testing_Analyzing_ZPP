"""
10_ppt_finding1.py — presentation-ready figures for Finding 1.

Finding 1: shared lipid activation program, 48 to 72 hpf.
Writes a four-figure set to chapter_assets/ppt_figures/finding1/, each 16:9
at 300 DPI, PDF + PNG.

    1a_table.{pdf,png}        — subsystem × regime effect table
    1b_bars.{pdf,png}         — grouped bar chart of mean log2FC (48→72)
    1c_trajectory.{pdf,png}   — per-subsystem × hpf transition line plot
    1d_density.{pdf,png}      — Leukotriene + Arachidonic densities at 5 hpf, 3 regimes

Sign convention in F2 CSV: earlier hpf is A, so negative log2FC = flux rose
across the transition. All figure axes report "72 vs 48" with positive =
rise, i.e. we flip the sign of the raw F2 log2FC for display.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import POOLED, SUBSETS, COMPARISON_TABLES, FIGURES

OUT = FIGURES.parent / "ppt_figures" / "finding1"
OUT.mkdir(parents=True, exist_ok=True)

# ---------- Finding 1 scope ----------
F1_SUBSYSTEMS = [
    "Leukotriene metabolism",
    "Arachidonic acid metabolism",
    "Eicosanoid metabolism",
    "Fatty acid oxidation",
    "Fatty acid activation (cytosolic)",
]
F1_SHORT = {
    "Leukotriene metabolism": "Leukotriene",
    "Arachidonic acid metabolism": "Arachidonic acid",
    "Eicosanoid metabolism": "Eicosanoid",
    "Fatty acid oxidation": "FA oxidation",
    "Fatty acid activation (cytosolic)": "FA activation (cyto)",
}
REGIMES = ("BL", "D", "LD")
REGIME_COLORS = {"BL": "#1976D2", "D": "#555555", "LD": "#F6B52C"}

# all four developmental transitions (midpoint x used for the line plot)
TRANSITIONS = [("24", "48"), ("48", "72"), ("72", "96"), ("96", "120")]
TRANSITION_MID = {("24","48"): 36, ("48","72"): 60, ("72","96"): 84, ("96","120"): 108}

# ---------- PPT style ----------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 14,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 13,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def _save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png")
    plt.close(fig)
    print(f"  wrote {name}.{{pdf,png}}")


def _rise_logfc(row):
    """F2 CSV stores earlier-is-A sign. Flip so positive = rose across transition."""
    return -row["mean_log2FC_centered"]


# ---------- load F2 subsystem summary ----------
print("[load] F2 subsystem summary")
f2 = pd.read_csv(COMPARISON_TABLES / "subsystem_summary_F2_centered.csv")
# Parse the "contrast" column to pull earlier/later hpf + regime
# contrasts look like "BL_24_vs_BL_48" or similar
def _parse_contrast(c):
    parts = c.split("_vs_")
    a = parts[0].split("_"); b = parts[1].split("_")
    return a[0], a[1], b[1]  # regime, earlier_hpf, later_hpf

f2["regime"], f2["hpf_a"], f2["hpf_b"] = zip(*f2.contrast.map(_parse_contrast))
f2["rise_log2FC"] = -f2["mean_log2FC_centered"]

# Finding 1 slice: 48→72 across 3 regimes × 5 subsystems
f1_core = f2[(f2.hpf_a == "48") & (f2.hpf_b == "72") &
             (f2.subsystem.isin(F1_SUBSYSTEMS))].copy()
print(f"[load] {len(f1_core)} Finding 1 rows at 48→72")


# ============================================================
# 1a — table card
# ============================================================
print("[fig] 1a table")
fig, ax = plt.subplots(figsize=(14, 7))
ax.axis("off")
rows = []
for sub in F1_SUBSYSTEMS:
    for reg in REGIMES:
        r = f1_core[(f1_core.subsystem == sub) & (f1_core.regime == reg)]
        if not len(r):
            rows.append([F1_SHORT[sub], reg, "—","—","—","—","—"]); continue
        r = r.iloc[0]
        q = r["enrichment_q"]
        q_str = f"{q:.2g}" if q < 1 else "ns"
        # F2 sign: n_reps_down = direction -1 = later > earlier = "rises"
        if r.rise_log2FC > 0:
            direction_str = f"rises ({int(r.n_reps_down)}/{int(r.n_reps_tested)})"
        else:
            direction_str = f"falls ({int(r.n_reps_up)}/{int(r.n_reps_tested)})"
        rows.append([
            F1_SHORT[sub], reg,
            f"{int(r.n_reps_tested)}",
            direction_str,
            f"{r.rise_log2FC:+.2f}",
            f"{r.coherence:.2f}",
            q_str,
        ])
cols = ["Subsystem", "Regime", "n", "Direction (72 vs 48)", "log2FC", "Coherence", "q"]
tbl = ax.table(cellText=rows, colLabels=cols, loc="center", cellLoc="center",
               colWidths=[0.22,0.07,0.06,0.21,0.11,0.12,0.11])
tbl.auto_set_font_size(False)
tbl.set_fontsize(13)
tbl.scale(1, 1.6)

# shade q column by significance
for i, row in enumerate(rows, start=1):
    try:
        qval = float(row[6]) if row[6] != "ns" else 1.0
    except Exception:
        qval = 1.0
    cell = tbl[(i, 6)]
    if qval < 0.001:
        cell.set_facecolor("#1a9850")
        cell.get_text().set_color("white")
    elif qval < 0.05:
        cell.set_facecolor("#91cf60")
    else:
        cell.set_facecolor("#f0f0f0")

# bold header + alternate row shading on subsystem groups
for j in range(len(cols)):
    tbl[(0, j)].set_facecolor("#333333")
    tbl[(0, j)].get_text().set_color("white")
    tbl[(0, j)].get_text().set_weight("bold")
for i, row in enumerate(rows, start=1):
    sub_idx = (i-1) // 3
    if sub_idx % 2 == 1:
        for j in range(6):  # not q (6 is the q column)
            if tbl[(i,j)].get_facecolor() != (0.0,0.0,0.0,0.0):
                tbl[(i,j)].set_facecolor("#fafafa")

ax.set_title("Finding 1 — shared lipid activation program, 48 → 72 hpf",
             loc="left", pad=18, fontsize=20, weight="bold")
_save(fig, "1a_table")


# ============================================================
# 1b — grouped bar chart (mean log2FC per subsystem × regime at 48→72)
# ============================================================
print("[fig] 1b bars")
fig, ax = plt.subplots(figsize=(14, 7))
x = np.arange(len(F1_SUBSYSTEMS))
w = 0.26
for i, reg in enumerate(REGIMES):
    vals, qvals = [], []
    for sub in F1_SUBSYSTEMS:
        r = f1_core[(f1_core.subsystem == sub) & (f1_core.regime == reg)]
        if len(r):
            vals.append(r.iloc[0].rise_log2FC); qvals.append(r.iloc[0].enrichment_q)
        else:
            vals.append(0); qvals.append(1)
    offset = (i - 1) * w
    bars = ax.bar(x + offset, vals, w, color=REGIME_COLORS[reg], label=reg,
                  edgecolor="black", linewidth=1)
    # hatch non-significant bars
    for b, q in zip(bars, qvals):
        if q >= 0.05:
            b.set_hatch("///")
            b.set_edgecolor("#555555")
        # annotate q-value above each bar
        h = b.get_height()
        if abs(h) > 0.05:
            q_str = f"q={q:.1g}" if q < 1 else "ns"
            ax.text(b.get_x() + b.get_width()/2, h + (0.08 if h > 0 else -0.14),
                    q_str, ha="center", fontsize=10, color="#333333")

ax.axhline(0, color="black", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels([F1_SHORT[s] for s in F1_SUBSYSTEMS], rotation=0)
ax.set_ylabel("Mean centered log$_2$FC\n(72 hpf vs 48 hpf)")
ax.set_title("Mean centered log$_2$FC per subsystem × regime, 48 → 72 hpf transition",
             loc="left", pad=10)
ax.legend(title="Regime", loc="upper right", frameon=True)
# add legend element for hatching
leg2 = ax.legend(handles=[
    Patch(facecolor=REGIME_COLORS["BL"], edgecolor="black", label="BL (12 h blue, 12 h blue)"),
    Patch(facecolor=REGIME_COLORS["D"],  edgecolor="black", label="D (constant dark)"),
    Patch(facecolor=REGIME_COLORS["LD"], edgecolor="black", label="LD (12:12 light/dark)"),
    Patch(facecolor="white", edgecolor="#555555", hatch="///", label="q ≥ 0.05"),
], loc="upper right", frameon=True, fontsize=12)
ax.add_artist(leg2)
ax.grid(axis="y", linestyle=":", alpha=0.4)
_save(fig, "1b_bars")


# ============================================================
# 1c — trajectory line plot (per subsystem across 4 transitions)
# ============================================================
print("[fig] 1c trajectory")
fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=True)
axes = axes.flatten()
for idx, sub in enumerate(F1_SUBSYSTEMS):
    ax = axes[idx]
    for reg in REGIMES:
        xs, ys, qs = [], [], []
        for (a, b) in TRANSITIONS:
            r = f2[(f2.subsystem == sub) & (f2.regime == reg) &
                   (f2.hpf_a == a) & (f2.hpf_b == b)]
            if len(r):
                xs.append(TRANSITION_MID[(a,b)])
                ys.append(r.iloc[0].rise_log2FC)
                qs.append(r.iloc[0].enrichment_q)
        if xs:
            ax.plot(xs, ys, "-o", color=REGIME_COLORS[reg], label=reg,
                    linewidth=2.5, markersize=9, markeredgecolor="white", markeredgewidth=1.2)
            # bold marker if significant
            for xi, yi, qi in zip(xs, ys, qs):
                if qi < 0.05:
                    ax.plot(xi, yi, "o", color=REGIME_COLORS[reg],
                            markersize=14, markeredgecolor="black", markeredgewidth=1.5,
                            zorder=10)
    ax.axhline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.5)
    ax.axvspan(48, 72, alpha=0.08, color="#d94801")  # highlight 48-72 band
    ax.set_title(F1_SHORT[sub], loc="left", fontsize=15)
    ax.set_xticks([36, 60, 84, 108])
    ax.set_xticklabels(["24→48", "48→72", "72→96", "96→120"], fontsize=11)
    ax.grid(True, linestyle=":", alpha=0.3)
    if idx in (0, 3):
        ax.set_ylabel("Mean centered log$_2$FC\n(later vs earlier hpf)")

axes[-1].axis("off")  # unused 6th slot
# unified legend in the empty subplot slot
legend_ax = axes[-1]
legend_elements = [
    plt.Line2D([0],[0], color=REGIME_COLORS["BL"], marker="o", markersize=10, linewidth=2.5, label="BL"),
    plt.Line2D([0],[0], color=REGIME_COLORS["D"],  marker="o", markersize=10, linewidth=2.5, label="D"),
    plt.Line2D([0],[0], color=REGIME_COLORS["LD"], marker="o", markersize=10, linewidth=2.5, label="LD"),
    plt.Line2D([0],[0], color="black", marker="o", markersize=12, linewidth=0, label="q < 0.05"),
    Patch(facecolor="#d94801", alpha=0.2, label="48→72 band (Finding 1)"),
]
legend_ax.legend(handles=legend_elements, loc="center", frameon=True, fontsize=14, title="Regime / gate")
fig.suptitle("Developmental trajectories across the five Finding 1 subsystems",
             fontsize=20, weight="bold", x=0.02, ha="left", y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.97])
_save(fig, "1c_trajectory")


# ============================================================
# 1d — density ridges for Leukotriene + Arachidonic across 5 hpf × 3 regimes
# ============================================================
print("[fig] 1d density (loading pooled samples)")
ZERO_TOL = 1e-10
HPFS = (24, 48, 72, 96, 120)

# load membership (any cell will do; subsystem is S-based)
mem = pd.read_csv(SUBSETS / "BL_24_subset_membership.csv")
# per-reaction subsystem lookup (use representative column + subsystem)
rxn_to_sub = dict(zip(mem.rxn, mem.subsystem))

def _logflux_for_subsystem(cell, subsystem):
    df = pd.read_csv(POOLED / f"{cell}_pooled.csv", index_col=0)
    sub_rxns = [r for r in df.index if rxn_to_sub.get(r) == subsystem]
    if not sub_rxns:
        return np.array([])
    vals = np.abs(df.loc[sub_rxns].values).astype(float).ravel()
    vals = vals[vals > ZERO_TOL]
    return np.log2(vals)

for sub, title in [
    ("Leukotriene metabolism", "Leukotriene metabolism"),
    ("Arachidonic acid metabolism", "Arachidonic acid metabolism"),
]:
    fig, axes = plt.subplots(1, 3, figsize=(16, 7), sharey=True, sharex=True)
    for ax, reg in zip(axes, REGIMES):
        offsets = []
        for i, h in enumerate(HPFS):
            cell = f"{reg}_{h}"
            vals = _logflux_for_subsystem(cell, sub)
            if len(vals) < 3:
                continue
            kde = gaussian_kde(vals)
            xr = np.linspace(vals.min()-0.5, vals.max()+0.5, 400)
            yr = kde(xr)
            yr = yr / yr.max() * 0.8  # normalize height
            baseline = i * 1.0
            ax.fill_between(xr, baseline, baseline + yr,
                            color=plt.cm.viridis(i / (len(HPFS)-1)),
                            alpha=0.65, edgecolor="black", linewidth=0.8)
            ax.text(xr[0] - 0.3, baseline + 0.3, f"{h} hpf",
                    fontsize=12, ha="right", va="center")
        ax.set_title(reg, fontsize=16)
        ax.set_xlabel("log$_2$ |flux|  (raw)")
        ax.spines["left"].set_visible(False)
        ax.yaxis.set_visible(False)
        ax.axvspan(-30, 30, alpha=0)  # ensure spans draw

    fig.suptitle(f"{title} — density of sampled flux across stages, by regime",
                 fontsize=20, weight="bold", x=0.02, ha="left", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    name = "1d_density_leukotriene" if "Leuk" in sub else "1d_density_arachidonic"
    _save(fig, name)

print("\n[done]  outputs in", OUT)
