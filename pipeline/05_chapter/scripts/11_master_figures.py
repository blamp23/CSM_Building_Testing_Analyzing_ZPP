"""
11_master_figures.py — Figure A, Figure B, Figure B′ and their supplementary masters.

Design principle (per editorial brief, 2026-10-27):
  Build the MASTER first (every contrast / every subsystem) and derive the
  chapter figure as a selection.  SB → B, SB' → B', SA → A1 → A.  Two figures
  can therefore never disagree.

Outputs under chapter_assets/ppt_figures/:
  A/
    SA_contrast_class_bars.{pdf,png}           — 27-panel master, one per contrast
    A_where_models_differ.{pdf,png}            — 2-panel chapter figure (pooled bars + PCA)
    class_map.csv                              — supplementary table: subsystem → class
  B/
    SB_all_subsystems_heatmap.{pdf,png}        — master heatmap, every qualifying subsystem
    B_findings_heatmap.{pdf,png}               — 14 chapter rows
  B_trajectories/
    SB_prime_all_subsystems_lines.{pdf,png}    — master small-multiples, every subsystem
    B_prime_findings_lines.{pdf,png}           — 14 chapter rows
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyArrow
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import COMPARISON_TABLES, FIGURES, CELLS

# For heatmap display: group columns by stage (BL/D/LD) within each hpf block.
# config.CELLS is grouped by regime, which is the wrong grouping for Figure B.
CELLS_DISPLAY = tuple(f"{c}_{h}" for h in (24, 48, 72, 96, 120) for c in ("BL", "D", "LD"))

OUT_A = FIGURES.parent / "ppt_figures" / "A"
OUT_B = FIGURES.parent / "ppt_figures" / "B"
OUT_BP = FIGURES.parent / "ppt_figures" / "B_trajectories"
for d in (OUT_A, OUT_B, OUT_BP):
    d.mkdir(parents=True, exist_ok=True)

# ---------- constants ----------
REGIMES = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
REGIME_COLORS = {"BL": "#1976D2", "D": "#333333", "LD": "#F6B52C"}
# PCA uses shape for timepoint, color for condition (per editorial brief).
HPF_MARKERS = {24: "o", 48: "X", 72: "P", 96: "s", 120: "D"}

# 14 chapter rows in the order they are discussed in §3.1
CHAPTER_14 = [
    "Leukotriene metabolism",
    "Arachidonic acid metabolism",
    "Eicosanoid metabolism",
    "Fatty acid oxidation",
    "Fatty acid activation (cytosolic)",
    "Beta oxidation of even-chain fatty acids (peroxisomal)",
    "Beta oxidation of unsaturated fatty acids (n-9) (peroxisomal)",
    "Fatty acid desaturation (even-chain)",
    "Omega-3 fatty acid metabolism",
    "Bile acid biosynthesis",
    "Steroid metabolism",
    "Serotonin and melatonin biosynthesis",
    "Phenylalanine, tyrosine and tryptophan biosynthesis",
    "Tyrosine metabolism",
]
# thread brackets for the right-margin label on Figure B
CHAPTER_THREADS = {
    "Lipid backbone":   ["Leukotriene metabolism", "Arachidonic acid metabolism",
                         "Eicosanoid metabolism", "Fatty acid oxidation",
                         "Fatty acid activation (cytosolic)"],
    "Peroxisomal":      ["Beta oxidation of even-chain fatty acids (peroxisomal)",
                         "Beta oxidation of unsaturated fatty acids (n-9) (peroxisomal)"],
    "Eicosanoid stall": ["Fatty acid desaturation (even-chain)",
                         "Omega-3 fatty acid metabolism"],
    "Sterol & bile":    ["Bile acid biosynthesis", "Steroid metabolism"],
    "Signaling / AA":   ["Serotonin and melatonin biosynthesis",
                         "Phenylalanine, tyrosine and tryptophan biosynthesis",
                         "Tyrosine metabolism"],
}

# ---------- class map (subsystem → metabolic class) ----------
# Classes used by Figure A stacked bars. Keyword lookup, deterministic.
CLASS_COLORS = {
    "Lipid":              "#2CA02C",   # green
    "Amino acid":         "#FF7F0E",   # orange
    "Sterol / bile":      "#D62728",   # red
    "Nucleotide":         "#9467BD",   # purple
    "Carbohydrate":       "#8C564B",   # brown
    "Cofactor / vitamin": "#E377C2",   # pink
    "Transport / other":  "#7F7F7F",   # grey
}
CLASS_ORDER = list(CLASS_COLORS.keys())


def classify(name):
    if not isinstance(name, str): return "Transport / other"
    n = name.lower()
    # sterol/bile first — "Formation and hydrolysis of cholesterol esters" is sterol, not lipid
    if any(k in n for k in ["bile acid", "steroid", "cholesterol", "androgen",
                             "estrogen", "hormone biosynthesis", "vitamin d",
                             "terpenoid"]):
        return "Sterol / bile"
    if any(k in n for k in ["fatty acid", "beta oxidation", "carnitine shuttle",
                             "glycerolipid", "glycerophospholipid", "sphingolipid",
                             "glycosphingolipid", "ether lipid", "acyl",
                             "linoleate", "leukotriene", "arachidonic",
                             "eicosanoid", "prostaglandin", "omega-", "lipoic",
                             "acylglyceride", "phosphatidylinositol phosphate"]):
        return "Lipid"
    if any(k in n for k in ["amino acid", "aminoacyl-trna", "alanine",
                             "arginine", "aspartate", "cysteine", "glycine",
                             "histidine", "lysine", "methionine", "phenyl",
                             "serine", "threonine", "tryptophan", "tyrosine",
                             "valine", "beta-alanine", "metabolism of other amino acids",
                             "protein assembly", "protein degradation",
                             "protein modification"]):
        return "Amino acid"
    if any(k in n for k in ["purine", "pyrimidine", "nucleotide"]):
        return "Nucleotide"
    if any(k in n for k in ["glycolysis", "gluconeogenesis", "tricarboxylic",
                             "oxidative phosphorylation", "pentose", "pyruvate",
                             "fructose", "galactose", "starch", "sucrose",
                             "mannose", "inositol phosphate", "propanoate",
                             "butanoate", "c5-branched", "amino sugar",
                             "glycan", "chondroitin", "heparan", "keratan",
                             "ascorbate", "blood group"]):
        return "Carbohydrate"
    if any(k in n for k in ["folate", "glutathione", "nicotinate", "pantothenate",
                             "riboflavin", "thiamine", "ubiquinone", "vitamin a",
                             "vitamin b", "vitamin e", "biopterin", "porphyrin",
                             "heme", "ros detox", "retinol",
                             "serotonin", "sulfur", "drug", "xenobiotics"]):
        return "Cofactor / vitamin"
    return "Transport / other"


# ---------- PPT style ----------
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11,
    "axes.titlesize": 14, "axes.labelsize": 12,
    "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})


def _save(fig, out_dir, name):
    fig.savefig(out_dir / f"{name}.pdf")
    fig.savefig(out_dir / f"{name}.png", dpi=200)
    plt.close(fig)
    print(f"  wrote {out_dir.name}/{name}.{{pdf,png}}")


# ============================================================
# LOAD everything once
# ============================================================
print("[load] per_reaction_flux_F1_centered")
f1 = pd.read_csv(COMPARISON_TABLES / "per_reaction_flux_F1_centered.csv")
print("[load] per_reaction_flux_F2_centered")
f2 = pd.read_csv(COMPARISON_TABLES / "per_reaction_flux_F2_centered.csv")
print("[load] subsystem_summary_F1_centered")
s1 = pd.read_csv(COMPARISON_TABLES / "subsystem_summary_F1_centered.csv")
print("[load] subsystem_summary_F2_centered")
s2 = pd.read_csv(COMPARISON_TABLES / "subsystem_summary_F2_centered.csv")
print("[load] reaction_inventory_F1")
inv = pd.read_csv(COMPARISON_TABLES / "reaction_inventory_F1.csv")
print("[load] pca_cells_centered")
pca_df = pd.read_csv(COMPARISON_TABLES / "pca_cells_centered.csv").rename(columns={"Unnamed: 0": "cell"})

# apply class map to the summary table
s1["class"] = s1.subsystem.map(classify)
s2["class"] = s2.subsystem.map(classify)

# write the class map as a supplementary table (one row per subsystem)
class_map = pd.DataFrame({"subsystem": sorted(s1.subsystem.dropna().unique())})
class_map["class"] = class_map.subsystem.map(classify)
class_map.to_csv(OUT_A / "class_map.csv", index=False)
print(f"[load] wrote class_map.csv ({len(class_map)} subsystems)")


# ============================================================
# FIGURE A — Where the models differ
# ============================================================

# ---- A1 helper: per-contrast notable counts by class ----
def notable_counts_by_class(df_rxn, contrast):
    sub = df_rxn[df_rxn.contrast == contrast].copy()
    sub["class"] = sub.subsystem.map(classify)
    notable = sub[sub.notable_centered & ~sub.loop_flagged]
    counts = notable.groupby("class").size().reindex(CLASS_ORDER, fill_value=0)
    total_tested = ((sub.carry_both) & (~sub.loop_flagged)).sum()
    return counts, total_tested


def _contrast_display(c):
    """BL_24_vs_D_24 → 'BL vs D @ 24' ; BL_24_vs_BL_48 → 'BL  24→48'."""
    m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_(\d+)", c)
    if not m: return c
    ra, ha, rb, hb = m.groups()
    if ha == hb:          # F1
        return f"{ra} vs {rb}  @ {ha} hpf"
    if ra == rb:          # F2
        return f"{ra}   {ha}→{hb} hpf"
    return c


# ---- SA: 27-panel master ----
print("\n[SA] building 27-contrast master")
# order F1 contrasts by hpf then regime pair (24 BL-D, 24 BL-LD, 24 D-LD, 48 ..., ...)
def _f1_sort_key(c):
    m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_\d+", c)
    reg_order = {"BL": 0, "D": 1, "LD": 2}
    return (int(m.group(2)), reg_order[m.group(1)], reg_order[m.group(3)]) if m else (999, 9, 9)
def _f2_sort_key(c):
    m = re.match(r"([A-Z]+)_(\d+)_vs_[A-Z]+_(\d+)", c)
    reg_order = {"BL": 0, "D": 1, "LD": 2}
    return (int(m.group(2)), reg_order[m.group(1)]) if m else (999, 9)
all_contrasts_f1 = sorted(s1.contrast.unique(), key=_f1_sort_key)  # 15
all_contrasts_f2 = sorted(s2.contrast.unique(), key=_f2_sort_key)  # 12

ncols, nrows = 5, 6  # 30 slots, use first 27
fig_sa, axes_sa = plt.subplots(nrows, ncols, figsize=(18, 18), sharey=False)
axes_sa = axes_sa.flatten()

all_contrasts = [(c, "F1", f1) for c in all_contrasts_f1] + \
                [(c, "F2", f2) for c in all_contrasts_f2]

sa_all_counts = []  # for pooling into A1 later
for i, (c, fam, src_df) in enumerate(all_contrasts):
    ax = axes_sa[i]
    counts, total = notable_counts_by_class(src_df, c)
    pct = 100 * counts.sum() / total if total else 0
    bottom = 0
    for cls in CLASS_ORDER:
        ax.bar([0], counts[cls], bottom=bottom, color=CLASS_COLORS[cls],
               edgecolor="white", linewidth=0.5)
        bottom += counts[cls]
    ax.set_xticks([])
    ax.set_title(_contrast_display(c), fontsize=9, loc="center")
    ax.text(0, bottom + max(1, bottom*0.02), f"n={counts.sum()} ({pct:.1f}%)",
            ha="center", fontsize=8, color="#444444")
    ax.set_ylim(0, max(bottom*1.15, 10))
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="y", labelsize=8)
    sa_all_counts.append({"contrast": c, "family": fam, "total_tested": total,
                           **{cls: counts[cls] for cls in CLASS_ORDER}})

# hide unused axes
for j in range(len(all_contrasts), nrows*ncols):
    axes_sa[j].axis("off")

# legend
leg_patches = [Patch(facecolor=CLASS_COLORS[c], edgecolor="white", label=c)
               for c in CLASS_ORDER]
fig_sa.legend(handles=leg_patches, loc="lower center", ncol=len(CLASS_ORDER),
               frameon=False, fontsize=11, bbox_to_anchor=(0.5, -0.01))
fig_sa.suptitle("Supplementary Figure SA — Notable reactions per contrast, by metabolic class",
                 fontsize=16, weight="bold", y=1.00)
fig_sa.text(0.5, 0.97,
             "Each panel: one contrast. Stack = count of representatives with |centered log2FC| ≥ 1 (notable), carried by both cells, loop-flagged excluded.",
             ha="center", fontsize=10, color="#555555")
fig_sa.tight_layout(rect=[0, 0.02, 1, 0.96])
_save(fig_sa, OUT_A, "SA_contrast_class_bars")

# save counts CSV so A1 reads from SA
sa_counts_df = pd.DataFrame(sa_all_counts)
sa_counts_df.to_csv(OUT_A / "SA_counts.csv", index=False)


# ---- A1: pooled version (sums over the 3 contrasts per stage/transition) ----
print("[A1] pooling SA counts into chapter view")
# F1 stages: 24, 48, 72, 96, 120 — pool the 3 contrasts at each
# F2 transitions: 24-48, 48-72, 72-96, 96-120 — pool the 3 regimes at each
def stage_from_f1_contrast(c):
    m = re.match(r"[A-Z]+_(\d+)_vs_[A-Z]+_\d+", c); return int(m.group(1)) if m else None

def trans_from_f2_contrast(c):
    m = re.match(r"[A-Z]+_(\d+)_vs_[A-Z]+_(\d+)", c)
    return (int(m.group(1)), int(m.group(2))) if m else None

f1_rows = sa_counts_df[sa_counts_df.family == "F1"].copy()
f1_rows["stage"] = f1_rows.contrast.map(stage_from_f1_contrast)
f1_pooled = f1_rows.groupby("stage")[CLASS_ORDER].sum()

f2_rows = sa_counts_df[sa_counts_df.family == "F2"].copy()
f2_rows["trans"] = f2_rows.contrast.map(trans_from_f2_contrast)
f2_pooled = f2_rows.groupby("trans")[CLASS_ORDER].sum()
f2_pooled.index = [f"{a}→{b}" for (a, b) in f2_pooled.index]

# ---- A2 PCA ----
# already loaded as pca_df (15 rows)

# compose Figure A (2 panels)
print("[A] composing chapter 2-panel figure")
fig_a = plt.figure(figsize=(18, 7))
gs = fig_a.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.25)
ax_a1 = fig_a.add_subplot(gs[0, 0])
ax_a2 = fig_a.add_subplot(gs[0, 1])

# A1 bars: F1 stages + F2 transitions side-by-side
f1_x = np.arange(len(f1_pooled))
f2_x = np.arange(len(f2_pooled)) + len(f1_pooled) + 1  # gap
all_x = np.concatenate([f1_x, f2_x])
all_labels = [f"{h} hpf" for h in f1_pooled.index] + list(f2_pooled.index)
pool_df = pd.concat([f1_pooled, f2_pooled]).reset_index(drop=True)

bottom = np.zeros(len(pool_df))
for cls in CLASS_ORDER:
    heights = pool_df[cls].values
    ax_a1.bar(all_x, heights, bottom=bottom, color=CLASS_COLORS[cls],
              edgecolor="white", linewidth=0.5, label=cls, width=0.75)
    bottom += heights
for xi, total in zip(all_x, pool_df.sum(axis=1)):
    ax_a1.text(xi, total + max(bottom)*0.015, f"{int(total)}",
               ha="center", fontsize=10, color="#333333")
ax_a1.set_xticks(all_x)
ax_a1.set_xticklabels(all_labels, rotation=0, fontsize=11)
# vertical divider between F1 and F2 blocks
ax_a1.axvline((f1_x[-1] + f2_x[0]) / 2, color="#aaaaaa", linestyle="--", alpha=0.6)
ax_a1.text(f1_x.mean(), -max(bottom)*0.10, "F1 — cross-regime at each stage\n(pooled over 3 contrasts)",
           ha="center", fontsize=11, color="#333333")
ax_a1.text(f2_x.mean(), -max(bottom)*0.10, "F2 — adjacent-stage within regime\n(pooled over 3 regimes)",
           ha="center", fontsize=11, color="#333333")
ax_a1.set_ylabel("Count of notable representatives\n(|centered log$_2$FC| ≥ 1)")
ax_a1.set_title("A1  Where the differences live",
                 loc="left", fontsize=15, pad=12, weight="bold")
ax_a1.legend(loc="upper left", ncol=1, fontsize=10, frameon=True, bbox_to_anchor=(0.01, 0.98))
ax_a1.grid(axis="y", linestyle=":", alpha=0.3)
ax_a1.set_ylim(bottom=0)

# A2 PCA — shape = timepoint, color = condition; manual label offsets to avoid clutter
PCA_LABEL_OFFSETS = {  # (dx, dy) in points, per cell
    "BL_24": (-18, -8), "D_24": (8, -14), "LD_24": (-14, 10),
    "BL_48": (-18, 8),  "D_48": (10, 8),   "LD_48": (-18, -14),
    "BL_72": (-20, -4), "D_72": (10, 10),  "LD_72": (10, -14),
    "BL_96": (10, -4),  "D_96": (-22, 8),  "LD_96": (10, 10),
    "BL_120": (-24, 2), "D_120": (10, -14),"LD_120": (10, 8),
}
for _, r in pca_df.iterrows():
    ax_a2.scatter(r.PC1, r.PC2, s=200, c=REGIME_COLORS[r.condition],
                   marker=HPF_MARKERS[r.hpf],
                   alpha=0.95,
                   edgecolor="black", linewidth=1.2, zorder=3)
    off = PCA_LABEL_OFFSETS.get(f"{r.condition}_{r.hpf}", (6, 6))
    ax_a2.annotate(f"{r.condition}_{r.hpf}", (r.PC1, r.PC2),
                    xytext=off, textcoords="offset points", fontsize=9)

ax_a2.axhline(0, color="black", linewidth=0.4, alpha=0.5)
ax_a2.axvline(0, color="black", linewidth=0.4, alpha=0.5)
ax_a2.set_xlabel("PC1")
ax_a2.set_ylabel("PC2")
ax_a2.set_title("A2  Centered PCA (per-cell median flux shift removed)",
                 loc="left", fontsize=15, pad=12, weight="bold")

# A2 legend: color = condition, shape = timepoint. Placed outside the axes.
leg_condition = [Line2D([0],[0], marker="o", color="w",
                         markerfacecolor=REGIME_COLORS[r], markersize=12,
                         label=r, markeredgecolor="black") for r in REGIMES]
leg_hpf = [Line2D([0],[0], marker=HPF_MARKERS[h], color="w",
                   markerfacecolor="#888888", markersize=12,
                   label=f"{h} hpf", markeredgecolor="black") for h in HPFS]
leg1 = ax_a2.legend(handles=leg_condition, loc="upper left",
                     bbox_to_anchor=(1.02, 1.0), fontsize=10,
                     frameon=True, title="Condition (colour)", title_fontsize=11)
leg2 = ax_a2.legend(handles=leg_hpf, loc="upper left",
                     bbox_to_anchor=(1.02, 0.70), fontsize=10,
                     frameon=True, title="Timepoint (shape)", title_fontsize=11)
ax_a2.add_artist(leg1)
ax_a2.grid(True, linestyle=":", alpha=0.3)

fig_a.suptitle("Figure A — Where the models differ",
                fontsize=20, weight="bold", y=1.03, x=0.01, ha="left")
_save(fig_a, OUT_A, "A_where_models_differ")


# ============================================================
# FIGURE B — findings heatmap (chapter) + SB master
# ============================================================
print("\n[B] computing per-cell per-subsystem centered log flux matrix")

# Melt per_reaction_flux into (cell, rxn, subsystem, mean_nz, median_log2)
# Each (cell, rxn) appears in 2 F1 contrasts; take first.
a = f1[["cell_A", "rxn_id", "subsystem", "mean_nz_A", "median_log2_A", "loop_flagged"]].rename(
    columns={"cell_A": "cell", "mean_nz_A": "mean_nz", "median_log2_A": "median_log2"})
b = f1[["cell_B", "rxn_id", "subsystem", "mean_nz_B", "median_log2_B", "loop_flagged"]].rename(
    columns={"cell_B": "cell", "mean_nz_B": "mean_nz", "median_log2_B": "median_log2"})
long = pd.concat([a, b]).drop_duplicates(subset=["cell", "rxn_id"]).reset_index(drop=True)
# drop loop-flagged, drop zero-flux
long = long[(~long.loop_flagged) & (long.mean_nz > 1e-10)].copy()
long["centered_log2"] = np.log2(long["mean_nz"]) - long["median_log2"]

# subsystem × cell matrix of mean centered log flux (stage-grouped column order)
mat = (long.groupby(["subsystem", "cell"])
            .centered_log2.mean()
            .unstack("cell")
            .reindex(columns=list(CELLS_DISPLAY)))
print(f"[B] subsystem×cell matrix: {mat.shape}")

# count representatives per cell per subsystem (for the "absent → grey" mask)
n_reps = (long.groupby(["subsystem", "cell"])
               .rxn_id.nunique()
               .unstack("cell")
               .reindex(columns=list(CELLS_DISPLAY))
               .fillna(0))

# Qualifying filter: subsystem must have ≥ 5 reps in at least one cell
qual_subs = n_reps[(n_reps >= 5).any(axis=1)].index.tolist()
print(f"[B] subsystems qualifying (≥5 reps in ≥1 cell): {len(qual_subs)}")

mat_q = mat.loc[qual_subs].copy()
n_reps_q = n_reps.loc[qual_subs]

# row z-score (ignoring NaNs). if row has zero variance → leave as-is.
row_mean = mat_q.mean(axis=1)
row_std = mat_q.std(axis=1).replace(0, np.nan)
mat_z = mat_q.sub(row_mean, axis=0).div(row_std, axis=0)

# order rows by class, then by hpf of peak abs(z)
peak_hpf = mat_z.abs().idxmax(axis=1).map(lambda c: int(c.split("_")[1]) if isinstance(c, str) else 999)
peak_class = pd.Series({s: classify(s) for s in mat_z.index})
class_rank = {c: i for i, c in enumerate(CLASS_ORDER)}
order_df = pd.DataFrame({"subsystem": mat_z.index,
                          "class": peak_class.values,
                          "class_rank": peak_class.map(class_rank).values,
                          "peak_hpf": peak_hpf.values})
order_df = order_df.sort_values(["class_rank", "peak_hpf", "subsystem"])
mat_z = mat_z.loc[order_df.subsystem]
n_reps_q = n_reps_q.loc[order_df.subsystem]


def _draw_heatmap(ax, data, n_reps_row_subset, row_labels_right=None,
                   title=None, show_class_bracket=True, bold_rows=None,
                   fontsize_rows=10):
    """Draw a heatmap with grey-for-absent cells and class brackets."""
    # mask: absent cells get grey (no value)
    absent_mask = (n_reps_row_subset == 0).values
    display = data.copy().values
    display[absent_mask] = np.nan

    vmax = np.nanmax(np.abs(display))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.cm.RdBu_r
    cmap.set_bad("#dddddd")

    im = ax.imshow(display, aspect="auto", cmap=cmap, norm=norm,
                     interpolation="nearest")

    # x ticks: stage blocks with regime columns (data columns are already stage-grouped)
    ncols = len(data.columns)
    ax.set_xticks(range(ncols))
    labels = []
    for c in data.columns:
        reg, hpf = c.split("_")
        labels.append(f"{reg}\n{hpf}")
    ax.set_xticklabels(labels, fontsize=9)
    # stage block separators (every 3 cols — one per hpf block)
    for i in range(1, len(HPFS)):
        ax.axvline(i*3 - 0.5, color="black", linewidth=1.2)
    # within-stage separators
    for i in range(ncols):
        if i % 3 != 0:
            ax.axvline(i - 0.5, color="#dddddd", linewidth=0.5)

    # y ticks (subsystem names)
    ax.set_yticks(range(len(data)))
    yticklabels = []
    for s in data.index:
        n_max = int(n_reps_row_subset.loc[s].max())
        yticklabels.append(f"{s}  (n={n_max})")
    ax.set_yticklabels(yticklabels, fontsize=fontsize_rows)
    # emphasise chapter rows if bold_rows given
    if bold_rows:
        for lab in ax.get_yticklabels():
            txt = lab.get_text().split("  (n=")[0]
            if txt in bold_rows:
                lab.set_weight("bold")

    ax.set_xlabel("Cell  (regime over hpf)")
    if title:
        ax.set_title(title, loc="left", fontsize=15, weight="bold", pad=10)
    cbar = plt.colorbar(im, ax=ax, pad=0.015, fraction=0.025, shrink=0.9)
    cbar.set_label("Row z-score of mean centered log$_2$|flux|", fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    return im


def _overlay_markers(ax, data, s1_df, s2_df):
    """Dots for F1 q<0.05 between regime cols; arrows for F2 q<0.05 between stage blocks."""
    cell_to_idx = {c: i for i, c in enumerate(data.columns)}
    row_to_idx = {s: i for i, s in enumerate(data.index)}

    # F1: dot between the two regime cols at that stage, in that row
    sig_f1 = s1_df[(s1_df.enrichment_q < 0.05) & (s1_df.n_reps_tested >= 5)]
    for _, r in sig_f1.iterrows():
        if r.subsystem not in row_to_idx: continue
        y = row_to_idx[r.subsystem]
        if r.cell_A not in cell_to_idx or r.cell_B not in cell_to_idx: continue
        xa = cell_to_idx[r.cell_A]; xb = cell_to_idx[r.cell_B]
        xm = (xa + xb) / 2
        ax.plot(xm, y, "o", color="black", markersize=4, markeredgecolor="white",
                markeredgewidth=0.8, zorder=10)

    # F2: arrow at the stage boundary, in the regime's row position
    sig_f2 = s2_df[(s2_df.enrichment_q < 0.05) & (s2_df.n_reps_tested >= 5)]
    for _, r in sig_f2.iterrows():
        if r.subsystem not in row_to_idx: continue
        y = row_to_idx[r.subsystem]
        if r.cell_A not in cell_to_idx or r.cell_B not in cell_to_idx: continue
        xa = cell_to_idx[r.cell_A]; xb = cell_to_idx[r.cell_B]
        # arrow direction: horizontal, pointing from earlier to later stage block
        xm = (xa + xb) / 2
        dx = (xb - xa) * 0.18
        ax.annotate("", xy=(xm + dx, y - 0.25), xytext=(xm - dx, y - 0.25),
                    arrowprops=dict(arrowstyle="->", color="black",
                                     lw=0.8, shrinkA=0, shrinkB=0), zorder=10)


# ---- SB master ----
print(f"[SB] building master heatmap ({len(mat_z)} subsystems)")
fig_sb, ax_sb = plt.subplots(figsize=(14, max(10, len(mat_z) * 0.14)))
_draw_heatmap(ax_sb, mat_z, n_reps_q, bold_rows=set(CHAPTER_14),
               title="Supplementary Figure SB — every qualifying subsystem (≥5 reps in ≥1 cell)",
               fontsize_rows=7.5)
_overlay_markers(ax_sb, mat_z, s1, s2)
# legend
ax_sb.text(0.01, 1.02,
            "F1 q<0.05 (dot between the two regime columns)   ·   F2 q<0.05 (arrow between stage blocks, in the regime's row)   ·   Chapter rows in bold   ·   Grey = subsystem absent from that cell",
            transform=ax_sb.transAxes, fontsize=8, color="#555555")
fig_sb.tight_layout()
_save(fig_sb, OUT_B, "SB_all_subsystems_heatmap")


# ---- B chapter: select the 14 rows ----
# Chapter rows are always included even if they don't meet the SB ≥5-reps gate
# (Serotonin and melatonin biosynthesis has 4 reps — explicit chapter exception).
print("[B] selecting 14 chapter rows (always-included for chapter)")
chapter_rows = []
for s in CHAPTER_14:
    if s not in mat.index:
        print(f"  WARN: chapter subsystem not in data at all: {s}")
        continue
    chapter_rows.append(s)
mat_b_raw = mat.loc[chapter_rows]
# z-score within each row across the 15 cells
rm = mat_b_raw.mean(axis=1); rs = mat_b_raw.std(axis=1).replace(0, np.nan)
mat_b = mat_b_raw.sub(rm, axis=0).div(rs, axis=0)
n_reps_b = n_reps.loc[chapter_rows]

fig_b, ax_b = plt.subplots(figsize=(14, 7))
_draw_heatmap(ax_b, mat_b, n_reps_b,
               title="Figure B — findings heatmap (14 chapter subsystems)",
               fontsize_rows=11)
_overlay_markers(ax_b, mat_b, s1, s2)

# thread brackets on the right (positioned outside the colorbar at x=1.14+)
row_idx = {s: i for i, s in enumerate(mat_b.index)}
for thread, members in CHAPTER_THREADS.items():
    ids = [row_idx[s] for s in members if s in row_idx]
    if not ids: continue
    y0, y1 = min(ids) - 0.4, max(ids) + 0.4
    ax_b.annotate("", xy=(1.14, y1), xytext=(1.14, y0),
                   xycoords=("axes fraction", "data"),
                   arrowprops=dict(arrowstyle="-", color="black", lw=2),
                   annotation_clip=False)
    ax_b.text(1.16, (y0 + y1)/2, thread,
               transform=ax_b.get_yaxis_transform(),
               fontsize=10, va="center", rotation=0, clip_on=False)

ax_b.text(0.01, 1.02,
          "F1 q<0.05 → dot between regime columns   ·   F2 q<0.05 → arrow at stage boundary   ·   Grey = subsystem absent   ·   Values are row-z-scored",
          transform=ax_b.transAxes, fontsize=9, color="#555555")
fig_b.tight_layout()
_save(fig_b, OUT_B, "B_findings_heatmap")


# ============================================================
# FIGURE B' — trajectories (chapter) + SB' master
# ============================================================
def _trajectory_panel(ax, subsystem, mat_values_raw, n_reps_subs):
    """One subsystem; three regime lines across 5 hpf; raw centered-log flux (not z-scored)."""
    for reg in REGIMES:
        xs, ys, present = [], [], []
        for h in HPFS:
            cell = f"{reg}_{h}"
            v = mat_values_raw.get(cell)
            n_r = n_reps_subs.get(cell, 0)
            if n_r == 0 or pd.isna(v):
                xs.append(h); ys.append(np.nan); present.append(False)
            else:
                xs.append(h); ys.append(v); present.append(True)
        # solid line where present
        ax.plot(xs, ys, "-o", color=REGIME_COLORS[reg], label=reg,
                 linewidth=1.8, markersize=6, markeredgecolor="white",
                 markeredgewidth=0.8)
        # open circles for absent
        for xi, yi, p in zip(xs, ys, present):
            if not p:
                ax.plot(xi, 0, "o", color=REGIME_COLORS[reg], mfc="white",
                         markersize=5, markeredgewidth=1.5)

    ax.axhline(0, color="black", linewidth=0.4, alpha=0.5, linestyle="--")
    ax.set_xticks(list(HPFS))
    ax.set_xticklabels([str(h) for h in HPFS], fontsize=9)
    n_max = int(n_reps_subs.max())
    ax.set_title(f"{subsystem}\n(n={n_max})", fontsize=9, loc="left")
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(True, linestyle=":", alpha=0.3)


# ---- B' chapter (14 panels, 3 cols × 5 rows) ----
print("[B'] building chapter trajectories (14 panels)")
rows_bp = 5; cols_bp = 3
fig_bp, axes_bp = plt.subplots(rows_bp, cols_bp, figsize=(14, 16),
                                 sharex=True)
axes_bp = axes_bp.flatten()
for i, sub in enumerate(CHAPTER_14):
    ax = axes_bp[i]
    if sub in mat.index:  # use unfiltered matrix so melatonin (n=4) is included
        _trajectory_panel(ax, sub, mat.loc[sub], n_reps.loc[sub])
    else:
        ax.axis("off")
# hide unused
for j in range(len(CHAPTER_14), rows_bp*cols_bp):
    axes_bp[j].axis("off")
# shared y-label on far left column
for i in range(0, rows_bp*cols_bp, cols_bp):
    axes_bp[i].set_ylabel("centered log$_2$|flux|", fontsize=10)
for i in range((rows_bp-1)*cols_bp, rows_bp*cols_bp):
    axes_bp[i].set_xlabel("hpf")

# legend in bottom-right empty slot (if any)
if len(CHAPTER_14) < rows_bp*cols_bp:
    axes_bp[-1].axis("off")
    axes_bp[-1].legend(handles=[
        Line2D([0],[0], color=REGIME_COLORS[r], marker="o", linewidth=1.8,
                markersize=8, label=r) for r in REGIMES
    ] + [Line2D([0],[0], color="black", marker="o", mfc="white",
                 markersize=8, linewidth=0, label="subsystem absent")],
         loc="center", fontsize=11, frameon=True, title="Regime")

fig_bp.suptitle("Figure B′ — findings trajectories (14 chapter subsystems)",
                 fontsize=16, weight="bold", y=1.00)
fig_bp.tight_layout(rect=[0, 0, 1, 0.98])
_save(fig_bp, OUT_BP, "B_prime_findings_lines")


# ---- SB' master (all qualifying subsystems) ----
print(f"[SB'] building master trajectories ({len(qual_subs)} subsystems)")
subs_sorted = order_df.subsystem.tolist()
n_sub = len(subs_sorted)
ncols = 5
nrows = (n_sub + ncols - 1) // ncols
fig_sbp, axes_sbp = plt.subplots(nrows, ncols, figsize=(18, nrows * 2.0),
                                   sharex=True)
axes_sbp = np.array(axes_sbp).flatten()
for i, sub in enumerate(subs_sorted):
    ax = axes_sbp[i]
    _trajectory_panel(ax, sub, mat_q.loc[sub], n_reps.loc[sub])
for j in range(n_sub, len(axes_sbp)):
    axes_sbp[j].axis("off")
fig_sbp.suptitle("Supplementary Figure SB′ — trajectories for every qualifying subsystem",
                  fontsize=18, weight="bold", y=1.00)
fig_sbp.tight_layout(rect=[0, 0, 1, 0.99])
_save(fig_sbp, OUT_BP, "SB_prime_all_subsystems_lines")


print("\n[done]")
print(f"  Figure A outputs:          {OUT_A}")
print(f"  Figure B outputs:          {OUT_B}")
print(f"  Figure B' outputs:         {OUT_BP}")
