"""
s180 -- Figures 3, 4, 5, 6 for the chapter.

Fig 3 — orientation heatmap: F1 between-condition + F2 adjacent-hpf,
       cells coloured by signed W1, Tier A filled, Tier B hatched.
Fig 4 — developmental lipid program (§4.2): 4 subsystems × 5 hpf × 3 regimes,
       ridgeline with a zero-fraction bar above each ridge.
Fig 5 — amino-acid separation at 72–96 hpf (§4.3): 3 subsystems × 2 ridges.
Fig 6 — convergence at 120 hpf (§4.3): 4 density panels + panel E dot plot
       of signed W1 for BL-vs-LD and D-vs-LD.

Pooled flux samples from grace_v4/results/pooled_v7/{cell}_pooled.csv.
Subsystem membership from grace_v4/results/subsets_v7/{cell}_subset_membership.csv.
Tier A / B tables from v7_final/.

All densities use log2(|v|/μ) with μ = pool-mean of |v| over non-zero flux
(method of record, §11.14). Zero fraction plotted as a thin bar above each
ridge so a reader sees the two-tier test's zero contribution.
"""
import sys, warnings
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
OUT = R / "chapter_figures"
OUT.mkdir(parents=True, exist_ok=True)

COND_COLOR = {"BL": "#1F77B4", "D": "#111111", "LD": "#F0A800"}
CONDS = ("BL", "D", "LD")
HPFS = [24, 48, 72, 96, 120]
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

tierA = pd.read_csv(R / "v7_final/tierA_final.csv")
tierB = pd.read_csv(R / "v7_final/tierB_final.csv")


# ============================================================
# Flux data helpers
# ============================================================
_pool_cache = {}
_sub_cache = {}

def _pooled(cell):
    if cell not in _pool_cache:
        _pool_cache[cell] = pd.read_csv(R / f"pooled_v7/{cell}_pooled.csv").set_index("rxn")
    return _pool_cache[cell]

def _subsets(cell):
    if cell not in _sub_cache:
        _sub_cache[cell] = pd.read_csv(R / f"subsets_v7/{cell}_subset_membership.csv")
    return _sub_cache[cell]


def _flux_pool(cell, subsystem, max_points=12000, seed=0):
    """Return (log2 normalised flux array of non-zeros, zero_fraction).
    Follows the method of record: μ = pool mean of |v| across non-zeros.
    log2(|v|/μ) is undefined at zero, so zeros are reported as a separate
    fraction (two-tier test logic).
    """
    sub = _subsets(cell)
    reps = sub.loc[sub.subsystem == subsystem, "representative"].unique()
    if len(reps) == 0:
        return np.array([]), 0.0
    pool = _pooled(cell)
    rep_in = [r for r in reps if r in pool.index]
    if not rep_in:
        return np.array([]), 0.0
    flat = pool.loc[rep_in].values.flatten()
    n_total = flat.size
    nonzero = flat[flat != 0]
    zero_frac = 1 - nonzero.size / n_total if n_total else 0.0
    if nonzero.size == 0:
        return np.array([]), zero_frac
    mu = np.abs(nonzero).mean()
    log2v = np.log2(np.abs(nonzero) / mu)
    rng = np.random.default_rng(seed)
    if log2v.size > max_points:
        log2v = rng.choice(log2v, max_points, replace=False)
    return log2v, zero_frac


# ============================================================
# FIGURE 3 — orientation heatmap
# ============================================================
# F1 contrasts: {BL vs D, BL vs LD, D vs LD} × 5 hpf = 15 columns.
# F2 contrasts: {24→48, 48→72, 72→96, 96→120} × 3 regimes = 12 columns.
# Rows = subsystems that appear in any Tier A or B of either family.

F1_CONTRASTS = [f"{a}_{h}_vs_{b}_{h}" for h in HPFS
                for a, b in [("BL","D"), ("BL","LD"), ("D","LD")]]
F2_CONTRASTS = [f"{c}_{h1}_vs_{c}_{h2}" for c in CONDS
                for h1, h2 in [(24,48),(48,72),(72,96),(96,120)]]

# Build long table of Tier A + B rows
tierA2 = tierA.assign(tier="A")
tierB2 = tierB.assign(tier="B")
all_sig = pd.concat([tierA2, tierB2], ignore_index=True)
all_sig["signed_W1"] = all_sig.W1_obs * all_sig.direction
# One row per (contrast, subsystem) — the row's family tags the panel
f1_subs = sorted(all_sig[all_sig.family == "F1"].subsystem.unique())
f2_subs = sorted(all_sig[all_sig.family == "F2"].subsystem.unique())

def _build_grid(contrasts, subsystems, family):
    arr = np.full((len(subsystems), len(contrasts)), np.nan)
    tier_mat = np.full((len(subsystems), len(contrasts)), "", dtype=object)
    sub_idx = {s: i for i, s in enumerate(subsystems)}
    con_idx = {c: j for j, c in enumerate(contrasts)}
    sub = all_sig[all_sig.family == family]
    for _, r in sub.iterrows():
        if r.contrast not in con_idx or r.subsystem not in sub_idx:
            continue
        i, j = sub_idx[r.subsystem], con_idx[r.contrast]
        arr[i, j] = r.signed_W1
        tier_mat[i, j] = r.tier
    return arr, tier_mat


fig3 = plt.figure(figsize=(16, 7.6))
gs = fig3.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.5)

for panel_i, (family, contrasts, subsystems) in enumerate([
    ("F1", F1_CONTRASTS, f1_subs),
    ("F2", F2_CONTRASTS, f2_subs),
]):
    ax = fig3.add_subplot(gs[0, panel_i])
    arr, tiermat = _build_grid(contrasts, subsystems, family)
    vmax = np.nanmax(np.abs(arr)) if np.isfinite(arr).any() else 1
    im = ax.imshow(arr, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_yticks(range(len(subsystems)))
    ax.set_yticklabels(subsystems, fontsize=7.5)
    ax.set_xticks(range(len(contrasts)))
    short_contrasts = [c.replace("_vs_", "–").replace("_", " ") for c in contrasts]
    ax.set_xticklabels(short_contrasts, rotation=55, ha="right", fontsize=7)
    # Hatch Tier B
    for i in range(len(subsystems)):
        for j in range(len(contrasts)):
            if tiermat[i, j] == "B":
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       fill=False, hatch="///",
                                       edgecolor="black", linewidth=0.3))
            if tiermat[i, j] in ("A", "B"):
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                                       fill=False, edgecolor="black",
                                       linewidth=0.6))
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("signed W1  (direction × W1)", fontsize=8)
    label = ("F1 — between-condition at each hpf"
             if family == "F1"
             else "F2 — adjacent-hpf within each regime")
    ax.set_title(f"{'A' if panel_i == 0 else 'B'} — {label}",
                 fontsize=11, loc="left")

# Legend
leg_handles = [
    Patch(facecolor="white", edgecolor="black", linewidth=0.8, label="Tier A"),
    Patch(facecolor="white", edgecolor="black", hatch="///", label="Tier B"),
]
fig3.legend(handles=leg_handles, loc="lower center", ncol=2, fontsize=9,
            frameon=False, bbox_to_anchor=(0.5, -0.02))
fig3.suptitle("Figure 3 — orientation heatmap: where in the (contrast × subsystem) grid things moved",
              fontsize=11, y=1.02)
fig3.text(0.5, -0.06,
          "Cells are coloured by signed W1 (direction × pooled Wasserstein distance). "
          "Tier A cells are solid-bordered; Tier B are hatched; empty cells carried no "
          "finding at q ≤ 0.10. This panel is an orientation device — it does not support "
          "any particular claim; the density panels (Figures 4–6) carry those.",
          ha="center", fontsize=8, style="italic", wrap=True)

fig3.savefig(OUT / "fig3_heatmap.pdf", bbox_inches="tight")
fig3.savefig(OUT / "fig3_heatmap.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig3_heatmap.pdf'} / .png")
plt.close(fig3)


# ============================================================
# Ridgeline helper
# ============================================================
from scipy.stats import gaussian_kde

def _ridge(ax, y_base, log2v, zero_frac, color, height=0.8, grid_x=None):
    """Draw a filled density ridge at vertical position y_base. Also draws
    a thin zero-fraction bar above the ridge (0..height_zero), coloured as
    the ridge, so the two-tier test's zero contribution is visible."""
    if grid_x is None:
        grid_x = np.linspace(-10, 10, 300)
    if log2v.size >= 3 and np.std(log2v) > 0:
        kde = gaussian_kde(log2v)
        y = kde(grid_x)
        y = y / y.max() * height
        ax.fill_between(grid_x, y_base, y_base + y, color=color, alpha=0.45,
                        edgecolor=color, linewidth=0.8)
    # Zero-fraction bar (max height = height * 0.25)
    zf_h = height * 0.25
    zf_y = y_base + height + 0.05
    ax.add_patch(Rectangle((-10, zf_y), 20 * zero_frac, zf_h,
                            facecolor=color, alpha=0.6, edgecolor="none"))
    ax.text(-10 + 20 * zero_frac + 0.15, zf_y + zf_h / 2,
            f"{zero_frac*100:.0f}% zero",
            fontsize=6, va="center", color=color)


def _panel_ridgelines(ax, subsystem, hpf_list, title, xlim=(-10, 10)):
    """One panel: subsystem × hpf_list × 3 regimes.  Returns (list of yticks, list of ylabels)."""
    yticks, ylabels = [], []
    grid_x = np.linspace(xlim[0], xlim[1], 300)
    for row_i, hpf in enumerate(hpf_list):
        y_base = row_i * 2.5  # spacing between ridges
        for cond in CONDS:
            cell = f"{cond}_{hpf}"
            log2v, zf = _flux_pool(cell, subsystem)
            if log2v.size == 0 and zf == 0:
                continue
            _ridge(ax, y_base, log2v, zf, COND_COLOR[cond], height=0.9, grid_x=grid_x)
        yticks.append(y_base + 0.6)
        ylabels.append(f"{hpf} hpf")
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlim(xlim)
    ax.set_xlabel("log2(|v| / pool-mean |v|)", fontsize=8)
    ax.set_title(title, fontsize=10, loc="left")
    ax.axvline(0, color="grey", lw=0.4, alpha=0.6)
    ax.grid(axis="x", alpha=0.25)


# ============================================================
# FIGURE 4 — developmental lipid program
# ============================================================
fig4 = plt.figure(figsize=(14, 9))
gs = fig4.add_gridspec(2, 2, hspace=0.45, wspace=0.3)
subs_fig4 = [
    ("Fatty acid activation (cytosolic)",  "A — Fatty acid activation (cytosolic)"),
    ("Fatty acid oxidation",               "B — Fatty acid oxidation"),
    ("Arachidonic acid metabolism",        "C — Arachidonic acid metabolism"),
    ("Leukotriene metabolism",             "D — Leukotriene metabolism"),
]
for i, (subsystem, title) in enumerate(subs_fig4):
    ax = fig4.add_subplot(gs[i // 2, i % 2])
    _panel_ridgelines(ax, subsystem, HPFS, title)

# Legend
legend_elements = [plt.Line2D([0], [0], marker="s", color=c, linestyle="",
                              markersize=10, label=k) for k, c in COND_COLOR.items()]
fig4.legend(handles=legend_elements, loc="lower center", ncol=3,
            fontsize=10, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig4.suptitle("Figure 4 — developmental lipid program  (§4.2; 24→120 hpf)",
              fontsize=12, y=0.995)
fig4.text(0.5, -0.03,
          "Each ridge is the pooled flux distribution for one subsystem at one hpf, "
          "overlaid for all three regimes (BL / D / LD). Thin bar above each ridge: "
          "zero-fraction (Fisher-test input of the two-tier method). "
          "The three regimes sit on top of each other within each hpf: light does not "
          "reverse the shared developmental transition — Fatty acid activation and "
          "oxidation both intensify from 24 to 72 hpf and stay elevated, Arachidonic "
          "acid and Leukotriene shift at 48→72 — but the shared trajectory is the "
          "signal; between-regime overlap at each hpf is the point.",
          ha="center", fontsize=8, style="italic", wrap=True)
fig4.savefig(OUT / "fig4_lipid_program.pdf", bbox_inches="tight")
fig4.savefig(OUT / "fig4_lipid_program.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig4_lipid_program.pdf'} / .png")
plt.close(fig4)


# ============================================================
# FIGURE 5 — amino-acid separation at 72–96 hpf
# ============================================================
fig5 = plt.figure(figsize=(13, 5.5))
gs = fig5.add_gridspec(1, 3, wspace=0.4)
subs_fig5 = [
    ("Arginine and proline metabolism",                    "A — Arginine and proline"),
    ("Phenylalanine, tyrosine and tryptophan biosynthesis","B — Phe/Tyr/Trp biosynthesis"),
    ("Valine, leucine, and isoleucine metabolism",         "C — Val/Leu/Ile metabolism"),
]
for i, (subsystem, title) in enumerate(subs_fig5):
    ax = fig5.add_subplot(gs[0, i])
    _panel_ridgelines(ax, subsystem, [72, 96], title)

fig5.legend(handles=legend_elements, loc="lower center", ncol=3,
            fontsize=10, frameon=False, bbox_to_anchor=(0.5, -0.03))
fig5.suptitle("Figure 5 — amino-acid separation at 72–96 hpf  (§4.3)",
              fontsize=12, y=1.02)
fig5.text(0.5, -0.08,
          "Three amino-acid subsystems with Tier A F1 signals at 72 and/or 96 hpf. "
          "Arginine and proline leads: BL_96 vs LD_96 W1 = 1.75 is the largest "
          "effect in the frozen F1 set. The three-regime overlap in Figure 4 is "
          "visibly absent here — at 72 and 96 hpf the three condition densities "
          "separate on the x-axis.",
          ha="center", fontsize=8, style="italic", wrap=True)
fig5.savefig(OUT / "fig5_amino_acid_72_96.pdf", bbox_inches="tight")
fig5.savefig(OUT / "fig5_amino_acid_72_96.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig5_amino_acid_72_96.pdf'} / .png")
plt.close(fig5)


# ============================================================
# FIGURE 6 — convergence at 120 hpf
# ============================================================
fig6 = plt.figure(figsize=(16, 8.5))
gs = fig6.add_gridspec(2, 3, hspace=0.5, wspace=0.4,
                       width_ratios=[1, 1, 1.1])
subs_fig6 = [
    ("Bile acid biosynthesis",      "A — Bile acid biosynthesis",       0),
    ("Steroid metabolism",          "B — Steroid metabolism",           1),
    ("Leukotriene metabolism",      "C — Leukotriene metabolism",       2),
    ("Acylglycerides metabolism",   "D — Acylglycerides metabolism",    3),
]
for subsystem, title, pos in subs_fig6:
    ax = fig6.add_subplot(gs[pos // 3, pos % 3])
    _panel_ridgelines(ax, subsystem, [120], title)

# Panel E — signed W1 dot plot: BL-vs-LD and D-vs-LD at 120 hpf
axE = fig6.add_subplot(gs[1, 1:])
e_subs = ["Bile acid biosynthesis", "Steroid metabolism",
          "Leukotriene metabolism", "Acylglycerides metabolism",
          "Valine, leucine, and isoleucine metabolism",
          "Phenylalanine, tyrosine and tryptophan biosynthesis"]
ys = np.arange(len(e_subs))
for subsystem in e_subs:
    pass  # placeholder; actually iterate below
for i, subsystem in enumerate(e_subs):
    for pair_lbl, color, offset in [("BL_120_vs_LD_120", COND_COLOR["BL"], -0.15),
                                     ("D_120_vs_LD_120",  COND_COLOR["D"],  +0.15)]:
        r = tierA[(tierA.contrast == pair_lbl) &
                  (tierA.subsystem == subsystem)]
        if len(r):
            signed_W1 = r.iloc[0].W1_obs * r.iloc[0].direction
            axE.scatter(signed_W1, ys[i] + offset, s=110, color=color,
                        edgecolor="black", linewidths=0.7, zorder=5,
                        label=pair_lbl if i == 0 else "")
            axE.plot([0, signed_W1], [ys[i] + offset, ys[i] + offset],
                     color=color, lw=1.2, alpha=0.5)
axE.axvline(0, color="grey", linewidth=0.6)
axE.set_yticks(ys)
axE.set_yticklabels(e_subs, fontsize=8)
axE.invert_yaxis()
axE.set_xlabel("signed W1  (direction × pooled W1)", fontsize=8)
axE.set_title("E — BL-vs-LD and D-vs-LD at 120 hpf: both non-standard regimes "
              "deviate from LD in the same direction",
              fontsize=10, loc="left")
axE.grid(axis="x", alpha=0.3)
leg_e = [plt.Line2D([0], [0], marker="o", color=COND_COLOR["BL"], linestyle="",
                     markersize=9, markeredgecolor="black",
                     label="BL vs LD at 120 hpf"),
         plt.Line2D([0], [0], marker="o", color=COND_COLOR["D"], linestyle="",
                     markersize=9, markeredgecolor="black",
                     label="D vs LD at 120 hpf")]
axE.legend(handles=leg_e, loc="lower right", fontsize=8, frameon=False)

# Density-panel legend
fig6.legend(handles=legend_elements, loc="lower center", ncol=3,
            fontsize=10, frameon=False, bbox_to_anchor=(0.5, -0.03))
fig6.suptitle("Figure 6 — convergence at 120 hpf: BL and D depart from LD in the same direction  (§4.3)",
              fontsize=12, y=1.0)
fig6.text(0.5, -0.055,
          "A–D: pooled flux distributions at 120 hpf for the four subsystems with "
          "the largest BL-vs-LD and D-vs-LD effects. E: signed W1 dot plot — points "
          "to the right of zero mean the regime's distribution is shifted upward in "
          "log2(|v|/μ) vs LD; points to the left are shifted downward. Both "
          "non-standard regimes (BL and D) deviate from LD in the same direction on "
          "every one of six amphipathic-metabolism subsystems — the central "
          "photoperiod convergence result.",
          ha="center", fontsize=8, style="italic", wrap=True)
fig6.savefig(OUT / "fig6_convergence_120.pdf", bbox_inches="tight")
fig6.savefig(OUT / "fig6_convergence_120.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig6_convergence_120.pdf'} / .png")
plt.close(fig6)

print("\n--- all four figures written ---")
print(f"  {OUT/'fig3_heatmap.pdf'}")
print(f"  {OUT/'fig4_lipid_program.pdf'}")
print(f"  {OUT/'fig5_amino_acid_72_96.pdf'}")
print(f"  {OUT/'fig6_convergence_120.pdf'}")
