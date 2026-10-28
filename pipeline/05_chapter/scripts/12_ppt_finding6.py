"""
12_ppt_finding6.py — presentation-ready figures for Finding 6.

Finding 6: at 120 hpf, LD alone steps up sterol and bile acid metabolism;
the same stage separates the core clock transcripts of LD from both
constant regimes.

Writes to chapter_assets/ppt_figures/finding6/:
    6a_sterol_bile_120.{pdf,png}      — 120 hpf regime-pair heatmap, 3 subsystems
    6b_journey.{pdf,png}              — trajectories, 120 hpf band highlighted
    6c_clock_panel.{pdf,png}          — 29-gene clock panel heatmap at 120 hpf
    6d_combined.{pdf,png}             — combined slide (sterol/bile on top, clock on bottom)
"""
import re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import gaussian_kde

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import POOLED, SUBSETS, COMPARISON_TABLES, TABLES, FIGURES

OUT = FIGURES.parent / "ppt_figures" / "finding6"
OUT.mkdir(parents=True, exist_ok=True)

REGIMES = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
REGIME_COLORS = {"BL": "#1976D2", "D": "#333333", "LD": "#F6B52C"}
MAROON = "#800020"

STEROL_BILE = [
    "Bile acid biosynthesis",
    "Steroid metabolism",
    "Cholesterol biosynthesis 1 (Bloch pathway)",
]
SHORT_NAMES = {
    "Bile acid biosynthesis": "Bile acid",
    "Steroid metabolism": "Steroid",
    "Cholesterol biosynthesis 1 (Bloch pathway)": "Cholesterol (Bloch)",
}
# 9 clock panel genes from §3.4 — LD exceeds both constants by ≥ 0.5 zFPKM
PANEL_GENES_BLD = ["cry1b", "cry2", "per1a", "per1b", "per3"]        # BL ≤ D < LD
PANEL_GENES_DBL = ["bhlhe40", "cry1a", "nr1d1", "per2"]              # D < BL < LD

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 12,
    "axes.titlesize": 15, "axes.labelsize": 13,
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


# ---------- Load data ----------
print("[load] subsystem_summary_F1, F2, clock panels")
s1 = pd.read_csv(COMPARISON_TABLES / "subsystem_summary_F1_centered.csv")
s2 = pd.read_csv(COMPARISON_TABLES / "subsystem_summary_F2_centered.csv")

def _parse_f1(c):
    m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_\d+", c)
    return (m.group(1), m.group(3), int(m.group(2))) if m else (None, None, None)
def _parse_f2(c):
    m = re.match(r"([A-Z]+)_(\d+)_vs_[A-Z]+_(\d+)", c)
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else (None, None, None)

s1[["pair_a", "pair_b", "hpf"]] = s1.contrast.apply(lambda c: pd.Series(_parse_f1(c)))
s2[["regime", "hpf_a", "hpf_b"]] = s2.contrast.apply(lambda c: pd.Series(_parse_f2(c)))
s1["hpf"] = s1.hpf.astype(int)
s2["hpf_a"] = s2.hpf_a.astype(int); s2["hpf_b"] = s2.hpf_b.astype(int)

clock_full = pd.read_csv(TABLES / "table_S8_clock_panel_full.csv")
print(f"[load] clock panel full: {len(clock_full)} rows")


# ============================================================
# 6a — sterol/bile at 120 hpf (regime-pair x subsystem)
# ============================================================
print("\n[6a] sterol/bile 120 hpf")

F1_PAIRS = [("BL", "D"), ("BL", "LD"), ("D", "LD")]
pair_titles = [f"{a} vs {b}" for a, b in F1_PAIRS]

mats = []
for a, b in F1_PAIRS:
    sub = s1[(s1.pair_a == a) & (s1.pair_b == b) & (s1.hpf == 120)]
    val = pd.Series({s: sub[sub.subsystem == s].mean_log2FC_centered.values[0]
                      if ((sub.subsystem == s).sum()) else np.nan
                      for s in STEROL_BILE})
    q = pd.Series({s: sub[sub.subsystem == s].enrichment_q.values[0]
                    if ((sub.subsystem == s).sum()) else np.nan
                    for s in STEROL_BILE})
    n = pd.Series({s: int(sub[sub.subsystem == s].n_reps_tested.values[0])
                    if ((sub.subsystem == s).sum()) else 0
                    for s in STEROL_BILE})
    mats.append((val, q, n))

VMAX = 3.0
norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)
cmap = plt.cm.RdBu_r
cmap.set_bad("#e5e5e5")

fig_6a, axes_6a = plt.subplots(1, 3, figsize=(12, 4.5), sharey=True,
                                 gridspec_kw={"wspace": 0.05})
for ax, (val, q, n), title in zip(axes_6a, mats, pair_titles):
    display = val.values.astype(float).reshape(-1, 1).copy()
    display[(n == 0).values.reshape(-1, 1)] = np.nan
    im = ax.imshow(display, aspect="auto", cmap=cmap, norm=norm,
                     interpolation="nearest")
    for i, (qv, nv) in enumerate(zip(q.values, n.values)):
        if pd.notna(qv) and qv < 0.05:
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), 1, 1, fill=False,
                                         edgecolor=MAROON, linewidth=2.5, zorder=10))
    ax.set_xticks([0]); ax.set_xticklabels(["120 hpf"])
    ax.set_title(title, fontsize=13, weight="bold")
    ax.set_yticks(range(len(STEROL_BILE)))
    row_labels = [f"{SHORT_NAMES.get(s, s)}\n(n={int(n.loc[s])})" for s in STEROL_BILE]
    if ax is axes_6a[0]:
        ax.set_yticklabels(row_labels, fontsize=11)
    else:
        ax.tick_params(labelleft=False)
cbar = fig_6a.colorbar(im, ax=axes_6a, pad=0.015, fraction=0.03, shrink=0.9,
                        extend="both")
cbar.set_label("mean centered log$_2$FC\n(red = alphabetically earlier regime higher)",
                fontsize=10)
fig_6a.suptitle("Figure 6a — Sterol & bile at 120 hpf: LD separates from BL and D",
                 fontsize=14, weight="bold", y=0.995)
fig_6a.text(0.5, 0.90,
             "Maroon box = q < 0.05   ·   Cholesterol is flat, so the step is downstream of sterol synthesis",
             ha="center", fontsize=10, color="#555555")
fig_6a.tight_layout(rect=[0, 0, 1, 0.86])
_save(fig_6a, "6a_sterol_bile_120")


# ============================================================
# 6b — journey (trajectory) with 120 hpf band highlighted
# ============================================================
print("\n[6b] journey trajectory, 120 hpf band highlighted")

# per-cell mean centered log flux for the 2 flux subsystems — reuse per_reaction_flux_F1_centered
f1 = pd.read_csv(COMPARISON_TABLES / "per_reaction_flux_F1_centered.csv")
a = f1[["cell_A", "rxn_id", "subsystem", "mean_nz_A", "median_log2_A", "loop_flagged"]].rename(
    columns={"cell_A": "cell", "mean_nz_A": "mean_nz", "median_log2_A": "median_log2"})
b = f1[["cell_B", "rxn_id", "subsystem", "mean_nz_B", "median_log2_B", "loop_flagged"]].rename(
    columns={"cell_B": "cell", "mean_nz_B": "mean_nz", "median_log2_B": "median_log2"})
long = pd.concat([a, b]).drop_duplicates(subset=["cell", "rxn_id"])
long = long[(~long.loop_flagged) & (long.mean_nz > 1e-10)].copy()
long["centered_log2"] = np.log2(long["mean_nz"]) - long["median_log2"]

CELLS_DISPLAY = tuple(f"{c}_{h}" for h in HPFS for c in REGIMES)
mat = (long.groupby(["subsystem", "cell"]).centered_log2.mean()
            .unstack("cell").reindex(columns=list(CELLS_DISPLAY)))
n_reps = (long.groupby(["subsystem", "cell"]).rxn_id.nunique()
               .unstack("cell").reindex(columns=list(CELLS_DISPLAY)).fillna(0))

sterol_bile_in_mat = [s for s in ["Bile acid biosynthesis", "Steroid metabolism"] if s in mat.index]

fig_6b, axes_6b = plt.subplots(1, 2, figsize=(12, 5.2), sharex=True)
for ax, sub in zip(axes_6b, sterol_bile_in_mat):
    for reg in REGIMES:
        xs, ys = [], []
        for h in HPFS:
            cell = f"{reg}_{h}"
            if n_reps.loc[sub, cell] == 0 or pd.isna(mat.loc[sub, cell]):
                xs.append(h); ys.append(np.nan)
            else:
                xs.append(h); ys.append(mat.loc[sub, cell])
        ax.plot(xs, ys, "-o", color=REGIME_COLORS[reg], label=reg,
                 linewidth=2.5, markersize=10, markeredgecolor="white",
                 markeredgewidth=1.2)
    ax.axvspan(108, 120, alpha=0.15, color="#d94801")  # 120 hpf band
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.set_xticks(list(HPFS))
    n_max = int(n_reps.loc[sub].max())
    ax.set_title(f"{sub}  (n={n_max})", fontsize=14, loc="left", weight="bold")
    ax.set_xlabel("hpf")
    ax.grid(True, linestyle=":", alpha=0.3)
axes_6b[0].set_ylabel("centered log$_2$|flux|")
axes_6b[-1].legend(title="Regime", loc="upper left", fontsize=11, frameon=True)
fig_6b.suptitle("Figure 6b — Sterol & bile journey: flat through 96 hpf, LD climbs at 120",
                 fontsize=14, weight="bold", y=0.99)
fig_6b.tight_layout(rect=[0, 0, 1, 0.95])
_save(fig_6b, "6b_journey")


# ============================================================
# 6c — clock gene panel heatmap at 120 hpf
# ============================================================
print("\n[6c] clock gene panel heatmap")

# order: panel genes first (grouped by ordering), then rest
ordered_genes = (
    PANEL_GENES_BLD + PANEL_GENES_DBL +
    [g for g in clock_full.symbol if g not in PANEL_GENES_BLD + PANEL_GENES_DBL]
)
clock = clock_full.set_index("symbol").reindex(ordered_genes)
# drop rows with no data
has_data = clock[["BL_120", "D_120", "LD_120"]].notna().all(axis=1)
clock_plot = clock[has_data].copy()
print(f"[6c] {len(clock_plot)} / {len(clock)} genes with zFPKM at 120 hpf")

# compute "LD margin over higher constant" to decide maroon box
clock_plot["ld_margin"] = clock_plot["LD_120"] - clock_plot[["BL_120", "D_120"]].max(axis=1)

# diverging colormap, zFPKM typically ranges -3 to 3
vmax = 3.0
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
cmap = plt.cm.RdBu_r
cmap.set_bad("#e5e5e5")

fig_6c, ax_6c = plt.subplots(figsize=(7, 10))
display = clock_plot[["BL_120", "D_120", "LD_120"]].values.astype(float)
im = ax_6c.imshow(display, aspect="auto", cmap=cmap, norm=norm,
                   interpolation="nearest")

# maroon box for the 9 panel genes (LD margin >= 0.5)
for i, (sym, row) in enumerate(clock_plot.iterrows()):
    if row.ld_margin >= 0.5:
        ax_6c.add_patch(plt.Rectangle((-0.5, i - 0.5), 3, 1, fill=False,
                                        edgecolor=MAROON, linewidth=2.5, zorder=10))

ax_6c.set_xticks(range(3))
ax_6c.set_xticklabels(["BL", "D", "LD"], fontsize=12, weight="bold")
ax_6c.set_yticks(range(len(clock_plot)))
# add thread markers in label
y_labels = []
for g in clock_plot.index:
    marker = ""
    if g in PANEL_GENES_BLD:
        marker = "  ★"
    elif g in PANEL_GENES_DBL:
        marker = "  ◆"
    y_labels.append(f"{g}{marker}")
ax_6c.set_yticklabels(y_labels, fontsize=10)
cbar = fig_6c.colorbar(im, ax=ax_6c, pad=0.02, fraction=0.05, shrink=0.9, extend="both")
cbar.set_label("zFPKM at 120 hpf", fontsize=11)

ax_6c.set_title("Figure 6c — Core clock panel at 120 hpf\n"
                 "Maroon box = LD exceeds both constants by ≥ 0.5 zFPKM (9/23 genes).\n"
                 "★ = BL ≤ D < LD ordering    ◆ = D < BL < LD ordering",
                 fontsize=12, weight="bold", loc="left", pad=14)
fig_6c.tight_layout()
_save(fig_6c, "6c_clock_panel")


# ============================================================
# 6d — combined slide
# ============================================================
print("\n[6d] combined slide")

fig_6d = plt.figure(figsize=(18, 11))
gs = fig_6d.add_gridspec(2, 2, height_ratios=[1, 1.1], hspace=0.4, wspace=0.3,
                           left=0.07, right=0.95, top=0.90, bottom=0.07)

# top-left: sterol/bile 120 hpf heatmap (compact 3 pairs)
gs_top_l = gs[0, 0].subgridspec(1, 3, wspace=0.08)
axes_top_l = [fig_6d.add_subplot(gs_top_l[0, i]) for i in range(3)]
for ax, (val, q, n), title in zip(axes_top_l, mats, pair_titles):
    display = val.values.astype(float).reshape(-1, 1).copy()
    display[(n == 0).values.reshape(-1, 1)] = np.nan
    im1 = ax.imshow(display, aspect="auto", cmap=cmap, norm=norm,
                      interpolation="nearest")
    for i, qv in enumerate(q.values):
        if pd.notna(qv) and qv < 0.05:
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), 1, 1, fill=False,
                                         edgecolor=MAROON, linewidth=2.5, zorder=10))
    ax.set_xticks([0]); ax.set_xticklabels(["120 hpf"], fontsize=10)
    ax.set_title(title, fontsize=12, weight="bold")
    ax.set_yticks(range(len(STEROL_BILE)))
    if ax is axes_top_l[0]:
        ax.set_yticklabels([f"{SHORT_NAMES.get(s, s)}\n(n={int(n.loc[s])})"
                             for s in STEROL_BILE], fontsize=10)
    else:
        ax.tick_params(labelleft=False)
cbar = fig_6d.colorbar(im1, ax=axes_top_l, pad=0.015, fraction=0.04, shrink=0.9, extend="both")
cbar.set_label("centered log$_2$FC", fontsize=10)
axes_top_l[1].text(0.5, 1.35, "A  Metabolic step at 120 hpf (sterol & bile)",
                     transform=axes_top_l[1].transAxes, fontsize=14, weight="bold",
                     ha="center")

# top-right: journey trajectories
gs_top_r = gs[0, 1].subgridspec(1, 2, wspace=0.3)
axes_top_r = [fig_6d.add_subplot(gs_top_r[0, i]) for i in range(2)]
for ax, sub in zip(axes_top_r, sterol_bile_in_mat):
    for reg in REGIMES:
        xs, ys = [], []
        for h in HPFS:
            cell = f"{reg}_{h}"
            if n_reps.loc[sub, cell] == 0 or pd.isna(mat.loc[sub, cell]):
                xs.append(h); ys.append(np.nan)
            else:
                xs.append(h); ys.append(mat.loc[sub, cell])
        ax.plot(xs, ys, "-o", color=REGIME_COLORS[reg], label=reg,
                 linewidth=2.0, markersize=7, markeredgecolor="white",
                 markeredgewidth=0.8)
    ax.axvspan(108, 120, alpha=0.15, color="#d94801")
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.set_xticks(list(HPFS))
    ax.set_title(sub, fontsize=11, loc="left")
    ax.set_xlabel("hpf", fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.3)
axes_top_r[0].set_ylabel("centered log$_2$|flux|", fontsize=10)
axes_top_r[-1].legend(title="Regime", loc="upper left", fontsize=9, frameon=True)
axes_top_r[0].text(1.08, 1.1, "B  Journey across hpf",
                     transform=axes_top_r[0].transAxes, fontsize=14, weight="bold",
                     ha="center")

# bottom spans full width: clock panel
gs_bot = gs[1, :].subgridspec(1, 2, width_ratios=[0.35, 1], wspace=0.08)
ax_6d_clock = fig_6d.add_subplot(gs_bot[0, 0])
im3 = ax_6d_clock.imshow(display, aspect="auto")  # placeholder
ax_6d_clock.clear()

display_cl = clock_plot[["BL_120", "D_120", "LD_120"]].values.astype(float)
im3 = ax_6d_clock.imshow(display_cl, aspect="auto", cmap=cmap, norm=norm,
                           interpolation="nearest")
for i, (sym, row) in enumerate(clock_plot.iterrows()):
    if row.ld_margin >= 0.5:
        ax_6d_clock.add_patch(plt.Rectangle((-0.5, i - 0.5), 3, 1, fill=False,
                                              edgecolor=MAROON, linewidth=2.5,
                                              zorder=10))
ax_6d_clock.set_xticks(range(3))
ax_6d_clock.set_xticklabels(["BL", "D", "LD"], fontsize=11, weight="bold")
ax_6d_clock.set_yticks(range(len(clock_plot)))
ax_6d_clock.set_yticklabels(y_labels, fontsize=9)
cbar3 = fig_6d.colorbar(im3, ax=ax_6d_clock, pad=0.015, fraction=0.06, shrink=0.85,
                          extend="both")
cbar3.set_label("zFPKM", fontsize=10)
ax_6d_clock.set_title("C  Core clock panel at 120 hpf",
                         fontsize=14, weight="bold", loc="left")

# bottom-right: annotation text + summary
ax_6d_text = fig_6d.add_subplot(gs_bot[0, 1])
ax_6d_text.axis("off")
ax_6d_text.text(0.0, 0.98,
                  "D  The two halves of Finding 6",
                  transform=ax_6d_text.transAxes, fontsize=14, weight="bold",
                  va="top")
summary = (
    "Metabolism (panels A, B):\n"
    "  • At 120 hpf, LD alone steps up bile acid biosynthesis (23/23 reactions rise, q = 1.2 × 10⁻³)\n"
    "    and steroid metabolism (18/18 reactions, q = 2.2 × 10⁻⁵)\n"
    "  • Flat in BL and D (96→120); cholesterol biosynthesis unchanged in all three\n"
    "    → the step sits downstream of sterol synthesis, in sterol modification and bile acid formation\n"
    "  • At 120 hpf: LD > BL on bile acid (20/24 reactions, q = 4.9 × 10⁻¹¹)\n"
    "               LD > D on bile acid (19/21, q = 1.7 × 10⁻⁶) and on steroid (17/18, q = 1.9 × 10⁻³)\n"
    "               BL and D do not differ on either subsystem\n"
    "\n"
    "Clock transcripts (panel C):\n"
    "  • LD_120 highest on 17 of 23 genes with data\n"
    "  • On 9 of them LD exceeds both constants by ≥ 0.5 zFPKM (★ 5 BL ≤ D < LD;  ◆ 4 D < BL < LD)\n"
    "  • No gene has D above BL by that margin\n"
    "\n"
    "The metabolic step and the clock separation group the two constant regimes together against\n"
    "the cycling one — the sterol/bile step is the candidate metabolic consequence of entrained\n"
    "clock output."
)
ax_6d_text.text(0.0, 0.92, summary, transform=ax_6d_text.transAxes,
                  fontsize=11, va="top", family="sans-serif", linespacing=1.4)

fig_6d.suptitle("Figure 6 — LD alone steps up sterol / bile at 120 hpf, when the clock separates it from BL and D",
                 fontsize=16, weight="bold", y=0.97)
_save(fig_6d, "6d_combined")

print("\n[done]  outputs in", OUT)
