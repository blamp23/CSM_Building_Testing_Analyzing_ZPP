"""
13_ppt_finding3_7.py — presentation-ready figures for Finding 3 (melatonin) and
Finding 7 (aromatic amino acids). Both live in the signaling / AA thread.

Finding 3: melatonin synthesis runs through a 48→96 hpf window that only the
cycling regime holds steady; BL loses it at 72, D at 96.
Finding 7: Phe/Tyr/Trp biosynthesis peaks at 72 hpf and falls in D and LD; BL
never makes the peak or the fall, so by 120 hpf BL sits above both.

Writes to chapter_assets/ppt_figures/finding3_7/:
    3a_melatonin_window.{pdf,png}  — trajectory + 48-96 band + 72/96 hpf heatmap
    7a_aromatic_flip.{pdf,png}     — Phe/Tyr/Trp + Tyrosine trajectories + 72 vs 120 heatmap
    37_combined.{pdf,png}          — one-slide summary
"""
import re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import COMPARISON_TABLES, FIGURES

OUT = FIGURES.parent / "ppt_figures" / "finding3_7"
OUT.mkdir(parents=True, exist_ok=True)

REGIMES = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
REGIME_COLORS = {"BL": "#1976D2", "D": "#333333", "LD": "#F6B52C"}
MAROON = "#800020"
VMAX = 3.0

MELATONIN = "Serotonin and melatonin biosynthesis"
AROMATIC = "Phenylalanine, tyrosine and tryptophan biosynthesis"
TYROSINE = "Tyrosine metabolism"

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


# ---------- load data ----------
print("[load] data")
f1 = pd.read_csv(COMPARISON_TABLES / "per_reaction_flux_F1_centered.csv")
s1 = pd.read_csv(COMPARISON_TABLES / "subsystem_summary_F1_centered.csv")

def _parse_f1(c):
    m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_\d+", c)
    return (m.group(1), m.group(3), int(m.group(2))) if m else (None, None, None)
s1[["pair_a", "pair_b", "hpf"]] = s1.contrast.apply(lambda c: pd.Series(_parse_f1(c)))
s1["hpf"] = s1.hpf.astype(int)

# per-cell matrix
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


# ---------- trajectory helper ----------
def _trajectory(ax, subsystem, highlight_bands=None, title=None):
    for reg in REGIMES:
        xs, ys = [], []
        for h in HPFS:
            cell = f"{reg}_{h}"
            if n_reps.loc[subsystem, cell] == 0 or pd.isna(mat.loc[subsystem, cell]):
                xs.append(h); ys.append(np.nan)
            else:
                xs.append(h); ys.append(mat.loc[subsystem, cell])
        ax.plot(xs, ys, "-o", color=REGIME_COLORS[reg], label=reg,
                 linewidth=2.5, markersize=10, markeredgecolor="white",
                 markeredgewidth=1.2, zorder=5)
    for band in (highlight_bands or []):
        ax.axvspan(band[0], band[1], alpha=0.12, color="#d94801", zorder=1)
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.set_xticks(list(HPFS))
    n_max = int(n_reps.loc[subsystem].max())
    if title:
        ax.set_title(f"{title}  (n={n_max})", fontsize=13, loc="left", weight="bold")
    ax.set_xlabel("hpf")
    ax.grid(True, linestyle=":", alpha=0.3)


# ---------- F1 single-stage heatmap helper ----------
def _f1_at_stage(subs, hpf):
    """Return (3 pairs × subs) matrix of mean_log2FC + q for a given hpf."""
    pairs = [("BL", "D"), ("BL", "LD"), ("D", "LD")]
    val = pd.DataFrame(index=subs, columns=[f"{a} vs {b}" for a, b in pairs], dtype=float)
    q = val.copy()
    n = val.copy()
    for a, b in pairs:
        sub = s1[(s1.pair_a == a) & (s1.pair_b == b) & (s1.hpf == hpf)]
        for s in subs:
            r = sub[sub.subsystem == s]
            if len(r):
                val.loc[s, f"{a} vs {b}"] = r.mean_log2FC_centered.values[0]
                q.loc[s, f"{a} vs {b}"] = r.enrichment_q.values[0]
                n.loc[s, f"{a} vs {b}"] = int(r.n_reps_tested.values[0])
            else:
                n.loc[s, f"{a} vs {b}"] = 0
    return val, q, n


def _draw_stage_heatmap(ax, val, q, n, title=None, show_ylabels=True):
    norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)
    cmap = plt.cm.RdBu_r
    cmap.set_bad("#e5e5e5")
    display = val.values.astype(float).copy()
    display[(n == 0).values] = np.nan
    im = ax.imshow(display, aspect="auto", cmap=cmap, norm=norm,
                     interpolation="nearest")
    for i in range(val.shape[0]):
        for j in range(val.shape[1]):
            qv = q.iloc[i, j]
            if pd.notna(qv) and qv < 0.05:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                         edgecolor=MAROON, linewidth=2.5, zorder=10))
    ax.set_xticks(range(len(val.columns)))
    ax.set_xticklabels(val.columns, fontsize=10)
    ax.set_yticks(range(len(val.index)))
    if show_ylabels:
        ax.set_yticklabels([f"{s}\n(n={int(n.iloc[i].max())})"
                             for i, s in enumerate(val.index)], fontsize=10)
    else:
        ax.tick_params(labelleft=False)
    if title:
        ax.set_title(title, fontsize=12, weight="bold", pad=6)
    return im


# ============================================================
# 3a — melatonin window
# ============================================================
print("\n[3a] melatonin 48-96 window")

fig_3a = plt.figure(figsize=(16, 6))
gs = fig_3a.add_gridspec(1, 2, width_ratios=[1.4, 1], wspace=0.3)
ax_traj = fig_3a.add_subplot(gs[0, 0])
_trajectory(ax_traj, MELATONIN, highlight_bands=[(48, 96)], title=None)
ax_traj.set_title("Melatonin synthesis trajectory  (n=4)\n"
                   "48–96 hpf band: LD holds, BL drops at 72, D collapses at 96",
                   fontsize=13, loc="left", weight="bold")
ax_traj.set_ylabel("centered log$_2$|flux|")
ax_traj.legend(title="Regime", loc="upper right", fontsize=11, frameon=True)

# two stacked heatmaps: 72 hpf and 96 hpf
gs_r = gs[0, 1].subgridspec(2, 1, hspace=0.6)
ax_72 = fig_3a.add_subplot(gs_r[0, 0])
ax_96 = fig_3a.add_subplot(gs_r[1, 0])
val72, q72, n72 = _f1_at_stage([MELATONIN], 72)
val96, q96, n96 = _f1_at_stage([MELATONIN], 96)
_draw_stage_heatmap(ax_72, val72, q72, n72, title="72 hpf  —  BL drops below D and LD")
_draw_stage_heatmap(ax_96, val96, q96, n96, title="96 hpf  —  D drops below LD")
# shared colorbar
from matplotlib.cm import ScalarMappable
sm = ScalarMappable(norm=TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX), cmap=plt.cm.RdBu_r)
cbar = fig_3a.colorbar(sm, ax=[ax_72, ax_96], pad=0.02, fraction=0.05, shrink=0.8,
                        extend="both")
cbar.set_label("mean centered log$_2$FC\n(red = earlier regime higher)", fontsize=10)

fig_3a.suptitle("Figure 3a — Melatonin biosynthesis 48–96 hpf window: only LD holds it steady",
                 fontsize=14, weight="bold", y=0.99)
_save(fig_3a, "3a_melatonin_window")


# ============================================================
# 7a — aromatic AA flip
# ============================================================
print("\n[7a] aromatic 72 vs 120 flip")

fig_7a = plt.figure(figsize=(16, 6))
gs = fig_7a.add_gridspec(1, 2, width_ratios=[1.4, 1], wspace=0.3)

# two trajectories stacked on the left
gs_l = gs[0, 0].subgridspec(2, 1, hspace=0.5)
ax_aromatic = fig_7a.add_subplot(gs_l[0, 0])
ax_tyrosine = fig_7a.add_subplot(gs_l[1, 0], sharex=ax_aromatic)
_trajectory(ax_aromatic, AROMATIC, highlight_bands=[(68, 76), (116, 124)],
             title="Phe, Tyr and Trp biosynthesis")
_trajectory(ax_tyrosine, TYROSINE, highlight_bands=[(68, 76), (116, 124)],
             title="Tyrosine metabolism")
ax_tyrosine.set_xlabel("hpf")
ax_aromatic.set_ylabel("centered log$_2$|flux|")
ax_tyrosine.set_ylabel("centered log$_2$|flux|")
ax_aromatic.legend(title="Regime", loc="upper right", fontsize=10, frameon=True)

# two stacked heatmaps: 72 hpf and 120 hpf, showing the sign flip
gs_r = gs[0, 1].subgridspec(2, 1, hspace=0.6)
ax_72 = fig_7a.add_subplot(gs_r[0, 0])
ax_120 = fig_7a.add_subplot(gs_r[1, 0])
val72, q72, n72 = _f1_at_stage([AROMATIC, TYROSINE], 72)
val120, q120, n120 = _f1_at_stage([AROMATIC, TYROSINE], 120)
_draw_stage_heatmap(ax_72, val72, q72, n72, title="72 hpf  —  peak; D > BL")
_draw_stage_heatmap(ax_120, val120, q120, n120, title="120 hpf  —  flip; BL > D")
cbar = fig_7a.colorbar(sm, ax=[ax_72, ax_120], pad=0.02, fraction=0.05, shrink=0.8,
                        extend="both")
cbar.set_label("mean centered log$_2$FC\n(red = earlier regime higher)", fontsize=10)

fig_7a.suptitle("Figure 7a — Aromatic amino acid flip: D peaks at 72, BL never does → BL ends highest at 120",
                 fontsize=14, weight="bold", y=0.99)
_save(fig_7a, "7a_aromatic_flip")


# ============================================================
# 37 combined slide
# ============================================================
print("\n[37] combined slide")

fig_37 = plt.figure(figsize=(18, 11))
gs = fig_37.add_gridspec(2, 2, height_ratios=[1, 1.1], width_ratios=[1.4, 1],
                           hspace=0.45, wspace=0.3,
                           left=0.07, right=0.95, top=0.90, bottom=0.07)

# top-left: melatonin trajectory
ax_tl = fig_37.add_subplot(gs[0, 0])
_trajectory(ax_tl, MELATONIN, highlight_bands=[(48, 96)])
ax_tl.set_title("A  Melatonin biosynthesis trajectory  (n=4)",
                 fontsize=13, loc="left", weight="bold")
ax_tl.set_ylabel("centered log$_2$|flux|")
ax_tl.legend(title="Regime", loc="upper right", fontsize=10, frameon=True)

# top-right: melatonin 72/96 heatmap stack
gs_tr = gs[0, 1].subgridspec(2, 1, hspace=0.7)
ax_tr_72 = fig_37.add_subplot(gs_tr[0, 0])
ax_tr_96 = fig_37.add_subplot(gs_tr[1, 0])
val72, q72, n72 = _f1_at_stage([MELATONIN], 72)
val96, q96, n96 = _f1_at_stage([MELATONIN], 96)
_draw_stage_heatmap(ax_tr_72, val72, q72, n72, title="B  Melatonin at 72 hpf")
_draw_stage_heatmap(ax_tr_96, val96, q96, n96, title="    Melatonin at 96 hpf")

# bottom-left: aromatic trajectories (both)
gs_bl = gs[1, 0].subgridspec(2, 1, hspace=0.6)
ax_ar = fig_37.add_subplot(gs_bl[0, 0])
ax_ty = fig_37.add_subplot(gs_bl[1, 0], sharex=ax_ar)
_trajectory(ax_ar, AROMATIC, highlight_bands=[(68, 76), (116, 124)])
_trajectory(ax_ty, TYROSINE, highlight_bands=[(68, 76), (116, 124)])
ax_ar.set_title("C  Phe/Tyr/Trp biosynthesis", fontsize=12, loc="left", weight="bold")
ax_ty.set_title("    Tyrosine metabolism", fontsize=12, loc="left", weight="bold")
ax_ar.set_ylabel("log$_2$|flux|")
ax_ty.set_ylabel("log$_2$|flux|")
ax_ty.set_xlabel("hpf")
ax_ar.legend(title="Regime", loc="upper right", fontsize=9, frameon=True)

# bottom-right: aromatic 72 and 120 hpf heatmap stack
gs_br = gs[1, 1].subgridspec(2, 1, hspace=0.6)
ax_br_72 = fig_37.add_subplot(gs_br[0, 0])
ax_br_120 = fig_37.add_subplot(gs_br[1, 0])
val72, q72, n72 = _f1_at_stage([AROMATIC, TYROSINE], 72)
val120, q120, n120 = _f1_at_stage([AROMATIC, TYROSINE], 120)
_draw_stage_heatmap(ax_br_72, val72, q72, n72, title="D  72 hpf — peak;  D > BL")
_draw_stage_heatmap(ax_br_120, val120, q120, n120, title="    120 hpf — flipped;  BL > D")

cbar = fig_37.colorbar(sm, ax=[ax_tr_72, ax_tr_96, ax_br_72, ax_br_120],
                         pad=0.015, fraction=0.03, shrink=0.65, extend="both")
cbar.set_label("mean centered log$_2$FC", fontsize=10)

fig_37.suptitle("Figure 3+7 — Melatonin held by cycling only; aromatic AAs peak-and-fall in D/LD, flat in BL",
                 fontsize=15, weight="bold", y=0.97)
_save(fig_37, "37_combined")


print("\n[done]  outputs in", OUT)
