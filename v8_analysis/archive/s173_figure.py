"""
s173 -- 3-panel §11.21 figure + 2 supplementary panels.

Panels:
  A. NADPH matched-floor capacity, hpf × condition (pairwise), attribution glyphs
  B. PAPS matched-floor capacity, hpf × condition (pairwise), attribution glyphs
  C. Xenobiotic burden slope per hpf (null consistency)

Supplement:
  S1. Subsystem presence heatmap (Probe 4) — 15 cells × 4 subsystems summary
  S2. gss zFPKM across 15 cells with exp/inex threshold bands
       + retained/pruned MAR04326 per cell

Colours: BL #1F77B4, D #111111, LD #F0A800.
"""
import sys
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "grace_v4/results/v7_susceptibility"

COND_COLOR = {"BL": "#1F77B4", "D": "#111111", "LD": "#F0A800"}
HPFS = [24, 48, 72, 96, 120]

p1 = pd.read_csv(OUT / "probe1_phase2_capacity.csv")
mf = p1[p1.regime == "matched_floor"].copy()

p3 = pd.read_csv(OUT / "probe3_xeno_burden.csv")

p4 = pd.read_csv(OUT / "probe4_presence.csv")

attr = pd.read_csv(OUT / "probe1_attribution_2x2.csv")
# A pair-key for attribution lookup
attr["pair_key"] = attr.apply(lambda r: tuple(sorted([r.cell_A, r.cell_B])), axis=1)


def _attr_glyph(hpf, cof, cell_A, cell_B):
    key = tuple(sorted([cell_A, cell_B]))
    sub = attr[(attr.hpf == hpf) & (attr.cofactor == cof) & (attr.pair_key == key)]
    if len(sub) == 0: return ""
    e = abs(sub.iloc[0].env_frac); t = abs(sub.iloc[0].topo_frac)
    if e > 0.7: return "E"       # envelope-dominant
    if t > 0.7: return "T"       # topology (flagged limitation)
    return "M"                   # mixed


# =================== FIGURE 1 (A + B + C) ===================
fig = plt.figure(figsize=(14, 4.8))
gs = fig.add_gridspec(1, 3, wspace=0.3)

def _panel_capacity(ax, cof, title):
    """Per-cell capacity, matched-floor pairwise. Three condition markers per hpf.
    For each hpf and each cell, plot ALL its pairwise matched-floor sinks as small
    dots, and the mean as a filled circle, in that cell's condition colour."""
    sub = mf[mf.cofactor == cof]
    for cond, color in COND_COLOR.items():
        for hpf in HPFS:
            cells_cond = sub[(sub.hpf == hpf) & (sub.condition == cond)]
            # Jitter x
            xs = [hpf + {"BL": -3, "D": 0, "LD": +3}[cond]] * len(cells_cond)
            ax.scatter(xs, cells_cond.sink_max, color=color, alpha=0.4, s=18, edgecolors='none')
            # Mean
            if len(cells_cond):
                ax.scatter([hpf + {"BL": -3, "D": 0, "LD": +3}[cond]],
                           [cells_cond.sink_max.mean()],
                           color=color, s=70, edgecolors='white', linewidths=1.2, zorder=5)
    # Attribution glyphs: pairs with >10% gap
    piv = sub.copy()
    piv["pair_key"] = piv.apply(lambda r: tuple(sorted([r.cell, r.partner])), axis=1)
    by = piv.groupby(["hpf", "pair_key"])["sink_max"].agg(["min", "max"]).reset_index()
    by["pct"] = (by["max"] - by["min"]) / by["min"].replace(0, 1e-9) * 100
    for _, r in by[by.pct > 10].iterrows():
        a, b = r.pair_key
        g = _attr_glyph(r.hpf, cof, a, b)
        if g:
            ax.text(r.hpf, r["max"] * 1.06, g, ha="center", va="bottom", fontsize=9,
                    fontweight="bold",
                    color={"E": "#2E7D32", "M": "#B07000", "T": "#B00020"}.get(g, "black"))
    ax.set_xticks(HPFS); ax.set_xlabel("hpf")
    ax.set_ylabel(f"{cof} sink capacity (rel. units)")
    ax.set_title(title, fontsize=11)
    ax.grid(axis="y", alpha=0.3)

# Panel A — NADPH
_panel_capacity(fig.add_subplot(gs[0, 0]), "NADPH",
                "A — NADPH supply, matched-floor (envelope E = green)")
# Panel B — PAPS
_panel_capacity(fig.add_subplot(gs[0, 1]), "PAPS",
                "B — PAPS supply, matched-floor (M = mixed)")

# Panel C — Probe 3 slope
axC = fig.add_subplot(gs[0, 2])
slope_df = p3[p3.u.isna()][["cell", "condition", "hpf", "variant", "slope"]]
for cond, color in COND_COLOR.items():
    sub = slope_df[(slope_df.condition == cond) & (slope_df.variant == "Phase2_glucuronide")]
    agg = sub.groupby("hpf").slope.mean().reindex(HPFS)
    axC.plot(HPFS, agg.values, marker="o", color=color, label=cond, linewidth=2, markersize=7)
axC.set_xticks(HPFS); axC.set_xlabel("hpf")
axC.set_ylabel("biomass cost per unit X")
axC.set_title("C — xeno burden slope (null at 96/120; nitrogen-bound)", fontsize=11)
axC.axhline(0, color="gray", linewidth=0.6)
axC.legend(fontsize=9, loc="upper right")
axC.grid(axis="y", alpha=0.3)

# Legend for A/B markers
legend_markers = [plt.Line2D([0], [0], marker="o", color=c, label=k,
                             markerfacecolor=c, markersize=9, linestyle="")
                  for k, c in COND_COLOR.items()]
fig.legend(legend_markers, list(COND_COLOR.keys()), loc="lower center",
           ncol=3, bbox_to_anchor=(0.5, -0.02), fontsize=10, frameon=False)

fig.suptitle("v7 susceptibility probes — matched-floor pairwise (relative units); "
             "glyphs: E=envelope-driven, M=mixed",
             fontsize=11, y=1.02)
fig.text(0.5, -0.06, "Relative units. Only within-hpf between-condition comparisons are interpretable.",
         ha="center", fontsize=9, style="italic")

fig.savefig(OUT / "fig_susceptibility.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig_susceptibility.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig_susceptibility.pdf'} / .png")
plt.close(fig)


# =================== SUPPLEMENT S1 — presence heatmap ===================
figS1, ax = plt.subplots(figsize=(10, 5.5))
# Count forced/optional/blocked per subsystem per cell
subs = ["Drug metabolism", "Xenobiotics metabolism",
        "Glutathione metabolism", "ROS detoxification"]
cells = [f"{c}_{h}" for c in ("BL", "D", "LD") for h in HPFS]
mat = np.zeros((len(cells), len(subs)))
for i, c in enumerate(cells):
    for j, s in enumerate(subs):
        sub = p4[(p4.cell == c) & (p4.subsystem == s)]
        mat[i, j] = len(sub)
im = ax.imshow(mat, cmap="viridis", aspect="auto")
ax.set_xticks(range(len(subs))); ax.set_xticklabels(subs, rotation=25, ha="right")
ax.set_yticks(range(len(cells))); ax.set_yticklabels(cells)
for i in range(len(cells)):
    for j in range(len(subs)):
        ax.text(j, i, int(mat[i, j]), ha="center", va="center",
                color="white" if mat[i, j] < mat.max()*0.6 else "black", fontsize=8)
plt.colorbar(im, ax=ax, label="n reactions carried into cell (all optional)")
ax.set_title("S1 — stress-response subsystem presence (reactions per cell × subsystem)\n"
             "all 1217 (rxn × cell) are OPTIONAL at Frame A — no stress reaction is biomass-required",
             fontsize=10)
figS1.savefig(OUT / "fig_S1_presence.pdf", bbox_inches="tight")
figS1.savefig(OUT / "fig_S1_presence.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig_S1_presence.pdf'} / .png")
plt.close(figS1)


# =================== SUPPLEMENT S2 — gss unknown-zone schematic ===================
# Load zFPKM + thresholds for gss (ENSDARG00000037706)
ZFPKM = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv").set_index("ensembl_gene_id")
THRESH = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv").set_index("cell")

gss_ens = "ENSDARG00000037706"
cells_order = [f"{c}_{h}" for c in ("BL", "D", "LD") for h in HPFS]
# Which cells retained MAR04326 in extraction? Reload from v7_models rxns lists.
import h5py
retained = {}
V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
for cell in cells_order:
    fp = V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat"
    with h5py.File(fp, "r") as f:
        rxns_raw = np.asarray(f["parent/rxns"]).ravel()
        rxns = ["".join(chr(int(c)) for c in np.asarray(f[r]).ravel()) if isinstance(r, h5py.Reference) else str(r) for r in rxns_raw]
    retained[cell] = ("MAR04326" in rxns)

figS2, ax = plt.subplots(figsize=(11, 5))
xs = np.arange(len(cells_order))
zs = [ZFPKM.loc[gss_ens, c] if c in ZFPKM.columns else np.nan for c in cells_order]
exp_thr = [THRESH.loc[c, "exp_threshold_zFPKM"] for c in cells_order]
inex_thr = [THRESH.loc[c, "inex_threshold_zFPKM"] for c in cells_order]
# Shade unknown zone per cell
for i, c in enumerate(cells_order):
    ax.add_patch(Rectangle((i-0.4, inex_thr[i]), 0.8, exp_thr[i]-inex_thr[i],
                           alpha=0.15, color="orange", edgecolor="none"))
# Plot gss values coloured by condition
for i, c in enumerate(cells_order):
    cond = c.split("_")[0]
    marker = "o" if retained[c] else "x"
    ax.scatter(i, zs[i], color=COND_COLOR[cond], s=120, marker=marker,
               edgecolor="black" if retained[c] else "none", linewidth=1.5, zorder=5)
# Threshold lines
ax.plot(xs, exp_thr, color="green", linewidth=1, linestyle="--", label="exp threshold")
ax.plot(xs, inex_thr, color="red", linewidth=1, linestyle="--", label="inex threshold")
ax.axhline(0, color="gray", linewidth=0.4)

ax.set_xticks(xs); ax.set_xticklabels(cells_order, rotation=40, ha="right")
ax.set_ylabel("zFPKM score")
ax.set_title("S2 — gss (glutathione synthetase, MAR04326) is near-floor in every v7 cell\n"
             "orange band = rFASTCORMICS 'unknown' zone (inex..exp); "
             "circles = reaction retained; ×= pruned", fontsize=10)
# Legend
import matplotlib.patches as mpatches
legend_handles = [plt.Line2D([0], [0], marker="o", color=COND_COLOR[k], linestyle="",
                             markersize=9, label=k) for k in ("BL","D","LD")]
legend_handles += [
    plt.Line2D([0], [0], marker="o", color="gray", linestyle="", markersize=9,
               markeredgecolor="black", markeredgewidth=1.5, label="MAR04326 retained"),
    plt.Line2D([0], [0], marker="x", color="gray", linestyle="", markersize=9,
               label="MAR04326 pruned"),
    mpatches.Patch(color="orange", alpha=0.3, label="unknown zone"),
]
ax.legend(handles=legend_handles, loc="lower right", fontsize=9, ncol=2)
ax.grid(axis="y", alpha=0.3)

figS2.text(0.02, -0.03,
           "gss sits in the unknown zone at hpf 24 in all three conditions and below the inex threshold at hpf 48–120.\n"
           "BL_24 zFPKM (−1.75) is LOWER than LD_24 (−1.17) and D_24 (−1.20), yet rFASTCORMICS retained MAR04326 in BL_24 and pruned it in D_24 / LD_24.\n"
           "The sign is wrong for a transcriptional story. The downstream 6× GSH capacity gap is unknown-zone amplification, not biology.",
           fontsize=9, style="italic", wrap=True)

figS2.savefig(OUT / "fig_S2_unknown_zone.pdf", bbox_inches="tight")
figS2.savefig(OUT / "fig_S2_unknown_zone.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig_S2_unknown_zone.pdf'} / .png")
plt.close(figS2)
