"""
s175 -- rebuilt §11.21 figure, promoting the unknown-zone finding to Panel A.

Main figure (3 panels):
  A. The unknown-zone amplification case (gss) — left: zFPKM across 15 cells
     with exp/inex threshold bands + retention markers; right: resulting GSH
     capacity gap (56 vs ~10 at BL_24 vs D_24 / LD_24).
  B. NADPH matched-floor capacity, pairwise, env_frac glyphs (clean envelope).
  C. Xenobiotic burden slope per hpf (null consistency with nitrogen-bound).

Supplement:
  S1. Subsystem presence heatmap (unchanged).
  S2. PAPS matched-floor capacity — only the 72 hpf column is a clean envelope
      finding; 96/120 BL-vs-LD mixed attribution is retracted (unknown-zone +
      geneless confound per s174).
"""
import sys
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import h5py

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "grace_v4/results/v7_susceptibility"
V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"

COND_COLOR = {"BL": "#1F77B4", "D": "#111111", "LD": "#F0A800"}
HPFS = [24, 48, 72, 96, 120]

p1 = pd.read_csv(OUT / "probe1_phase2_capacity.csv")
mf = p1[p1.regime == "matched_floor"].copy()
p3 = pd.read_csv(OUT / "probe3_xeno_burden.csv")
p4 = pd.read_csv(OUT / "probe4_presence.csv")
attr = pd.read_csv(OUT / "probe1_attribution_2x2.csv")
attr["pair_key"] = attr.apply(lambda r: tuple(sorted([r.cell_A, r.cell_B])), axis=1)

ZFPKM = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv").set_index("ensembl_gene_id")
THRESH = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv").set_index("cell")

CELLS_ALL = [f"{c}_{h}" for c in ("BL", "D", "LD") for h in HPFS]
GSS_ENS = "ENSDARG00000037706"


def _attr_glyph(hpf, cof, cell_A, cell_B):
    key = tuple(sorted([cell_A, cell_B]))
    sub = attr[(attr.hpf == hpf) & (attr.cofactor == cof) & (attr.pair_key == key)]
    if len(sub) == 0: return ""
    e = abs(sub.iloc[0].env_frac); t = abs(sub.iloc[0].topo_frac)
    if e > 0.7: return "E"
    if t > 0.7: return "T"
    return "M"


# which cells retained MAR04326 (GSS) ?
retained_gss = {}
for cell in CELLS_ALL:
    with h5py.File(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "r") as f:
        rxns_raw = np.asarray(f["parent/rxns"]).ravel()
        rxns = ["".join(chr(int(c)) for c in np.asarray(f[r]).ravel())
                if isinstance(r, h5py.Reference) else str(r) for r in rxns_raw]
    retained_gss[cell] = ("MAR04326" in rxns)


# =================== MAIN FIGURE ===================
fig = plt.figure(figsize=(16, 5.2))
gs = fig.add_gridspec(1, 4, width_ratios=[1.4, 1.0, 1.0, 1.0], wspace=0.35)

# ---------- Panel A (left sub-panel): gss zFPKM across cells ----------
axA1 = fig.add_subplot(gs[0, 0])
xs = np.arange(len(CELLS_ALL))
zs = [ZFPKM.loc[GSS_ENS, c] if c in ZFPKM.columns else np.nan for c in CELLS_ALL]
exp_thr = [THRESH.loc[c, "exp_threshold_zFPKM"] for c in CELLS_ALL]
inex_thr = [THRESH.loc[c, "inex_threshold_zFPKM"] for c in CELLS_ALL]
for i, c in enumerate(CELLS_ALL):
    axA1.add_patch(Rectangle((i - 0.42, inex_thr[i]), 0.84, exp_thr[i] - inex_thr[i],
                              facecolor="orange", alpha=0.17, edgecolor="none"))
for i, c in enumerate(CELLS_ALL):
    cond = c.split("_")[0]
    if retained_gss[c]:
        axA1.scatter(i, zs[i], color=COND_COLOR[cond], s=120, marker="o",
                     edgecolor="black", linewidth=1.4, zorder=5)
    else:
        axA1.scatter(i, zs[i], facecolor="white", edgecolor=COND_COLOR[cond],
                     s=120, marker="o", linewidth=1.4, zorder=5)
axA1.plot(xs, exp_thr, color="green", linewidth=1, linestyle="--", alpha=0.8)
axA1.plot(xs, inex_thr, color="red", linewidth=1, linestyle="--", alpha=0.8)
axA1.axhline(0, color="gray", linewidth=0.3)
axA1.set_xticks(xs); axA1.set_xticklabels(CELLS_ALL, rotation=45, ha="right", fontsize=8)
axA1.set_ylabel("gss zFPKM score", fontsize=10)
axA1.set_title("A — unknown-zone amplification\n"
               "gss (MAR04326) near-floor in every cell;\n"
               "band = unknown zone; filled = reaction retained",
               fontsize=10, loc="left")
axA1.grid(axis="y", alpha=0.3)
# Legend sub-panel
leg_handles = [
    plt.Line2D([0], [0], marker="o", color=COND_COLOR["BL"], linestyle="", markersize=9,
               markeredgecolor="black", markeredgewidth=1.4, label="BL retained"),
    plt.Line2D([0], [0], marker="o", color="white", linestyle="", markersize=9,
               markeredgecolor=COND_COLOR["BL"], markeredgewidth=1.4, label="BL pruned"),
    plt.Line2D([0], [0], linestyle="--", color="green", label="exp"),
    plt.Line2D([0], [0], linestyle="--", color="red", label="inex"),
    mpatches.Patch(color="orange", alpha=0.3, label="unknown zone"),
]
axA1.legend(handles=leg_handles, loc="lower right", fontsize=7, ncol=2, frameon=False)

# ---------- Panel A (right sub-panel): resulting GSH capacity gap at 24 hpf ----------
axA2 = fig.add_subplot(gs[0, 1])
# GSH matched_floor capacity for BL_24, D_24, LD_24 (mean over partners)
gsh24 = mf[(mf.cofactor == "GSH") & (mf.hpf == 24)]
means = gsh24.groupby("condition").sink_max.mean()
positions = [0, 1, 2]
conds = ["BL", "D", "LD"]
bars = axA2.bar(positions, [means.get(c, 0) for c in conds],
                color=[COND_COLOR[c] for c in conds], width=0.65,
                edgecolor="black", linewidth=0.8)
for i, c in enumerate(conds):
    v = means.get(c, 0)
    axA2.text(i, v + 1, f"{v:.1f}", ha="center", va="bottom", fontsize=10,
              fontweight="bold")
    # retained/pruned annotation
    cell_name = f"{c}_24"
    label = "RETAINED" if retained_gss[cell_name] else "PRUNED"
    axA2.text(i, 2, label, ha="center", va="center", fontsize=8,
              color="white" if c != "LD" else "black", rotation=90)
axA2.set_xticks(positions); axA2.set_xticklabels(conds)
axA2.set_ylabel("GSH sink capacity at 24 hpf\n(matched floor; relative units)", fontsize=9)
axA2.set_title("same confound, 6× capacity gap\n"
               "driven by MAR04326 retention,\nnot by gss expression",
               fontsize=10)
axA2.set_ylim(0, max(60, means.max() * 1.2))
axA2.grid(axis="y", alpha=0.3)

# ---------- Panel B: NADPH envelope ----------
axB = fig.add_subplot(gs[0, 2])
sub = mf[mf.cofactor == "NADPH"]
for cond, color in COND_COLOR.items():
    for hpf in HPFS:
        cells_cond = sub[(sub.hpf == hpf) & (sub.condition == cond)]
        offset = {"BL": -3, "D": 0, "LD": +3}[cond]
        xs = [hpf + offset] * len(cells_cond)
        axB.scatter(xs, cells_cond.sink_max, color=color, alpha=0.4, s=18, edgecolors='none')
        if len(cells_cond):
            axB.scatter([hpf + offset], [cells_cond.sink_max.mean()],
                        color=color, s=70, edgecolors='white', linewidths=1.2, zorder=5)
piv = sub.copy()
piv["pair_key"] = piv.apply(lambda r: tuple(sorted([r.cell, r.partner])), axis=1)
by = piv.groupby(["hpf", "pair_key"])["sink_max"].agg(["min", "max"]).reset_index()
by["pct"] = (by["max"] - by["min"]) / by["min"].replace(0, 1e-9) * 100
for _, r in by[by.pct > 10].iterrows():
    a, b = r.pair_key
    g = _attr_glyph(r.hpf, "NADPH", a, b)
    if g:
        axB.text(r.hpf, r["max"] * 1.06, g, ha="center", va="bottom", fontsize=9,
                 fontweight="bold",
                 color={"E": "#2E7D32", "M": "#B07000", "T": "#B00020"}.get(g, "black"))
axB.set_xticks(HPFS); axB.set_xlabel("hpf")
axB.set_ylabel("NADPH sink capacity\n(matched floor; relative units)", fontsize=9)
axB.set_title("B — NADPH, envelope-driven\n"
              "(E = env_frac ≥ 0.7 on >10% gap)", fontsize=10)
axB.grid(axis="y", alpha=0.3)

# ---------- Panel C: xeno burden slope ----------
axC = fig.add_subplot(gs[0, 3])
slope_df = p3[p3.u.isna()][["condition", "hpf", "variant", "slope"]]
for cond, color in COND_COLOR.items():
    sub2 = slope_df[(slope_df.condition == cond) & (slope_df.variant == "Phase2_glucuronide")]
    agg = sub2.groupby("hpf").slope.mean().reindex(HPFS)
    axC.plot(HPFS, agg.values, marker="o", color=color, label=cond,
             linewidth=2, markersize=7)
axC.set_xticks(HPFS); axC.set_xlabel("hpf")
axC.set_ylabel("biomass cost per unit X", fontsize=9)
axC.set_title("C — xenobiotic burden slope\n"
              "(null at 96/120: nitrogen-bound)", fontsize=10)
axC.axhline(0, color="gray", linewidth=0.6)
axC.legend(fontsize=8, loc="upper right")
axC.grid(axis="y", alpha=0.3)

fig.suptitle("Fig. v7 susceptibility probes — the methodological finding (Panel A), one clean envelope result (Panel B), one null consistency check (Panel C)",
             fontsize=10.5, y=1.02)
fig.text(0.5, -0.04,
         "Relative units only. Within-hpf between-condition comparisons are interpretable. "
         "GSH @ 24 and 120, UDP-GlcA @ 72–120, PAPS @ 96–120 BL-vs-LD attributions had topology components; "
         "GPR audits (gsh_pathway_audit.csv, udp_glca_pathway_audit.csv, paps_pathway_audit.csv) classified those as "
         "unknown-zone amplification (same class as the gss case in A), so they are retracted as biology. See §11.21.5.",
         ha="center", fontsize=8, style="italic", wrap=True)

fig.savefig(OUT / "fig_susceptibility.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig_susceptibility.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig_susceptibility.pdf'} / .png")
plt.close(fig)


# =================== SUPPLEMENT S1 unchanged — presence heatmap ===================
figS1, ax = plt.subplots(figsize=(10, 5.5))
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
                color="white" if mat[i, j] < mat.max() * 0.6 else "black", fontsize=8)
plt.colorbar(im, ax=ax, label="n reactions carried into cell (all optional)")
ax.set_title("S1 — stress-response subsystem presence (reactions per cell × subsystem)\n"
             "all 1217 (rxn × cell) are OPTIONAL at Frame A — no stress reaction is biomass-required",
             fontsize=10)
figS1.savefig(OUT / "fig_S1_presence.pdf", bbox_inches="tight")
figS1.savefig(OUT / "fig_S1_presence.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig_S1_presence.pdf'} / .png")
plt.close(figS1)


# =================== SUPPLEMENT S2 — PAPS with proper confound annotation ===================
figS2, ax = plt.subplots(figsize=(10, 4.5))
sub = mf[mf.cofactor == "PAPS"]
for cond, color in COND_COLOR.items():
    for hpf in HPFS:
        cells_cond = sub[(sub.hpf == hpf) & (sub.condition == cond)]
        offset = {"BL": -3, "D": 0, "LD": +3}[cond]
        xs = [hpf + offset] * len(cells_cond)
        ax.scatter(xs, cells_cond.sink_max, color=color, alpha=0.4, s=18, edgecolors='none')
        if len(cells_cond):
            ax.scatter([hpf + offset], [cells_cond.sink_max.mean()],
                       color=color, s=70, edgecolors='white', linewidths=1.2, zorder=5)
piv = sub.copy()
piv["pair_key"] = piv.apply(lambda r: tuple(sorted([r.cell, r.partner])), axis=1)
by = piv.groupby(["hpf", "pair_key"])["sink_max"].agg(["min", "max"]).reset_index()
by["pct"] = (by["max"] - by["min"]) / by["min"].replace(0, 1e-9) * 100
# audited unknown-zone-confound pairs (BL vs LD at 96 and 120)
confound_pairs = {(96, tuple(sorted(("BL_96", "LD_96")))),
                  (120, tuple(sorted(("BL_120", "LD_120"))))}
for _, r in by[by.pct > 10].iterrows():
    a, b = r.pair_key
    g = _attr_glyph(r.hpf, "PAPS", a, b)
    if g:
        is_confound = (r.hpf, r.pair_key) in confound_pairs
        label = f"{g}✗" if is_confound else g
        color = "#B00020" if is_confound else {"E": "#2E7D32", "M": "#B07000"}.get(g, "black")
        ax.text(r.hpf, r["max"] * 1.06, label, ha="center", va="bottom", fontsize=9,
                fontweight="bold", color=color)
ax.set_xticks(HPFS); ax.set_xlabel("hpf")
ax.set_ylabel("PAPS sink capacity (matched floor; relative units)")
ax.set_title("S2 — PAPS: clean envelope at 72 hpf (D ≈ 1.6× BL/LD); "
             "BL-vs-LD at 96/120 (M✗) is unknown-zone confound per audit, retracted as biology",
             fontsize=9.5)
ax.grid(axis="y", alpha=0.3)
figS2.savefig(OUT / "fig_S2_paps.pdf", bbox_inches="tight")
figS2.savefig(OUT / "fig_S2_paps.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig_S2_paps.pdf'} / .png")
plt.close(figS2)
