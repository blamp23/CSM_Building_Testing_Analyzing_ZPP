"""
s181 -- optional Figure 7: unknown-zone amplification confound.

The chapter's fourth methodological contribution has no figure anywhere.
This panel carries §1.3 and §5 visually.

Two panels:
  A — three-level zFPKM call per cell (expressed / unknown / off) as a
      stacked bar. 39–53 % of genes sit in the unknown zone per cell —
      the pool from which rFASTCORMICS makes a FASTCORE-connectivity
      call. gss shown as a per-cell point, since it is the headline
      example of a near-floor gene that resolves differently per cell.
  B — Zebrafish1 reactions: 60.9 % carry a GPR, 39.1 % are gene-less
      and extraction-pruned at FASTCORE's discretion. Companion to A.
      The two bars together show the pool out of which unknown-zone
      amplification arises.
"""
import sys
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
OUT = R / "chapter_figures"

COND_COLOR = {"BL": "#1F77B4", "D": "#111111", "LD": "#F0A800"}
HPFS = [24, 48, 72, 96, 120]
CELLS = [f"{c}_{h}" for c in ("BL", "D", "LD") for h in HPFS]

z = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv").set_index("ensembl_gene_id")
th = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv").set_index("cell")

GSS_ENS = "ENSDARG00000037706"

# compute three-level counts and gss state per cell
rows = []
for c in CELLS:
    if c not in z.columns or c not in th.index: continue
    vals = z[c].dropna()
    exp_t = th.loc[c, "exp_threshold_zFPKM"]
    inex_t = th.loc[c, "inex_threshold_zFPKM"]
    n_exp = int((vals >= exp_t).sum())
    n_off = int((vals <= inex_t).sum())
    n_unk = int(len(vals) - n_exp - n_off)
    gss_z = z.loc[GSS_ENS, c] if GSS_ENS in z.index else np.nan
    if pd.notna(gss_z) and gss_z >= exp_t:   gss_state = "expressed"
    elif pd.notna(gss_z) and gss_z <= inex_t: gss_state = "off"
    else:                                    gss_state = "unknown"
    rows.append(dict(cell=c, exp=n_exp, unk=n_unk, off=n_off,
                     gss_z=gss_z, gss_state=gss_state,
                     exp_t=exp_t, inex_t=inex_t))
df = pd.DataFrame(rows)
tot = df["exp"] + df["unk"] + df["off"]
df["exp_pct"] = df["exp"] / tot * 100
df["unk_pct"] = df["unk"] / tot * 100
df["off_pct"] = df["off"] / tot * 100

# ---------- Figure ----------
fig = plt.figure(figsize=(15, 5.5))
gs = fig.add_gridspec(1, 2, width_ratios=[2.0, 1.0], wspace=0.3)

# Panel A — three-level stack per cell
axA = fig.add_subplot(gs[0, 0])
xs = np.arange(len(df))
axA.bar(xs, df.exp_pct, color="#2E7D32", edgecolor="black", linewidth=0.4,
        label="expressed  (zFPKM ≥ exp threshold)")
axA.bar(xs, df.unk_pct, bottom=df.exp_pct, color="#F0A800",
        edgecolor="black", linewidth=0.4,
        label="unknown-zone  (inex < zFPKM < exp)")
axA.bar(xs, df.off_pct, bottom=df.exp_pct + df.unk_pct, color="#B00020",
        edgecolor="black", linewidth=0.4,
        label="off  (zFPKM ≤ inex threshold)")
# gss marker per cell: a small open circle at its state's segment
for i, (_, r) in enumerate(df.iterrows()):
    if r.gss_state == "expressed":
        y_gss = r.exp_pct / 2
    elif r.gss_state == "unknown":
        y_gss = r.exp_pct + r.unk_pct / 2
    else:
        y_gss = r.exp_pct + r.unk_pct + r.off_pct / 2
    axA.scatter(xs[i], y_gss, s=70, marker="o", facecolor="white",
                edgecolor="black", linewidths=1.2, zorder=5)
axA.set_xticks(xs); axA.set_xticklabels(df.cell, rotation=45, ha="right", fontsize=8)
axA.set_ylabel("% of genes with zFPKM in each tier", fontsize=9)
axA.set_ylim(0, 100)
axA.set_title("A — rFASTCORMICS three-level zFPKM call per cell\n"
              "white dot = gss (ENSDARG00000037706); unknown-zone bars are the\n"
              "pool from which unknown-zone amplification arises",
              fontsize=10, loc="left")
axA.legend(fontsize=8, loc="upper right", bbox_to_anchor=(1.0, 1.15), ncol=3,
           frameon=False)
# Annotation for 24 hpf unknown-zone bulge
axA.annotate("53 % unknown\nat 24 hpf", xy=(0, 48), xytext=(1.5, 70),
             fontsize=8, fontweight="bold", color="#B07000",
             arrowprops=dict(arrowstyle="->", color="#B07000", lw=0.8))

# Panel B — reaction-side gene-less fraction
axB = fig.add_subplot(gs[0, 1])
counts = {"with GPR\n(60.9 %)\n7,861 rxns": 60.9,
          "gene-less\n(39.1 %)\n5,048 rxns": 39.1}
labels = list(counts.keys())
values = list(counts.values())
colors = ["#2E7D32", "#808080"]
axB.bar([0, 1], values, color=colors, edgecolor="black", linewidth=0.6,
        width=0.6)
for i, (lab, v) in enumerate(zip(labels, values)):
    axB.text(i, v + 1.5, lab, ha="center", va="bottom", fontsize=9,
             fontweight="bold")
axB.set_xticks([0, 1]); axB.set_xticklabels([])
axB.set_ylabel("% of Zebrafish1 reactions", fontsize=9)
axB.set_ylim(0, 85)
axB.set_title("B — Zebrafish1 reaction GPR status\n"
              "(reactions prunable by FASTCORE without any gene call)",
              fontsize=10, loc="left")

fig.suptitle("Figure 7 — unknown-zone amplification: the pool and the mechanism  (§1.3 / §5)",
             fontsize=11, y=1.03)
fig.text(0.5, -0.08,
         "A: rFASTCORMICS discretises each gene's zFPKM per cell into one of three "
         "tiers; the unknown-zone fraction is ~39–53 % per cell and ~53 % at 24 hpf. "
         "Reactions whose GPR evaluates to 'unknown' are included or excluded by "
         "FASTCORE connectivity rather than by condition expression. gss (white dot) "
         "sits in the unknown band at hpf 24 for all three regimes and in the 'off' "
         "band at later hpf — the gss case documented in §5.1 is one draw from this "
         "pool. B: on the reaction side, 39.1 % of Zebrafish1 reactions have no gene "
         "annotation at all, so FASTCORE's include/exclude decision for them has no "
         "expression input — the second route by which capacity differences can arise "
         "from structural choices rather than from biology.",
         ha="center", fontsize=8, style="italic", wrap=True)
fig.savefig(OUT / "fig7_unknown_zone.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig7_unknown_zone.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig7_unknown_zone.pdf'} / .png")
