"""
s191 -- figures for the reaction-set footprint.

Figure 8 — stage comparison. 2 rows × 5 columns.
  Row A: 3-set Venn per hpf (BL / D / LD). Each unique-to-one region
         coloured by flux-carrying fraction (viridis), shared-by-two in
         grey, shared-by-all in light grey.
  Row B: stacked GPR-class bar beneath each Venn, for ALL non-shared
         reactions at that hpf.

Figure 9 — regime developmental sequence. 3 rows × 4 columns.
  Each cell: 2-set Venn of adjacent-stage reactions within one regime.
  Shaded by direction (a_only / b_only). Companion stacked bar of GPR
  class for the non-shared reactions below each row.
"""
import sys
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn3, venn2
import matplotlib.cm as cm

REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
IN = R / "v7_final_addenda/footprint"
OUT = R / "chapter_figures"
OUT.mkdir(parents=True, exist_ok=True)

HPFS = [24, 48, 72, 96, 120]
CONDS = ("BL", "D", "LD")

GPR_CLASS_COLOR = {
    "GENELESS": "#808080",
    "UNMAPPABLE": "#C0C0C0",
    "EXPRESSED": "#2E7D32",
    "UNKNOWN_ZONE": "#F0A800",
    "UNKNOWN_ZONE_mixed": "#CC7700",
    "OFF": "#B00020",
    "": "#DDDDDD",
}

# ============== FIGURE 8 — stage comparison =====================
fig8 = plt.figure(figsize=(21, 7))
gs = fig8.add_gridspec(2, 5, hspace=0.3, wspace=0.25,
                       height_ratios=[2.3, 1.0])

for i, hpf in enumerate(HPFS):
    axV = fig8.add_subplot(gs[0, i])
    df = pd.read_csv(IN / f"stage_{hpf}_three_way.csv")
    # Rebuild the three sets from the holders column
    sets = {c: set(df[df.holders.str.contains(f"{c}_{hpf}", na=False)].rxn)
            for c in CONDS}
    v = venn3(subsets=[sets["BL"], sets["D"], sets["LD"]],
              set_labels=(f"BL_{hpf}", f"D_{hpf}", f"LD_{hpf}"),
              ax=axV,
              layout_algorithm=None,
              set_colors=("#1F77B4", "#111111", "#F0A800"),
              alpha=0.45)
    # Colour the 3 unique-to-one regions by flux-carrying fraction
    unique_ids = {"BL": "100", "D": "010", "LD": "001"}
    for cond, region_id in unique_ids.items():
        patch = v.get_patch_by_id(region_id)
        if patch is None: continue
        subdf = df[(df.tier == "unique_to_one") &
                   df.holders.str.contains(f"{cond}_{hpf}", na=False)]
        if len(subdf):
            frac = subdf.flux_carrying.mean()
            # Overlay with intensity-shaded patch (viridis)
            patch.set_facecolor(cm.viridis(frac))
            patch.set_alpha(0.55)
        lbl = v.get_label_by_id(region_id)
        if lbl:
            n = (df.tier == "unique_to_one") & df.holders.str.contains(f"{cond}_{hpf}", na=False)
            lbl.set_text(f"{n.sum()}")
    axV.set_title(f"hpf {hpf}", fontsize=12, loc="center")

    # --- Stacked GPR bar underneath this Venn ---
    axB = fig8.add_subplot(gs[1, i])
    non_shared = df[df.tier != "shared_by_all"]
    counts = non_shared.gpr_class.fillna("").value_counts()
    total = counts.sum()
    bottom = 0
    for cls in ("EXPRESSED", "UNKNOWN_ZONE", "UNKNOWN_ZONE_mixed", "OFF", "GENELESS", "UNMAPPABLE"):
        n = counts.get(cls, 0)
        if n == 0: continue
        pct = 100 * n / total if total else 0
        color = GPR_CLASS_COLOR.get(cls, "#DDDDDD")
        axB.bar(0, pct, bottom=bottom, color=color, edgecolor="black",
                linewidth=0.3, width=0.9)
        if pct > 6:
            axB.text(0, bottom + pct / 2, f"{cls[:3]}\n{n} ({pct:.0f}%)",
                     ha="center", va="center", fontsize=7,
                     color="white" if color in ("#808080","#2E7D32","#B00020","#111111") else "black")
        bottom += pct
    axB.set_xlim(-0.6, 0.6); axB.set_ylim(0, 100)
    axB.set_xticks([]); axB.set_yticks([0, 50, 100])
    axB.set_ylabel("non-shared rxns (%)", fontsize=8)
    axB.set_title(f"GPR class  n={total}", fontsize=8)
    for spine in ("top", "right"): axB.spines[spine].set_visible(False)

# Global legend
leg_handles = [
    mpatches.Patch(facecolor=cm.viridis(0.0), label="0% flux-carrying"),
    mpatches.Patch(facecolor=cm.viridis(0.5), label="50% flux-carrying"),
    mpatches.Patch(facecolor=cm.viridis(1.0), label="100% flux-carrying"),
]
leg_handles += [mpatches.Patch(color=c, label=k)
                for k, c in GPR_CLASS_COLOR.items() if k
                and k in ("EXPRESSED","UNKNOWN_ZONE","UNKNOWN_ZONE_mixed","OFF","GENELESS","UNMAPPABLE")]
fig8.legend(handles=leg_handles, loc="lower center", ncol=9, fontsize=8,
            bbox_to_anchor=(0.5, -0.04), frameon=False)
fig8.suptitle("Figure 8 — reaction-set footprint across BL / D / LD at each stage",
              fontsize=12, y=1.02)
fig8.text(0.5, -0.09,
          "Top: three-set Venn per hpf. Unique-to-one regions coloured by the "
          "fraction of those reactions that carry non-zero flux in the owning cell "
          "(viridis). Shared-by-two regions and the central shared-by-all region "
          "are left at default colour. "
          "Bottom: stacked bar of GPR class for all non-shared reactions at that hpf. "
          "Expressed = ≥1 mapped gene with zFPKM ≥ exp threshold in the owning cell. "
          "Non-shared reaction counts range 346–571 per stage; GENELESS fraction is "
          "~30 % everywhere.",
          ha="center", fontsize=8, style="italic", wrap=True)
fig8.savefig(OUT / "fig8_footprint_stage.pdf", bbox_inches="tight")
fig8.savefig(OUT / "fig8_footprint_stage.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig8_footprint_stage.pdf'}")
plt.close(fig8)


# ============== FIGURE 9 — regime adjacent-stage sequence =========
fig9 = plt.figure(figsize=(17, 10))
gs = fig9.add_gridspec(6, 4, hspace=0.4, wspace=0.25,
                       height_ratios=[2.0, 0.8, 2.0, 0.8, 2.0, 0.8])

for row_i, cond in enumerate(CONDS):
    df = pd.read_csv(IN / f"regime_{cond}_adjacent.csv")
    df_rxn = df[df.rxn.notna()]
    for col_i, (h1, h2) in enumerate([(24,48),(48,72),(72,96),(96,120)]):
        axV = fig9.add_subplot(gs[row_i * 2, col_i])
        pair = f"{h1}→{h2}"
        pair_rxns = df_rxn[df_rxn.pair == pair]
        a_only = set(pair_rxns[pair_rxns.direction=="a_only"].rxn)
        b_only = set(pair_rxns[pair_rxns.direction=="b_only"].rxn)
        # Compute shared count from summary row (no rxn-level record)
        header = df[(df.pair==pair) & df.rxn.isna()].iloc[0]
        shared_n = int(header.shared)
        # venn2 needs pretend sets
        v = venn2(subsets=(len(a_only) + shared_n,
                            len(b_only) + shared_n,
                            shared_n),
                   set_labels=(f"{cond}_{h1}", f"{cond}_{h2}"),
                   ax=axV, layout_algorithm=None,
                   set_colors=("#AAAAAA", "#AAAAAA"))
        # Flux-carrying overlays
        for region_id, rxn_set in [("10", a_only), ("01", b_only)]:
            patch = v.get_patch_by_id(region_id)
            if patch is None: continue
            if not rxn_set:
                continue
            frac_flux = pair_rxns[pair_rxns.rxn.isin(rxn_set)].flux_carrying.mean()
            patch.set_facecolor(cm.viridis(frac_flux))
            patch.set_alpha(0.55)
        axV.set_title(f"{cond}  {h1}→{h2}", fontsize=10)

        # Stacked GPR bar beneath (shared across row via subplot grid)
        axB = fig9.add_subplot(gs[row_i * 2 + 1, col_i])
        non_shared = pair_rxns
        counts = non_shared.gpr_class.fillna("").value_counts()
        total = counts.sum()
        bottom = 0
        for cls in ("EXPRESSED", "UNKNOWN_ZONE", "UNKNOWN_ZONE_mixed", "OFF", "GENELESS", "UNMAPPABLE"):
            n = counts.get(cls, 0)
            if n == 0: continue
            pct = 100 * n / total if total else 0
            color = GPR_CLASS_COLOR.get(cls, "#DDDDDD")
            axB.bar(0, pct, bottom=bottom, color=color, edgecolor="black",
                    linewidth=0.3, width=0.9)
            if pct > 10:
                axB.text(0, bottom + pct / 2, f"{cls[:3]}",
                         ha="center", va="center", fontsize=6.5,
                         color="white" if color in ("#808080","#2E7D32","#B00020","#111111") else "black")
            bottom += pct
        axB.set_xlim(-0.6, 0.6); axB.set_ylim(0, 100)
        axB.set_xticks([]); axB.set_yticks([0, 50, 100])
        axB.tick_params(axis="y", labelsize=6)
        axB.set_title(f"n={total}", fontsize=7)
        for spine in ("top", "right"): axB.spines[spine].set_visible(False)
    # Row header
    fig9.text(0.01, 1.0 - (row_i * 2 + 1) / 6.5 + 0.03, cond,
              fontsize=14, fontweight="bold", va="center",
              color={"BL": "#1F77B4", "D": "#111111", "LD": "#F0A800"}[cond])

leg_handles = [
    mpatches.Patch(facecolor=cm.viridis(0.0), label="0% flux-carrying"),
    mpatches.Patch(facecolor=cm.viridis(0.5), label="50% flux-carrying"),
    mpatches.Patch(facecolor=cm.viridis(1.0), label="100% flux-carrying"),
]
leg_handles += [mpatches.Patch(color=GPR_CLASS_COLOR[k], label=k)
                for k in ("EXPRESSED","UNKNOWN_ZONE","UNKNOWN_ZONE_mixed","OFF","GENELESS","UNMAPPABLE")]
fig9.legend(handles=leg_handles, loc="lower center", ncol=9, fontsize=8,
            bbox_to_anchor=(0.5, -0.02), frameon=False)
fig9.suptitle("Figure 9 — reaction-set turnover across adjacent stages within each regime",
              fontsize=12, y=1.02)
fig9.text(0.5, -0.055,
          "Each panel is a two-set Venn of extracted reactions at consecutive hpf "
          "within one regime. a_only and b_only regions coloured by flux-carrying "
          "fraction. Companion stacked bar below each Venn gives GPR class breakdown "
          "of the non-shared reactions.",
          ha="center", fontsize=8, style="italic", wrap=True)
fig9.savefig(OUT / "fig9_footprint_regime.pdf", bbox_inches="tight")
fig9.savefig(OUT / "fig9_footprint_regime.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig9_footprint_regime.pdf'}")
plt.close(fig9)
