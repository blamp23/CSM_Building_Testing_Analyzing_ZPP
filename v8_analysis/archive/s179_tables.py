"""
s179 -- build the six tables for the chapter.

Writes grace_v4/results/chapter_figures/tables.md with:

  MAIN TEXT
    Table 1 — Zebrafish1 as used (12,909 rxns / 8,344 mets / 9 compartments)
    Table 5 — Tier 0 2×2 attribution (6 columns incl. interaction)

  SUPPLEMENTARY
    Table S2 — 50 Tier A rows (compact)
    Table S3 — parent-vs-child growth capacity per cell
    Table S4 — Tier 0 robustness (10% uptake relaxation)
    Table S6 — Tier A rolled up by subsystem

Also writes each table as a standalone CSV so pandocable.
"""
from pathlib import Path
import pandas as pd, numpy as np

REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
OUT = R / "chapter_figures"
OUT.mkdir(parents=True, exist_ok=True)
TABLES_DIR = OUT / "tables_csv"
TABLES_DIR.mkdir(parents=True, exist_ok=True)


def _md_table(df, headers=None, aligns=None):
    """Render a pandas DataFrame as a GitHub-flavoured markdown table."""
    cols = list(df.columns) if headers is None else headers
    if aligns is None:
        aligns = ["---"] * len(cols)
    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join(aligns) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            if isinstance(v, float):
                if np.isnan(v): cells.append("—")
                elif abs(v) >= 1e4 or (0 < abs(v) < 0.01):
                    cells.append(f"{v:.2e}")
                else: cells.append(f"{v:g}")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


# -------------- TABLE 1 -------------------------------------------------------
t1 = pd.DataFrame([
    ("Reactions",                "12,909",       "`len(m.rxns)` on baked_ocr_24_v7.mat"),
    ("Metabolites",              "8,344",        "`len(m.mets)`"),
    ("Compartments",             "9",            "e, x, m, c, l, r, g, n, i"),
    ("Reactions with ≥1 gene",   "7,861  (60.9 %)", "regex over m.grRules"),
    ("Gene-less reactions",      "5,048  (39.1 %)", "transporters, exchanges, templates"),
    ("Unique gene symbols",      "2,697",        "across all GPRs"),
    ("Symbol → Ensembl map rows", "5,268",       "`results/qc/gene_map_ensdarg_to_model.csv`"),
    ("Model symbols mapping",    "2,592 / 2,697 (96.1 %)", "to ≥1 Ensembl ID"),
], columns=["quantity", "value", "source"])


# -------------- TABLE 5 (expanded) --------------------------------------------
# 6 columns: pair, total gap (raw AA−BB), topology, envelope, interaction, verdict
# From hard-coded Table 5 ground truth (same values Fig 5B shows); plus the 9
# "no material gap" pairs for completeness.
# The 9 non-material pairs we read from the C3 CSV with interaction = raw-gap - topo - env.
c3 = pd.read_csv(R / "atp_sweep_v7/tier0_C3_2x2_before_after.csv")
# Build a complete 15-pair Table 5
t5_rows = []
# Values from the chapter table (ground truth) for the 6 material-gap pairs:
material_shares = {
    "BL_96_vs_LD_96":  (-4.39, -0.89, +0.08, -3.59, "not separable"),
    "BL_96_vs_D_96":   (-4.04, -0.89, -0.72, -2.43, "not separable"),
    "BL_120_vs_LD_120":(-0.84, -0.04, -0.80, +0.00, "envelope"),
    "D_120_vs_LD_120": (-0.73,  0.00, -0.73, +0.00, "envelope"),
    "BL_72_vs_D_72":   (+1.57,  0.00, +1.57, +0.00, "envelope"),
    "D_72_vs_LD_72":   (-1.35,  0.00, -1.35, +0.00, "envelope"),
}
for _, r in c3.iterrows():
    key = r.pair
    pair_display = r.pair.replace("_vs_", " – ")
    if key in material_shares:
        gap, topo, env, inter, verdict = material_shares[key]
    else:
        gap = r.total_gap_new
        topo = r.topology_new
        env = r.envelope_new
        # interaction consistent with AA-BB = topo+env+interaction:
        inter = r.total_gap_new - topo - env
        if abs(gap) < 0.1:   # "no material gap" threshold
            verdict = "no material gap"
        else:
            verdict = "mixed/small"
    t5_rows.append(dict(pair=pair_display,
                        gap_raw=f"{gap:+.2f}",
                        topology=f"{topo:+.2f}",
                        envelope=f"{env:+.2f}",
                        interaction=f"{inter:+.2f}" if abs(inter) >= 0.005 else "0.00",
                        verdict=verdict))
t5 = pd.DataFrame(t5_rows)


# -------------- TABLE S2 — 50 Tier A rows -------------------------------------
tA = pd.read_csv(R / "v7_final/tierA_final.csv")
t2 = tA[["family","contrast","subsystem","N_rep","W1_obs","snr","direction",
         "q_bh","p_matched_floor","n_sig_contrasts"]].copy()
t2.columns = ["family","contrast","subsystem","N","W1","snr","dir","q_BH","p_floor","n_rep_sig"]
t2["dir"] = t2["dir"].map(lambda x: "+" if x > 0 else ("−" if x < 0 else "0") if pd.notna(x) else "—")
t2 = t2.sort_values(["family","contrast","subsystem"]).reset_index(drop=True)
t2["N"] = t2["N"].astype(int)
t2["W1"] = t2["W1"].round(2)
t2["snr"] = t2["snr"].round(1)
t2["q_BH"] = t2["q_BH"].apply(lambda v: f"{v:.1e}" if pd.notna(v) else "—")
t2["p_floor"] = t2["p_floor"].apply(lambda v: f"{v:.3f}" if pd.notna(v) else "—")


# -------------- TABLE S3 — parent vs child biomass -----------------------------
pc = pd.read_csv(R / "tier0_v7/tier0_parent_child_v7.csv")
t3 = pc[["cell","condition","hpf","parent_biomax","child_biomax","parent_over_child"]].copy()
t3.columns = ["cell","condition","hpf","parent biomax","child biomax","parent / child"]
for col in ("parent biomax","child biomax","parent / child"):
    t3[col] = t3[col].round(2)


# -------------- TABLE S4 — Tier 0 robustness -----------------------------------
tr = pd.read_csv(R / "v7_final_addenda/tier0_robustness.csv")
t4 = tr[["cell","condition","hpf","biomax_baseline",
         "pct_valine","pct_histidine","pct_isoleucine","pct_glucose"]].copy()
t4.columns = ["cell","condition","hpf","biomax", "Δ% valine","Δ% histidine","Δ% isoleucine","Δ% glucose"]
for col in ("biomax", "Δ% valine","Δ% histidine","Δ% isoleucine","Δ% glucose"):
    t4[col] = t4[col].round(2)


# -------------- TABLE S6 — Tier A by subsystem ---------------------------------
ts = pd.read_csv(R / "v7_final/tierA_by_subsystem.csv")
t6 = ts[["subsystem","n_contrasts_F1","n_contrasts_F2","n_contrasts",
         "families","hpf_list","directions","W1_min","W1_max","snr_min","snr_max"]].copy()
t6.columns = ["subsystem","F1","F2","total","families","hpf","directions","W1_min","W1_max","snr_min","snr_max"]
for col in ("W1_min","W1_max","snr_min","snr_max"):
    t6[col] = t6[col].round(2)
t6 = t6.sort_values(["total","subsystem"], ascending=[False, True]).reset_index(drop=True)


# -------------- WRITE CSVs -----------------------------------------------------
t1.to_csv(TABLES_DIR / "table1_model_sizes.csv", index=False)
t2.to_csv(TABLES_DIR / "tableS2_tierA_50rows.csv", index=False)
t3.to_csv(TABLES_DIR / "tableS3_growth_capacity.csv", index=False)
t4.to_csv(TABLES_DIR / "tableS4_tier0_robustness.csv", index=False)
t5.to_csv(TABLES_DIR / "table5_attribution.csv", index=False)
t6.to_csv(TABLES_DIR / "tableS6_tierA_by_subsystem.csv", index=False)


# -------------- WRITE MARKDOWN -------------------------------------------------
md = ["# Chapter tables", "",
      "Six tables, keyed to the chapter's reading order.",
      "Main-text set: Table 1 (Zebrafish1 counts) and Table 5 (Tier 0 attribution).",
      "Supplementary set: Tables S2–S4 and S6 (log tables; prose carries the finding).",
      "", "---", "", "## MAIN TEXT", "", "### Table 1 — Zebrafish1 as used", ""]
md.append(_md_table(t1))
md += ["", "Source: direct counts from `ocr_anchored_extraction/models/baked_ocr_24_v7.mat` ",
       "(structurally identical to pristine Zebrafish1; OCR anchoring sets a single ",
       "lower bound on MAR09048 and does not add or remove reactions, metabolites, ",
       "genes, or compartments). See `grace_v4/results/v7_final_addenda/model_direct_counts.md` ",
       "for the generating code.",
       "",
       "### Table 5 — Tier 0 2×2 topology × envelope attribution (all 15 pairs)",
       ""]
md.append(_md_table(t5))
md += ["",
       "**Columns.** *gap_raw* is AA − BB on the extracted-model biomass_max in ",
       "the shared v7 base model with per-cell reaction masks (post-§11.19 corrected ",
       "values). *topology* and *envelope* are the two single-swap effects. *interaction* ",
       "= gap_raw − topology − envelope (sign convention: AA − BB = topology + envelope ",
       "+ interaction). *verdict* uses the three-share classification of Figure 5B.",
       "",
       "**Key reading.** For the two BL_96 pairs the interaction term carries 60–79 % of ",
       "the raw deficit (|−2.43| / |−4.04| = 60 %; |−3.59| / |−4.39| = 82 %); the gap is ",
       "**not separable** into additive topology + envelope effects. Figure 5B shows the ",
       "three-way share directly.",
       "",
       "---", "", "## SUPPLEMENTARY", "", "### Table S2 — 50 Tier A rows (compact)", ""]
md.append(_md_table(t2))
md += ["",
       "**Columns.** *N* is the number of representative-collapsed reactions in the ",
       "contrast's subsystem (replicates > 1 are two-tier eligible). *W1* is the pooled ",
       "Wasserstein distance on log₂(|v|/μ) after dropping zeros. *snr* = W1 / 95th-pct ",
       "of split-chain null. *dir* is the sign of the across-cell median shift. *q_BH* is ",
       "Benjamini-Hochberg within family (F1 or F2). *p_floor* is the matched-floor ",
       "re-sample empirical p when the row was flagged floor_suspect (NaN otherwise; see ",
       "§11.10 / §11.20). *n_rep_sig* is the number of replicating contrasts that also ",
       "cross q ≤ 0.10. Source: `grace_v4/results/v7_final/tierA_final.csv`.",
       "",
       "### Table S3 — Growth capacity per cell (parent base vs extracted child)", ""]
md.append(_md_table(t3))
md += ["",
       "**Reading.** *parent biomax* is biomass_max for the v7 base model under the ",
       "cell's (cond, hpf) envelope; *child biomax* uses the extracted model. Ratio ",
       "quantifies capability lost to extraction (BL_96 the extreme at 3.12). Plotted ",
       "as Figure 2A (dumbbell). Source: `grace_v4/results/tier0_v7/tier0_parent_child_v7.csv`.",
       "",
       "### Table S4 — Tier 0 robustness: 10 % uptake relaxation", ""]
md.append(_md_table(t4))
md += ["",
       "**Reading.** For each cell, lb of each exchange is relaxed by 10 % toward more ",
       "uptake; Δ% is the resulting change in biomass_max. His/Ile/Val non-zero rows ",
       "identify amino-acid binding under the stock Zebrafish1 BOF; glucose is near-zero ",
       "everywhere (biomass is not carbon-bound at these hpf). The row-level numbers were ",
       "pre-§11.21 the backbone of the §4.1 nitrogen-boundedness claim; see §4.1 for ",
       "the parameterisation-sensitivity restatement. Source: ",
       "`grace_v4/results/v7_final_addenda/tier0_robustness.csv`.",
       "",
       "### Table S6 — Tier A rolled up by subsystem", ""]
md.append(_md_table(t6))
md += ["",
       "**Reading.** One row per subsystem that contributes any Tier A finding. ",
       "*F1 / F2 / total* count contrasts within each family. *hpf* and *directions* ",
       "list the specific developmental transitions and signed shifts. ",
       "*W1_min* / *W1_max* / *snr_min* / *snr_max* bound the pooled Wasserstein ",
       "effect size and the signal-to-noise over the subsystem's Tier A contrasts. ",
       "Figure 3 renders this as the F1 × F2 orientation heatmap. ",
       "Source: `grace_v4/results/v7_final/tierA_by_subsystem.csv`.",
       ""]

(OUT / "tables.md").write_text("\n".join(md))
print(f"[write] {OUT/'tables.md'}")
print(f"[write] {TABLES_DIR}/ (six CSVs)")
print()
print("Tables summary:")
print(f"  Table 1 — {len(t1)} rows   (main)")
print(f"  Table 5 — {len(t5)} rows   (main; 6 material-gap + 9 no-gap)")
print(f"  Table S2 — {len(t2)} rows (supp)")
print(f"  Table S3 — {len(t3)} rows (supp)")
print(f"  Table S4 — {len(t4)} rows (supp)")
print(f"  Table S6 — {len(t6)} rows (supp)")
