"""
s272 -- PCA on cells + subsystem enrichment via hypergeometric.

PART A — PCA
  Reduce each of the 15 cells to a vector of per-reaction log mean non-zero
  flux. Standardize columns, PCA to 3 components. Plot PC1-PC2 colored by
  condition, shaped by hpf. Export loadings and explained variance.

PART B — Subsystem enrichment
  For each (subsystem, contrast), test whether notable reactions
  (|log2FC_mean| >= 1 AND carry_both AND not loop_flagged) are
  over-represented in the subsystem relative to the full model via
  Fisher exact. BH-correct within contrast. Writes two new columns to
  the subsystem summary CSVs: enrichment_p, enrichment_q.

Outputs:
  grace_v4/results/chapter_assets/comparison_tables/
    pca_cells.csv                     cell × PC coordinates
    pca_loadings.csv                  reaction × PC loadings
    pca_explained_variance.csv        variance per PC
    pca_cells.pdf / .png              scatter plot
    subsystem_summary_F1.csv          REWRITTEN with enrichment cols
    subsystem_summary_F2.csv          REWRITTEN with enrichment cols
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import fisher_exact, false_discovery_control

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import POOLED as POOL, SUBSETS as SUBS, COMPARISON_TABLES as OUT
OUT.mkdir(parents=True, exist_ok=True)

CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]
ZERO_TOL = 1e-10
LOG2FC_NOTABLE = 1.0

# ================================================================
# PART A — PCA
# ================================================================
print("[PCA] building cell × reaction matrix of log mean non-zero flux")

# Collect per-cell log mean non-zero flux vectors
cell_vectors = {}
for name in CELLS:
    df = pd.read_csv(POOL / f"{name}_pooled.csv", index_col=0)
    vals = np.abs(df.values.astype(float))
    # mean of non-zero samples per reaction (nan if all zero)
    nz_mask = vals > ZERO_TOL
    masked = np.where(nz_mask, vals, np.nan)
    with np.errstate(all="ignore"):
        mean_nz = np.nanmean(masked, axis=1)
    cell_vectors[name] = pd.Series(mean_nz, index=df.index)
    print(f"  {name:7s}  {len(mean_nz):5d} rxns  "
          f"{np.isfinite(mean_nz).sum()} with non-zero flux")

# Intersection of reactions carrying in all 15 cells
all_rxns = set.intersection(*[set(v[v.notna()].index) for v in cell_vectors.values()])
all_rxns = sorted(all_rxns)
print(f"\n[PCA] intersection: {len(all_rxns)} reactions present in all 15 cells")

# Build matrix
mat = np.full((len(CELLS), len(all_rxns)), np.nan)
for i, name in enumerate(CELLS):
    v = cell_vectors[name]
    mat[i, :] = v.reindex(all_rxns).values
# log2 transform
mat_log = np.log2(mat + 1e-6)
mat_log = pd.DataFrame(mat_log, index=CELLS, columns=all_rxns)

# Standardize columns (so no single high-flux reaction dominates)
scaler = StandardScaler()
X = scaler.fit_transform(mat_log.values)

# PCA
pca = PCA(n_components=3)
pcs = pca.fit_transform(X)
var_ratio = pca.explained_variance_ratio_
print(f"\n[PCA] variance explained: PC1={var_ratio[0]*100:.1f}%, "
      f"PC2={var_ratio[1]*100:.1f}%, PC3={var_ratio[2]*100:.1f}%")

# Save PC coordinates
pc_df = pd.DataFrame(pcs, index=CELLS, columns=["PC1","PC2","PC3"])
pc_df["condition"] = [c.split("_")[0] for c in CELLS]
pc_df["hpf"] = [int(c.split("_")[1]) for c in CELLS]
pc_df.to_csv(OUT / "pca_cells.csv")

# Save loadings
loadings = pd.DataFrame(pca.components_.T, index=all_rxns,
                        columns=["PC1","PC2","PC3"])
# Add subsystem for each rxn
mem = pd.read_csv(SUBS / "BL_24_subset_membership.csv")
rep_to_sub = (mem.drop_duplicates("representative")
                 .set_index("representative")["subsystem"].to_dict())
loadings["subsystem"] = loadings.index.map(lambda r: rep_to_sub.get(r, ""))
loadings.to_csv(OUT / "pca_loadings.csv")

# Save variance
pd.DataFrame({
    "PC": ["PC1","PC2","PC3"],
    "variance_ratio": var_ratio,
    "cumulative": np.cumsum(var_ratio),
}).to_csv(OUT / "pca_explained_variance.csv", index=False)

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
cond_color = {"BL": "#2E7D32", "D": "#333333", "LD": "#E65100"}
hpf_marker = {24: "o", 48: "s", 72: "^", 96: "D", 120: "P"}
for ax_idx, (xlab, ylab) in enumerate([("PC1","PC2"), ("PC1","PC3")]):
    ax = axes[ax_idx]
    for cell in CELLS:
        c, h = cell.split("_"); hpf = int(h)
        row = pc_df.loc[cell]
        ax.scatter(row[xlab], row[ylab],
                   c=cond_color[c], marker=hpf_marker[hpf],
                   s=180, edgecolors="black", linewidths=1.0, alpha=0.9)
        ax.annotate(cell, (row[xlab], row[ylab]),
                    xytext=(5,5), textcoords="offset points", fontsize=8)
    xi = int(xlab[-1])-1; yi = int(ylab[-1])-1
    ax.set_xlabel(f"{xlab} ({var_ratio[xi]*100:.1f}%)", fontsize=11)
    ax.set_ylabel(f"{ylab} ({var_ratio[yi]*100:.1f}%)", fontsize=11)
    ax.axhline(0, color="gray", linewidth=0.5, alpha=0.3)
    ax.axvline(0, color="gray", linewidth=0.5, alpha=0.3)
    ax.grid(alpha=0.3)
# Legend
from matplotlib.lines import Line2D
leg_cond = [Line2D([0],[0], marker="o", linestyle="none",
                    markerfacecolor=cond_color[c], markeredgecolor="black",
                    markersize=10, label=c) for c in CONDS]
leg_hpf  = [Line2D([0],[0], marker=hpf_marker[h], linestyle="none",
                    markerfacecolor="#AAAAAA", markeredgecolor="black",
                    markersize=10, label=f"{h} hpf") for h in HPFS]
axes[0].legend(handles=leg_cond + leg_hpf, loc="best", fontsize=9,
               ncol=2, framealpha=0.9)
fig.suptitle(f"PCA on cell × log(mean non-zero flux) across "
             f"{len(all_rxns)} shared reactions", fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "pca_cells.pdf", bbox_inches="tight")
fig.savefig(OUT / "pca_cells.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"[write] {OUT/'pca_cells.pdf'}")

# Top 10 reactions on each PC
print("\n[loadings] top 10 reactions on each PC (absolute loading):")
for pc in ["PC1","PC2","PC3"]:
    top = loadings.reindex(loadings[pc].abs().sort_values(ascending=False).index).head(10)
    print(f"\n  {pc}:")
    for rxn, row in top.iterrows():
        print(f"    {rxn:12s} {row[pc]:+.3f}  {row.subsystem}")


# ================================================================
# PART B — Hypergeometric subsystem enrichment
# ================================================================
print("\n[enrichment] subsystem enrichment via Fisher exact")

for fam in ("F1", "F2"):
    rxn = pd.read_csv(OUT / f"per_reaction_flux_{fam}.csv")
    roll = pd.read_csv(OUT / f"subsystem_summary_{fam}.csv")
    # notable_rxn flag
    rxn["notable"] = (rxn.carry_both & (~rxn.loop_flagged) &
                      (rxn.log2FC_mean.abs() >= LOG2FC_NOTABLE))

    # Per-contrast totals
    rows_out = []
    for contrast, g in rxn.groupby("contrast"):
        total_tested = int((g.carry_both & (~g.loop_flagged)).sum())
        total_notable = int(g.notable.sum())
        for sub, gs in g.groupby("subsystem"):
            if not sub: continue
            sub_tested = int((gs.carry_both & (~gs.loop_flagged)).sum())
            sub_notable = int(gs.notable.sum())
            # 2x2:
            #        notable  not-notable
            # in-sub a         b
            # elsew  c         d
            a = sub_notable
            b = sub_tested - sub_notable
            c = total_notable - sub_notable
            d = total_tested - total_notable - b
            if a + b < 1 or a + c < 1:
                p = np.nan; odds = np.nan
            else:
                try:
                    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
                except Exception:
                    odds, p = np.nan, np.nan
            rows_out.append(dict(contrast=contrast, subsystem=sub,
                                 n_notable_in_sub=a,
                                 n_tested_in_sub=sub_tested,
                                 n_notable_elsewhere=c,
                                 n_tested_elsewhere=total_tested - sub_tested,
                                 odds_ratio=odds, enrichment_p=p))
    enr = pd.DataFrame(rows_out)
    # BH correct within contrast
    enr["enrichment_q"] = np.nan
    for contrast, idx in enr.groupby("contrast").groups.items():
        ps = enr.loc[idx, "enrichment_p"].values
        good = np.isfinite(ps)
        if good.sum() > 0:
            q = np.full_like(ps, np.nan, dtype=float)
            q[good] = false_discovery_control(ps[good], method="bh")
            enr.loc[idx, "enrichment_q"] = q

    # Merge into roll
    roll_new = roll.merge(enr[["contrast","subsystem","odds_ratio",
                               "enrichment_p","enrichment_q"]],
                          on=["contrast","subsystem"], how="left")
    roll_new.to_csv(OUT / f"subsystem_summary_{fam}.csv", index=False)
    n_sig = (roll_new.enrichment_q < 0.05).sum()
    print(f"  {fam}: wrote {len(roll_new)} rows, "
          f"{n_sig} subsystems with enrichment_q < 0.05")

# Show the top enriched subsystems per family
print("\n=== top 15 enriched subsystems (F2) ===")
r2 = pd.read_csv(OUT / "subsystem_summary_F2.csv")
top2 = r2.loc[r2.enrichment_q.notna()].nsmallest(15, "enrichment_q")
cols = ["contrast","subsystem","n_notable_in_sub","n_tested_in_sub",
        "odds_ratio","enrichment_q","verdict","mean_log2FC","coherence"]
print(top2[cols].to_string(index=False))

print("\n=== top 15 enriched subsystems (F1) ===")
r1 = pd.read_csv(OUT / "subsystem_summary_F1.csv")
top1 = r1.loc[r1.enrichment_q.notna()].nsmallest(15, "enrichment_q")
print(top1[cols].to_string(index=False))

print("\nDone.")
