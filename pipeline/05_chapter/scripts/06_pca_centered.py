"""
s273 -- Per-cell-centered PCA.

Same input as s272 (log mean non-zero flux per reaction per cell, 15 cells
on 3,905 shared reactions), with one additional preprocessing step:
subtract each cell's median log flux across reactions BEFORE column
standardization and PCA.

Removes per-cell magnitude scaling (OCR-with-hpf drift, biomass-budget
asymmetry at BL_96) and leaves pathway composition differences.

Writes alongside the existing PCA files:
    pca_cells_centered.csv
    pca_loadings_centered.csv
    pca_explained_variance_centered.csv
    pca_cells_centered.pdf / .png
    pca_cells_comparison.pdf / .png     (raw vs centered side-by-side)
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import POOLED as POOL, SUBSETS as SUBS, COMPARISON_TABLES as OUT
OUT.mkdir(parents=True, exist_ok=True)

CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]
ZERO_TOL = 1e-10

cond_color = {"BL": "#2E7D32", "D": "#333333", "LD": "#E65100"}
hpf_marker = {24: "o", 48: "s", 72: "^", 96: "D", 120: "P"}

# ---------- Build cell × reaction log flux matrix (same as s272) ----------
print("[PCA] rebuilding cell × reaction log(mean non-zero flux) matrix")
cell_vectors = {}
for name in CELLS:
    df = pd.read_csv(POOL / f"{name}_pooled.csv", index_col=0)
    vals = np.abs(df.values.astype(float))
    nz_mask = vals > ZERO_TOL
    masked = np.where(nz_mask, vals, np.nan)
    with np.errstate(all="ignore"):
        mean_nz = np.nanmean(masked, axis=1)
    cell_vectors[name] = pd.Series(mean_nz, index=df.index)

all_rxns = set.intersection(*[set(v[v.notna()].index) for v in cell_vectors.values()])
all_rxns = sorted(all_rxns)
print(f"[PCA] {len(all_rxns)} reactions present in all 15 cells")

mat = np.full((len(CELLS), len(all_rxns)), np.nan)
for i, name in enumerate(CELLS):
    v = cell_vectors[name]
    mat[i, :] = v.reindex(all_rxns).values

mat_log = np.log2(mat + 1e-6)

# ---------- Report per-cell median log flux (what we're removing) ----------
cell_medians = np.nanmedian(mat_log, axis=1)
print("\n[centering] per-cell median log flux (will be subtracted):")
for i, name in enumerate(CELLS):
    print(f"  {name:7s}  median log flux = {cell_medians[i]:+.3f}")

# ---------- Center each row ----------
mat_centered = mat_log - cell_medians[:, None]

# Standardize columns
scaler_centered = StandardScaler()
X_centered = scaler_centered.fit_transform(mat_centered)

# PCA
pca_c = PCA(n_components=3)
pcs_c = pca_c.fit_transform(X_centered)
var_c = pca_c.explained_variance_ratio_
print(f"\n[PCA centered] variance explained: "
      f"PC1={var_c[0]*100:.1f}%  PC2={var_c[1]*100:.1f}%  PC3={var_c[2]*100:.1f}%")

# Save coordinates
pc_df = pd.DataFrame(pcs_c, index=CELLS, columns=["PC1","PC2","PC3"])
pc_df["condition"] = [c.split("_")[0] for c in CELLS]
pc_df["hpf"] = [int(c.split("_")[1]) for c in CELLS]
pc_df["cell_median_log_flux"] = cell_medians
pc_df.to_csv(OUT / "pca_cells_centered.csv")

# Loadings
mem = pd.read_csv(SUBS / "BL_24_subset_membership.csv")
rep_to_sub = (mem.drop_duplicates("representative")
                 .set_index("representative")["subsystem"].to_dict())
loadings_c = pd.DataFrame(pca_c.components_.T, index=all_rxns,
                           columns=["PC1","PC2","PC3"])
loadings_c["subsystem"] = loadings_c.index.map(lambda r: rep_to_sub.get(r, ""))
loadings_c.to_csv(OUT / "pca_loadings_centered.csv")

pd.DataFrame({
    "PC": ["PC1","PC2","PC3"],
    "variance_ratio": var_c,
    "cumulative": np.cumsum(var_c),
}).to_csv(OUT / "pca_explained_variance_centered.csv", index=False)

# ---------- Standalone centered PCA figure ----------
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
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
    ax.set_xlabel(f"{xlab} ({var_c[xi]*100:.1f}%)", fontsize=11)
    ax.set_ylabel(f"{ylab} ({var_c[yi]*100:.1f}%)", fontsize=11)
    ax.axhline(0, color="gray", linewidth=0.5, alpha=0.3)
    ax.axvline(0, color="gray", linewidth=0.5, alpha=0.3)
    ax.grid(alpha=0.3)
leg_cond = [Line2D([0],[0], marker="o", linestyle="none",
                    markerfacecolor=cond_color[c], markeredgecolor="black",
                    markersize=10, label=c) for c in CONDS]
leg_hpf  = [Line2D([0],[0], marker=hpf_marker[h], linestyle="none",
                    markerfacecolor="#AAAAAA", markeredgecolor="black",
                    markersize=10, label=f"{h} hpf") for h in HPFS]
axes[0].legend(handles=leg_cond + leg_hpf, loc="best", fontsize=9,
               ncol=2, framealpha=0.9)
fig.suptitle(f"Per-cell-centered PCA (magnitude removed, pathway composition only)\n"
             f"{len(all_rxns)} shared reactions", fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "pca_cells_centered.pdf", bbox_inches="tight")
fig.savefig(OUT / "pca_cells_centered.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"[write] pca_cells_centered.pdf")

# ---------- Comparison figure: raw vs centered PCA ----------
# Load raw PCA
pc_raw = pd.read_csv(OUT / "pca_cells.csv", index_col=0)
var_raw = pd.read_csv(OUT / "pca_explained_variance.csv")
vr = var_raw.variance_ratio.values

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax_idx, (label, pc_data, var_vals) in enumerate(
    [("A — Raw (magnitude + composition)", pc_raw, vr),
     ("B — Per-cell centered (composition only)", pc_df, var_c)]):
    ax = axes[ax_idx]
    for cell in CELLS:
        c, h = cell.split("_"); hpf = int(h)
        row = pc_data.loc[cell]
        ax.scatter(row.PC1, row.PC2,
                   c=cond_color[c], marker=hpf_marker[hpf],
                   s=180, edgecolors="black", linewidths=1.0, alpha=0.9)
        ax.annotate(cell, (row.PC1, row.PC2),
                    xytext=(5,5), textcoords="offset points", fontsize=8)
    ax.set_xlabel(f"PC1 ({var_vals[0]*100:.1f}%)", fontsize=11)
    ax.set_ylabel(f"PC2 ({var_vals[1]*100:.1f}%)", fontsize=11)
    ax.set_title(label, fontsize=11)
    ax.axhline(0, color="gray", linewidth=0.5, alpha=0.3)
    ax.axvline(0, color="gray", linewidth=0.5, alpha=0.3)
    ax.grid(alpha=0.3)

axes[0].legend(handles=leg_cond + leg_hpf, loc="best", fontsize=9,
               ncol=2, framealpha=0.9)
fig.suptitle("PCA before and after per-cell magnitude normalization",
             fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "pca_cells_comparison.pdf", bbox_inches="tight")
fig.savefig(OUT / "pca_cells_comparison.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"[write] pca_cells_comparison.pdf")

# ---------- Top 10 reactions per PC ----------
print("\n[loadings centered] top 10 reactions per PC (absolute loading):")
for pc in ["PC1","PC2","PC3"]:
    top = loadings_c.reindex(loadings_c[pc].abs().sort_values(ascending=False).index).head(10)
    print(f"\n  {pc}:")
    for rxn, row in top.iterrows():
        print(f"    {rxn:12s} {row[pc]:+.3f}  {row.subsystem}")

# ---------- Summary of what moved ----------
print("\n[comparison] Position shift from raw to centered PCA:")
print(f"  {'cell':<8} {'raw PC1':>10} {'cen PC1':>10} {'ΔPC1':>8} "
      f"{'raw PC2':>10} {'cen PC2':>10} {'ΔPC2':>8}  cell_median")
for cell in CELLS:
    rp = pc_raw.loc[cell]
    cp = pc_df.loc[cell]
    print(f"  {cell:<8} {rp.PC1:>10.2f} {cp.PC1:>10.2f} {cp.PC1-rp.PC1:>+8.2f} "
          f"{rp.PC2:>10.2f} {cp.PC2:>10.2f} {cp.PC2-rp.PC2:>+8.2f}  "
          f"{cell_medians[CELLS.index(cell)]:+.2f}")
