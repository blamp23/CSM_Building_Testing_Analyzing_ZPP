"""
s200 — easy figures from the v7 CSV deliverables.

Scope: Fig 1 (1a schematic + 1b funnel), Fig 2 (F2 heatmap), Fig 3 (F1 heatmap),
Fig 5 (3 panels), figS4 (Rhat/ESS violin), figS5 (mean/median concordance).

Deferred (need pooled samples + null regen): Fig 1c, Fig 4, figS1, S2, S3.
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

ROOT = Path(__file__).resolve().parents[2]
MOR  = ROOT / "grace_v4/results/findings_v7_mean_dropzero"
MED  = ROOT / "grace_v4/results/findings_v7_median"
T0   = ROOT / "grace_v4/results/tier0_v7"
S146 = ROOT / "reviewer_packet_v4_slim/s146_tier0_decomposition.csv"
BMAX = ROOT / "reviewer_packet_v4_slim/biomass_max_tier0.csv"
RHAT = ROOT / "grace_v4/results/hopsy_v7_local/rank_rhat_summary_hopsy_v7_local.csv"
OUT  = ROOT / "grace_v4/results/figures_v7"
OUT.mkdir(parents=True, exist_ok=True)

CCOL = {"BL": "#1F77B4", "D": "#111111", "LD": "#F0A800"}
HPFS = [24, 48, 72, 96, 120]
HPF_RAMP = {h: plt.get_cmap("viridis")(i/4) for i, h in enumerate(HPFS)}

# ----- short names ---------------------------------------------------------
SHORT = {
    "Phenylalanine, tyrosine and tryptophan biosynthesis": "Aromatic AA",
    "Valine, leucine, and isoleucine metabolism": "BCAA",
    "Arginine and proline metabolism": "Arg/Pro",
    "Glycolysis / Gluconeogenesis": "Glycolysis",
    "Fatty acid activation (cytosolic)": "FA activation (cyt)",
    "Fatty acid activation (endoplasmic reticular)": "FA activation (ER)",
}
def short(s):
    s = str(s)
    if s in SHORT: return SHORT[s]
    return re.sub(r"\s+metabolism$", "", s)

# ----- classes -------------------------------------------------------------
CLASS_ORDER = ["Amino acid", "Nucleotide",
               "Eicosanoid / FA oxidation", "Sterol / lipid / glyceride",
               "Early lipid synthesis", "Central carbon / other"]
AA_SET = {"Aromatic AA","BCAA","Arg/Pro","Histidine","Tyrosine","Tryptophan","Lysine",
          "Cysteine and methionine","Glycine, serine and threonine",
          "Alanine, aspartate and glutamate","Beta-alanine"}
NUC_SET = {"Nucleotide","Pyrimidine","Amino sugar and nucleotide sugar"}
EICO_PRE = ("Arachidonic acid","Leukotriene","Fatty acid oxidation",
            "FA activation (cyt)","FA activation (ER)","Acyl-CoA hydrolysis")
STEROL_PRE = ("Steroid","Cholesterol","Bile acid biosynthesis","Acylglycerides",
              "Glycerophospholipid","Sphingolipid","Ganglio")
EARLY_PRE = ("FA elongation","FA biosynthesis (unsaturated)","FA desaturation","Omega-3")
def classify(sh):
    s_low = sh.lower()
    if sh in AA_SET: return "Amino acid"
    if sh in NUC_SET: return "Nucleotide"
    if "beta-oxidation" in s_low or any(sh.startswith(p) for p in EICO_PRE): return "Eicosanoid / FA oxidation"
    if any(p.lower() in s_low for p in STEROL_PRE): return "Sterol / lipid / glyceride"
    if any(sh.startswith(p) for p in EARLY_PRE): return "Early lipid synthesis"
    return "Central carbon / other"

# ----- load ----------------------------------------------------------------
tA = pd.read_csv(MOR/"tierA_final.csv")
tB = pd.read_csv(MOR/"tierB_final.csv")
tC = pd.read_csv(MOR/"tierC_final.csv")
allf = pd.read_csv(MOR/"all_rows_in_family.csv")
negc = pd.read_csv(MOR/"negcontrol_per_subsystem.csv")
s146 = pd.read_csv(S146)
bmax_t0 = pd.read_csv(BMAX)
tt2x2 = pd.read_csv(T0/"tier0_2x2_v7.csv")
greedy = pd.read_csv(T0/"bl96_greedy_addback.csv")
rhat = pd.read_csv(RHAT)

allf["short"] = allf.subsystem.apply(short)
allf["class"] = allf["short"].apply(classify)
tA["short"] = tA.subsystem.apply(short); tA["class"] = tA["short"].apply(classify)
tB["short"] = tB.subsystem.apply(short); tB["class"] = tB["short"].apply(classify)
tC["short"] = tC.subsystem.apply(short); tC["class"] = tC["short"].apply(classify)

bm_cell = s146.set_index("cell")["child_biomax"].to_dict()
parent_cell = s146.set_index("cell")["parent_biomax"].to_dict()

# =========================================================================
# FIG 1 — 1a schematic + 1b funnel
# =========================================================================
fig = plt.figure(figsize=(11, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.4], wspace=0.3)

# ----- 1a: 3x5 grid -------
ax1a = fig.add_subplot(gs[0, 0])
for i, c in enumerate(("BL","D","LD")):
    for j, h in enumerate(HPFS):
        ax1a.add_patch(plt.Rectangle((j, 2-i), 0.9, 0.9, facecolor=CCOL[c], edgecolor="k", lw=0.5))
        ax1a.text(j+0.45, 2-i+0.45, f"{c}_{h}", ha="center", va="center", fontsize=7,
                  color="white" if c=="D" else "black")
ax1a.set_xlim(-0.2, 6); ax1a.set_ylim(-0.5, 3.5)
ax1a.set_xticks([j+0.45 for j in range(5)]); ax1a.set_xticklabels([f"{h} hpf" for h in HPFS], fontsize=7)
ax1a.set_yticks([2+0.45,1+0.45,0+0.45]); ax1a.set_yticklabels(["BL","D","LD"], fontsize=7)
ax1a.set_aspect("equal"); ax1a.set_frame_on(False)
# inputs text
for yy, txt in zip([3.15, 2.95, 2.75],
                   ["RNA-seq → network (rFASTCORMICS)",
                    "metabolomics → envelope (105 exchanges)",
                    "OCR → O₂ cap"]):
    ax1a.text(0, yy, "→ "+txt, fontsize=7, va="center")
ax1a.set_title("(a) design: 15 cell-contexts", fontsize=9, loc="left")

# ----- 1b: funnel -------
ax1b = fig.add_subplot(gs[0, 1])
stages = ["in family","q ≤ 0.10","snr+zero+env","floor-safe","replicated","Tier A"]
real_n = [1687, 146, 125, 115, 68, 62]
nc_n   = [976, 5, 0, 0, 0, 0]
y = np.arange(len(stages))[::-1]
h = 0.35
ax1b.barh(y+h/2, real_n, h, color="#808080", label=f"real (1687→62)")
ax1b.barh(y-h/2, nc_n,   h, color="#F3A6A6", label=f"neg-ctrl (976→0)")
for yi, rn, nn in zip(y, real_n, nc_n):
    ax1b.text(rn, yi+h/2, f"  {rn}", va="center", fontsize=7)
    ax1b.text(nn, yi-h/2, f"  {nn}", va="center", fontsize=7)
ax1b.set_yticks(y); ax1b.set_yticklabels(stages, fontsize=7)
ax1b.set_xscale("log"); ax1b.set_xlabel("row count"); ax1b.set_xlim(0.5, 3000)
ax1b.set_title("(b) funnel", fontsize=9, loc="left")
ax1b.legend(fontsize=7, loc="lower right", frameon=False)
for sp in ("top","right"): ax1b.spines[sp].set_visible(False)

plt.tight_layout()
fig.savefig(OUT/"fig1.pdf", bbox_inches="tight")
fig.savefig(OUT/"fig1.png", dpi=300, bbox_inches="tight")
plt.close(fig)
pd.DataFrame({"stage":stages,"real":real_n,"neg_ctrl":nc_n}).to_csv(OUT/"fig1_data.csv", index=False)
print("[write] fig1.pdf + png + _data.csv")


# =========================================================================
# Helpers for heatmaps
# =========================================================================
def best_signed_W1(sub_df, flip_sign_for_f2):
    """One row per contrast, pick signed W1. Returns dict (contrast)->signed_W1."""
    out = {}
    for _, r in sub_df.iterrows():
        s = -r.direction*r.W1_obs if flip_sign_for_f2 else r.direction*r.W1_obs
        out[r.contrast] = s
    return out

def row_order_in_class(df_class, max_abs_colname="max_abs"):
    return df_class.sort_values(max_abs_colname, ascending=False)

# =========================================================================
# FIG 2 — F2 heatmap with growth-capacity strip
# =========================================================================
f2 = allf[allf.family=="F2"].copy()
# keep subsystems having any A/B/C row
subs_with_sig = sorted(set(f2.loc[f2.tier.isin(["A","B","C"]), "short"]))
f2 = f2[f2["short"].isin(subs_with_sig)].copy()
f2["signed_W1"] = -f2.direction * f2.W1_obs

# column scheme: 12 base + 4 consensus (one after each hpf-transition group)
TRANSITIONS = [(24,48),(48,72),(72,96),(96,120)]
CONDS = ["BL","D","LD"]
col_specs = []
for a,b in TRANSITIONS:
    for c in CONDS:
        col_specs.append(("data", c, a, b, f"{a}→{b}", c))
    col_specs.append(("consensus", None, a, b, f"{a}→{b}", "consensus"))

# subsystem order: by class, within class by max|signed_W1|
rows_df = (f2.groupby("short")["signed_W1"].apply(lambda s: s.abs().max())
           .reset_index(name="max_abs"))
rows_df["class"] = rows_df["short"].apply(classify)
rows_df["class_order"] = rows_df["class"].apply(lambda c: CLASS_ORDER.index(c))
rows_df = rows_df.sort_values(["class_order","max_abs"], ascending=[True,False]).reset_index(drop=True)
ordered_subs = rows_df["short"].tolist()

# Signed_W1 matrix
n_rows = len(ordered_subs); n_cols = len(col_specs)
M = np.full((n_rows, n_cols), np.nan)
tier_mat = np.full((n_rows, n_cols), "", dtype=object)
floor_mat = np.zeros((n_rows, n_cols), dtype=bool)
in_family = np.zeros((n_rows, n_cols), dtype=bool)
sub_lookup = {s: i for i, s in enumerate(ordered_subs)}
# Build cell-lookup for F2 contrasts
f2_lookup = {}
for _, r in f2.iterrows():
    f2_lookup[(r["short"], r.cell_A, r.cell_B)] = r

for j, (kind, cond, a, b, tr_lbl, c_lbl) in enumerate(col_specs):
    if kind == "consensus": continue
    cA, cB = f"{cond}_{a}", f"{cond}_{b}"
    for i, sh in enumerate(ordered_subs):
        row = f2_lookup.get((sh, cA, cB))
        if row is None: continue
        in_family[i, j] = True
        M[i, j] = -row.direction * row.W1_obs
        if row.tier in ("A","B","C"):
            tier_mat[i, j] = row.tier
        # floor flag
        bmA, bmB = bm_cell.get(cA), bm_cell.get(cB)
        if bmA and bmB:
            if max(bmA,bmB)/min(bmA,bmB) > 1.15:
                floor_mat[i, j] = True

# consensus columns
for j, (kind, cond, a, b, tr_lbl, c_lbl) in enumerate(col_specs):
    if kind != "consensus": continue
    data_cols = [jj for jj,(k,c,aa,bb,_,_) in enumerate(col_specs) if k=="data" and aa==a and bb==b]
    tiers = tier_mat[:, data_cols]
    signs = np.sign(M[:, data_cols])
    for i in range(n_rows):
        if all(t in ("A","B") for t in tiers[i]):
            s0 = signs[i,0]
            if s0 != 0 and (signs[i]==s0).all():
                M[i, j] = s0 * 1.0  # use ±1 for consensus flag; just needs sign for cmap
                in_family[i, j] = True
                tier_mat[i, j] = "consensus"

vmax = max(np.nanmax(np.abs(M[:,:12])), 0.01) if np.isfinite(np.nanmax(M[:,:12])) else 1.0

# plot
fig = plt.figure(figsize=(9, 10))
gs = fig.add_gridspec(2, 1, height_ratios=[1.5, 10], hspace=0.08)
ax_top = fig.add_subplot(gs[0,0])
for c in CONDS:
    ys = [bm_cell[f"{c}_{h}"] for h in HPFS]
    ax_top.plot(HPFS, ys, "-o", color=CCOL[c], markersize=4, lw=1.5, label=c)
ax_top.set_xticks(HPFS); ax_top.set_xticklabels(HPFS, fontsize=7)
ax_top.set_ylabel("growth capacity\n(biomass_max)", fontsize=7)
ax_top.legend(fontsize=7, frameon=False, loc="lower right", ncol=3)
ax_top.set_title("Fig 2 — F2 developmental program (signed_W1; red = rises with development)",
                 fontsize=9, loc="left")
ax_top.tick_params(labelsize=7)
for sp in ("top","right"): ax_top.spines[sp].set_visible(False)

ax = fig.add_subplot(gs[1,0])
# background: fill "not in family" with hatch
for i in range(n_rows):
    for j in range(n_cols):
        if not in_family[i,j]:
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="#F5F5F5",
                                       hatch="////", edgecolor="#D0D0D0", lw=0))
# fill tier C cells grey + c glyph
Mplot = M.copy()
im = ax.imshow(Mplot, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
               aspect="auto", extent=(0, n_cols, n_rows, 0), interpolation="none")
# tier overlays
for i in range(n_rows):
    for j in range(n_cols):
        t = tier_mat[i,j]
        if t == "A":
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="none", edgecolor="black", lw=1.5))
        elif t == "B":
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="none", edgecolor="#555", lw=0.8))
        elif t == "C":
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="#D9D9D9", edgecolor="none"))
            ax.text(j+0.5, i+0.5, "c", fontsize=6, ha="center", va="center", color="black")
        if floor_mat[i,j]:
            ax.text(j+0.5, i+0.7, "×", fontsize=6, ha="center", va="center", color="black")

# column labels
col_labels_top = [cs[4] for cs in col_specs]  # transition
col_labels_bot = [cs[5] for cs in col_specs]  # cond or consensus
ax.set_xticks([j+0.5 for j in range(n_cols)])
ax.set_xticklabels([f"{a}\n{b}" for a,b in zip(col_labels_top, col_labels_bot)],
                   fontsize=6, rotation=0)
# row labels + class rules
ax.set_yticks([i+0.5 for i in range(n_rows)])
ax.set_yticklabels(ordered_subs, fontsize=6)
# horizontal class rules
class_bounds = []
prev_class = None
for i, s in enumerate(ordered_subs):
    cls = classify(s)
    if cls != prev_class:
        if prev_class is not None:
            class_bounds.append(i)
        prev_class = cls
for b in class_bounds:
    ax.axhline(b, color="black", lw=0.6)

ax.set_xlim(0, n_cols); ax.set_ylim(n_rows, 0)
cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
cbar.set_label("signed W1 (red = rises with development)", fontsize=7)
cbar.ax.tick_params(labelsize=6)
# legend patches
legend_patches = [
    mpatches.Patch(facecolor="none", edgecolor="black", lw=1.5, label="Tier A"),
    mpatches.Patch(facecolor="none", edgecolor="#555", lw=0.8, label="Tier B"),
    mpatches.Patch(facecolor="#D9D9D9", edgecolor="none", label="Tier C"),
    mpatches.Patch(facecolor="#F5F5F5", hatch="////", edgecolor="#D0D0D0", label="not in family"),
    mpatches.Patch(facecolor="white", edgecolor="black", lw=0.3, label="floor-exposed ×"),
]
ax.legend(handles=legend_patches, loc="lower right", bbox_to_anchor=(1.3, -0.15),
          fontsize=6, frameon=False)

fig.savefig(OUT/"fig2.pdf", bbox_inches="tight")
fig.savefig(OUT/"fig2.png", dpi=300, bbox_inches="tight")
plt.close(fig)
# data csv
rows_out = []
for i, sh in enumerate(ordered_subs):
    for j, (kind, cond, a, b, tr_lbl, c_lbl) in enumerate(col_specs):
        rows_out.append(dict(subsystem=sh, cls=classify(sh),
                             col_kind=kind, condition=cond or "",
                             transition=f"{a}_{b}",
                             signed_W1=M[i,j], tier=tier_mat[i,j],
                             floor_exposed=bool(floor_mat[i,j]),
                             in_family=bool(in_family[i,j])))
pd.DataFrame(rows_out).to_csv(OUT/"fig2_data.csv", index=False)
print(f"[write] fig2.pdf + png + _data.csv  ({n_rows} rows × {n_cols} cols)")


# =========================================================================
# FIG 3 — F1 heatmap
# =========================================================================
f1 = allf[allf.family=="F1"].copy()
f1_subs = sorted(set(f1.loc[f1.tier.isin(["A","B","C"]), "short"]))
f1 = f1[f1["short"].isin(f1_subs)].copy()
f1["signed_W1"] = f1.direction * f1.W1_obs

# column order: hpf × (BL vs LD, D vs LD, BL vs D)
F1_PAIRS_PER_HPF = [("BL","LD"), ("D","LD"), ("BL","D")]
col_specs_f1 = []
for h in HPFS:
    for a,b in F1_PAIRS_PER_HPF:
        col_specs_f1.append((h, a, b, f"{a}vs{b}"))

rows_df_f1 = (f1.groupby("short")["signed_W1"].apply(lambda s: s.abs().max())
              .reset_index(name="max_abs"))
rows_df_f1["class"] = rows_df_f1["short"].apply(classify)
rows_df_f1["class_order"] = rows_df_f1["class"].apply(lambda c: CLASS_ORDER.index(c))
rows_df_f1 = rows_df_f1.sort_values(["class_order","max_abs"], ascending=[True,False]).reset_index(drop=True)
ordered_subs_f1 = rows_df_f1["short"].tolist()

f1_lookup = {}
for _, r in f1.iterrows():
    f1_lookup[(r["short"], r.cell_A, r.cell_B)] = r

n_rows1 = len(ordered_subs_f1); n_cols1 = len(col_specs_f1)
M1 = np.full((n_rows1, n_cols1), np.nan)
tier_mat1 = np.full((n_rows1, n_cols1), "", dtype=object)
in_fam1 = np.zeros((n_rows1, n_cols1), dtype=bool)
is_bl96_col = np.zeros(n_cols1, dtype=bool)
for j, (h, a, b, lbl) in enumerate(col_specs_f1):
    cA = f"{a}_{h}"; cB = f"{b}_{h}"
    if h == 96:
        if a == "BL" or b == "BL":
            is_bl96_col[j] = True
    for i, sh in enumerate(ordered_subs_f1):
        row = f1_lookup.get((sh, cA, cB))
        if row is None: continue
        in_fam1[i, j] = True
        M1[i, j] = row.direction * row.W1_obs
        if row.tier in ("A","B","C"):
            tier_mat1[i, j] = row.tier

vmax1 = max(np.nanmax(np.abs(M1)), 0.01) if np.isfinite(np.nanmax(M1)) else 1.0

# Compute annotation glyph per row
glyph_per_row = []
for sh in ordered_subs_f1:
    glyph = ""
    for h in HPFS:
        rBL_LD = f1_lookup.get((sh, f"BL_{h}", f"LD_{h}"))
        rD_LD  = f1_lookup.get((sh, f"D_{h}",  f"LD_{h}"))
        rBL_D  = f1_lookup.get((sh, f"BL_{h}", f"D_{h}"))
        AB = lambda r: (r is not None and r.tier in ("A","B"))
        if AB(rBL_LD) and AB(rD_LD) and np.sign(rBL_LD.direction*rBL_LD.W1_obs) == np.sign(rD_LD.direction*rD_LD.W1_obs):
            glyph = "◆"; break
    if not glyph:
        # BL-specific: BL_vs_LD and BL_vs_D A/B same side for BL, D_vs_LD blank
        for h in HPFS:
            rBL_LD = f1_lookup.get((sh, f"BL_{h}", f"LD_{h}"))
            rBL_D  = f1_lookup.get((sh, f"BL_{h}", f"D_{h}"))
            rD_LD  = f1_lookup.get((sh, f"D_{h}",  f"LD_{h}"))
            AB = lambda r: (r is not None and r.tier in ("A","B"))
            NS = lambda r: (r is None or r.tier == "")
            if AB(rBL_LD) and AB(rBL_D) and NS(rD_LD):
                if np.sign(rBL_LD.direction) == np.sign(rBL_D.direction):
                    glyph = "▲"; break
    if not glyph:
        for h in HPFS:
            rD_LD  = f1_lookup.get((sh, f"D_{h}",  f"LD_{h}"))
            rBL_D  = f1_lookup.get((sh, f"BL_{h}", f"D_{h}"))
            rBL_LD = f1_lookup.get((sh, f"BL_{h}", f"LD_{h}"))
            AB = lambda r: (r is not None and r.tier in ("A","B"))
            NS = lambda r: (r is None or r.tier == "")
            # D-specific: D_vs_LD and BL_vs_D A/B with D on same side; BL_vs_LD blank
            # In BL_vs_D, cell_A=BL, cell_B=D; direction=+1 means BL higher i.e. D LOWER.
            # In D_vs_LD, cell_A=D, cell_B=LD; direction=+1 means D higher.
            # "D on same side" means D is consistently higher or consistently lower across both.
            if AB(rD_LD) and AB(rBL_D) and NS(rBL_LD):
                d_in_D_LD = np.sign(rD_LD.direction)      # +1 means D higher than LD
                d_in_BL_D = -np.sign(rBL_D.direction)     # +1 means D higher than BL
                if d_in_D_LD == d_in_BL_D and d_in_D_LD != 0:
                    glyph = "▼"; break
    glyph_per_row.append(glyph)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
# background hatched
for i in range(n_rows1):
    for j in range(n_cols1):
        if not in_fam1[i,j]:
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="#F5F5F5",
                                       hatch="////", edgecolor="#D0D0D0", lw=0))
im = ax.imshow(M1, cmap="RdBu_r", vmin=-vmax1, vmax=vmax1,
               aspect="auto", extent=(0, n_cols1, n_rows1, 0), interpolation="none")
for i in range(n_rows1):
    for j in range(n_cols1):
        t = tier_mat1[i,j]
        if t == "A":
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="none", edgecolor="black", lw=1.5))
        elif t == "B":
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="none", edgecolor="#555", lw=0.8))
        elif t == "C":
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="#D9D9D9", edgecolor="none"))
            ax.text(j+0.5, i+0.5, "c", fontsize=6, ha="center", va="center")
# BL_96 hatch overlay
for j in range(n_cols1):
    if is_bl96_col[j]:
        ax.add_patch(plt.Rectangle((j, 0), 1, n_rows1, facecolor="none",
                                   hatch=r"///", edgecolor="#555", alpha=0.25, lw=0))
# right annotation column of glyphs
for i, g in enumerate(glyph_per_row):
    if g:
        ax.text(n_cols1 + 0.3, i+0.5, g, fontsize=9, ha="center", va="center")
ax.set_xticks([j+0.5 for j in range(n_cols1)])
ax.set_xticklabels([f"{h}\n{lbl}" for (h,_,_,lbl) in col_specs_f1], fontsize=6)
ax.set_yticks([i+0.5 for i in range(n_rows1)])
ax.set_yticklabels(ordered_subs_f1, fontsize=6)
# class rules
prev_class = None
for i, s in enumerate(ordered_subs_f1):
    cls = classify(s)
    if cls != prev_class:
        if prev_class is not None: ax.axhline(i, color="black", lw=0.6)
        prev_class = cls
ax.set_xlim(0, n_cols1 + 1.2); ax.set_ylim(n_rows1, 0)
cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.08)
cbar.set_label("signed W1 (+ = first-named condition higher)", fontsize=7)
cbar.ax.tick_params(labelsize=6)
ax.set_title("Fig 3 — F1 photoperiod axis", fontsize=9, loc="left")

legend_patches = [
    mpatches.Patch(facecolor="none", edgecolor="black", lw=1.5, label="Tier A"),
    mpatches.Patch(facecolor="none", edgecolor="#555", lw=0.8, label="Tier B"),
    mpatches.Patch(facecolor="#D9D9D9", edgecolor="none", label="Tier C"),
    mpatches.Patch(facecolor="#F5F5F5", hatch="////", edgecolor="#D0D0D0", label="not in family"),
    mpatches.Patch(facecolor="white", hatch=r"///", edgecolor="#555", alpha=0.25,
                   label="BL_96 extraction flag"),
]
ax.legend(handles=legend_patches, loc="lower center", bbox_to_anchor=(0.5, -0.12),
          fontsize=6, frameon=False, ncol=5)
# glyph legend below
fig.text(0.5, 0.015,
         "◆ shared BL-vs-LD and D-vs-LD   ▲ BL-specific   ▼ D-specific",
         fontsize=7, ha="center")

fig.savefig(OUT/"fig3.pdf", bbox_inches="tight")
fig.savefig(OUT/"fig3.png", dpi=300, bbox_inches="tight")
plt.close(fig)
rows_out_f1 = []
for i, sh in enumerate(ordered_subs_f1):
    for j, (h, a, b, lbl) in enumerate(col_specs_f1):
        rows_out_f1.append(dict(subsystem=sh, cls=classify(sh), hpf=h, contrast=lbl,
                                signed_W1=M1[i,j], tier=tier_mat1[i,j],
                                in_family=bool(in_fam1[i,j])))
pd.DataFrame(rows_out_f1).to_csv(OUT/"fig3_data.csv", index=False)
print(f"[write] fig3.pdf + png + _data.csv  ({n_rows1} rows × {n_cols1} cols)")


# =========================================================================
# FIG 5 — growth capacity attribution
# =========================================================================
fig = plt.figure(figsize=(13, 4))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.2, 1], wspace=0.3)

# ----- 5a: child vs parent bars -----
ax = fig.add_subplot(gs[0,0])
width = 0.25
xs = np.arange(len(HPFS))
for k, c in enumerate(CONDS):
    children = [bm_cell[f"{c}_{h}"] for h in HPFS]
    parents  = [parent_cell[f"{c}_{h}"] for h in HPFS]
    ax.bar(xs + (k-1)*width, parents, width, facecolor="none", edgecolor=CCOL[c], lw=0.8, alpha=0.5,
           label=f"{c} parent" if k==0 else None)
    ax.bar(xs + (k-1)*width, children, width, color=CCOL[c], label=c)
ax.set_xticks(xs); ax.set_xticklabels(HPFS)
ax.set_ylabel("biomass_max"); ax.set_xlabel("hpf")
ax.legend(fontsize=7, ncol=2, frameon=False)
ax.set_title("(a) child vs parent biomass_max", fontsize=9, loc="left")
for sp in ("top","right"): ax.spines[sp].set_visible(False)

# ----- 5b: 2x2 attribution stacked bars -----
ax = fig.add_subplot(gs[0,1])
tt = tt2x2.sort_values(["hpf","pair"]).reset_index(drop=True)
y = np.arange(len(tt))[::-1]
abs_top = tt.topology_effect.abs().values
abs_env = tt.envelope_effect.abs().values
abs_int = tt.interaction.abs().values
no_mat = (tt.gap_pct_of_smaller_corner < 5).values
hatch_pat = ["////" if nm else None for nm in no_mat]
# stacked bars
ax.barh(y, abs_top, color="#CC4444", edgecolor="black", lw=0.3, label="|topology|")
ax.barh(y, abs_env, left=abs_top, color="#4477CC", edgecolor="black", lw=0.3, label="|envelope|")
ax.barh(y, abs_int, left=abs_top+abs_env, color="#999999", edgecolor="black", lw=0.3, label="|interaction|")
for i, nm in enumerate(no_mat):
    yi = y[i]; total = abs_top[i]+abs_env[i]+abs_int[i]
    if nm:
        ax.barh(yi, max(total,0.05), color="none", edgecolor="black", hatch="////", lw=0.3)
        ax.text(max(total,0.05)+0.05, yi, "no material gap", va="center", fontsize=6)
    else:
        verd = tt.iloc[i].dominance
        ax.text(total+0.05, yi, verd, va="center", fontsize=6)
ax.set_yticks(y); ax.set_yticklabels(tt.pair, fontsize=6)
ax.set_xlabel("|effect| (biomass units)")
ax.legend(fontsize=7, frameon=False, loc="lower right")
ax.set_title("(b) 2×2 topology × envelope decomposition", fontsize=9, loc="left")
for sp in ("top","right"): ax.spines[sp].set_visible(False)

# ----- 5c: BL_96 greedy curve -----
ax = fig.add_subplot(gs[0,2])
steps = [0] + greedy.step.tolist()
biomax = [8.535] + greedy.biomax_after.tolist()
ax.plot(steps, biomax, "-o", color="#1F77B4", markersize=5, lw=1.5)
ax.axhline(12.928, ls="--", color="#F0A800", lw=1, label="LD_96 = 12.93")
ax.axhline(12.255, ls="--", color="#555", lw=1, label="ceiling = 12.26")
ax.set_xlabel("step"); ax.set_ylabel("BL_96 biomass_max")
ax.set_xticks(steps)
# annotate step 1
ax.annotate("MAR00466\nPS cytosol→ER\n(no GPR)",
            xy=(1, biomax[1]), xytext=(1.3, biomax[1]-0.9),
            fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.5))
ax.legend(fontsize=7, frameon=False, loc="lower right")
ax.set_title("(c) BL_96 greedy add-back", fontsize=9, loc="left")
for sp in ("top","right"): ax.spines[sp].set_visible(False)

fig.savefig(OUT/"fig5.pdf", bbox_inches="tight")
fig.savefig(OUT/"fig5.png", dpi=300, bbox_inches="tight")
plt.close(fig)
# data csv: concatenate three sub-tables
pd5a = pd.DataFrame([dict(panel="5a", cell=cell, condition=cell.split("_")[0],
                          hpf=int(cell.split("_")[1]),
                          parent_biomax=parent_cell[cell], child_biomax=bm_cell[cell])
                     for cell in [f"{c}_{h}" for c in CONDS for h in HPFS]])
pd5b = tt[["pair","hpf","topology_effect","envelope_effect","interaction",
           "gap_pct_of_smaller_corner","dominance"]].assign(panel="5b")
pd5c = pd.DataFrame({"panel":"5c","step":steps,"biomax_after":biomax})
pd.concat([pd5a, pd5b, pd5c], ignore_index=True).to_csv(OUT/"fig5_data.csv", index=False)
print("[write] fig5.pdf + png + _data.csv")


# =========================================================================
# figS4 — Rhat/ESS violins
# =========================================================================
fig = plt.figure(figsize=(10, 4))
gs = fig.add_gridspec(1, 2, wspace=0.25)
# Rhat: one point per cell; show "distribution" via loading rhat_rank.csv per cell
rhat_dir = ROOT/"grace_v4/results/hopsy_v7_local"
cells = sorted(rhat.cell.tolist())
rhat_data = []
ess_data = []
labels = []
for c in cells:
    fp = rhat_dir/f"{c}_rhat_rank.csv"
    if fp.exists():
        d = pd.read_csv(fp)
        if "rhat_rank" in d.columns:
            rh = d.rhat_rank.dropna().values
        elif "rhat" in d.columns:
            rh = d.rhat.dropna().values
        else:
            rh = np.array([])
        if "ess_bulk" in d.columns:
            es = d.ess_bulk.dropna().values
        else:
            es = np.array([])
        rhat_data.append(rh); ess_data.append(es); labels.append(c)

ax = fig.add_subplot(gs[0,0])
parts = ax.violinplot(rhat_data, showmedians=True)
for i, pc in enumerate(parts["bodies"]):
    cond = labels[i].split("_")[0]
    pc.set_facecolor(CCOL[cond]); pc.set_alpha(0.6)
ax.set_xticks(range(1, len(labels)+1)); ax.set_xticklabels(labels, rotation=90, fontsize=6)
ax.set_ylabel("rank-R̂"); ax.set_ylim(1.0, 1.3)
ax.axhline(1.05, color="#555", ls="--", lw=0.5)
ax.axhline(1.01, color="#AAA", ls=":", lw=0.5)
ax.set_title("S4 (a) rank-R̂ per cell", fontsize=9, loc="left")
for sp in ("top","right"): ax.spines[sp].set_visible(False)

ax = fig.add_subplot(gs[0,1])
parts = ax.violinplot(ess_data, showmedians=True)
for i, pc in enumerate(parts["bodies"]):
    cond = labels[i].split("_")[0]
    pc.set_facecolor(CCOL[cond]); pc.set_alpha(0.6)
ax.set_xticks(range(1, len(labels)+1)); ax.set_xticklabels(labels, rotation=90, fontsize=6)
ax.set_ylabel("bulk ESS")
ax.axhline(400, color="#555", ls="--", lw=0.5)
ax.axhline(100, color="#AAA", ls=":", lw=0.5)
ax.set_title("S4 (b) bulk ESS per cell", fontsize=9, loc="left")
for sp in ("top","right"): ax.spines[sp].set_visible(False)

fig.savefig(OUT/"figS4.pdf", bbox_inches="tight")
fig.savefig(OUT/"figS4.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("[write] figS4.pdf + png")


# =========================================================================
# figS5 — mean vs median concordance stacked bar
# =========================================================================
mean_hits = set(zip(tA.contrast, tA.subsystem))
med_f1 = pd.read_csv(MED/"findings_family1.csv")
med_f2 = pd.read_csv(MED/"findings_family2.csv")
med = pd.concat([med_f1, med_f2], ignore_index=True)
med_hits = set(zip(med[med.q_bh<=0.10].contrast, med[med.q_bh<=0.10].subsystem))

only_mean = len(mean_hits - med_hits)
only_med  = len(med_hits - mean_hits)
both      = len(mean_hits & med_hits)

fig, ax = plt.subplots(figsize=(5, 3.5))
cats = ["mean only", "median only", "both"]
vals = [only_mean, only_med, both]
cols = ["#CC4444", "#999999", "#1F77B4"]
ax.bar(cats, vals, color=cols)
for xi, v in zip(cats, vals):
    ax.text(xi, v+1, str(v), ha="center", fontsize=8)
ax.set_ylabel("hit count")
ax.set_title("S5 mean (method-of-record) vs median concordance", fontsize=9, loc="left")
for sp in ("top","right"): ax.spines[sp].set_visible(False)
fig.savefig(OUT/"figS5.pdf", bbox_inches="tight")
fig.savefig(OUT/"figS5.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("[write] figS5.pdf + png")

# =========================================================================
# assertions
# =========================================================================
# Every Tier A row (62) is in Fig 2 (F2) or Fig 3 (F1)
tA_subs_f2 = set(tA[tA.family=="F2"]["short"])
tA_subs_f1 = set(tA[tA.family=="F1"]["short"])
missing_f2 = tA_subs_f2 - set(ordered_subs)
missing_f1 = tA_subs_f1 - set(ordered_subs_f1)
if missing_f2 or missing_f1:
    print(f"ASSERT WARN: F2 Tier A subs missing: {missing_f2}; F1: {missing_f1}")
else:
    print("[assert] every Tier A row appears in Fig 2 or Fig 3 — OK")

# no value < -30 in log2 panels — we don't have log2 panels in this scope
print("[assert] no CDF log2 panels in scope")

print("\nDONE. All deliverables in:", OUT)
for f in sorted(OUT.iterdir()): print("  ", f.name)
