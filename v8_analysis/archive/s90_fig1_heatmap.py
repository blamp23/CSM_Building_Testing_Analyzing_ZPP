"""
s90 -- Fig 1: F1 + F2 overview heatmaps for Tier A/B.

For each family:
  rows    = subsystems appearing in Tier A or B (dropping Tier-C and non-sig rows)
  columns = contrasts in the family (hpf × condition-pair for F1; condition × transition for F2)
  cell    = signed W1 (direction × W1_obs), diverging colormap
            grey        = not in §1 family
            white       = in family but not significant (q > 0.10)
  outline = thick red  = Tier A
            thin red   = Tier B
  marks   = ●  on cells with floor_resolved=survives
            ╱╲ hatch on cells with floor_exposed=yes (pending matched floor)

Output: results/findings_v7_mean_dropzero/plots/fig1_F1_heatmap.pdf, fig1_F2_heatmap.pdf
"""
from pathlib import Path
import numpy as np, pandas as pd, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

ROOT = Path(__file__).resolve().parents[1]
MOR  = ROOT/"results/findings_v7_mean_dropzero"
OUT  = MOR/"plots"
OUT.mkdir(exist_ok=True, parents=True)

f1 = pd.read_csv(MOR/"findings_family1.csv"); f1['family']='F1'
f2 = pd.read_csv(MOR/"findings_family2.csv"); f2['family']='F2'
all_rows = pd.concat([f1,f2], ignore_index=True)

tA = pd.read_csv(MOR/"tierA_final.csv")
tB = pd.read_csv(MOR/"tierB_final.csv")
tier_rows = pd.concat([tA.assign(tier='A'), tB.assign(tier='B')], ignore_index=True)

HPFS = [24, 48, 72, 96, 120]
CONDS = ["BL", "D", "LD"]

def row_order(df_tier):
    """Order subsystems by n Tier A desc then Tier B desc."""
    counts = df_tier.groupby(['subsystem','tier']).size().unstack(fill_value=0)
    counts['total'] = counts.sum(axis=1)
    counts = counts.sort_values(['A' if 'A' in counts.columns else 'B', 'total'],
                                 ascending=False)
    return list(counts.index)

def build_F1_matrix(all_rows_fam, tier_rows_fam):
    pairs = [('BL','D'), ('BL','LD'), ('D','LD')]
    col_labels = [f'{a}vs{b}\n{h}hpf' for h in HPFS for a,b in pairs]
    col_keys   = [(h, a, b) for h in HPFS for a,b in pairs]
    subs = row_order(tier_rows_fam)
    W   = np.full((len(subs), len(col_keys)), np.nan)          # signed W1
    T   = np.full((len(subs), len(col_keys)), "", dtype=object) # tier A/B/''
    FR  = np.zeros((len(subs), len(col_keys)), dtype=bool)     # floor-resolved mark
    FX  = np.zeros((len(subs), len(col_keys)), dtype=bool)     # floor-exposed hatch
    IN  = np.zeros((len(subs), len(col_keys)), dtype=bool)     # row in §1 family
    for i, sub in enumerate(subs):
        for j, (h, a, b) in enumerate(col_keys):
            c = f'{a}_{h}_vs_{b}_{h}'
            r = all_rows_fam[(all_rows_fam.contrast==c) & (all_rows_fam.subsystem==sub)]
            if len(r) == 0: continue
            r = r.iloc[0]
            IN[i,j] = True
            if r.q_bh <= 0.10:
                W[i,j] = r.direction * r.W1_obs
            # Tier marks come from tier_rows
            tr = tier_rows_fam[(tier_rows_fam.contrast==c) & (tier_rows_fam.subsystem==sub)]
            if len(tr): T[i,j] = tr.iloc[0].tier
            if 'floor_resolved' in r and 'survives' in str(r.get('floor_resolved','')): FR[i,j] = True
            if 'floor_exposed' in r.index and str(r.get('floor_exposed','')).startswith('yes'): FX[i,j] = True
    return subs, col_labels, W, T, FR, FX, IN

def build_F2_matrix(all_rows_fam, tier_rows_fam):
    transitions = [(HPFS[i], HPFS[i+1]) for i in range(4)]
    col_labels = [f'{c}\n{a}→{b}' for c in CONDS for a,b in transitions]
    col_keys   = [(c, a, b) for c in CONDS for a,b in transitions]
    subs = row_order(tier_rows_fam)
    W   = np.full((len(subs), len(col_keys)), np.nan)
    T   = np.full((len(subs), len(col_keys)), "", dtype=object)
    FR  = np.zeros((len(subs), len(col_keys)), dtype=bool)
    FX  = np.zeros((len(subs), len(col_keys)), dtype=bool)
    IN  = np.zeros((len(subs), len(col_keys)), dtype=bool)
    for i, sub in enumerate(subs):
        for j, (cond, a, b) in enumerate(col_keys):
            c = f'{cond}_{a}_vs_{cond}_{b}'
            r = all_rows_fam[(all_rows_fam.contrast==c) & (all_rows_fam.subsystem==sub)]
            if len(r) == 0: continue
            r = r.iloc[0]
            IN[i,j] = True
            if r.q_bh <= 0.10:
                W[i,j] = r.direction * r.W1_obs
            tr = tier_rows_fam[(tier_rows_fam.contrast==c) & (tier_rows_fam.subsystem==sub)]
            if len(tr): T[i,j] = tr.iloc[0].tier
            if 'floor_resolved' in r and 'survives' in str(r.get('floor_resolved','')): FR[i,j] = True
            if 'floor_exposed' in r.index and str(r.get('floor_exposed','')).startswith('yes'): FX[i,j] = True
    return subs, col_labels, W, T, FR, FX, IN

def draw_heatmap(subs, col_labels, W, T, FR, FX, IN, title, out_fp):
    n_rows, n_cols = len(subs), len(col_labels)
    fig_h = max(4, 0.28 * n_rows + 1.5)
    fig_w = max(6, 0.55 * n_cols + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Background: grey for not-in-family, white for in-family-not-sig, signed W1 colormap otherwise
    vmax = np.nanmax(np.abs(W)) if np.isfinite(W).any() else 1.0
    vmax = max(vmax, 0.5)
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")

    for i in range(n_rows):
        for j in range(n_cols):
            x0, y0 = j, n_rows - i - 1
            if not IN[i,j]:
                ax.add_patch(mpatches.Rectangle((x0, y0), 1, 1, facecolor="#dddddd",
                                                 edgecolor='#999', linewidth=0.3))
                continue
            if np.isnan(W[i,j]):
                ax.add_patch(mpatches.Rectangle((x0, y0), 1, 1, facecolor="white",
                                                 edgecolor='#999', linewidth=0.3))
            else:
                face = cmap(norm(W[i,j]))
                ax.add_patch(mpatches.Rectangle((x0, y0), 1, 1, facecolor=face,
                                                 edgecolor='#999', linewidth=0.3))
            if FX[i,j]:
                ax.add_patch(mpatches.Rectangle((x0, y0), 1, 1, facecolor='none',
                                                 edgecolor='black', hatch='////',
                                                 linewidth=0.0))
            if T[i,j] == 'A':
                ax.add_patch(mpatches.Rectangle((x0+.04, y0+.04), 0.92, 0.92,
                                                 facecolor='none', edgecolor='red',
                                                 linewidth=1.8))
            elif T[i,j] == 'B':
                ax.add_patch(mpatches.Rectangle((x0+.06, y0+.06), 0.88, 0.88,
                                                 facecolor='none', edgecolor='red',
                                                 linewidth=0.7, linestyle='--'))
            if FR[i,j]:
                ax.plot(x0 + 0.5, y0 + 0.5, 'o', color='black', ms=4,
                        markerfacecolor='white', markeredgewidth=1.2)

    ax.set_xlim(0, n_cols); ax.set_ylim(0, n_rows)
    ax.set_xticks(np.arange(n_cols) + 0.5); ax.set_xticklabels(col_labels, fontsize=7, rotation=45, ha='right')
    short = [s if len(s) <= 42 else s[:39]+'…' for s in subs]
    ax.set_yticks(np.arange(n_rows)[::-1] + 0.5); ax.set_yticklabels(short, fontsize=7)
    ax.set_title(title, fontsize=11, pad=10)

    cax = fig.add_axes([0.92, 0.3, 0.015, 0.4])
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label('signed W1\n(direction × W1_obs)', fontsize=7)
    cb.ax.tick_params(labelsize=6)

    legend_elems = [
        mpatches.Patch(facecolor='#dddddd', edgecolor='#999', label='not in §1 family'),
        mpatches.Patch(facecolor='white', edgecolor='#999', label='in family, q > 0.10'),
        mpatches.Patch(facecolor='none', edgecolor='red', linewidth=1.6, label='Tier A outline'),
        mpatches.Patch(facecolor='none', edgecolor='red', linewidth=0.8, linestyle='--', label='Tier B outline'),
        mpatches.Patch(facecolor='none', edgecolor='black', hatch='////', label='floor-exposed (pending)'),
        plt.Line2D([0], [0], marker='o', color='black', mfc='white', mec='black', ms=5, lw=0, label='floor-resolved: survives'),
    ]
    ax.legend(handles=legend_elems, loc='upper left', bbox_to_anchor=(1.015, 0.98),
              fontsize=6, frameon=False)

    plt.tight_layout(rect=[0, 0, 0.9, 1])
    plt.savefig(out_fp, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f"[write] {out_fp}", flush=True)

# ============ Family 1 ============
all_f1 = all_rows[all_rows.family == 'F1']
tier_f1 = tier_rows[tier_rows.family == 'F1']
subs, labels, W, T, FR, FX, IN = build_F1_matrix(all_f1, tier_f1)
print(f"[F1] {len(subs)} subsystems × {len(labels)} contrasts", flush=True)
draw_heatmap(subs, labels, W, T, FR, FX, IN,
             "Fig 1a — Family 1 (between-condition at hpf): signed W1 heatmap",
             OUT/"fig1_F1_heatmap.pdf")

# ============ Family 2 ============
all_f2 = all_rows[all_rows.family == 'F2']
tier_f2 = tier_rows[tier_rows.family == 'F2']
subs, labels, W, T, FR, FX, IN = build_F2_matrix(all_f2, tier_f2)
print(f"[F2] {len(subs)} subsystems × {len(labels)} contrasts", flush=True)
draw_heatmap(subs, labels, W, T, FR, FX, IN,
             "Fig 1b — Family 2 (within-condition developmental): signed W1 heatmap",
             OUT/"fig1_F2_heatmap.pdf")

print("DONE")
