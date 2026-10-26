"""
s90b -- Fig 1a (faceted): three separate F1 heatmaps, one per condition pair.
  BL vs D   — subsystems × 5 hpf
  BL vs LD  — subsystems × 5 hpf
  D  vs LD  — subsystems × 5 hpf
Same encoding as s90: signed W1 fill, Tier outlines, floor-resolved dots,
floor-exposed hatch, grey=out-of-family, white=in-family-not-sig.
Combined into one PDF with 3 pages (and one combined 3-panel figure for inspection).
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

# Condition palette (matches Figs 2 & 3)
COND_COLORS = {
    "BL": "#1f77b4",   # blue
    "D":  "#000000",   # black
    "LD": "#f4c430",   # yellow
}
# Neutral centre — pure white gives maximum contrast against both ends
CENTRE = "#ffffff"

def pair_cmap(a_cond, b_cond):
    """Diverging cmap: positive signed_W1 (A higher) → A's color,
    zero → CENTRE, negative (B higher) → B's color.
    So the two ends of the colorbar directly correspond to the two
    contrasted conditions."""
    neg_color = COND_COLORS[b_cond]   # A < B
    pos_color = COND_COLORS[a_cond]   # A > B
    return mcolors.LinearSegmentedColormap.from_list(
        f"{a_cond}_v_{b_cond}", [neg_color, CENTRE, pos_color], N=256)

ROOT = Path(__file__).resolve().parents[1]
MOR  = ROOT/"results/findings_v7_mean_dropzero"
OUT  = MOR/"plots"
OUT.mkdir(exist_ok=True, parents=True)

f1 = pd.read_csv(MOR/"findings_family1.csv").assign(family='F1')
tA = pd.read_csv(MOR/"tierA_final.csv")
tB = pd.read_csv(MOR/"tierB_final.csv")
tier_rows = pd.concat([tA.assign(tier='A'), tB.assign(tier='B')], ignore_index=True)
tier_f1 = tier_rows[tier_rows.family == 'F1']

HPFS = [24, 48, 72, 96, 120]
PAIRS = [('BL','D'), ('BL','LD'), ('D','LD')]

def row_order(sub_set):
    """Order subsystems by # Tier A desc then alphabetical."""
    counts = tier_f1[tier_f1.subsystem.isin(sub_set)].groupby(['subsystem','tier']).size().unstack(fill_value=0)
    for col in ('A', 'B'):
        if col not in counts.columns: counts[col] = 0
    counts = counts.reindex(sorted(sub_set))
    counts['order'] = counts['A']*100 + counts['B']
    return list(counts.sort_values('order', ascending=False).index)

def build_matrix_pair(a_cond, b_cond):
    """subsystems × 5 hpf matrix for one condition pair."""
    rows_pair = f1[f1.contrast.str.startswith(f"{a_cond}_") & f1.contrast.str.contains(f"_vs_{b_cond}_")]
    # restrict to Tier A/B subsystems in this family
    sub_set = set(tier_f1.subsystem) & set(rows_pair.subsystem)
    subs = row_order(sub_set)
    n_rows = len(subs)
    W   = np.full((n_rows, 5), np.nan)
    T   = np.full((n_rows, 5), "", dtype=object)
    FR  = np.zeros((n_rows, 5), dtype=bool)
    FX  = np.zeros((n_rows, 5), dtype=bool)
    IN  = np.zeros((n_rows, 5), dtype=bool)
    for i, sub in enumerate(subs):
        for j, h in enumerate(HPFS):
            c = f'{a_cond}_{h}_vs_{b_cond}_{h}'
            r = rows_pair[(rows_pair.contrast==c) & (rows_pair.subsystem==sub)]
            if len(r) == 0: continue
            r = r.iloc[0]
            IN[i,j] = True
            # Color every in-family cell by direction × W1_obs; the red outline
            # (Tier A/B) signals which ones clear BH.
            if pd.notna(r.W1_obs) and pd.notna(r.direction):
                W[i,j] = r.direction * r.W1_obs
            tr = tier_f1[(tier_f1.contrast==c) & (tier_f1.subsystem==sub)]
            if len(tr): T[i,j] = tr.iloc[0].tier
            if 'floor_resolved' in r and 'survives' in str(r.get('floor_resolved','')): FR[i,j] = True
            if 'floor_exposed' in r.index and str(r.get('floor_exposed','')).startswith('yes'): FX[i,j] = True
    return subs, W, T, FR, FX, IN

def draw(ax, subs, W, T, FR, FX, IN, title, vmax, cmap):
    n_rows, n_cols = len(subs), 5
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    for i in range(n_rows):
        for j in range(n_cols):
            x0, y0 = j, n_rows - i - 1
            if not IN[i,j]:
                ax.add_patch(mpatches.Rectangle((x0,y0),1,1, facecolor="#dddddd",
                                                 edgecolor='#999', linewidth=0.3))
                continue
            if np.isnan(W[i,j]):
                ax.add_patch(mpatches.Rectangle((x0,y0),1,1, facecolor="white",
                                                 edgecolor='#999', linewidth=0.3))
            else:
                face = cmap(norm(W[i,j]))
                ax.add_patch(mpatches.Rectangle((x0,y0),1,1, facecolor=face,
                                                 edgecolor='#999', linewidth=0.3))
            if FX[i,j]:
                ax.add_patch(mpatches.Rectangle((x0,y0),1,1, facecolor='none',
                                                 edgecolor='black', hatch='////',
                                                 linewidth=0.0))
            if T[i,j] == 'A':
                ax.add_patch(mpatches.Rectangle((x0+.04,y0+.04), 0.92, 0.92,
                                                 facecolor='none', edgecolor='red', linewidth=1.8))
            elif T[i,j] == 'B':
                ax.add_patch(mpatches.Rectangle((x0+.06,y0+.06), 0.88, 0.88,
                                                 facecolor='none', edgecolor='red', linewidth=0.7, linestyle='--'))
            if FR[i,j]:
                ax.plot(x0+0.5, y0+0.5, 'o', color='black', ms=4,
                        markerfacecolor='white', markeredgewidth=1.2)
    ax.set_xlim(0, n_cols); ax.set_ylim(0, n_rows)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels([f'{h} hpf' for h in HPFS], fontsize=7)
    short = [s if len(s) <= 42 else s[:39]+'…' for s in subs]
    ax.set_yticks(np.arange(n_rows)[::-1] + 0.5); ax.set_yticklabels(short, fontsize=7)
    ax.set_title(title, fontsize=10, pad=8)
    return norm

# compute shared vmax across all 3 pairs for consistent color scale
vmax_global = 0
for a, b in PAIRS:
    _, W, *_ = build_matrix_pair(a, b)
    if np.isfinite(W).any():
        vmax_global = max(vmax_global, np.nanmax(np.abs(W)))
vmax_global = max(vmax_global, 0.5)
print(f"[color scale] shared vmax = {vmax_global:.2f}", flush=True)

# ---- individual single-pair PDFs
for a, b in PAIRS:
    subs, W, T, FR, FX, IN = build_matrix_pair(a, b)
    cmap = pair_cmap(a, b)
    fig_h = max(4, 0.3*len(subs) + 1.5)
    fig, ax = plt.subplots(figsize=(6, fig_h))
    norm = draw(ax, subs, W, T, FR, FX, IN, f"Fig 1a — F1: {a} vs {b}", vmax_global, cmap)
    cax = fig.add_axes([0.92, 0.3, 0.015, 0.4])
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(f'signed W1\n← {b} higher     {a} higher →', fontsize=7, labelpad=6)
    cb.ax.tick_params(labelsize=6)
    plt.tight_layout(rect=[0, 0, 0.9, 1])
    out_fp = OUT/f"fig1_F1_{a}_vs_{b}.pdf"
    plt.savefig(out_fp, dpi=140, bbox_inches='tight'); plt.close(fig)
    print(f"[write] {out_fp}  rows={len(subs)}", flush=True)

# ---- combined side-by-side 3-panel figure
# use the max row count for a shared row axis
all_subs = sorted(set.union(*[set(build_matrix_pair(a,b)[0]) for a,b in PAIRS]))
all_subs = row_order(set(all_subs))
# build matrices aligned on all_subs
mats = []
for a, b in PAIRS:
    _, W_p, T_p, FR_p, FX_p, IN_p = build_matrix_pair(a, b)
    sub_p, _, _, _, _, _ = build_matrix_pair(a, b)
    # realign to all_subs
    W_a = np.full((len(all_subs), 5), np.nan)
    T_a = np.full((len(all_subs), 5), "", dtype=object)
    FR_a = np.zeros((len(all_subs), 5), dtype=bool)
    FX_a = np.zeros((len(all_subs), 5), dtype=bool)
    IN_a = np.zeros((len(all_subs), 5), dtype=bool)
    for k, s in enumerate(all_subs):
        if s in sub_p:
            ki = sub_p.index(s)
            W_a[k] = W_p[ki]; T_a[k] = T_p[ki]
            FR_a[k] = FR_p[ki]; FX_a[k] = FX_p[ki]; IN_a[k] = IN_p[ki]
    mats.append((a, b, W_a, T_a, FR_a, FX_a, IN_a))

fig_h = max(4, 0.3*len(all_subs) + 1.5)
fig = plt.figure(figsize=(16, fig_h))
# Three heatmap axes + 3 narrow colorbar axes beside them + legend on far right
gs = fig.add_gridspec(1, 7, width_ratios=[1.0, 0.07, 1.0, 0.07, 1.0, 0.07, 0.65],
                       wspace=0.3)
norms = []
for col, (a, b, W, T, FR, FX, IN) in enumerate(mats):
    ax = fig.add_subplot(gs[0, col*2])
    cmap = pair_cmap(a, b)
    norm = draw(ax, all_subs, W, T, FR, FX, IN, f"{a} vs {b}", vmax_global, cmap)
    norms.append(norm)
    if col > 0:
        ax.set_yticklabels([])
    else:
        ax.set_yticks(np.arange(len(all_subs))[::-1] + 0.5)
        ax.set_yticklabels([s if len(s) <= 42 else s[:39]+'…' for s in all_subs], fontsize=7)
    # per-panel colorbar immediately to the right
    cax = fig.add_subplot(gs[0, col*2+1])
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(f'← {b} higher    {a} higher →', fontsize=7, labelpad=4)
    cb.ax.tick_params(labelsize=6)

legend_ax = fig.add_subplot(gs[0, -1]); legend_ax.axis('off')
legend_elems = [
    mpatches.Patch(facecolor='#dddddd', edgecolor='#999', label='not in §1 family'),
    mpatches.Patch(facecolor='white', edgecolor='#999', label='in family, no W1 computed'),
    mpatches.Patch(facecolor=CENTRE, edgecolor='#999',
                   label='in family, W1 ≈ 0 (BH may not pass)'),
    mpatches.Patch(facecolor='none', edgecolor='red', linewidth=1.6,
                   label='Tier A — clears BH + full funnel'),
    mpatches.Patch(facecolor='none', edgecolor='red', linewidth=0.8, linestyle='--',
                   label='Tier B — clears BH, partial funnel'),
    mpatches.Patch(facecolor='none', edgecolor='black', hatch='////',
                   label='floor-exposed (pending)'),
    plt.Line2D([0],[0], marker='o', color='black', mfc='white', mec='black', ms=5, lw=0,
               label='floor-resolved: survives'),
]
legend_ax.legend(handles=legend_elems, loc='upper left', fontsize=7, frameon=False,
                  bbox_to_anchor=(-0.1, 0.95))

fig.suptitle("Fig 1a — Family 1 heatmap faceted by condition pair\n"
             "(BL=blue, D=black, LD=yellow — colorbar ends match the two cells of each contrast)",
             fontsize=11, y=0.98)
combo_fp = OUT/"fig1_F1_faceted.pdf"
plt.savefig(combo_fp, dpi=140, bbox_inches='tight'); plt.close(fig)
print(f"[write] {combo_fp}  rows={len(all_subs)}", flush=True)
