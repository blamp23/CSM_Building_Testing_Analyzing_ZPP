"""
s92 -- Fig 3: developmental axis, one row per Tier A subsystem,
3 columns BL / D / LD.
  Each column: 5-hpf CDFs on one axis, hpf as sequential ramp.
    Normaliser μ_r = mean of non-zero values pooled over the 5 cells of
    that column's condition (within-condition panel-pooled).
  Bracket significant adjacent transitions above the panel with their W1.
    Hatched bracket for floor-exposed, pending matched-floor.
  Below each column: centroid trajectory — median log₂ relative flux
    with IQR band vs hpf.

Palette: hpf sequential ramp (light grey → dark) per condition.
Reps only, loop-flagged excluded, zeros dropped.

Output: results/findings_v7_mean_dropzero/plots/fig3_developmental_axis.pdf
"""
from pathlib import Path
import numpy as np, pandas as pd, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
POOL = ROOT/"results/pooled_v7"
SUBS = ROOT/"results/subsets_v7"
MOR  = ROOT/"results/findings_v7_mean_dropzero"
OUT  = MOR/"plots"
OUT.mkdir(exist_ok=True, parents=True)

CONDS = ["BL", "D", "LD"]
HPFS = [24, 48, 72, 96, 120]
ZERO_TOL = 1e-10
LOOP_MAX = 100.0

# Base hues for each condition (hpf ramp interpolates to this)
BASE_HUE = {"BL": "#1f77b4", "D": "#000000", "LD": "#ff7f0e"}
def hpf_ramp(cond, hpf):
    """Light → dark ramp per condition."""
    base = np.array(mcolors.to_rgb(BASE_HUE[cond]))
    t = (hpf - 24) / (120 - 24)   # 0 at 24 hpf, 1 at 120 hpf
    # start at 70% white blend, end at pure base
    light = np.array([0.9, 0.9, 0.9])
    return tuple((1 - 0.3 - 0.7*t) * light + (0.3 + 0.7*t) * base)

tA = pd.read_csv(MOR/"tierA_final.csv")
sig_subs = sorted(tA.subsystem.unique())
print(f"[plan] {len(sig_subs)} Tier A subsystems", flush=True)

# ---------- caches
pooled_cache, mem_cache = {}, {}
def get_pooled(cell):
    if cell not in pooled_cache:
        pooled_cache[cell] = pd.read_csv(POOL/f"{cell}_pooled.csv", index_col=0)
    return pooled_cache[cell]
def get_mem(cell):
    if cell not in mem_cache:
        m = pd.read_csv(SUBS/f"{cell}_subset_membership.csv")
        mem_cache[cell] = m.drop_duplicates("representative").set_index("representative")
    return mem_cache[cell]

def panel_data_within_condition(subsystem, cond):
    """CDF data per hpf for one condition, with μ_r pooled across the 5 cells."""
    cells = [f"{cond}_{h}" for h in HPFS]
    # Intersect reps across all 5 hpfs
    rep_sets = []
    for cell in cells:
        m = get_mem(cell)
        rep_sets.append(set(m.index[m.subsystem == subsystem].intersection(get_pooled(cell).index)))
    common = sorted(set.intersection(*rep_sets)) if rep_sets else []
    if len(common) < 3: return None
    Xs = [np.abs(get_pooled(cell).loc[common].values) for cell in cells]
    max_abs = np.stack([x.max(1) for x in Xs]).max(0)
    keep = max_abs <= LOOP_MAX
    if keep.sum() < 3: return None
    Xs = [x[keep] for x in Xs]
    all_pooled = np.concatenate(Xs, axis=1)
    masked = np.where(all_pooled >= ZERO_TOL, all_pooled, np.nan)
    with np.errstate(all="ignore"):
        mu = np.nanmean(masked, axis=1)
    mu = np.where(np.isfinite(mu) & (mu > 0), mu, 1.0)
    out = dict(hpfs=HPFS, cells=cells, mu_rep=mu, n_reps=keep.sum(),
                log_ratio={}, median_logv={}, iqr_lo={}, iqr_hi={})
    for h, X in zip(HPFS, Xs):
        mask = X >= ZERO_TOL
        rep_idx, samp_idx = np.where(mask)
        lr = np.log2(X[rep_idx, samp_idx] / mu[rep_idx])
        out['log_ratio'][h] = lr
        # per-rep median → centroid trajectory
        medians = np.array([np.log2(np.median(X[r][mask[r]]) / mu[r]) if mask[r].any() else np.nan
                             for r in range(X.shape[0])])
        out['median_logv'][h] = medians
    return out

# Tier A hits per (subsystem, condition, transition) — for brackets
tA_dict = {}  # (sub, cond, (hpf_a, hpf_b)) → row
for _, r in tA.iterrows():
    m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_(\d+)", r.contrast)
    cA, hA, cB, hB = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
    if r.family == 'F2' and cA == cB:
        tA_dict[(r.subsystem, cA, (hA, hB))] = r

pdf_fp = OUT/"fig3_developmental_axis.pdf"
with PdfPages(pdf_fp) as pdf:
    pages = 0
    ROWS_PER_PAGE = 3
    rows_this_page = []
    for sub in sig_subs:
        rows_this_page.append(sub)
        if len(rows_this_page) < ROWS_PER_PAGE and sub != sig_subs[-1]:
            continue
        n = len(rows_this_page)
        fig, axes = plt.subplots(n*2, 3, figsize=(13, 2.4*n*2),
                                  gridspec_kw={'height_ratios':[2.4, 1.0]*n},
                                  squeeze=False)
        for i, s in enumerate(rows_this_page):
            for j, cond in enumerate(CONDS):
                ax_cdf = axes[i*2,     j]
                ax_tr  = axes[i*2 + 1, j]
                pd_ = panel_data_within_condition(s, cond)
                if pd_ is None:
                    ax_cdf.text(0.5,0.5,"no reps",ha='center',va='center',transform=ax_cdf.transAxes,fontsize=7)
                    ax_tr.axis('off')
                    continue
                # (a) CDFs by hpf
                for h in HPFS:
                    v = pd_['log_ratio'][h]; v = v[~np.isnan(v)]
                    if len(v) == 0: continue
                    vs = np.sort(v)
                    ys = np.arange(1, len(vs)+1)/len(vs)
                    ax_cdf.plot(vs, ys, color=hpf_ramp(cond, h), lw=1.1, drawstyle='steps-post',
                                 label=f'{h}' if j==0 else None)
                ax_cdf.set_xlim(-6, 6); ax_cdf.set_ylim(0, 1)
                ax_cdf.axvline(0, color='#ccc', lw=0.3, ls=':')
                ax_cdf.tick_params(labelsize=6)
                if j == 0:
                    ax_cdf.set_ylabel('CDF', fontsize=7)
                    ax_cdf.legend(title='hpf', fontsize=5, title_fontsize=6, loc='lower right', frameon=False)
                ax_cdf.set_title(f"{cond}  (n_reps={pd_['n_reps']})", fontsize=8)

                # Bracket significant transitions
                y_bracket = 1.06
                for ia in range(len(HPFS)-1):
                    hA, hB = HPFS[ia], HPFS[ia+1]
                    key = (s, cond, (hA, hB))
                    if key not in tA_dict: continue
                    r = tA_dict[key]
                    # position bracket at the midpoint between the two hpf CDFs' x-centroids
                    xa = np.nanmedian(pd_['log_ratio'][hA])
                    xb = np.nanmedian(pd_['log_ratio'][hB])
                    if np.isnan(xa) or np.isnan(xb): continue
                    xc = (xa + xb) / 2
                    ax_cdf.annotate(f"W1={r.W1_obs:.2f}",
                                     xy=(xc, y_bracket), xycoords=('data','axes fraction'),
                                     ha='center', va='bottom', fontsize=5.5, color='red',
                                     bbox=dict(facecolor='white', alpha=0.75, edgecolor='red', pad=1))
                    # draw a bracket
                    hatch = '/' if str(r.get('floor_exposed','')).startswith('yes') else None
                    ax_cdf.axvspan(min(xa,xb), max(xa,xb), ymin=1.01, ymax=1.05,
                                    color='red', alpha=0.3, hatch=hatch, clip_on=False)

                # (b) centroid trajectory
                mids = [np.nanmedian(pd_['median_logv'][h]) for h in HPFS]
                q25  = [np.nanquantile(pd_['median_logv'][h], 0.25) for h in HPFS]
                q75  = [np.nanquantile(pd_['median_logv'][h], 0.75) for h in HPFS]
                ax_tr.fill_between(HPFS, q25, q75, color=BASE_HUE[cond], alpha=0.2, linewidth=0)
                ax_tr.plot(HPFS, mids, 'o-', color=BASE_HUE[cond], lw=1.2, ms=4)
                ax_tr.axhline(0, color='k', lw=0.3)
                ax_tr.set_xticks(HPFS); ax_tr.tick_params(labelsize=6)
                ax_tr.set_xlabel('hpf', fontsize=7)
                if j == 0:
                    ax_tr.set_ylabel('median log₂(|v|/μ_r)\n± IQR', fontsize=6)
                ax_tr.set_xlim(20, 125)
                ax_tr.set_ylim(-5, 5)
            # row label on left
            axes[i*2, 0].annotate(s[:50],
                                   xy=(-0.3, 0.5), xycoords=axes[i*2,0].transAxes,
                                   ha='right', va='center', fontsize=8, color='#333',
                                   weight='bold')
        pages += 1
        plt.tight_layout(rect=[0.04, 0.02, 0.99, 0.97])
        fig.suptitle(f"Fig 3 — Developmental axis (page {pages})", fontsize=10, y=0.995)
        pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)
        print(f"  page {pages}  ({n} subsystems)", flush=True)
        rows_this_page = []
print(f"[write] {pdf_fp}", flush=True)
