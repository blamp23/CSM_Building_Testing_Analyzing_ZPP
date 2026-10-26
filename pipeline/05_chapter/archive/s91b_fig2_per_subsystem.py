"""
s91 -- Fig 2: condition-axis, one row of 3 panels per Tier A contrast.
  (a) CDFs of all 3 conditions at that hpf (F1) or at hpf_A/hpf_B (F2).
      Shade area between the two contrasted curves (= W1_obs). W1, snr, q
      printed in corner.
  (b) Shift function Δq(p) = q_A(p) − q_B(p) for p ∈ [0,1], with grey band
      for the 5–95th percentile of Δq(p) under random null draws.
  (c) Per-rep dot plot: one dot per rep per condition at log₂(median|v|/μ_r)
      sorted by effect size, lines connect same rep across conditions,
      thick line = subsystem median; small grey bar below shows per-rep
      zero fraction.

Palette: BL = blue (#1f77b4), D = black (#000000), LD = orange (#ff7f0e).
Zeros DROPPED; μ_r = mean of non-zero values pooled over the 3 cells shown
in the panel. Reps only; loop-flagged excluded.

Output: results/findings_v7_mean_dropzero/plots/fig2_condition_axis.pdf
"""
from pathlib import Path
import numpy as np, pandas as pd, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parents[1]
POOL = ROOT/"results/pooled_v7"
SUBS = ROOT/"results/subsets_v7"
MOR  = ROOT/"results/findings_v7_mean_dropzero"
OUT  = MOR/"plots"
OUT.mkdir(exist_ok=True, parents=True)

CONDS = ["BL", "D", "LD"]
COLORS = {"BL": "#1f77b4", "D": "#000000", "LD": "#f4c430"}
HPFS = [24, 48, 72, 96, 120]
ZERO_TOL = 1e-10
LOOP_MAX = 100.0
NDRAW = 500   # for null-shift band (smaller than 2000 since computed per panel)

tA_only = pd.read_csv(MOR/"tierA_final.csv").assign(_tier='A')
tB_only = pd.read_csv(MOR/"tierB_final.csv").assign(_tier='B')
tAB_keys = set(zip(
    pd.concat([tA_only, tB_only]).contrast,
    pd.concat([tA_only, tB_only]).subsystem,
))

f1 = pd.read_csv(MOR/"findings_family1.csv").assign(family='F1')
f2 = pd.read_csv(MOR/"findings_family2.csv").assign(family='F2')
tA = pd.concat([f1, f2], ignore_index=True)
tA = tA[tA.N_tier2.notna() & (tA.N_tier2 >= 3)].copy()

# GAG-only scope: keratan / chondroitin / dermatan / heparan / glycosaminoglycan
# / sulfate (any GAG kind) / O-glycan / N-glycan / Glycosphingolipid.
# Case-insensitive match.
import re as _re
GAG_PAT = _re.compile(
    r"(keratan|chondroitin|dermatan|heparan|heparin|glycosaminoglycan|"
    r"sulfate|o-glycan|o_glycan|n-glycan|n_glycan|glycosphingolipid|"
    r"glycan metabolism)",
    _re.IGNORECASE,
)
tA = tA[tA.subsystem.astype(str).apply(lambda s: bool(GAG_PAT.search(s)))].copy()
print(f"[GAG scope] {len(tA)} (contrast × subsystem) rows across "
      f"{tA.subsystem.nunique()} GAG-related subsystems: "
      f"{sorted(tA.subsystem.unique())}", flush=True)

def _tier_of(row):
    if (row.contrast, row.subsystem) in tAB_keys:
        rA = tA_only[(tA_only.contrast==row.contrast) & (tA_only.subsystem==row.subsystem)]
        rB = tB_only[(tB_only.contrast==row.contrast) & (tB_only.subsystem==row.subsystem)]
        if len(rA) > 0: return 'A'
        if len(rB) > 0: return 'B'
    return ''
tA['tier_label'] = tA.apply(_tier_of, axis=1)
tA = tA.sort_values(['family','subsystem','contrast'])

# ----- caches
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

def panel_data(subsystem, cells_dict):
    """cells_dict = {display_label: cell_name}. Return log-ratios per label."""
    rep_sets = []
    for lab, cell in cells_dict.items():
        m = get_mem(cell)
        rs = m.index[m.subsystem == subsystem].intersection(get_pooled(cell).index)
        rep_sets.append(set(rs))
    common = sorted(set.intersection(*rep_sets)) if rep_sets else []
    if len(common) < 3: return None
    Xs = {lab: np.abs(get_pooled(cell).loc[common].values) for lab, cell in cells_dict.items()}
    # loop filter per rep
    max_abs = np.stack([Xs[lab].max(1) for lab in cells_dict]).max(0)
    keep = max_abs <= LOOP_MAX
    if keep.sum() < 3: return None
    for lab in cells_dict: Xs[lab] = Xs[lab][keep]
    reps_kept = [common[i] for i in np.where(keep)[0]]
    # μ_r = mean of nonzero values pooled over cells in the panel
    all_pooled = np.concatenate([Xs[lab] for lab in cells_dict], axis=1)
    masked = np.where(all_pooled >= ZERO_TOL, all_pooled, np.nan)
    with np.errstate(all="ignore"):
        mu_rep = np.nanmean(masked, axis=1)
    mu_rep = np.where(np.isfinite(mu_rep) & (mu_rep > 0), mu_rep, 1.0)
    out = dict(reps=reps_kept, mu_rep=mu_rep)
    for lab in cells_dict:
        X = Xs[lab]
        mask = X >= ZERO_TOL
        rep_idx, samp_idx = np.where(mask)
        v = X[rep_idx, samp_idx]
        mu_broad = mu_rep[rep_idx]
        lr = np.log2(v / mu_broad)
        zf = 1.0 - mask.mean(axis=1)
        median_logv = np.array([np.log2(np.median(X[r][mask[r]])/mu_rep[r])
                                 if mask[r].any() else np.nan
                                 for r in range(X.shape[0])])
        out[lab] = dict(log_ratio=lr, zero_frac_per_rep=zf, median_logv=median_logv)
    return out


def hpf_ramp(cond, hpf):
    """Lighter at 24, darker at 120 for the same base condition color."""
    base = np.array(matplotlib.colors.to_rgb(COLORS[cond]))
    t = (hpf - 24) / (120 - 24)
    light = np.array([0.85, 0.85, 0.85])
    return tuple((1 - 0.25 - 0.75*t) * light + (0.25 + 0.75*t) * base)
import matplotlib

def panel_W1_shift(subsystem, cell_a, cell_b, mu_rep_source=None):
    """Return (q_A, q_B, delta_q, null_band_5, null_band_95) over 201 quantiles."""
    rep_sets = [set(get_mem(c).index[get_mem(c).subsystem == subsystem]
                                     .intersection(get_pooled(c).index))
                for c in (cell_a, cell_b)]
    common = sorted(set.intersection(*rep_sets))
    if len(common) < 3: return None
    AA = np.abs(get_pooled(cell_a).loc[common].values)
    AB = np.abs(get_pooled(cell_b).loc[common].values)
    max_abs = np.maximum(AA.max(1), AB.max(1))
    keep = max_abs <= LOOP_MAX
    AA, AB = AA[keep], AB[keep]
    AB_concat = np.concatenate([AA, AB], axis=1)
    masked = np.where(AB_concat >= ZERO_TOL, AB_concat, np.nan)
    with np.errstate(all="ignore"):
        mu = np.nanmean(masked, axis=1)
    mu = np.where(np.isfinite(mu) & (mu > 0), mu, 1.0)
    log_AA = np.where(AA >= ZERO_TOL, np.log2(np.maximum(AA, ZERO_TOL)/mu[:,None]), np.nan)
    log_AB = np.where(AB >= ZERO_TOL, np.log2(np.maximum(AB, ZERO_TOL)/mu[:,None]), np.nan)
    q = np.linspace(0, 1, 201)
    qa = np.nanquantile(log_AA.ravel(), q)
    qb = np.nanquantile(log_AB.ravel(), q)
    dq = qa - qb

    # null band: random subsets of same size from the subsystem-rep background
    all_mem = get_mem(cell_a).subsystem.reindex(common[:len(AA)]).fillna("UNMAPPED")
    # background: non-exchange, non-loop reps in both pooled sets with carriage both
    carry_both = (np.median(AA, axis=1) > 1e-4) & (np.median(AB, axis=1) > 1e-4)
    not_loop = max_abs[keep] <= LOOP_MAX
    bg_mask = carry_both & not_loop
    bg_idx = np.where(bg_mask)[0]
    if len(bg_idx) <= len(AA):
        return dict(q=q, qa=qa, qb=qb, dq=dq, band_lo=np.full_like(q, np.nan), band_hi=np.full_like(q, np.nan))
    rng = np.random.default_rng(0)
    dq_null = np.zeros((NDRAW, 201))
    for i in range(NDRAW):
        pick = rng.choice(bg_idx, len(AA), replace=False)
        la = log_AA[pick].ravel()
        lb = log_AB[pick].ravel()
        dq_null[i] = np.nanquantile(la, q) - np.nanquantile(lb, q)
    band_lo = np.nanpercentile(dq_null, 5, axis=0)
    band_hi = np.nanpercentile(dq_null, 95, axis=0)
    return dict(q=q, qa=qa, qb=qb, dq=dq, band_lo=band_lo, band_hi=band_hi)

def parse_contrast(c):
    m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_(\d+)", c)
    return m.group(1), int(m.group(2)), m.group(3), int(m.group(4))

# Pages: one Tier A row per page, full-size (density + CDF + shift + per-rep dots).
# F1 first (sorted by hpf then subsystem), then F2 (sorted by condition, hpf_A, subsystem).

def parse_cols(df):
    rows = []
    for _, r in df.iterrows():
        m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_(\d+)", r.contrast)
        rr = r.copy()
        rr['cA'], rr['hA'], rr['cB'], rr['hB'] = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
        rows.append(rr)
    return pd.DataFrame(rows).reset_index(drop=True)

tAp = parse_cols(tA)
tAp = pd.concat([
    tAp[tAp.family=='F1'].sort_values(['hA','subsystem','cA','cB']),
    tAp[tAp.family=='F2'].sort_values(['cA','hA','subsystem']),
], ignore_index=True)
print(f"[plan] {len(tAp)} Tier A rows, one per page", flush=True)

pdf_fp = OUT/"fig2_per_subsystem_GAG.pdf"
with PdfPages(pdf_fp) as pdf:
    for page_num, (_, r) in enumerate(tAp.iterrows(), 1):
        n = 1
        rows_this_page = [r]
        # Single Tier A row per page, laid out generously (letter-ish landscape).
        fig = plt.figure(figsize=(16, 7.5))
        gs = fig.add_gridspec(2, 3, height_ratios=[1, 4], hspace=0.08, wspace=0.28)
        axes = np.empty((n, 3), dtype=object)
        axes_dens = np.empty(n, dtype=object)
        axes_dens[0] = fig.add_subplot(gs[0, 0])
        axes[0, 0]   = fig.add_subplot(gs[1, 0], sharex=axes_dens[0])
        axes[0, 1]   = fig.add_subplot(gs[1, 1])
        axes[0, 2]   = fig.add_subplot(gs[1, 2])
        for i, r in enumerate(rows_this_page):
            cA, hA, cB, hB = parse_contrast(r.contrast)
            # Build cells_dict and per-label colors.
            # F1 (between-condition): 3 condition CDFs at the single hpf; contrast shade = (cA, cB).
            # F2 (within-condition, adjacent hpf): 2 CDFs for the same condition at the two hpfs,
            #   using hpf-ramp of the condition's base color; shade between them.
            if r.family == 'F1':
                hpf = hA
                cells_dict = {c: f"{c}_{hpf}" for c in CONDS}
                label_colors = {c: COLORS[c] for c in CONDS}
                legend_labels = {c: c for c in CONDS}
                contrast_pair = (cA, cB)
            else:
                lab_A = f"{cA} {hA} hpf"
                lab_B = f"{cB} {hB} hpf"
                cells_dict = {lab_A: f"{cA}_{hA}", lab_B: f"{cB}_{hB}"}
                label_colors = {lab_A: hpf_ramp(cA, hA), lab_B: hpf_ramp(cB, hB)}
                legend_labels = {lab_A: lab_A, lab_B: lab_B}
                contrast_pair = (lab_A, lab_B)
            pd_ = panel_data(r.subsystem, cells_dict)
            # Density strip (top, small)
            ax_d = axes_dens[i]
            # Panel (a): CDF (below density)
            ax = axes[i,0]
            if pd_ is None:
                ax_d.axis('off')
                ax.text(0.5,0.5,"no eligible reps",ha='center',va='center',transform=ax.transAxes,fontsize=8)
            else:
                # Density KDEs
                x_grid = np.linspace(-6, 6, 240)
                for lab in cells_dict:
                    v = pd_[lab]['log_ratio']; v = v[~np.isnan(v)]
                    if len(v) < 5 or np.ptp(v) < 1e-6: continue
                    if len(v) > 20000:
                        v = np.random.default_rng(0).choice(v, 20000, replace=False)
                    try:
                        kde = gaussian_kde(v, bw_method=0.3)
                        ax_d.plot(x_grid, kde(x_grid), color=label_colors[lab], lw=0.9, alpha=0.85)
                    except Exception:
                        pass
                ax_d.axvline(0, color='#ccc', lw=0.3, ls=':')
                ax_d.set_xlim(-6, 6)
                ax_d.set_yticks([]); ax_d.tick_params(labelbottom=False, length=0)
                for sp in ('top','right','left'): ax_d.spines[sp].set_visible(False)
                ax_d.spines['bottom'].set_linewidth(0.3)
                if i == 0:
                    ax_d.set_title(f"(a) density + CDF of log₂(|v|/μ_r)", fontsize=8)
                # CDFs
                for lab in cells_dict:
                    v = pd_[lab]['log_ratio']; v = v[~np.isnan(v)]
                    if len(v) == 0: continue
                    vs = np.sort(v)
                    ys = np.arange(1, len(vs)+1)/len(vs)
                    ax.plot(vs, ys, color=label_colors[lab], lw=1.2, drawstyle='steps-post', label=legend_labels[lab])
                # shade between the two contrasted curves (= W1_obs)
                lab_a, lab_b = contrast_pair
                if lab_a in pd_ and lab_b in pd_:
                    va = np.sort(pd_[lab_a]['log_ratio']); va = va[~np.isnan(va)]
                    vb = np.sort(pd_[lab_b]['log_ratio']); vb = vb[~np.isnan(vb)]
                    if len(va) and len(vb):
                        q = np.linspace(0,1,201)
                        qa = np.nanquantile(va, q); qb = np.nanquantile(vb, q)
                        ax.fill_betweenx(q, np.minimum(qa,qb), np.maximum(qa,qb),
                                          color='red', alpha=0.08, linewidth=0)
            ax.set_xlim(-6, 6); ax.set_ylim(0, 1)
            ax.axvline(0, color='#ccc', lw=0.3, ls=':')
            ax.set_xlabel('log₂(|v|/μ_r)', fontsize=7)
            ax.set_ylabel('CDF', fontsize=7)
            ax.legend(fontsize=6, loc='lower right', frameon=False)
            ax.tick_params(labelsize=6)
            ax.text(0.02, 0.98,
                    f"W1={r.W1_obs:.2f}  snr={r.snr:.1f}  q={r.q_bh:.3f}\n"
                    f"dir={'↑' if r.direction>0 else '↓'}  Tier {r.tier_label or r.tier}",
                    ha='left', va='top', fontsize=6, color='#333',
                    transform=ax.transAxes,
                    bbox=dict(facecolor='white', alpha=0.75, edgecolor='none', pad=1.5))

            # Panel (b): shift function
            ax = axes[i,1]
            sh = panel_W1_shift(r.subsystem, f"{cA}_{hA}", f"{cB}_{hB}")
            if sh is not None:
                ax.fill_between(sh['q'], sh['band_lo'], sh['band_hi'],
                                 color='#999', alpha=0.25, linewidth=0, label='null 5–95%')
                ax.plot(sh['q'], sh['dq'], color='red', lw=1.1, label=f'{cA}−{cB}')
                ax.axhline(0, color='k', lw=0.3)
            ax.set_xlabel('quantile p', fontsize=7); ax.set_ylabel('Δq(p) = q_A − q_B', fontsize=7)
            ax.set_ylim(-6, 6); ax.tick_params(labelsize=6)
            ax.legend(fontsize=6, loc='upper left', frameon=False)
            if i == 0:
                ax.set_title("(b) Shift function (null band = 5–95%)", fontsize=8)

            # Panel (c): per-rep dot plot
            ax = axes[i,2]
            if pd_ is not None:
                labs = list(cells_dict.keys())
                sort_lab = contrast_pair[0]
                effect = pd_[sort_lab]['median_logv']
                order = np.argsort(np.nan_to_num(effect, nan=0))
                xs = np.arange(len(order))
                for lab in labs:
                    ys = pd_[lab]['median_logv'][order]
                    ax.plot(xs, ys, '.', color=label_colors[lab], ms=3, alpha=0.85,
                            label=legend_labels[lab])
                # thin lines between consecutive labs
                for a_i, b_i in zip(range(len(labs)-1), range(1, len(labs))):
                    for k, idx in enumerate(order):
                        ya = pd_[labs[a_i]]['median_logv'][idx]
                        yb = pd_[labs[b_i]]['median_logv'][idx]
                        if np.isnan(ya) or np.isnan(yb): continue
                        ax.plot([k, k], [ya, yb], color='#aaa', lw=0.2, alpha=0.4)
                for lab in labs:
                    m = np.nanmedian(pd_[lab]['median_logv'])
                    ax.axhline(m, color=label_colors[lab], lw=0.8, alpha=0.6)
                # zero-fraction bars under
                ymin = -8; bar_h = 0.6
                for ci, lab in enumerate(labs):
                    zf = pd_[lab]['zero_frac_per_rep'][order]
                    ys = ymin + ci*bar_h
                    for k, z in enumerate(zf):
                        ax.add_patch(plt.Rectangle((k-.4, ys), 0.8, bar_h*z,
                                                    facecolor=label_colors[lab], alpha=0.3,
                                                    edgecolor='none'))
                ax.set_ylim(-10, 8)
                ax.axhline(0, color='k', lw=0.3, ls=':')
                ax.set_xlim(-0.5, len(xs)-0.5)
                ax.set_xlabel(f'rep (sorted by {sort_lab} median)', fontsize=6)
                ax.set_ylabel('log₂(median|v|/μ_r)', fontsize=7)
                ax.tick_params(labelsize=6)
                ax.legend(fontsize=6, loc='upper right', frameon=False)
            if i == 0:
                ax.set_title("(c) Per-rep medians (zero fraction bars below)", fontsize=8)

            # row label on the left (anchored to the density strip so it sits mid-row)
            axes_dens[i].annotate(f"{r.family}: {r.contrast}\n{r.subsystem[:50]}",
                                   xy=(-0.35, -1.5), xycoords=axes_dens[i].transAxes,
                                   ha='right', va='center', fontsize=7, color='#333')
        row = rows_this_page[0]
        fig.text(0.5, 0.985,
                 f"Fig 2 — {row.family}: {row.contrast} — {row.subsystem}  "
                 f"(page {page_num}/{len(tAp)})",
                 ha='center', va='top', fontsize=12, weight='bold')
        fig.text(0.5, 0.961,
                 "(a) density + CDF of log₂(|v|/μ_r) per cell — "
                 "red-shaded area = W1 between the two contrasted curves    │    "
                 "(b) shift function Δq(p) = q_A(p) − q_B(p) across quantiles — "
                 "grey band = 5–95% of 2000 null draws    │    "
                 "(c) per-rep medians sorted by cell_A — "
                 "dots per cell, thin lines connect same rep across cells, "
                 "horizontal lines = per-cell medians, grey bars below = per-rep zero fractions",
                 ha='center', va='top', fontsize=8, wrap=True, color='#444')
        pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)
        if page_num % 10 == 0:
            print(f"  page {page_num}/{len(tAp)} done", flush=True)
print(f"[write] {pdf_fp}", flush=True)
