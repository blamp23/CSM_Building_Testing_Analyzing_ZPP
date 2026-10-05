"""
s81b -- Density plots v2 for v7, Option 1 spec.

What the two-tier W1 test actually consumes:
    for each representative r (reaction) in the subsystem,
        mu_r = mean(|v_r|) pooled across the conditions plotted on this panel
    x = log2(|v_r,sample| / mu_r)   dropping (|v| < 1e-10)

That's the test's statistic. Not raw log|flux|.

Filters (match §1 + §11.9 of FAMILY_PREREG):
  - enzyme-subset representatives only (one row per `representative` from
    subset_v7)
  - drop loop-flagged reps (max |v| > 100 anywhere in the 3 conditions)
  - drop UNMAPPED, Exchange, Artificial, Transport, Isolated, Pool reactions
  - skip panel if < 3 reps survive

Each panel:
  - overlay KDEs for BL / D / LD at that hpf, same x axis
  - top-right text: `BL zf=0.XX  D zf=0.YY  LD zf=0.ZZ` (zero fractions
    before log — the mass we dropped)
  - red border if any Family 1 or Family 2 contrast involving that hpf is
    at q <= 0.10 for that subsystem

Outputs:
  results/findings_v7/plots/density_sig_subsystems_median.pdf     (14 sig subs)
  results/findings_v7/plots/density_all_subsystems_median.pdf     (paginated)
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
OUT  = ROOT/"results/findings_v7/plots"
OUT.mkdir(parents=True, exist_ok=True)

CONDS = ["BL", "D", "LD"]
HPFS  = [24, 48, 72, 96, 120]
COLORS = {"BL": "#1f77b4", "D": "#d62728", "LD": "#2ca02c"}
ZERO_TOL = 1e-10
LOOP_MAX = 100.0
EXCLUDE = {"Exchange reactions","Artificial reactions","Transport reactions",
           "Exchange","Artificial","Transport","Isolated","Pool reactions",
           "UNMAPPED",""}
X_LO, X_HI = -6, 6     # log2(|v|/mu) typically lives here; test uses finer scale
X_GRID = np.linspace(X_LO, X_HI, 240)

print("[load] pooled samples + memberships...", flush=True)
pooled = {}; mem = {}
for cond in CONDS:
    for h in HPFS:
        cell = f"{cond}_{h}"
        pooled[cell] = pd.read_csv(POOL/f"{cell}_pooled.csv", index_col=0)
        m = pd.read_csv(SUBS/f"{cell}_subset_membership.csv")
        mem[cell] = m.drop_duplicates("representative").set_index("representative")
print("[load] done.", flush=True)

# Family-eligible subsystems (union across cells, minus EXCLUDE)
all_subs = set()
for m in mem.values():
    all_subs.update(m.subsystem.dropna().unique())
all_subs = sorted(s for s in all_subs if s not in EXCLUDE)
print(f"[plan] {len(all_subs)} family-eligible subsystems", flush=True)

# Load findings for red borders
f1 = pd.read_csv(ROOT/"results/findings_v7/findings_family1.csv")
f2 = pd.read_csv(ROOT/"results/findings_v7/findings_family2.csv")
sig = pd.concat([f1[f1.q_bh <= 0.10], f2[f2.q_bh <= 0.10]], ignore_index=True)

def parse_contrast(c):
    m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_(\d+)", c)
    return m.group(1), int(m.group(2)), m.group(3), int(m.group(4))

sig_cells = set()
for _, r in sig.iterrows():
    cA, hA, cB, hB = parse_contrast(r.contrast)
    sig_cells.add((r.subsystem, hA))
    sig_cells.add((r.subsystem, hB))
sig_subs = sorted(set(sig.subsystem) - EXCLUDE)
print(f"[plan] sig subsystems after §11.9 Pool exclusion: {len(sig_subs)}", flush=True)

def compute_panel(subsystem, hpf):
    """Return dict with per-condition log-ratio arrays + zero fractions.
    Returns None if < 3 reps survive the filters."""
    # Intersect reactions in the subsystem across the 3 cells + find common
    # reps that aren't loop-flagged and have at least one non-zero flux
    cells = {c: f"{c}_{hpf}" for c in CONDS}
    rep_sets = []
    for c, cell in cells.items():
        m = mem[cell]
        rs = m.index[m.subsystem == subsystem]
        rs = rs.intersection(pooled[cell].index)
        rep_sets.append(set(rs))
    common = sorted(set.intersection(*rep_sets)) if rep_sets else []
    if len(common) < 3:
        return None

    # Load flux matrices for just these reps (n_reps x 2000)
    Xs = {c: np.abs(pooled[cells[c]].loc[common].values) for c in CONDS}

    # Loop-flag per rep: max |v| across all 3 cells
    max_abs = np.stack([Xs[c].max(axis=1) for c in CONDS], axis=0).max(axis=0)
    keep = max_abs <= LOOP_MAX
    if keep.sum() < 3:
        return None

    for c in CONDS:
        Xs[c] = Xs[c][keep]
    reps_kept = [common[i] for i in np.where(keep)[0]]

    # §11.12: median of non-zero values pooled across the 3 conditions
    all_pooled = np.concatenate([Xs[c] for c in CONDS], axis=1)   # n_rep x 6000
    masked = np.where(all_pooled >= ZERO_TOL, all_pooled, np.nan)
    with np.errstate(all="ignore"):
        mu_rep = np.nanmedian(masked, axis=1)                     # n_rep
    # A rep with mu_rep=0 would make log undefined; drop those too
    use_rep = mu_rep > 0
    if use_rep.sum() < 3:
        return None
    mu_rep = mu_rep[use_rep]
    for c in CONDS:
        Xs[c] = Xs[c][use_rep]

    out = {"reps_used": int(use_rep.sum())}
    for c in CONDS:
        X = Xs[c]                                                 # n_rep x 2000
        zero_mask = X < ZERO_TOL
        zero_frac = float(zero_mask.sum()) / X.size
        nz = X[~zero_mask]
        mu_broadcast = mu_rep[:, None]
        # For log we need to re-align: use the per-rep mu for each entry
        rep_idx, sample_idx = np.where(~zero_mask)
        v_nz = X[rep_idx, sample_idx]
        mu_nz = mu_rep[rep_idx]
        log_ratio = np.log2(v_nz / mu_nz)
        # subsample for KDE if huge
        if len(log_ratio) > 20000:
            log_ratio = np.random.default_rng(0).choice(log_ratio, 20000, replace=False)
        out[c] = {"log_ratio": log_ratio, "zero_frac": zero_frac}
    return out

def draw_row(axes_row, subsystem, box_sig):
    for j, h in enumerate(HPFS):
        ax = axes_row[j]
        panel = compute_panel(subsystem, h)
        if panel is None:
            ax.text(0.5, 0.5, "no eligible reps", ha="center", va="center",
                    fontsize=5, color="gray", transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_edgecolor("#ccc"); spine.set_linewidth(0.3)
            continue
        has_any = False
        for c in CONDS:
            v = panel[c]["log_ratio"]
            if len(v) < 5 or np.ptp(v) < 1e-6:
                continue
            try:
                kde = gaussian_kde(v, bw_method=0.3)
                ys = kde(X_GRID)
                ax.plot(X_GRID, ys, color=COLORS[c], lw=1.0, alpha=0.9)
                has_any = True
            except Exception:
                pass
        # zero fractions text
        zf_text = "  ".join(
            f"{c} zf={panel[c]['zero_frac']:.2f}" for c in CONDS
        )
        ax.text(0.02, 0.96, zf_text, fontsize=4.5, color="#444",
                transform=ax.transAxes, ha="left", va="top")
        if box_sig and (subsystem, h) in sig_cells:
            for spine in ax.spines.values():
                spine.set_edgecolor("red"); spine.set_linewidth(1.6)
        else:
            for spine in ax.spines.values():
                spine.set_edgecolor("#999"); spine.set_linewidth(0.4)
        if not has_any:
            ax.text(0.5, 0.4, "no density", ha="center", va="center",
                    fontsize=5, color="gray", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(X_LO, X_HI)
        ax.axvline(0, color="#bbb", lw=0.4, ls=":")

def plot_set(sub_list, outfp, title, box_sig=False, chunk=None):
    if chunk is None: chunk = len(sub_list)
    pages = [sub_list[i:i+chunk] for i in range(0, len(sub_list), chunk)]
    with PdfPages(outfp) as pdf:
        for pi, page in enumerate(pages, 1):
            n = len(page)
            fig, axes = plt.subplots(n, 5, figsize=(11, max(2, n*0.85)),
                                      squeeze=False)
            for i, sub in enumerate(page):
                draw_row(axes[i], sub, box_sig=box_sig)
                label = sub if len(sub) <= 36 else sub[:33] + "…"
                axes[i,0].set_ylabel(label, fontsize=6, rotation=0,
                                      ha="right", va="center", labelpad=4)
            for j, h in enumerate(HPFS):
                axes[0, j].set_title(f"{h} hpf", fontsize=9)
            if pi == 1:
                handles = [plt.Line2D([0],[0], color=COLORS[c], lw=1.4, label=c)
                           for c in CONDS]
                fig.legend(handles=handles, ncol=3, loc="upper right",
                           bbox_to_anchor=(0.995, 0.995),
                           fontsize=7, frameon=False)
            fig.suptitle(f"{title}  (page {pi}/{len(pages)})"
                          if len(pages) > 1 else title,
                          fontsize=10, y=0.995)
            fig.supxlabel("log₂(|v| / per-rep pooled median)   —   zeros dropped",
                           fontsize=8)
            plt.tight_layout(rect=[0.24, 0.03, 0.99, 0.97])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print(f"  page {pi}/{len(pages)}  ({len(page)} subsystems)", flush=True)
    print(f"[write] {outfp}", flush=True)

print("\n[plot 2] significant subsystems, boxes on sig cells...", flush=True)
plot_set(sig_subs, OUT/"density_sig_subsystems_median.pdf",
         "Significant subsystems (q ≤ 0.10) — W1 test statistic",
         box_sig=True)

print(f"\n[plot 1] all {len(all_subs)} family-eligible subsystems...", flush=True)
plot_set(all_subs, OUT/"density_all_subsystems_median.pdf",
         "All family-eligible subsystems — W1 test statistic",
         box_sig=False, chunk=18)

print("\nDONE")
