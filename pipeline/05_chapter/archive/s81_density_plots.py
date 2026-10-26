"""
s81 -- Density plots of per-subsystem log|flux| for v7.

Output 1 (density_all_subsystems.pdf):
    Multi-page PDF, 20 subsystems/page, 5 columns (hpf 24/48/72/96/120).
    Each subplot overlays BL/D/LD KDEs. Blank cell if the subsystem has
    no representatives in all three cells at that hpf.

Output 2 (density_sig_subsystems.pdf):
    Single page, 14 significant subsystems (q <= 0.10 anywhere in
    Family 1 or Family 2). Red box around any (subsystem, hpf) cell
    that was significant in a contrast involving that hpf.
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
EPS   = 1e-30
COLORS = {"BL": "#1f77b4", "D": "#d62728", "LD": "#2ca02c"}
X_LO, X_HI = -105, 15
X_GRID = np.linspace(X_LO, X_HI, 240)

print("[load] pooled samples + memberships for 15 cells...", flush=True)
pooled = {}; mem = {}
for cond in CONDS:
    for h in HPFS:
        cell = f"{cond}_{h}"
        pooled[cell] = pd.read_csv(POOL/f"{cell}_pooled.csv", index_col=0)
        m = pd.read_csv(SUBS/f"{cell}_subset_membership.csv")
        mem[cell] = m.drop_duplicates("representative").set_index("representative")
print("[load] done.", flush=True)

all_subs = set()
for cell, m in mem.items():
    all_subs.update(m.subsystem.dropna().unique())
all_subs = sorted(s for s in all_subs if s and s != "UNMAPPED")
print(f"[plan] {len(all_subs)} subsystems × 5 hpf.", flush=True)

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
sig_subs = sorted(set(sig.subsystem))
print(f"[plan] {len(sig_subs)} significant subsystems; "
      f"{len(sig_cells)} (sub, hpf) sig cells.", flush=True)

def get_logflux(subsystem, cell, max_n=20000):
    X = pooled[cell]; m = mem[cell]
    rxns = m.index[m.subsystem == subsystem]
    rxns = rxns.intersection(X.index)
    if len(rxns) == 0:
        return None
    A = np.abs(X.loc[rxns].values).ravel()
    v = np.log2(A + EPS)
    if len(v) > max_n:
        v = np.random.default_rng(0).choice(v, max_n, replace=False)
    return v

def draw_subsystem_row(axes_row, subsystem, box_sig=False):
    """Draw 5 subplots for one subsystem across hpfs on provided axes."""
    plotted_any = False
    for j, h in enumerate(HPFS):
        ax = axes_row[j]
        has_data = False
        for cond in CONDS:
            cell = f"{cond}_{h}"
            v = get_logflux(subsystem, cell)
            if v is None or len(v) < 5 or np.ptp(v) < 1e-6:
                continue
            try:
                kde = gaussian_kde(v, bw_method=0.25)
                ys = kde(X_GRID)
                ax.plot(X_GRID, ys, color=COLORS[cond], lw=1.0, alpha=0.9,
                        label=cond)
                has_data = True; plotted_any = True
            except Exception:
                pass
        if box_sig and (subsystem, h) in sig_cells:
            for spine in ax.spines.values():
                spine.set_edgecolor("red"); spine.set_linewidth(1.6)
        else:
            for spine in ax.spines.values():
                spine.set_edgecolor("#999"); spine.set_linewidth(0.4)
        if not has_data:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    fontsize=5, color="gray", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(X_LO, X_HI)
    return plotted_any

def plot_set(sub_list, outfp, title, box_sig=False, chunk=None):
    """Make 5-col grid of subsystem rows. If chunk given, paginate."""
    if chunk is None: chunk = len(sub_list)
    pages = [sub_list[i:i+chunk] for i in range(0, len(sub_list), chunk)]
    with PdfPages(outfp) as pdf:
        for pi, page in enumerate(pages, 1):
            n = len(page)
            fig, axes = plt.subplots(n, 5, figsize=(10, max(2, n*0.8)),
                                      squeeze=False)
            for i, sub in enumerate(page):
                draw_subsystem_row(axes[i], sub, box_sig=box_sig)
                label = sub if len(sub) <= 36 else sub[:33] + "…"
                axes[i,0].set_ylabel(label, fontsize=6, rotation=0,
                                      ha="right", va="center",
                                      labelpad=4)
            for j, h in enumerate(HPFS):
                axes[0, j].set_title(f"{h} hpf", fontsize=9)
            # legend only on page 1
            if pi == 1:
                handles = [plt.Line2D([0],[0], color=COLORS[c], lw=1.4, label=c)
                           for c in CONDS]
                fig.legend(handles=handles, ncol=3, loc="upper right",
                           bbox_to_anchor=(0.995, 0.995),
                           fontsize=7, frameon=False)
            fig.suptitle(f"{title}  (page {pi}/{len(pages)})"
                          if len(pages) > 1 else title,
                          fontsize=10, y=0.995)
            fig.supxlabel("log₂(|flux| + 1e-30)", fontsize=8)
            plt.tight_layout(rect=[0.22, 0.02, 0.99, 0.97])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print(f"  page {pi}/{len(pages)}  ({len(page)} subsystems)", flush=True)
    print(f"[write] {outfp}", flush=True)

print("\n[plot 2] 14 significant subsystems, boxes on sig cells...", flush=True)
plot_set(sig_subs, OUT/"density_sig_subsystems.pdf",
         "Significant subsystems (q ≤ 0.10 anywhere) — red = sig at that hpf",
         box_sig=True)

print(f"\n[plot 1] all {len(all_subs)} subsystems, 20 per page...", flush=True)
plot_set(all_subs, OUT/"density_all_subsystems.pdf",
         "All subsystems — log|flux| densities",
         box_sig=False, chunk=20)

print("\nDONE")
