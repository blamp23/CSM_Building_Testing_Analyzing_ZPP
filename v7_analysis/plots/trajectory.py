"""
s81e -- Per-reaction trajectory plots for significant subsystems.

For each subsystem that has any q <= 0.10 hit in Family 1 or Family 2:
  - One row per subsystem, 5 panels (hpf 24/48/72/96/120).
  - Each panel: y = log2(median |v| per reaction), x = 3 cells (BL / D / LD).
  - Each reaction = one thin line connecting its 3 points.
  - Thick overlay = subsystem-level mean of per-reaction log2(median |v|).
  - Red border = any contrast involving that hpf for that subsystem was
    at q <= 0.10 under the MEAN primary pipeline (findings_v7/).

This plot shows directly which reactions move and in which direction.
If the subsystem is "really" shifting, lines slope the same way and the
thick mean line slopes. If reactions scatter in different directions,
mean stays flat but line fan-out tells the story.
"""
from pathlib import Path
import numpy as np, pandas as pd, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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

print("[load] pooled samples + memberships...", flush=True)
pooled = {}; mem = {}
for cond in CONDS:
    for h in HPFS:
        cell = f"{cond}_{h}"
        pooled[cell] = pd.read_csv(POOL/f"{cell}_pooled.csv", index_col=0)
        m = pd.read_csv(SUBS/f"{cell}_subset_membership.csv")
        mem[cell] = m.drop_duplicates("representative").set_index("representative")
print("[load] done", flush=True)

f1 = pd.read_csv(ROOT/"results/findings_v7/findings_family1.csv")
f2 = pd.read_csv(ROOT/"results/findings_v7/findings_family2.csv")
sig = pd.concat([f1[f1.q_bh <= 0.10], f2[f2.q_bh <= 0.10]], ignore_index=True)

def parse_contrast(c):
    m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_(\d+)", c)
    return m.group(1), int(m.group(2)), m.group(3), int(m.group(4))

sig_cells = set()
for _, r in sig.iterrows():
    cA, hA, cB, hB = parse_contrast(r.contrast)
    sig_cells.add((r.subsystem, hA)); sig_cells.add((r.subsystem, hB))
EXCLUDE = {"Exchange reactions","Artificial reactions","Transport reactions",
           "Exchange","Artificial","Transport","Isolated","Pool reactions"}
sig_subs = sorted(set(sig.subsystem) - EXCLUDE)
print(f"[plan] {len(sig_subs)} significant subsystems", flush=True)

def panel_data(subsystem, hpf):
    """Return per-reaction log2(median |v|) per cell + loop-filtered rep count,
    or None if < 3 reps common across cells."""
    rep_sets = []
    cells = {c: f"{c}_{hpf}" for c in CONDS}
    for c, cell in cells.items():
        rs = mem[cell].index[mem[cell].subsystem == subsystem]
        rs = rs.intersection(pooled[cell].index)
        rep_sets.append(set(rs))
    common = sorted(set.intersection(*rep_sets)) if rep_sets else []
    if len(common) < 3:
        return None
    Xs = {c: np.abs(pooled[cells[c]].loc[common].values) for c in CONDS}
    max_abs = np.stack([Xs[c].max(1) for c in CONDS]).max(0)
    keep = max_abs <= LOOP_MAX
    if keep.sum() < 3:
        return None
    for c in CONDS: Xs[c] = Xs[c][keep]
    # median absolute flux per rep per cell; drop reps where all 3 cells have zero median
    med = {c: np.median(Xs[c], axis=1) for c in CONDS}
    all_zero = np.all(np.stack([med[c] <= ZERO_TOL for c in CONDS]).all(axis=0))
    # Convert to log2 scale, floor zeros at ZERO_TOL for display
    log_med = {c: np.log2(np.maximum(med[c], ZERO_TOL)) for c in CONDS}
    return log_med

def draw_row(axes_row, subsystem):
    # Determine y-axis bounds consistent across hpf for this subsystem
    all_vals = []
    panels = {}
    for h in HPFS:
        d = panel_data(subsystem, h)
        panels[h] = d
        if d is None: continue
        for c in CONDS:
            all_vals.append(d[c])
    if not all_vals:
        for ax in axes_row: ax.axis("off")
        return
    vals = np.concatenate(all_vals)
    y_lo, y_hi = np.quantile(vals, 0.02), np.quantile(vals, 0.98)
    pad = max(1.0, 0.1 * (y_hi - y_lo))
    y_lo -= pad; y_hi += pad
    for j, h in enumerate(HPFS):
        ax = axes_row[j]
        d = panels[h]
        if d is None:
            ax.text(0.5, 0.5, "no reps", ha="center", va="center", fontsize=6,
                    color="gray", transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        n_rep = len(d["BL"])
        xs = np.arange(3)
        # thin per-reaction lines
        for i in range(n_rep):
            ys = [d[c][i] for c in CONDS]
            ax.plot(xs, ys, color="#555", lw=0.4, alpha=0.25)
            for k, c in enumerate(CONDS):
                ax.plot(xs[k], ys[k], ".", color=COLORS[c], ms=2.0, alpha=0.6)
        # thick condition means
        means = [float(np.mean(d[c])) for c in CONDS]
        ax.plot(xs, means, color="black", lw=1.6, zorder=5)
        for k, c in enumerate(CONDS):
            ax.plot(xs[k], means[k], "o", color=COLORS[c], ms=6, mec="black",
                    mew=0.6, zorder=6)
        ax.set_xticks(xs); ax.set_xticklabels(CONDS, fontsize=6)
        ax.set_ylim(y_lo, y_hi)
        ax.tick_params(axis="y", labelsize=5, length=2, pad=1)
        ax.annotate(f"n={n_rep}", (0.98, 0.02), xycoords="axes fraction",
                    ha="right", va="bottom", fontsize=5, color="#444")
        # red border on sig cells
        if (subsystem, h) in sig_cells:
            for spine in ax.spines.values():
                spine.set_edgecolor("red"); spine.set_linewidth(1.6)
        else:
            for spine in ax.spines.values():
                spine.set_edgecolor("#999"); spine.set_linewidth(0.4)

print("[plot] drawing trajectory PDF (one row per sig subsystem)...", flush=True)
out_fp = OUT/"trajectory_sig_subsystems.pdf"
with PdfPages(out_fp) as pdf:
    CHUNK = 8   # subsystems per page; trajectory rows need more vertical space
    pages = [sig_subs[i:i+CHUNK] for i in range(0, len(sig_subs), CHUNK)]
    for pi, page in enumerate(pages, 1):
        n = len(page)
        fig, axes = plt.subplots(n, 5, figsize=(11, max(2, n*1.4)),
                                  squeeze=False)
        for i, sub in enumerate(page):
            draw_row(axes[i], sub)
            label = sub if len(sub) <= 36 else sub[:33] + "…"
            axes[i,0].set_ylabel(label, fontsize=6, rotation=0,
                                  ha="right", va="center", labelpad=4)
        for j, h in enumerate(HPFS):
            axes[0, j].set_title(f"{h} hpf", fontsize=9)
        if pi == 1:
            handles = [plt.Line2D([0],[0], color=COLORS[c], marker="o",
                                   lw=0, ms=6, label=c) for c in CONDS]
            handles.append(plt.Line2D([0],[0], color="black", lw=1.6,
                                       label="subsystem mean"))
            fig.legend(handles=handles, ncol=4, loc="upper right",
                       bbox_to_anchor=(0.995, 0.995),
                       fontsize=7, frameon=False)
        fig.suptitle(f"Per-reaction trajectory  (page {pi}/{len(pages)})"
                      if len(pages) > 1 else "Per-reaction trajectory",
                      fontsize=10, y=0.995)
        fig.supylabel("log₂( median |v| per reaction )", fontsize=8, x=0.02)
        plt.tight_layout(rect=[0.24, 0.03, 0.99, 0.97])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        print(f"  page {pi}/{len(pages)}  ({len(page)} subsystems)", flush=True)
print(f"[write] {out_fp}")
print("DONE")
