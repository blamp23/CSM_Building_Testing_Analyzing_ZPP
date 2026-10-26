"""
s276 -- Two sets of per-subsystem density plots.

Folder 1 (developmental/):
    One PDF per subsystem.
    Three panels per PDF: BL, D, LD.
    Each panel overlays 5 density curves (one per hpf, viridis gradient).
    Uses raw log flux (preserves developmental magnitude scaling).

Folder 2 (light_regime/):
    One PDF per subsystem that has a cross-regime finding.
    Five panels per PDF (24, 48, 72, 96, 120 hpf).
    Each panel overlays 3 density curves (BL blue, D grey, LD yellow).
    Uses per-cell-median-CENTERED log flux (matches the centered analysis).
    Sig box on the panel(s) where the cross-regime finding lives.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
import matplotlib.cm as cm

plt.style.use("ggplot")
plt.rcParams.update({
    "axes.facecolor": "#F0F0F0",
    "axes.edgecolor": "#999999",
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": "white",
    "grid.linewidth": 1.0,
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "axes.labelcolor": "#333333",
    "font.family": "sans-serif",
    "font.size": 10,
})

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import POOLED as POOL, SUBSETS as SUBS, FIGURES
OUT_DEV = FIGURES / "subsystem_density" / "developmental"
OUT_LR  = FIGURES / "subsystem_density" / "light_regime"
OUT_DEV.mkdir(parents=True, exist_ok=True)
OUT_LR.mkdir(parents=True, exist_ok=True)

CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
ZERO_TOL = 1e-10

LR_COLORS = {"BL": "#1976D2", "D": "#555555", "LD": "#F6B52C"}
LR_ALPHA = 0.5
# Viridis gradient for hpf (24 = light, 120 = dark)
HPF_COLORS = {h: cm.viridis(i / (len(HPFS) - 1))
              for i, h in enumerate(HPFS)}
HPF_ALPHA = 0.4

# Subsystems with a developmental story (plot in developmental/)
DEV_SUBSYSTEMS = [
    "Fatty acid activation (cytosolic)",
    "Fatty acid oxidation",
    "Arachidonic acid metabolism",
    "Leukotriene metabolism",
    "Eicosanoid metabolism",
    "Cholesterol biosynthesis 1 (Bloch pathway)",
    "Beta oxidation of di-unsaturated fatty acids (n-6) (mitochondrial)",
    "Fatty acid elongation (even-chain)",
    "Fatty acid biosynthesis (even-chain)",
    "Valine, leucine, and isoleucine metabolism",
    "Lysine metabolism",
    "Serotonin and melatonin biosynthesis",
    "Bile acid biosynthesis",
    "Steroid metabolism",
    "Arginine and proline metabolism",
    "Nucleotide metabolism",
    "Pyrimidine metabolism",
    "Phenylalanine, tyrosine and tryptophan biosynthesis",
    "Tyrosine metabolism",
    "Fatty acid desaturation (even-chain)",
    "Omega-3 fatty acid metabolism",
]

# Subsystems with cross-regime (F1) findings — sig_hpf = where the finding lives
LR_SIG = {
    "Fatty acid activation (cytosolic)": [48],              # BL high
    "Serotonin and melatonin biosynthesis": [72],           # BL low
    "Arachidonic acid metabolism": [96],                    # BL reversal (centered)
    "Leukotriene metabolism": [96],                         # BL reversal (centered)
    "Fatty acid desaturation (even-chain)": [96],           # BL low (centered)
    "Omega-3 fatty acid metabolism": [96],                  # BL low (centered)
    "Bile acid biosynthesis": [120],                        # LD high
    "Steroid metabolism": [120],                            # LD high
    "Phenylalanine, tyrosine and tryptophan biosynthesis": [120],  # BL high
    "Tyrosine metabolism": [120],                           # BL high
    "Arginine and proline metabolism": [72],                # D > LD
    "Pyrimidine metabolism": [48],                          # mixed
    "Valine, leucine, and isoleucine metabolism": [120],    # mixed
}

print(f"[prep] loading pooled samples for 15 cells...", flush=True)
cell_data = {}
for cond in CONDS:
    for h in HPFS:
        cell = f"{cond}_{h}"
        samples = pd.read_csv(POOL / f"{cell}_pooled.csv", index_col=0)
        mem = pd.read_csv(SUBS / f"{cell}_subset_membership.csv")
        rep_to_sub = (mem.drop_duplicates("representative")
                        .set_index("representative")["subsystem"].to_dict())
        cell_data[cell] = dict(samples=samples, rep_to_sub=rep_to_sub)
print("[prep] done")

all_subs_seen = set()
for d in cell_data.values():
    all_subs_seen.update(d["rep_to_sub"].values())

# Compute per-cell median log flux (for centering in LR plots)
# Based on shared reactions across all 15 cells
cell_vectors = {}
for name, d in cell_data.items():
    vals = np.abs(d["samples"].values.astype(float))
    nz_mask = vals > ZERO_TOL
    masked = np.where(nz_mask, vals, np.nan)
    with np.errstate(all="ignore"):
        mean_nz = np.nanmean(masked, axis=1)
    cell_vectors[name] = pd.Series(mean_nz, index=d["samples"].index)
shared_rxns = set.intersection(*[set(v[v.notna()].index) for v in cell_vectors.values()])
shared_rxns = sorted(shared_rxns)
cell_medians = {}
for name in cell_data:
    v = cell_vectors[name].reindex(shared_rxns).values
    cell_medians[name] = float(np.nanmedian(np.log2(v + 1e-6)))
print(f"[prep] per-cell median log flux computed from {len(shared_rxns)} shared reactions")


def _collect_log_fluxes(cell, subsystem, center=False):
    """Collect log2(|flux|) samples for one subsystem in one cell."""
    d = cell_data[cell]
    reps = [r for r, s in d["rep_to_sub"].items() if s == subsystem]
    reps_in = [r for r in reps if r in d["samples"].index]
    if not reps_in:
        return np.array([])
    vals = np.abs(d["samples"].loc[reps_in].values.astype(float)).ravel()
    nz = vals[vals > ZERO_TOL]
    if len(nz) == 0:
        return np.array([])
    log_vals = np.log2(nz)
    if center:
        log_vals = log_vals - cell_medians[cell]
    return log_vals


def _slug(name):
    s = (name.replace(" ", "_").replace("/", "_")
             .replace("(", "").replace(")", "").replace(",", ""))
    while "__" in s:
        s = s.replace("__", "_")
    return s


def _dev_plot(subsystem):
    """3 panels (BL, D, LD); each panel shows 5 overlaid hpf density curves."""
    all_fluxes = {}
    for cond in CONDS:
        for h in HPFS:
            cell = f"{cond}_{h}"
            all_fluxes[cell] = _collect_log_fluxes(cell, subsystem, center=False)
    if not any(len(v) > 0 for v in all_fluxes.values()):
        print(f"  [skip dev] {subsystem}: no data")
        return False
    concat = np.concatenate([v for v in all_fluxes.values() if len(v) > 0])
    xmin, xmax = np.percentile(concat, [0.5, 99.5])
    xmin -= 0.5; xmax += 0.5
    xgrid = np.linspace(xmin, xmax, 400)

    # Compute densities and find global ymax
    densities = {}
    ymax = 0
    for cond in CONDS:
        for h in HPFS:
            cell = f"{cond}_{h}"
            flux = all_fluxes[cell]
            if len(flux) < 10:
                densities[cell] = None; continue
            try:
                kde = gaussian_kde(flux, bw_method=0.3)
                d = kde(xgrid)
                densities[cell] = d
                ymax = max(ymax, d.max())
            except Exception:
                densities[cell] = None
    if ymax == 0:
        print(f"  [skip dev] {subsystem}: KDE failed"); return False

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for i, cond in enumerate(CONDS):
        ax = axes[i]
        for h in HPFS:
            cell = f"{cond}_{h}"
            d = densities[cell]
            if d is None: continue
            ax.fill_between(xgrid, 0, d, color=HPF_COLORS[h], alpha=HPF_ALPHA,
                            edgecolor=HPF_COLORS[h], linewidth=1.3)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, ymax * 1.08)
        ax.set_title(cond, fontsize=12, color="#333333", fontweight="bold")
        ax.set_xlabel("log₂(|flux|)", fontsize=10)
        if i == 0:
            ax.set_ylabel("density", fontsize=10)
    # Legend for hpf
    handles = [Patch(facecolor=HPF_COLORS[h], alpha=HPF_ALPHA,
                      edgecolor=HPF_COLORS[h], linewidth=1.5, label=f"{h} hpf")
               for h in HPFS]
    fig.legend(handles=handles, loc="upper center", ncol=5,
               bbox_to_anchor=(0.5, 1.02), fontsize=10, frameon=False)
    fig.suptitle(subsystem, fontsize=13, y=1.09, color="#111111")
    fig.tight_layout()
    fname = _slug(subsystem)
    fig.savefig(OUT_DEV / f"{fname}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DEV / f"{fname}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  [dev] {fname}.pdf")
    return True


def _lr_plot(subsystem, sig_hpfs):
    """5 panels (one per hpf); each panel shows 3 overlaid regime density curves.
    Uses CENTERED log flux (per-cell median subtracted)."""
    all_fluxes = {}
    for cond in CONDS:
        for h in HPFS:
            cell = f"{cond}_{h}"
            all_fluxes[cell] = _collect_log_fluxes(cell, subsystem, center=True)
    if not any(len(v) > 0 for v in all_fluxes.values()):
        print(f"  [skip lr] {subsystem}: no data")
        return False
    concat = np.concatenate([v for v in all_fluxes.values() if len(v) > 0])
    xmin, xmax = np.percentile(concat, [0.5, 99.5])
    xmin -= 0.5; xmax += 0.5
    xgrid = np.linspace(xmin, xmax, 400)

    densities = {}
    ymax = 0
    for cond in CONDS:
        for h in HPFS:
            cell = f"{cond}_{h}"
            flux = all_fluxes[cell]
            if len(flux) < 10:
                densities[cell] = None; continue
            try:
                kde = gaussian_kde(flux, bw_method=0.3)
                d = kde(xgrid)
                densities[cell] = d
                ymax = max(ymax, d.max())
            except Exception:
                densities[cell] = None
    if ymax == 0:
        print(f"  [skip lr] {subsystem}: KDE failed"); return False

    fig, axes = plt.subplots(1, len(HPFS), figsize=(14, 3.6), sharey=True)
    for i, h in enumerate(HPFS):
        ax = axes[i]
        for cond in CONDS:
            cell = f"{cond}_{h}"
            d = densities[cell]
            if d is None: continue
            ax.fill_between(xgrid, 0, d, color=LR_COLORS[cond], alpha=LR_ALPHA,
                            edgecolor=LR_COLORS[cond], linewidth=1.3)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, ymax * 1.08)
        ax.set_title(f"{h} hpf", fontsize=11, color="#333333")
        ax.set_xlabel("centered log₂(|flux|)", fontsize=10)
        if i == 0:
            ax.set_ylabel("density", fontsize=10)
        if h in sig_hpfs:
            for spine in ax.spines.values():
                spine.set_edgecolor("#2E2E2E")
                spine.set_linewidth(2.3)
            ax.set_title(f"{h} hpf", fontsize=11, color="#000000", fontweight="bold")
    handles = [Patch(facecolor=LR_COLORS[c], alpha=LR_ALPHA,
                      edgecolor=LR_COLORS[c], linewidth=1.5, label=c) for c in CONDS]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.02), fontsize=11, frameon=False)
    fig.suptitle(subsystem, fontsize=13, y=1.08, color="#111111")
    fig.tight_layout()
    fname = _slug(subsystem)
    fig.savefig(OUT_LR / f"{fname}.pdf", bbox_inches="tight")
    fig.savefig(OUT_LR / f"{fname}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  [lr]  {fname}.pdf  sig={sig_hpfs}")
    return True


def _resolve(name):
    if name in all_subs_seen:
        return name
    matches = [s for s in all_subs_seen
               if name.lower() in s.lower() or s.lower() in name.lower()]
    return matches[0] if matches else None


# ---------- Developmental folder ----------
print(f"\n[dev] developmental plots for {len(DEV_SUBSYSTEMS)} subsystems")
n_dev = 0
for subs in DEV_SUBSYSTEMS:
    resolved = _resolve(subs)
    if resolved is None:
        print(f"  [skip dev] '{subs}' not found")
        continue
    if resolved != subs:
        print(f"  [match dev] '{subs}' → '{resolved}'")
    if _dev_plot(resolved):
        n_dev += 1
print(f"[dev done] {n_dev} plots → {OUT_DEV}")

# ---------- Light-regime folder ----------
print(f"\n[lr] light-regime plots for {len(LR_SIG)} subsystems")
n_lr = 0
for subs, sig_hpfs in LR_SIG.items():
    resolved = _resolve(subs)
    if resolved is None:
        print(f"  [skip lr] '{subs}' not found")
        continue
    if resolved != subs:
        print(f"  [match lr] '{subs}' → '{resolved}'")
    if _lr_plot(resolved, sig_hpfs):
        n_lr += 1
print(f"[lr done] {n_lr} plots → {OUT_LR}")
