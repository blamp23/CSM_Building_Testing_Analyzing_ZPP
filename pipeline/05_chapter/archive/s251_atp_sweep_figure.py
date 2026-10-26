"""
ATP sweep heatmap figure: rows = models grouped base/BL/D/LD × hpf,
cols = substrates. Panel A: ATP per C aerobic (greyed where unusable).
Panel B: binary usable matrix.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mc

REPO = Path(__file__).resolve().parents[2]
DIR  = REPO / "grace_v4/results/atp_sweep_v7"
df   = pd.read_csv(DIR / "atp_sweep_v7.csv")
panel = pd.read_csv(DIR / "atp_sweep_substrates.csv")

HPFS = (24, 48, 72, 96, 120)
CONDS = ("BL", "D", "LD")
ROW_ORDER = [f"base_{h}" for h in HPFS] + [f"{c}_{h}" for c in CONDS for h in HPFS]

# column order: group substrates by rough class, then by carbon count
CLASS = {
    "glucose":"sugar","fructose":"sugar","galactose":"sugar","mannose":"sugar",
    "glycerol":"polyol","inositol":"polyol",
    "L-lactate":"short-C","pyruvate":"short-C","acetate":"short-C",
    "beta-hydroxybutyrate":"ketone","acetoacetate":"ketone",
    "palmitate":"FA","stearate":"FA","oleate":"FA","linoleate":"FA",
    "linolenate":"FA","DHA":"FA","EPA":"FA","choline":"other",
    "adenosine":"nucleoside","guanosine":"nucleoside","cytidine":"nucleoside",
    "uridine":"nucleoside","thymidine":"nucleoside",
}
AA = [s for s in panel.label if s.startswith("L-")]
for s in AA: CLASS[s] = "AA"
CLASS_ORDER = ["sugar","polyol","short-C","ketone","FA","AA","nucleoside","other"]

panel["cls"] = panel.label.map(CLASS).fillna("other")
panel = panel.sort_values(by=["cls","carbon","label"],
                          key=lambda s: s.map({c:i for i,c in enumerate(CLASS_ORDER)}) if s.name=="cls" else s)
SUBS = panel.label.tolist()

# Build matrices
atp_mat = np.full((len(ROW_ORDER), len(SUBS)), np.nan)
usable_mat = np.zeros_like(atp_mat, dtype=bool)
for i, m in enumerate(ROW_ORDER):
    sub = df[df.model == m].set_index("substrate")
    for j, s in enumerate(SUBS):
        if s in sub.index:
            r = sub.loc[s]
            atp_mat[i, j] = r.atp_per_C_aerobic
            usable_mat[i, j] = bool(r.usable) and np.isfinite(r.atp_per_C_aerobic)

# Mask unusable in panel A
A = np.where(usable_mat, atp_mat, np.nan)

# Figure
fig = plt.figure(figsize=(max(12, 0.32*len(SUBS) + 3), 10))
gs = fig.add_gridspec(2, 1, hspace=0.35, height_ratios=[1, 1])

ax1 = fig.add_subplot(gs[0])
vmin = np.nanpercentile(A, 2)
vmax = np.nanpercentile(A, 98)
cmap = plt.get_cmap("viridis").copy()
cmap.set_bad("#d9d9d9")  # grey for NaN/unusable
im1 = ax1.imshow(A, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
ax1.set_yticks(range(len(ROW_ORDER))); ax1.set_yticklabels(ROW_ORDER, fontsize=8)
ax1.set_xticks(range(len(SUBS))); ax1.set_xticklabels(SUBS, rotation=90, fontsize=7)
ax1.set_title("A. ATP per carbon (aerobic)  — grey = unusable (ATP_aer ≤ 0.5)")
# group separators
for sep in [5, 10, 15]:
    ax1.axhline(sep - 0.5, color="k", lw=0.6)
# class dividers
prev = panel.cls.iloc[0]; xpos = 0
for k, c in enumerate(panel.cls):
    if c != prev:
        ax1.axvline(k - 0.5, color="k", lw=0.6, alpha=0.4)
        prev = c
plt.colorbar(im1, ax=ax1, fraction=0.015, pad=0.01, label="ATP / C")

ax2 = fig.add_subplot(gs[1])
im2 = ax2.imshow(usable_mat.astype(int), aspect="auto",
                 cmap=mc.ListedColormap(["#ffffff", "#2c7fb8"]),
                 vmin=0, vmax=1, interpolation="nearest")
ax2.set_yticks(range(len(ROW_ORDER))); ax2.set_yticklabels(ROW_ORDER, fontsize=8)
ax2.set_xticks(range(len(SUBS))); ax2.set_xticklabels(SUBS, rotation=90, fontsize=7)
ax2.set_title("B. Usable (ATP_aer > 0.5)  — blue = yes")
for sep in [5, 10, 15]:
    ax2.axhline(sep - 0.5, color="k", lw=0.6)
prev = panel.cls.iloc[0]
for k, c in enumerate(panel.cls):
    if c != prev:
        ax2.axvline(k - 0.5, color="k", lw=0.6, alpha=0.4)
        prev = c

fig.suptitle("ATP yield capacity across v7 models (base + 15 extracted)", y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig(DIR / "atp_sweep_heatmap.pdf", bbox_inches="tight")
fig.savefig(DIR / "atp_sweep_heatmap.png", dpi=200, bbox_inches="tight")
print(f"[write] atp_sweep_heatmap.{{pdf,png}}  shape={atp_mat.shape}")
