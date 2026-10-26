"""
Figure for pFBA essentials:
  A. Per-probe stacked bar: n reactions by classification
  B. Recurrent-rxn table: reactions appearing in >=2 probes, flux × probe heatmap
  C. Spotlight on 'extraction_artefact_despite_expression' hits
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
DIR  = REPO / "grace_v4/results/atp_sweep_v7"
d = pd.read_csv(DIR / "atp_sweep_pfba_essentials.csv")

d["probe"] = d.substrate + " → " + d.failing_cell
probes = sorted(d.probe.unique(),
                key=lambda p: (int(d[d.probe == p].hpf.iloc[0]), p))

CLASS_ORDER = ["extraction_artefact_despite_expression", "transcriptional",
               "no_gene_data", "extraction_artefact_no_gpr"]
CLASS_COLOR = {
    "extraction_artefact_despite_expression": "#d7301f",
    "transcriptional":                        "#2b8cbe",
    "no_gene_data":                           "#bdbdbd",
    "extraction_artefact_no_gpr":             "#636363",
}

# ---------- Panel A: stacked bar --------------------------------------------
counts = (d.groupby(["probe","classification"]).size()
          .unstack(fill_value=0).reindex(probes))
for c in CLASS_ORDER:
    if c not in counts.columns: counts[c] = 0
counts = counts[CLASS_ORDER]

# ---------- Panel B: recurrent reactions ------------------------------------
rxn_probe_cnt = d.groupby("rxn").probe.nunique().sort_values(ascending=False)
recurrent = rxn_probe_cnt[rxn_probe_cnt >= 2].index.tolist()
flux_mat = (d[d.rxn.isin(recurrent)]
            .pivot_table(index="rxn", columns="probe",
                         values="pfba_flux", aggfunc="first")
            .reindex(index=recurrent, columns=probes))

# ---------- figure ----------------------------------------------------------
fig = plt.figure(figsize=(13.5, 10))
gs = fig.add_gridspec(2, 1, height_ratios=[1.1, 1.4], hspace=0.45)

# A
axA = fig.add_subplot(gs[0])
bottom = np.zeros(len(probes))
for c in CLASS_ORDER:
    vals = counts[c].values
    axA.barh(probes, vals, left=bottom, color=CLASS_COLOR[c],
             edgecolor="white", linewidth=0.5, label=c)
    bottom += vals
axA.invert_yaxis()
axA.set_xlabel("# pFBA-route reactions missing in failing cell")
axA.set_title("A.  Classification of pFBA catabolic-route reactions "
              "present in base, absent in failing cell")
axA.legend(loc="lower right", fontsize=8, frameon=True)
for i, p in enumerate(probes):
    axA.text(bottom[i] + 0.3, i, f"n={int(bottom[i])}",
             va="center", fontsize=8)

# B
axB = fig.add_subplot(gs[1])
if len(recurrent):
    import matplotlib.colors as mc
    v_abs_max = np.nanmax(np.abs(flux_mat.values))
    cmap = plt.get_cmap("RdBu_r")
    im = axB.imshow(flux_mat.values, aspect="auto", cmap=cmap,
                    vmin=-v_abs_max, vmax=v_abs_max, interpolation="nearest")
    axB.set_yticks(range(len(recurrent))); axB.set_yticklabels(recurrent, fontsize=8)
    axB.set_xticks(range(len(probes))); axB.set_xticklabels(probes, rotation=45, ha="right", fontsize=8)
    plt.colorbar(im, ax=axB, fraction=0.025, pad=0.01, label="pFBA flux")
    # annotate missing cells in grey
    for i in range(len(recurrent)):
        for j in range(len(probes)):
            if np.isnan(flux_mat.values[i, j]):
                axB.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           facecolor="#eeeeee", edgecolor="none"))
    axB.set_title(f"B.  Recurrent reactions ({len(recurrent)} appearing in ≥2 probes)  "
                  "— cell = pFBA flux in base; grey = not in that probe's route")

fig.suptitle("pFBA essentials: catabolic-route gaps in condition-specific failures",
             y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig(DIR / "atp_sweep_pfba_essentials.pdf", bbox_inches="tight")
fig.savefig(DIR / "atp_sweep_pfba_essentials.png", dpi=200, bbox_inches="tight")
print(f"[write] atp_sweep_pfba_essentials.{{pdf,png}}")
print(f"  n_probes={len(probes)}  n_recurrent={len(recurrent)}")

# ---------- Panel C: spotlight on the 2 'despite_expression' hits ----------
hits = d[d.classification == "extraction_artefact_despite_expression"].copy()
if len(hits):
    fig2, axes = plt.subplots(1, len(hits), figsize=(5 * len(hits), 4.5))
    if len(hits) == 1: axes = [axes]
    for ax, (_, r) in zip(axes, hits.iterrows()):
        parts = r.gene_zfpkms.split(";")
        names, zs = [], []
        for p in parts:
            if "=" not in p: continue
            n, z = p.split("=")
            try:
                zs.append(float(z)); names.append(n)
            except ValueError:
                continue
        # threshold line
        from s250_atp_sweep import o2_demand_for_complete_oxidation  # (unused but keeps import tested)
        thr_fp = REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv"
        thr_tab = pd.read_csv(thr_fp)
        thr_val = float(thr_tab.loc[thr_tab.iloc[:, 0] == r.failing_cell,
                                    "exp_threshold_zFPKM"].iloc[0])
        colors = ["#d7301f" if z >= thr_val else "#2b8cbe" for z in zs]
        ax.bar(range(len(names)), zs, color=colors)
        ax.axhline(thr_val, color="k", lw=1, ls="--",
                   label=f"exp threshold = {thr_val:.2f}")
        ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=0, fontsize=9)
        ax.set_ylabel("zFPKM")
        ax.set_title(f"{r.substrate} → {r.failing_cell}\n{r.rxn} ({r.subsys})\nGPR: {r.gpr}", fontsize=9)
        ax.legend(fontsize=8)
    fig2.suptitle("C.  'Extraction artefact despite expression' hits "
                  "(GPR had ≥1 gene above threshold)", y=1.02)
    fig2.tight_layout()
    fig2.savefig(DIR / "atp_sweep_pfba_despite_expression.pdf", bbox_inches="tight")
    fig2.savefig(DIR / "atp_sweep_pfba_despite_expression.png", dpi=200, bbox_inches="tight")
    print(f"[write] atp_sweep_pfba_despite_expression.{{pdf,png}}  n={len(hits)}")
