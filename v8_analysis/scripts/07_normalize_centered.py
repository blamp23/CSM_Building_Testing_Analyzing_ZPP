"""
s274 -- Apply per-cell-median centering to F1/F2 per-reaction tables and
regenerate subsystem summaries with normalized log2FC.

Shift: log2FC_centered = log2FC_raw - (median_logflux_A - median_logflux_B)
where median_logflux_X is the median of log2(mean non-zero flux) across
all 3,905 shared reactions in cell X (same quantity used in PCA centering).

Then recompute:
  - per-reaction `notable` flag using |log2FC_centered| >= 1
  - subsystem rollup counts (up/down/notable_up/notable_down, coherence)
  - enrichment Fisher exact + BH q within contrast

Writes:
  per_reaction_flux_F1_centered.csv
  per_reaction_flux_F2_centered.csv
  subsystem_summary_F1_centered.csv
  subsystem_summary_F2_centered.csv
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, false_discovery_control

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import POOLED as POOL, COMPARISON_TABLES as OUT
OUT.mkdir(parents=True, exist_ok=True)

CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]
ZERO_TOL = 1e-10
LOG2FC_NOTABLE = 1.0

# ---------- 1. Compute per-cell median log flux across shared reactions ----------
print("[prep] computing per-cell median log flux")
cell_vectors = {}
for name in CELLS:
    df = pd.read_csv(POOL / f"{name}_pooled.csv", index_col=0)
    vals = np.abs(df.values.astype(float))
    nz_mask = vals > ZERO_TOL
    masked = np.where(nz_mask, vals, np.nan)
    with np.errstate(all="ignore"):
        mean_nz = np.nanmean(masked, axis=1)
    cell_vectors[name] = pd.Series(mean_nz, index=df.index)

all_rxns = set.intersection(*[set(v[v.notna()].index) for v in cell_vectors.values()])
all_rxns = sorted(all_rxns)
cell_medians = {}
for name in CELLS:
    v = cell_vectors[name].reindex(all_rxns).values
    cell_medians[name] = float(np.nanmedian(np.log2(v + 1e-6)))
    print(f"  {name:7s}  median log2 flux (across {len(all_rxns)} shared rxns): {cell_medians[name]:+.3f}")


# ---------- 2. Apply centering to per-reaction tables ----------
def _coherence(up, dn):
    tot = up + dn
    return abs(up - dn) / tot if tot > 0 else np.nan


for fam in ("F1", "F2"):
    print(f"\n[normalize] {fam}")
    rxn = pd.read_csv(OUT / f"per_reaction_flux_{fam}.csv")
    # Compute median_shift per contrast
    rxn["median_log2_A"] = rxn.cell_A.map(cell_medians)
    rxn["median_log2_B"] = rxn.cell_B.map(cell_medians)
    rxn["median_shift"] = rxn["median_log2_A"] - rxn["median_log2_B"]
    # Centered log2FC
    rxn["log2FC_centered"] = rxn["log2FC_mean"] - rxn["median_shift"]
    rxn["direction_centered"] = np.sign(rxn["log2FC_centered"]).fillna(0).astype(int)
    rxn["notable_centered"] = (rxn.carry_both & (~rxn.loop_flagged) &
                                (rxn.log2FC_centered.abs() >= LOG2FC_NOTABLE))

    rxn.to_csv(OUT / f"per_reaction_flux_{fam}_centered.csv", index=False)
    print(f"  per_reaction_flux_{fam}_centered.csv: {len(rxn)} rows")
    print(f"    notable_raw: {int(((rxn.carry_both)&(~rxn.loop_flagged)&(rxn.log2FC_mean.abs()>=1)).sum())}")
    print(f"    notable_centered: {int(rxn.notable_centered.sum())}")

    # ---------- Subsystem rollup with centered log2FC ----------
    roll_rows = []
    for (contrast, subsystem), gs in rxn.groupby(["contrast", "subsystem"]):
        if not subsystem: continue
        cA = gs.cell_A.iloc[0]; cB = gs.cell_B.iloc[0]
        tested = gs[gs.carry_both & (~gs.loop_flagged)]
        n_reps_tested = len(tested)
        n_up = int((tested.direction_centered == 1).sum())
        n_down = int((tested.direction_centered == -1).sum())
        n_notable_up = int((tested[tested.log2FC_centered >= LOG2FC_NOTABLE]).shape[0])
        n_notable_down = int((tested[tested.log2FC_centered <= -LOG2FC_NOTABLE]).shape[0])
        frac_notable = (n_notable_up + n_notable_down) / n_reps_tested if n_reps_tested else np.nan
        coherence = _coherence(n_up, n_down)
        if len(tested) and tested.log2FC_centered.notna().any():
            vals = tested.log2FC_centered.dropna()
            mean_lfc = float(vals.mean()); median_lfc = float(vals.median())
            idx = tested.log2FC_centered.abs().idxmax()
            if pd.notna(idx):
                top = tested.loc[idx]
                max_abs_lfc = float(top.log2FC_centered)
                max_abs_rxn = str(top.rxn_id)
            else:
                max_abs_lfc = np.nan; max_abs_rxn = ""
        else:
            mean_lfc = median_lfc = max_abs_lfc = np.nan; max_abs_rxn = ""

        # Verdict logic (same decision tree as before, applied to centered values)
        n_rxns_shared = len(gs)
        n_rxns_only_A = int(gs.carrying_asymmetric.sum() if 'carrying_asymmetric' in gs else 0)
        verdict = "flat"
        if (n_reps_tested >= 5 and pd.notna(coherence) and coherence >= 0.7
            and pd.notna(mean_lfc) and abs(mean_lfc) >= 0.5):
            verdict = "flux_A_higher" if mean_lfc > 0 else "flux_B_higher"
        elif n_reps_tested >= 5 and pd.notna(coherence) and coherence < 0.5:
            verdict = "mixed_flux"

        roll_rows.append(dict(
            family=fam, contrast=contrast, cell_A=cA, cell_B=cB,
            subsystem=subsystem,
            n_reps_shared=n_rxns_shared, n_reps_tested=n_reps_tested,
            n_reps_up=n_up, n_reps_down=n_down,
            n_reps_notable_up=n_notable_up,
            n_reps_notable_down=n_notable_down,
            frac_notable=frac_notable, coherence=coherence,
            mean_log2FC_centered=mean_lfc,
            median_log2FC_centered=median_lfc,
            max_abs_log2FC_centered=max_abs_lfc,
            max_abs_rxn=max_abs_rxn,
            verdict=verdict,
        ))

    roll = pd.DataFrame(roll_rows)

    # ---------- Hypergeometric enrichment on centered-notable ----------
    enr_rows = []
    for contrast, g in rxn.groupby("contrast"):
        tested = g[g.carry_both & (~g.loop_flagged)]
        total_tested = len(tested)
        total_notable = int(tested.notable_centered.sum())
        for sub, gs in tested.groupby("subsystem"):
            if not sub: continue
            sub_tested = len(gs)
            sub_notable = int(gs.notable_centered.sum())
            a = sub_notable
            b = sub_tested - sub_notable
            c = total_notable - sub_notable
            d = total_tested - total_notable - b
            if a + b < 1 or a + c < 1:
                p = np.nan; odds = np.nan
            else:
                try:
                    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
                except Exception:
                    odds, p = np.nan, np.nan
            enr_rows.append(dict(contrast=contrast, subsystem=sub,
                                 odds_ratio=odds, enrichment_p=p))
    enr = pd.DataFrame(enr_rows)
    enr["enrichment_q"] = np.nan
    for contrast, idx in enr.groupby("contrast").groups.items():
        ps = enr.loc[idx, "enrichment_p"].values
        good = np.isfinite(ps)
        if good.sum() > 0:
            q = np.full_like(ps, np.nan, dtype=float)
            q[good] = false_discovery_control(ps[good], method="bh")
            enr.loc[idx, "enrichment_q"] = q

    roll_merged = roll.merge(enr[["contrast","subsystem","odds_ratio",
                                   "enrichment_p","enrichment_q"]],
                              on=["contrast","subsystem"], how="left")
    roll_merged.to_csv(OUT / f"subsystem_summary_{fam}_centered.csv", index=False)
    n_sig = (roll_merged.enrichment_q < 0.05).sum()
    print(f"  subsystem_summary_{fam}_centered.csv: {len(roll_merged)} rows, "
          f"{n_sig} with enrichment_q < 0.05 (centered)")

print("\n[done]")

# ---------- Comparison summary ----------
print("\n=== Comparison: raw vs centered subsystem findings ===")
for fam in ("F1", "F2"):
    raw = pd.read_csv(OUT / f"subsystem_summary_{fam}.csv")
    cen = pd.read_csv(OUT / f"subsystem_summary_{fam}_centered.csv")
    n_raw_sig = (raw.enrichment_q < 0.05).sum() if "enrichment_q" in raw else 0
    n_cen_sig = (cen.enrichment_q < 0.05).sum() if "enrichment_q" in cen else 0
    print(f"\n{fam}: raw q<0.05: {n_raw_sig}  |  centered q<0.05: {n_cen_sig}")

    # Verdict distribution
    print(f"  verdict distribution:")
    print(f"    raw:")
    for v, n in raw.verdict.value_counts().items():
        print(f"      {v}: {n}")
    print(f"    centered:")
    for v, n in cen.verdict.value_counts().items():
        print(f"      {v}: {n}")

print("\n=== Top enriched subsystems in CENTERED tables (q<0.05) ===")
for fam in ("F2", "F1"):
    df = pd.read_csv(OUT / f"subsystem_summary_{fam}_centered.csv")
    top = df.loc[df.enrichment_q < 0.05].nsmallest(10, "enrichment_q")
    print(f"\n{fam} top 10 (centered):")
    cols = ["contrast","subsystem","n_reps_tested","n_reps_notable_up",
            "n_reps_notable_down","mean_log2FC_centered","coherence",
            "enrichment_q","verdict"]
    print(top[cols].to_string(index=False))
