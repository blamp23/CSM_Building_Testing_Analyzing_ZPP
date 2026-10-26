"""
s270 -- Per-reaction DESeq-style analysis on hopsy-sampled, CFF-corrected
flux distributions.

For every (representative reaction, contrast) pair:
  - Mann-Whitney U test on non-zero |flux| samples
  - log2 fold change of non-zero mean (and median)
  - direction = sign(log2FC_mean)
  - zero-fraction Fisher test (asymmetric zero-inflation)
  - loop flag (max |v| > 100 in either cell)

Then BH-correct Mann-Whitney p per contrast (not per subsystem — standard
DESeq-shaped per-item multiple testing). Also roll up to subsystem by
aggregation (counts up/down/sig, mean |log2FC|, top mover, coherence).

Output (into grace_v4/results/chapter_assets/tables/):
    per_reaction_F1.csv
    per_reaction_F2.csv
    subsystem_rollup_F1.csv
    subsystem_rollup_F2.csv

Conventions:
  - Family 1 (same hpf, cross condition): A = alphabetically earlier
    (BL < D < LD). sign +1 means A > B.
  - Family 2 (same regime, adjacent hpf): A = earlier hpf. sign +1
    means A > B (flux falls from A to B over development).
"""
import sys, re, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact, false_discovery_control

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
R = REPO / "grace_v4/results"
POOL = R / "pooled_v7"
SUBS = R / "subsets_v7"
OUT = R / "chapter_assets/tables"

ZERO_TOL = 1e-10
CARRY_TOL = 1e-4
LOOP_MAG = 100.0
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

# ---------- Load all cell data up front ----------
print("[load] pooled samples + subset membership for 15 cells", flush=True)
cell_data = {}
for cell in CELLS:
    samples = pd.read_csv(POOL / f"{cell}_pooled.csv", index_col=0)
    membership = pd.read_csv(SUBS / f"{cell}_subset_membership.csv")
    # rep → subsystem lookup
    rep_to_sub = (membership.drop_duplicates("representative")
                            .set_index("representative")["subsystem"]
                            .to_dict())
    # rep → subset size lookup (how many physical reactions each rep stands for)
    rep_to_size = (membership.groupby("representative").size().to_dict())
    cell_data[cell] = dict(samples=samples, rep_to_sub=rep_to_sub,
                           rep_to_size=rep_to_size)
    print(f"  {cell:7s}  {samples.shape[0]:5d} reactions × {samples.shape[1]:5d} samples  "
          f"({len(rep_to_sub)} reps)", flush=True)

# ---------- Build contrast list ----------
# F1: same hpf, cross condition. Convention: A alphabetically earlier (BL<D<LD).
# F2: same condition, adjacent hpf. Convention: A earlier hpf.
contrasts = []
for hpf in HPFS:
    cells_at_hpf = [f"{c}_{hpf}" for c in CONDS]
    for i, A in enumerate(cells_at_hpf):
        for B in cells_at_hpf[i+1:]:
            contrasts.append(("F1", A, B))
for cond in CONDS:
    cells_in_regime = [f"{cond}_{h}" for h in HPFS]
    for A, B in zip(cells_in_regime[:-1], cells_in_regime[1:]):
        contrasts.append(("F2", A, B))
print(f"\n[contrasts] {sum(1 for c in contrasts if c[0]=='F1')} F1  "
      f"+ {sum(1 for c in contrasts if c[0]=='F2')} F2  = {len(contrasts)} total")


def _per_rxn_test(sA, sB):
    """One (rxn, contrast). sA, sB: 1D arrays of flux samples.
    Returns dict of statistics.
    """
    aA = np.abs(sA); aB = np.abs(sB)
    nz_A = aA[aA > ZERO_TOL]; nz_B = aB[aB > ZERO_TOL]
    n_nz_A, n_nz_B = len(nz_A), len(nz_B)
    nA, nB = len(aA), len(aB)
    zfA = 1 - n_nz_A / nA if nA else np.nan
    zfB = 1 - n_nz_B / nB if nB else np.nan
    loop = bool((aA.max() if nA else 0) > LOOP_MAG or (aB.max() if nB else 0) > LOOP_MAG)
    carry_both = bool(n_nz_A > 0 and n_nz_B > 0 and
                      np.median(aA) > CARRY_TOL and np.median(aB) > CARRY_TOL)
    carrying_asymmetric = bool((n_nz_A == 0) ^ (n_nz_B == 0))
    # Fisher on zero fractions
    try:
        fisher_p = fisher_exact([[n_nz_A, nA - n_nz_A], [n_nz_B, nB - n_nz_B]])[1]
    except Exception:
        fisher_p = np.nan
    # MW on non-zero |v|
    if n_nz_A >= 3 and n_nz_B >= 3:
        try:
            U, p = mannwhitneyu(nz_A, nz_B, alternative="two-sided")
        except Exception:
            U, p = np.nan, np.nan
    else:
        U, p = np.nan, np.nan
    # Effect sizes
    median_A = float(np.median(nz_A)) if n_nz_A else np.nan
    median_B = float(np.median(nz_B)) if n_nz_B else np.nan
    mean_A = float(np.mean(nz_A)) if n_nz_A else np.nan
    mean_B = float(np.mean(nz_B)) if n_nz_B else np.nan
    if n_nz_A and n_nz_B and mean_A > 0 and mean_B > 0:
        log2FC_mean = float(np.log2(mean_A / mean_B))
        log2FC_median = float(np.log2(median_A / median_B)) if (median_A > 0 and median_B > 0) else np.nan
    else:
        log2FC_mean = np.nan; log2FC_median = np.nan
    direction = int(np.sign(log2FC_mean)) if np.isfinite(log2FC_mean) else 0
    return dict(
        n_samples_A=nA, n_samples_B=nB,
        n_nz_A=n_nz_A, n_nz_B=n_nz_B,
        zero_frac_A=zfA, zero_frac_B=zfB,
        fisher_p_zeros=fisher_p,
        median_nz_A=median_A, median_nz_B=median_B,
        mean_nz_A=mean_A, mean_nz_B=mean_B,
        log2FC_mean=log2FC_mean, log2FC_median=log2FC_median,
        direction=direction,
        MW_U=U, MW_p=p,
        carrying_asymmetric=carrying_asymmetric,
        loop_flagged=loop,
        carry_both=carry_both,
    )


# ---------- Loop over contrasts, build per-reaction table ----------
print(f"\n[run] per-reaction MW + log2FC across {len(contrasts)} contrasts", flush=True)
all_rows = []
t0 = time.time()
for fam, A, B in contrasts:
    dA = cell_data[A]; dB = cell_data[B]
    reps_A = set(dA["rep_to_sub"]); reps_B = set(dB["rep_to_sub"])
    common_reps = sorted(reps_A & reps_B)
    sA_df = dA["samples"]; sB_df = dB["samples"]
    for rep in common_reps:
        if rep not in sA_df.index or rep not in sB_df.index:
            continue
        sA = sA_df.loc[rep].values.astype(float)
        sB = sB_df.loc[rep].values.astype(float)
        stats = _per_rxn_test(sA, sB)
        row = dict(
            family=fam, contrast=f"{A}_vs_{B}", cell_A=A, cell_B=B,
            rxn_id=rep,
            subsystem=dA["rep_to_sub"].get(rep, dB["rep_to_sub"].get(rep, "UNMAPPED")),
            subset_size_A=dA["rep_to_size"].get(rep, 1),
            subset_size_B=dB["rep_to_size"].get(rep, 1),
        )
        row.update(stats)
        all_rows.append(row)
    print(f"  {fam} {A:7s} vs {B:7s}  reps_tested = {len(common_reps):4d}  "
          f"elapsed = {(time.time()-t0)/60:.1f}m", flush=True)

df_rxn = pd.DataFrame(all_rows)
print(f"\n[df_rxn] {len(df_rxn)} rows")

# ---------- BH correct MW_p per contrast ----------
df_rxn["MW_q"] = np.nan
for contrast, idx in df_rxn.groupby("contrast").groups.items():
    sub = df_rxn.loc[idx, "MW_p"].values
    good = np.isfinite(sub)
    if good.sum() > 0:
        q = np.full_like(sub, np.nan, dtype=float)
        q[good] = false_discovery_control(sub[good], method="bh")
        df_rxn.loc[idx, "MW_q"] = q
print(f"[bh]  computed MW_q per contrast; "
      f"{(df_rxn.MW_q < 0.05).sum()} rxn×contrast rows at q<0.05 "
      f"(hopsy samples are NOT independent — MW_q is massively over-powered; "
      f"use |log2FC_mean| as primary filter)")

# ---------- Split by family and write ----------
col_order = [
    "family", "contrast", "cell_A", "cell_B", "rxn_id", "subsystem",
    "subset_size_A", "subset_size_B",
    "n_samples_A", "n_samples_B", "n_nz_A", "n_nz_B",
    "zero_frac_A", "zero_frac_B", "fisher_p_zeros",
    "median_nz_A", "median_nz_B", "mean_nz_A", "mean_nz_B",
    "log2FC_mean", "log2FC_median", "direction",
    "MW_U", "MW_p", "MW_q",
    "carrying_asymmetric", "loop_flagged", "carry_both",
]
df_rxn = df_rxn[col_order]

for fam in ("F1", "F2"):
    sub = df_rxn[df_rxn.family == fam].drop(columns="family").reset_index(drop=True)
    fp = OUT / f"per_reaction_{fam}.csv"
    sub.to_csv(fp, index=False)
    sig = (sub.MW_q < 0.05).sum()
    print(f"  wrote per_reaction_{fam}.csv  {len(sub)} rows  "
          f"({sig} at q<0.05)")

# ---------- Subsystem aggregation ----------
print(f"\n[rollup] aggregating to subsystem × contrast")

def _coherence(up, dn):
    tot = up + dn
    return abs(up - dn) / tot if tot > 0 else np.nan

rollup = []
for (fam, contrast), g in df_rxn.groupby(["family", "contrast"]):
    cA = g.cell_A.iloc[0]; cB = g.cell_B.iloc[0]
    for sub, gs in g.groupby("subsystem"):
        tested = gs[gs.carry_both & ~gs.loop_flagged]
        n_up = int((tested.direction == 1).sum())
        n_down = int((tested.direction == -1).sum())
        sig = tested[tested.MW_q < 0.05]
        n_sig_up = int((sig.direction == 1).sum())
        n_sig_down = int((sig.direction == -1).sum())
        if len(tested) and tested.log2FC_mean.notna().any():
            vals = tested.log2FC_mean.dropna()
            mean_lfc = float(vals.mean())
            median_lfc = float(vals.median())
            max_abs_idx = tested.log2FC_mean.abs().idxmax() if tested.log2FC_mean.abs().notna().any() else None
            if max_abs_idx is not None:
                top = tested.loc[max_abs_idx]
                max_abs_lfc = float(top.log2FC_mean)
                max_abs_rxn = top.rxn_id
            else:
                max_abs_lfc = np.nan; max_abs_rxn = ""
        else:
            mean_lfc = median_lfc = max_abs_lfc = np.nan; max_abs_rxn = ""
        if len(tested) and tested.MW_q.notna().any():
            min_q = float(tested.MW_q.min())
            min_q_rxn = str(tested.loc[tested.MW_q.idxmin(), "rxn_id"])
        else:
            min_q = np.nan; min_q_rxn = ""
        rollup.append(dict(
            family=fam, contrast=contrast, cell_A=cA, cell_B=cB,
            subsystem=sub,
            n_reps_total=len(gs),
            n_reps_tested=len(tested),
            n_reps_loop_flagged=int(gs.loop_flagged.sum()),
            n_reps_not_carrying=int((~gs.carry_both).sum()),
            n_reps_up=n_up, n_reps_down=n_down,
            n_reps_sig_q05=len(sig),
            n_reps_sig_up_q05=n_sig_up,
            n_reps_sig_down_q05=n_sig_down,
            mean_log2FC=mean_lfc, median_log2FC=median_lfc,
            max_abs_log2FC=max_abs_lfc, max_abs_rxn=max_abs_rxn,
            min_MW_q=min_q, min_MW_q_rxn=min_q_rxn,
            coherence=_coherence(n_up, n_down),
        ))

df_roll = pd.DataFrame(rollup)
for fam in ("F1", "F2"):
    sub = df_roll[df_roll.family == fam].drop(columns="family").reset_index(drop=True)
    fp = OUT / f"subsystem_rollup_{fam}.csv"
    sub.to_csv(fp, index=False)
    print(f"  wrote subsystem_rollup_{fam}.csv  {len(sub)} rows")

print("\nDone.")
