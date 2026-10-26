"""
s271 -- Final between-condition comparison tables.

Three layers per contrast:
  1. Reaction inventory — which reactions does each cell retain?
     (presence/absence of the physical reaction in each extracted model)
  2. Per-reaction flux comparison — for reactions in both cells,
     magnitude + shape + direction of the flux-distribution shift.
  3. Subsystem summary — aggregates inventory + flux into a single
     per-(subsystem, contrast) row with a categorical verdict.

Writes to grace_v4/results/chapter_assets/comparison_tables/
"""
import sys, time, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, wasserstein_distance

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import POOLED as POOL, SUBSETS as SUBS, COMPARISON_TABLES as OUT
OUT.mkdir(parents=True, exist_ok=True)

ZERO_TOL = 1e-10
CARRY_TOL = 1e-4
LOOP_MAG = 100.0
LOG2FC_NOTABLE = 1.0
COHERENCE_THRESH = 0.70
STRUCT_RATIO = 2.0
STRUCT_FRAC = 0.20

CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]


# ---------- load every cell ----------
print("[load] 15 cells (samples + membership)", flush=True)
cell = {}
for name in CELLS:
    samples = pd.read_csv(POOL / f"{name}_pooled.csv", index_col=0)
    mem = pd.read_csv(SUBS / f"{name}_subset_membership.csv")
    rep_to_sub = (mem.drop_duplicates("representative")
                     .set_index("representative")["subsystem"].to_dict())
    rxn_to_sub = dict(zip(mem.rxn, mem.subsystem))
    rxn_to_rep = dict(zip(mem.rxn, mem.representative))
    cell[name] = dict(samples=samples, mem=mem,
                      rxns=set(mem.rxn), reps=set(mem.representative),
                      rep_to_sub=rep_to_sub, rxn_to_sub=rxn_to_sub,
                      rxn_to_rep=rxn_to_rep)
    print(f"  {name:7s}  {len(cell[name]['rxns']):5d} rxns  "
          f"{len(cell[name]['reps']):5d} reps  "
          f"{len(set(rep_to_sub.values()))} subsystems", flush=True)


# ---------- contrast list ----------
contrasts = []
for hpf in HPFS:
    cs = [f"{c}_{hpf}" for c in CONDS]
    for i, A in enumerate(cs):
        for B in cs[i+1:]:
            contrasts.append(("F1", A, B))
for cond in CONDS:
    cs = [f"{cond}_{h}" for h in HPFS]
    for A, B in zip(cs[:-1], cs[1:]):
        contrasts.append(("F2", A, B))


# ---------- per-reaction flux comparison ----------
def _rxn_stats(sA, sB):
    aA = np.abs(sA); aB = np.abs(sB)
    nz_A = aA[aA > ZERO_TOL]; nz_B = aB[aB > ZERO_TOL]
    nA, nB = len(aA), len(aB)
    n_nz_A, n_nz_B = len(nz_A), len(nz_B)
    zfA = 1 - n_nz_A / nA if nA else np.nan
    zfB = 1 - n_nz_B / nB if nB else np.nan
    loop = bool((aA.max() if nA else 0) > LOOP_MAG or (aB.max() if nB else 0) > LOOP_MAG)
    carry_both = bool(n_nz_A > 0 and n_nz_B > 0 and
                      np.median(aA) > CARRY_TOL and np.median(aB) > CARRY_TOL)
    asymmetric = bool((n_nz_A == 0) ^ (n_nz_B == 0))
    try:
        fisher_p = fisher_exact([[n_nz_A, nA-n_nz_A], [n_nz_B, nB-n_nz_B]])[1]
    except Exception:
        fisher_p = np.nan

    median_A = float(np.median(nz_A)) if n_nz_A else np.nan
    median_B = float(np.median(nz_B)) if n_nz_B else np.nan
    mean_A = float(np.mean(nz_A)) if n_nz_A else np.nan
    mean_B = float(np.mean(nz_B)) if n_nz_B else np.nan

    log2FC_mean = log2FC_median = np.nan
    direction = 0
    if n_nz_A and n_nz_B and mean_A > 0 and mean_B > 0:
        log2FC_mean = float(np.log2(mean_A / mean_B))
        direction = int(np.sign(log2FC_mean))
        if median_A > 0 and median_B > 0:
            log2FC_median = float(np.log2(median_A / median_B))

    # Shape columns
    W1_rxn = np.nan
    log_var_ratio = np.nan
    if n_nz_A >= 3 and n_nz_B >= 3:
        W1_rxn = float(wasserstein_distance(nz_A, nz_B))
        vA = float(np.var(nz_A)); vB = float(np.var(nz_B))
        if vA > 1e-20 and vB > 1e-20:
            log_var_ratio = float(np.log2(vA / vB))

    return dict(
        n_samples_A=nA, n_samples_B=nB,
        n_nz_A=n_nz_A, n_nz_B=n_nz_B,
        zero_frac_A=zfA, zero_frac_B=zfB,
        fisher_p_zeros=fisher_p,
        median_nz_A=median_A, median_nz_B=median_B,
        mean_nz_A=mean_A, mean_nz_B=mean_B,
        log2FC_mean=log2FC_mean, log2FC_median=log2FC_median,
        direction=direction,
        W1_rxn=W1_rxn, log_var_ratio=log_var_ratio,
        carrying_asymmetric=asymmetric, loop_flagged=loop,
        carry_both=carry_both,
    )


# ---------- iterate contrasts, build three tables ----------
inv_rows = []
rxn_rows = []
roll_rows = []

print(f"\n[run] {len(contrasts)} contrasts", flush=True)
t0 = time.time()

for fam, A, B in contrasts:
    dA, dB = cell[A], cell[B]

    # ---- INVENTORY LAYER ----
    all_rxns = sorted(dA["rxns"] | dB["rxns"])
    for rxn in all_rxns:
        in_A = rxn in dA["rxns"]; in_B = rxn in dB["rxns"]
        sub_A = dA["rxn_to_sub"].get(rxn, "")
        sub_B = dB["rxn_to_sub"].get(rxn, "")
        subs = sub_A if in_A else sub_B
        if in_A and in_B: cat = "both"
        elif in_A: cat = "A_only"
        else: cat = "B_only"
        inv_rows.append(dict(
            family=fam, contrast=f"{A}_vs_{B}", cell_A=A, cell_B=B,
            rxn_id=rxn, subsystem=subs,
            in_A=in_A, in_B=in_B, category=cat,
            representative_A=dA["rxn_to_rep"].get(rxn, ""),
            representative_B=dB["rxn_to_rep"].get(rxn, ""),
        ))

    # ---- PER-REACTION FLUX LAYER (common representatives only) ----
    sA_df = dA["samples"]; sB_df = dB["samples"]
    common_reps = sorted(dA["reps"] & dB["reps"] & set(sA_df.index) & set(sB_df.index))
    for rep in common_reps:
        sA = sA_df.loc[rep].values.astype(float)
        sB = sB_df.loc[rep].values.astype(float)
        stats = _rxn_stats(sA, sB)
        sub = dA["rep_to_sub"].get(rep, dB["rep_to_sub"].get(rep, "UNMAPPED"))
        row = dict(
            family=fam, contrast=f"{A}_vs_{B}", cell_A=A, cell_B=B,
            rxn_id=rep, subsystem=sub,
            subset_size_A=int((dA["mem"].representative == rep).sum()),
            subset_size_B=int((dB["mem"].representative == rep).sum()),
        )
        row.update(stats)
        rxn_rows.append(row)
    print(f"  {fam} {A:7s} vs {B:7s}  "
          f"inv={len(all_rxns):4d}  flux_tested={len(common_reps):4d}  "
          f"({(time.time()-t0)/60:.1f}m)", flush=True)

print(f"\n[inv] {len(inv_rows):,} rows")
print(f"[flux] {len(rxn_rows):,} rows")

df_inv = pd.DataFrame(inv_rows)
df_rxn = pd.DataFrame(rxn_rows)

# ---- SUBSYSTEM SUMMARY LAYER ----
print("\n[rollup] subsystem summaries")

def _coherence(up, dn):
    tot = up + dn
    return abs(up - dn) / tot if tot > 0 else np.nan

for (fam, contrast), g_inv in df_inv.groupby(["family", "contrast"]):
    cA = g_inv.cell_A.iloc[0]; cB = g_inv.cell_B.iloc[0]
    g_flux = df_rxn[(df_rxn.family == fam) & (df_rxn.contrast == contrast)]
    for sub in sorted(g_inv.subsystem.dropna().unique()):
        if not sub: continue
        inv_sub = g_inv[g_inv.subsystem == sub]
        flux_sub = g_flux[g_flux.subsystem == sub]

        # Inventory counts (physical reactions)
        n_both = int((inv_sub.category == "both").sum())
        n_A_only = int((inv_sub.category == "A_only").sum())
        n_B_only = int((inv_sub.category == "B_only").sum())
        n_total_A = n_both + n_A_only
        n_total_B = n_both + n_B_only
        n_union = n_both + n_A_only + n_B_only

        # Flux comparison (representatives, only shared ones)
        tested = flux_sub[flux_sub.carry_both & ~flux_sub.loop_flagged]
        n_reps_shared = len(flux_sub)
        n_reps_tested = len(tested)
        n_up = int((tested.direction == 1).sum())
        n_down = int((tested.direction == -1).sum())
        n_notable_up = int((tested[tested.log2FC_mean >= LOG2FC_NOTABLE]).shape[0])
        n_notable_down = int((tested[tested.log2FC_mean <= -LOG2FC_NOTABLE]).shape[0])
        frac_notable = ((n_notable_up + n_notable_down) / n_reps_tested
                        if n_reps_tested else np.nan)
        coherence = _coherence(n_up, n_down)

        # Magnitude
        if len(tested) and tested.log2FC_mean.notna().any():
            vals = tested.log2FC_mean.dropna()
            mean_lfc = float(vals.mean())
            median_lfc = float(vals.median())
            idx = tested.log2FC_mean.abs().idxmax()
            if pd.notna(idx):
                top = tested.loc[idx]
                max_abs_lfc = float(top.log2FC_mean)
                max_abs_rxn = str(top.rxn_id)
            else:
                max_abs_lfc = np.nan; max_abs_rxn = ""
            median_W1 = float(tested.W1_rxn.median()) if tested.W1_rxn.notna().any() else np.nan
            median_log_var_ratio = float(tested.log_var_ratio.median()) if tested.log_var_ratio.notna().any() else np.nan
        else:
            mean_lfc = median_lfc = max_abs_lfc = median_W1 = median_log_var_ratio = np.nan
            max_abs_rxn = ""

        # Verdict categorical
        verdict = "flat"
        # Structural asymmetry wins first
        if n_union >= 5:
            if (n_A_only >= STRUCT_RATIO * max(n_B_only, 1) and
                n_A_only / n_union >= STRUCT_FRAC):
                verdict = "structural_A_rich"
            elif (n_B_only >= STRUCT_RATIO * max(n_A_only, 1) and
                  n_B_only / n_union >= STRUCT_FRAC):
                verdict = "structural_B_rich"
        if verdict == "flat":
            if (n_reps_tested >= 5
                and pd.notna(coherence) and coherence >= COHERENCE_THRESH
                and pd.notna(mean_lfc) and abs(mean_lfc) >= 0.5):
                verdict = "flux_A_higher" if mean_lfc > 0 else "flux_B_higher"
            elif (n_reps_tested >= 5 and pd.notna(coherence)
                  and coherence < 0.5):
                verdict = "mixed_flux"

        roll_rows.append(dict(
            family=fam, contrast=contrast, cell_A=cA, cell_B=cB, subsystem=sub,
            # Inventory
            n_rxns_total_A=n_total_A, n_rxns_total_B=n_total_B,
            n_rxns_shared=n_both,
            n_rxns_A_only=n_A_only, n_rxns_B_only=n_B_only,
            n_rxns_union=n_union,
            # Flux comparison
            n_reps_shared=n_reps_shared,
            n_reps_tested=n_reps_tested,
            n_reps_loop_flagged=int(flux_sub.loop_flagged.sum()),
            n_reps_not_carrying=int((~flux_sub.carry_both).sum()),
            n_reps_up=n_up, n_reps_down=n_down,
            n_reps_notable_up=n_notable_up, n_reps_notable_down=n_notable_down,
            frac_notable=frac_notable, coherence=coherence,
            # Magnitude
            mean_log2FC=mean_lfc, median_log2FC=median_lfc,
            max_abs_log2FC=max_abs_lfc, max_abs_rxn=max_abs_rxn,
            median_W1_rxn=median_W1,
            median_log_var_ratio=median_log_var_ratio,
            # Verdict
            verdict=verdict,
        ))

df_roll = pd.DataFrame(roll_rows)

# ---------- write per family ----------
print("\n[write]")
for fam in ("F1", "F2"):
    inv_f = df_inv[df_inv.family == fam].drop(columns="family").reset_index(drop=True)
    rxn_f = df_rxn[df_rxn.family == fam].drop(columns="family").reset_index(drop=True)
    roll_f = df_roll[df_roll.family == fam].drop(columns="family").reset_index(drop=True)
    inv_f.to_csv(OUT / f"reaction_inventory_{fam}.csv", index=False)
    rxn_f.to_csv(OUT / f"per_reaction_flux_{fam}.csv", index=False)
    roll_f.to_csv(OUT / f"subsystem_summary_{fam}.csv", index=False)
    print(f"  {fam}: inventory={len(inv_f):,}  flux={len(rxn_f):,}  rollup={len(roll_f):,}")

# ---------- verdict summary ----------
print("\n[verdict distribution]")
for fam in ("F1", "F2"):
    roll_f = df_roll[df_roll.family == fam]
    print(f"\n{fam}:")
    print(roll_f.verdict.value_counts().to_string())

print(f"\n[time] {(time.time()-t0)/60:.1f} min total")
