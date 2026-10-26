# Between-condition flux comparison tables

Primary output of `grace_v4/scripts/s271_final_comparison_tables.py`.

These tables describe the differences between pairs of extracted
rFASTCORMICS models (one per cell = condition × hpf) at three layers:

1. **Reaction inventory** — which physical reactions each cell retains.
2. **Per-reaction flux comparison** — for reactions in both cells, how
   their flux distributions differ.
3. **Subsystem summary** — aggregates of the above into one row per
   (subsystem, contrast) with a categorical verdict.

Each layer has two files, one per contrast family:

| contrast family | what it compares | file count |
|---|---|---|
| **F1** — same hpf, cross-condition | BL / D / LD at the same developmental stage | 3 files |
| **F2** — same regime, adjacent hpf | within a single light condition, next-stage change | 3 files |

Convention for every pairwise comparison in both families:
- **F1**: `A` is the alphabetically earlier condition (BL before D before LD).
- **F2**: `A` is the earlier hpf.
- `direction = +1` means A > B. `direction = −1` means B > A.
- `log2FC_mean = log2(mean_nz_A / mean_nz_B)` so positive = A higher, negative = B higher.

---

## Layer 1 — Reaction inventory

**Files:** `reaction_inventory_F1.csv`, `reaction_inventory_F2.csv`

One row per (physical reaction × contrast). Expected size: ~150,000 rows.

**Columns:**

| column | meaning |
|---|---|
| `contrast, cell_A, cell_B` | contrast identifier |
| `rxn_id` | reaction MAR ID |
| `subsystem` | Wang GEM subsystem annotation |
| `in_A, in_B` | boolean: is this reaction in each cell's extracted model? |
| `category` | one of `both` / `A_only` / `B_only` |
| `representative_A, representative_B` | which enzyme-subset representative this reaction maps to in each cell (may differ between cells) |

**What this layer answers:**
- "Which reactions are present in BL_96 but not in LD_96?"
- "How many reactions in Drug metabolism does cell A retain vs cell B?"
- "Is this reaction in both cells or only one?"

**Important:** This layer is at the **physical-reaction level**. Reactions
that are mathematically coupled (share a GPR, move in lockstep) have
different rows but the same `representative`. For the flux-level
comparison use Layer 2.

---

## Layer 2 — Per-reaction flux comparison

**Files:** `per_reaction_flux_F1.csv`, `per_reaction_flux_F2.csv`

One row per (shared representative × contrast). Expected size: ~65,000 rows.

Only representatives that exist in **both** cells are tested here — these
are the reactions where a flux comparison is well-defined. Presence-only
hits live in Layer 1 (`category == "A_only"` or `"B_only"`).

**Columns:**

**Identifiers**

| column | meaning |
|---|---|
| `contrast, cell_A, cell_B` | contrast |
| `rxn_id` | representative reaction |
| `subsystem` | Wang GEM subsystem |
| `subset_size_A, subset_size_B` | how many physical reactions this rep stands for in each cell |

**Status flags**

| column | meaning |
|---|---|
| `carry_both` | median non-zero flux > 1e-4 in both cells (yes = the reaction is "on" in both) |
| `loop_flagged` | max \|v\| > 100 in either cell (possible thermodynamic loop; interpret with care) |
| `carrying_asymmetric` | one cell is all-zero on this reaction (present but inactive in that cell) |

**Raw counts** — how robust each estimate is

| column | meaning |
|---|---|
| `n_samples_A/B` | number of hopsy samples per cell (2000 unless NaN) |
| `n_nz_A/B` | number of non-zero samples per cell |
| `zero_frac_A/B` | fraction of samples at zero |
| `fisher_p_zeros` | Fisher exact test on the zero-fraction contingency table (defensible — proportions test, not sample-independent) |

**Location (effect size — primary)**

| column | meaning |
|---|---|
| `mean_nz_A, mean_nz_B` | raw non-zero means — interpretable in flux units |
| `median_nz_A, median_nz_B` | robust central tendency |
| **`log2FC_mean`** | **log2(mean_nz_A / mean_nz_B) — THE primary effect size** |
| `log2FC_median` | log2 of median ratio — robust variant |
| **`direction`** | **sign(log2FC_mean): +1 if A higher, −1 if B higher, 0 if NaN or tied** |

**Shape (effect size — secondary)**

| column | meaning |
|---|---|
| `W1_rxn` | Wasserstein-1 distance between the two non-zero flux distributions. Captures full-distribution difference including shape change at matched means. |
| `log_var_ratio` | log2(var_nz_A / var_nz_B). Values near 0 = similar spread. Values far from 0 = one cell has much wider/narrower flux range. |

**What this layer answers:**
- "By what factor did this reaction's flux change between cells?"
- "Is the shift mostly in the mean, or does the distribution shape also change?"
- "Which direction is the shift?"

---

## Layer 3 — Subsystem summary

**Files:** `subsystem_summary_F1.csv`, `subsystem_summary_F2.csv`

One row per (subsystem × contrast). Expected size: ~3,200 rows across both families.

**Columns:**

**Identifiers**

| column | meaning |
|---|---|
| `contrast, cell_A, cell_B, subsystem` | — |

**Inventory counts** (from Layer 1)

| column | meaning |
|---|---|
| `n_rxns_total_A, n_rxns_total_B` | reaction counts in each cell's extracted model, within this subsystem |
| `n_rxns_shared` | reactions in both |
| `n_rxns_A_only, n_rxns_B_only` | asymmetric presence |
| `n_rxns_union` | sum: all reactions ever present in either cell within this subsystem |

**Flux counts** (from Layer 2)

| column | meaning |
|---|---|
| `n_reps_shared` | representatives present in both cells |
| `n_reps_tested` | reps with `carry_both=True` and `loop_flagged=False` — the usable test set |
| `n_reps_loop_flagged, n_reps_not_carrying` | excluded by status flags |
| `n_reps_up, n_reps_down` | direction counts among tested reps |
| `n_reps_notable_up, n_reps_notable_down` | with `\|log2FC\| ≥ 1` additionally applied |
| `frac_notable` | (notable_up + notable_down) / n_reps_tested |
| `coherence` | \|up − down\| / (up + down). 1.0 = unanimous direction, 0.0 = evenly split |

**Magnitude summary**

| column | meaning |
|---|---|
| `mean_log2FC` | central tendency of per-rep log2FC |
| `median_log2FC` | robust variant |
| `max_abs_log2FC, max_abs_rxn` | biggest single mover, its rxn_id |
| `median_W1_rxn` | typical per-rep distribution-level distance |
| `median_log_var_ratio` | typical per-rep shape asymmetry |

**Categorical verdict**

| column | meaning |
|---|---|
| **`verdict`** | one of six: `structural_A_rich`, `structural_B_rich`, `flux_A_higher`, `flux_B_higher`, `mixed_flux`, `flat` |

### Verdict decision tree

Applied in order. The first that matches wins.

1. **`structural_A_rich`** — Cell A retains at least 2× as many reactions as B in this subsystem, AND the A-only count is at least 20 % of the union. Reading: "This pathway is better represented in A."
2. **`structural_B_rich`** — Opposite.
3. **`flux_A_higher`** — At least 5 reps tested, coherence ≥ 0.7, mean_log2FC positive and ≥ 0.5. Reading: "Shared reactions carry more flux in A."
4. **`flux_B_higher`** — Same with mean_log2FC negative.
5. **`mixed_flux`** — At least 5 reps tested, coherence < 0.5. Reading: "Shared reactions move, but half up and half down — no coherent direction."
6. **`flat`** — No clear structural or flux difference by the above criteria.

---

## How to find "important hits"

### For a reaction-level hit

```python
import pandas as pd
rxn = pd.read_csv("per_reaction_flux_F2.csv")

# Classical "notable" filter
hits = rxn.query('carry_both & not loop_flagged & abs(log2FC_mean) >= 1')
```

This filter keeps reactions that:
- Carry non-zero flux in both cells (comparison is defined).
- Are not loop-flagged (not an artefact of a thermodynamic loop).
- Shifted by at least a factor of 2 (log2FC_mean ≥ 1 or ≤ −1).

### For a subsystem-level hit

```python
roll = pd.read_csv("subsystem_summary_F2.csv")

# Shared-reaction flux shift
flux_hits = roll.query('verdict in ["flux_A_higher", "flux_B_higher"]')

# Structural asymmetry
struct_hits = roll.query('verdict in ["structural_A_rich", "structural_B_rich"]')

# Bidirectional (mixed) subsystems — worth looking at individually
mixed = roll.query('verdict == "mixed_flux"')
```

### For presence/absence findings

```python
inv = pd.read_csv("reaction_inventory_F1.csv")

# What reactions are in BL_96 but not LD_96?
only_in_BL = inv.query('contrast == "BL_96_vs_LD_96" and category == "A_only"')
```

---

## What the tables deliberately do NOT include

- **p-values on per-reaction tests.** Hopsy samples are not independent
  (uniform draws from a convex polytope), so classical tests on 2000
  samples per cell are massively overpowered. Rather than report
  inflated q-values, the tables rely on magnitude thresholds
  (`log2FC_mean`, `coherence`) that are not sensitive to sample
  independence. See `s270_per_reaction_deseq.py` for a direct
  demonstration of the overpowering (56–74 % of raw MW tests significant).
- **Within-sample temporal dynamics.** Hopsy sampling is a steady-state
  flux probe; this analysis compares equilibria, not dynamics.
- **Reaction coupling / module structure.** Pairwise correlations
  between representatives are not computed. The `subset_size` column
  reports extraction-time coupling only.
- **Sign flips in reversible reactions.** All comparisons are on \|v\|,
  so forward-vs-backward flux reversals are not resolved.

---

## Suggested chapter language

Instead of "p < 0.05," report:

- "**MAR04326 (glutathione synthetase)** carried mean non-zero flux 0.1 in
  BL_24 and 2.4 in D_24, a 24-fold increase (log2FC +4.6). The reaction
  is carrying in both cells and is not loop-flagged."
- "In **Nucleotide metabolism** between D_24 and LD_24, 60 of 85 tested
  representatives carried higher flux in LD_24 (coherence 0.71, mean
  log2FC +1.4). The subsystem verdict is `flux_B_higher`; the biggest
  single mover is MAR04675 at log2FC +5.2."
- "**Retinol metabolism** at D_72 vs LD_72 shows 13 of 13 shared
  representatives with higher flux in LD_72 (coherence 1.0, mean log2FC
  +3.4). The subsystem is also structurally symmetric (18 vs 17
  reactions retained)."
- "**Drug metabolism** at BL_96 vs LD_96 is `structural_A_rich` — BL_96
  retains 27 reactions vs 14 in LD_96. The 11 shared representatives
  tested also show 11/11 higher flux in BL_96 (log2FC +1.4)."

Each sentence stands on effect size + direction + consistency — no
p-value required.

---

## Methods paragraph for the chapter

> Flux distributions were sampled using hopsy (billiard walk, 4 chains ×
> 500 draws, thinning 300, CycleFreeFlux-corrected) under the measured
> metabolomic envelope and hpf-specific OCR cap applied via
> `apply_context_bounds`. Each reaction's flux comparison is
> characterized by three statistics: `log2FC_mean` (location shift),
> `log_var_ratio` (spread shift), and `W1_rxn` (Wasserstein distance,
> full-distribution). Direction is `sign(log2FC_mean)`. Subsystem
> summaries aggregate per-representative statistics into counts (up,
> down, notable), magnitude summaries (mean, median, max), and a
> categorical verdict distinguishing structural asymmetry from flux
> asymmetry. Because flux samples describe the extracted model's
> feasible region rather than biological replicates, we report effect
> sizes rather than classical frequentist p-values; the shape columns
> are included to catch cases where the distribution differs in ways
> the mean alone would miss. Reactions pruned by extraction from one
> cell are reported as presence/absence findings in the reaction
> inventory layer.

---

## File list

```
reaction_inventory_F1.csv       Layer 1, F1 contrasts
reaction_inventory_F2.csv       Layer 1, F2 contrasts
per_reaction_flux_F1.csv        Layer 2, F1 contrasts
per_reaction_flux_F2.csv        Layer 2, F2 contrasts
subsystem_summary_F1.csv        Layer 3, F1 contrasts
subsystem_summary_F2.csv        Layer 3, F2 contrasts
README.md                       this file
```
