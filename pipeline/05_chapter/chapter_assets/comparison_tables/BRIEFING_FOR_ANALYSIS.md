# Briefing — Zebrafish embryo metabolic model comparisons

A self-contained description of the dataset for a downstream analysis.
Primary question: **how does light regime (BL / D / LD) affect embryonic
metabolism?** Secondary question: **how does metabolism change across
development (24 → 48 → 72 → 96 → 120 hpf)?**

---

## 1. Experimental setup

Zebrafish embryos were reared under three light regimes from fertilization:

- **BL** — constant blue light
- **D** — constant darkness
- **LD** — standard 14 h light / 10 h dark cycle (the "normal" condition)

Samples were taken at five developmental stages: **24, 48, 72, 96, 120 hpf**.
That gives **15 cell samples** (`BL_24`, `BL_48`, …, `LD_120`), each with:
- a bulk transcriptome (RNA-seq → zFPKM)
- a metabolomic exchange envelope (105 measured intake rates with fc values)
- a measured whole-embryo oxygen uptake rate for its hpf

No chemical exposure, no perturbation beyond lighting. These are baseline
embryos reared to each timepoint.

## 2. From data to metabolic model

**Base model.** Zebrafish1 (Wang et al.) — 12,909 reactions, 8,344
metabolites, 9 compartments. A whole-genome scale metabolic reconstruction.

**Per-cell extraction (rFASTCORMICS).** For each of the 15 cells:

1. zFPKM discretized to {expressed, unknown, off} at two thresholds.
2. FASTCORE builds a context-specific model containing only reactions
   either (a) with expressed GPR, (b) connectivity-required with
   unknown GPR, or (c) necessary for a biomass objective.
3. Result: 15 extracted models at `reviewer_packet_v4_slim/v7_models/
   trans_rfastcormics_{cell}.mat`, each with ~4,800-5,600 reactions
   depending on cell.

**Flux envelope.** For each cell, `apply_context_bounds()` sets the
exchange bounds:
- Measured metabolite intakes: `lb = -10 × fc` on 105 measured exchanges
- Hpf-specific OCR cap on `MAR09048` (O₂): 97/173/227/273/313 µmol O₂/g
  DW/hr at 24/48/72/96/120 hpf
- Medium (bulk inorganics) open
- Four vitamin exchanges opened with per-cell demand bounds (§11.19)
- All other exchanges closed to uptake

## 3. From model to flux samples

Each extracted model has ~5,000 reactions and infinite feasible flux
states consistent with the bounds. We sample them uniformly:

**hopsy billiard walk** on the simplified+rounded polytope. Four
independent chains × 500 draws × thinning 300 = 2,000 pooled samples
per reaction per cell. Each cell is a matrix: **(reactions × 2000)**.

**CycleFreeFlux post-correction (Desouki 2015).** Removes thermodynamic
loops in each sample while preserving exchange fluxes and biomass. The
samples on disk are post-CFF.

Output per cell: `grace_v4/results/pooled_v7/{cell}_pooled.csv`

**Representative collapse.** Many reactions are mathematically coupled
(share a GPR, lie in series — their fluxes move in lockstep). Each
coupled group is represented by one "representative" reaction. ~1,500
representatives per cell. The membership map is at
`grace_v4/results/subsets_v7/{cell}_subset_membership.csv`.

## 4. Statistical framing (important)

Hopsy's 2,000 samples per reaction per cell describe the **feasible
flux region of the extracted model**, not biological replicates. The
samples are correlated (uniform draws from a convex polytope), so
**classical p-values from Mann-Whitney, t-tests, or DESeq don't apply**
— they would be massively overpowered.

The analysis therefore uses:
- **log2 fold change on non-zero flux means** as the primary effect
  size per reaction
- **direction** as the sign of the fold change
- **coherence** at the subsystem level — fraction of reactions moving
  in the same direction
- **Fisher exact test** on zero fractions (defensible — proportions
  test, no sample-independence assumption)
- **Hypergeometric enrichment** on counts of notable reactions per
  subsystem (defensible — proportions test; gives real q-values)

The one q-value you can trust in these tables is `enrichment_q`.
Everything else is magnitude-based.

## 5. Direction conventions (critical — do NOT mix these up)

- **Family 1 (F1) — same hpf, cross-condition.** `A` is the
  alphabetically earlier condition (BL < D < LD). Example: `BL_24_vs_D_24`
  has `A = BL_24`, `B = D_24`. `direction = +1` means BL higher than D.
- **Family 2 (F2) — same regime, adjacent hpf.** `A` is the earlier
  hpf. Example: `D_48_vs_D_72` has `A = D_48`, `B = D_72`.
  `direction = +1` means flux falls from 48 → 72. `direction = −1`
  means flux rises from 48 → 72.

Positive log2FC = A > B. Negative = B > A.

For developmental questions use F2. For light-regime questions use F1.

## 6. The three table layers

Files in `grace_v4/results/chapter_assets/comparison_tables/`:

### Layer 1 — Reaction inventory (presence / absence)

Which reactions does each cell's extracted model retain?

- `reaction_inventory_F1.csv` and `reaction_inventory_F2.csv`
- One row per (physical reaction × contrast)
- Columns: `rxn_id, subsystem, in_A, in_B, category`
  (`category` ∈ {`both`, `A_only`, `B_only`})

Use this to answer "what reactions does BL_96 have that LD_96 doesn't?"

### Layer 2 — Per-reaction flux comparison

For reactions in both cells, how their flux distributions differ.

- `per_reaction_flux_F1.csv` and `per_reaction_flux_F2.csv`
- One row per (shared representative × contrast)
- Key columns:
  - **`log2FC_mean`** — primary effect size (log2 of mean non-zero
    flux ratio A/B)
  - **`direction`** — +1 (A higher), −1 (B higher), 0 (undefined/tied)
  - **`log2FC_median`** — robust variant
  - `W1_rxn` — Wasserstein distance between non-zero flux distributions
    (catches shape-only differences at matched means)
  - `log_var_ratio` — log2 of variance ratio (spread asymmetry)
  - `mean_nz_A`, `mean_nz_B`, `median_nz_A`, `median_nz_B` — raw units
  - `zero_frac_A/B`, `fisher_p_zeros` — presence/absence asymmetry
  - `carry_both` — bool, both cells show meaningful non-zero flux
  - `loop_flagged` — bool, max \|v\| > 100 in either cell
    (possible thermodynamic loop artefact)
  - `carrying_asymmetric` — bool, one cell is all-zero

**Primary notable filter:** `carry_both AND not loop_flagged AND |log2FC_mean| >= 1`

### Layer 3 — Subsystem summary

Aggregates layers 1 and 2 into one row per (subsystem, contrast).

- `subsystem_summary_F1.csv` and `subsystem_summary_F2.csv`
- Key columns:
  - **Inventory counts:** `n_rxns_total_A/B, n_rxns_shared,
    n_rxns_A_only, n_rxns_B_only, n_rxns_union`
  - **Flux counts:** `n_reps_tested, n_reps_up, n_reps_down,
    n_reps_notable_up, n_reps_notable_down`
  - **`coherence`** — \|up − down\| / (up + down). 1.0 = unanimous
    direction, 0.0 = evenly split.
  - **Magnitudes:** `mean_log2FC, median_log2FC, max_abs_log2FC,
    max_abs_rxn`
  - **`median_W1_rxn, median_log_var_ratio`** — subsystem-level
    distribution-shape summaries
  - **`verdict`** — categorical: `flux_A_higher`, `flux_B_higher`,
    `structural_A_rich`, `structural_B_rich`, `mixed_flux`, or `flat`
  - **`enrichment_q`** — BH-corrected Fisher exact q, tests whether
    notable reactions are over-represented in this subsystem within
    this contrast. This is the one defensible q-value.

### Layer 4 — PCA on cells

- `pca_cells.csv` — 15 cells × (PC1, PC2, PC3, condition, hpf)
- `pca_loadings.csv` — reactions × PC loadings with subsystem labels
- `pca_explained_variance.csv` — PC1 explains 45.5 %, PC2 20.9 %,
  PC3 11.8 % (first 3 = 78 %)
- `pca_cells.pdf` / `.png` — 2-panel scatter

**Known pattern:** PC1 is almost perfectly ordered by hpf
(developmental axis). Light regime is a second-order perturbation.
BL_96 is the biggest outlier.

## 7. What biologically sensible questions look like

### Primary (light regime)

Use **F1** files (same hpf, cross-condition).

- **At each hpf, which subsystems show large coherent flux differences
  between BL/D/LD?** Filter `subsystem_summary_F1.csv` on
  `verdict in ("flux_A_higher","flux_B_higher") AND enrichment_q < 0.05`.
- **At each hpf, which reactions are retained in one condition but
  absent in another?** Filter `reaction_inventory_F1.csv` on
  `category in ("A_only","B_only")`, group by subsystem.
- **Which hpf shows the strongest light-regime signal?** Count hits
  per hpf in each F1 contrast.
- **Known patterns worth checking:**
  - BL_48 vs D_48 Fatty acid activation (cytosolic): BL much higher.
  - BL_96 is a general outlier (lowest biomass_max, many flux-higher
    signatures vs LD_96).
  - 120 hpf: Bile acid biosynthesis + Steroid metabolism + Leukotriene
    metabolism tend to favor LD or D over BL.
  - 120 hpf is when circadian clock genes become fully expressed —
    look for metabolic signatures of this.

### Secondary (development)

Use **F2** files (same regime, adjacent-hpf transitions).

- **Which pathways are activated across development, regardless of
  light?** Filter `subsystem_summary_F2.csv` on
  `verdict == "flux_B_higher" AND enrichment_q < 0.05 AND coherence >= 0.7`,
  then look for subsystems that show up in all three regimes at the
  same hpf transition.
- **Which pathways are downregulated across development?** Same with
  `verdict == "flux_A_higher"` (A is earlier hpf, so A higher = flux
  falls over development).
- **Known patterns worth confirming:**
  - Fatty acid oxidation, Arachidonic acid metabolism, Leukotriene
    metabolism, and Fatty acid activation all rise across 48 → 72
    across all three regimes (the "lipid activation program").
  - Steroid metabolism, Bile acid biosynthesis rise at 96 → 120.
  - Nucleotide / pyrimidine metabolism rises generally with development.

### Cross-cutting questions

- **Do light-regime hits cluster at specific hpf?** Group F1 significant
  hits by hpf.
- **Does one condition "lead" or "lag" developmental transitions?**
  Compare F2 transitions across the three regimes at matched hpf —
  e.g., does the 48 → 72 lipid program complete earlier in LD than D?
- **Which subsystems have `structural_A_rich` / `structural_B_rich`
  verdicts?** These are presence/absence findings — the model retains
  more reactions in one cell than another, independent of flux. Often
  reflects transcriptional on/off patterns.

## 8. Things to NOT over-interpret

- **MW_q columns don't exist in these tables — don't invent them.**
  Use `enrichment_q` as the only q-value.
- **`loop_flagged = True` reactions** should be filtered out of
  interpretations unless the subsystem has unusually high flux across
  the board (TCA cycle, oxphos, etc. — those are expected high-magnitude).
- **The log2FC on reactions where one cell has very few non-zero samples
  (e.g. n_nz_A < 20)** can be unstable. Check `n_nz_A` and `n_nz_B` are
  both substantial before quoting a number.
- **Transport reactions as a subsystem** typically dominate significance
  because it's a huge subsystem — real but not usually biologically
  specific.
- **Subsystems with unknown-zone genes** (gss, gfus.1, some drug/retinol
  metabolism) can show false flux differences driven by rFASTCORMICS's
  include/exclude decision at near-floor-expressed genes, not by actual
  transcriptional differences. See `grace_v4/results/v7_susceptibility/
  s195_amendment_11p22.md` or ask about the "unknown-zone amplification
  confound" if a finding is in Drug metabolism or Retinol metabolism.
- **Hopsy samples are not biological replicates.** Treat log2FC values
  as the model's prediction for the fold-change, not as a frequentist
  claim with uncertainty.

## 9. Files to hand over

Minimum useful set (should all be in `grace_v4/results/chapter_assets/
comparison_tables/`):

```
BRIEFING_FOR_ANALYSIS.md          this file
README.md                         table-level column dictionary

per_reaction_flux_F1.csv          Layer 2, light-regime contrasts
per_reaction_flux_F2.csv          Layer 2, developmental contrasts
subsystem_summary_F1.csv          Layer 3, light-regime subsystems
subsystem_summary_F2.csv          Layer 3, developmental subsystems
reaction_inventory_F1.csv         Layer 1, presence/absence (optional)
reaction_inventory_F2.csv         Layer 1, presence/absence (optional)

pca_cells.csv                     cell PCA coordinates
pca_loadings.csv                  reaction loadings
pca_cells.pdf                     PCA figure
pca_explained_variance.csv        variance per PC
```

If the downstream analysis also wants context on reactions:
`reviewer_packet_v4_slim/baked_ocr_24_v7.mat` (base model with reaction
names + metabolite names + GPR). Not strictly needed if the AI can
look up MAR IDs in Human-GEM / Wang online.

## 10. Suggested analysis prompt

Something like:

> You have per-reaction and per-subsystem flux comparison tables for 15
> zebrafish embryo samples (3 light regimes × 5 developmental stages)
> from genome-scale metabolic models sampled via hopsy. Primary question:
> what metabolic subsystems differ between light regimes, and are the
> differences stage-specific? Secondary question: what metabolic programs
> unfold across development, and do the three regimes follow the same
> schedule?
>
> Use `subsystem_summary_F1.csv` for cross-condition questions and
> `subsystem_summary_F2.csv` for developmental questions. The primary
> significance column is `enrichment_q` (BH-corrected hypergeometric
> on notable reactions). The primary effect size is `mean_log2FC` with
> `coherence` as the directional consistency. `verdict` gives a
> categorical read: `flux_A_higher`, `flux_B_higher`, `structural_A_rich`,
> `structural_B_rich`, `mixed_flux`, or `flat`.
>
> Convention: in F1, `A = BL` or earlier of (BL,D,LD); in F2, `A =` earlier
> hpf. `+log2FC` means A higher. `direction = +1` means A higher.
>
> Flag any findings that make biological sense (consistent with known
> embryo development, circadian metabolism, lipid programs) and flag any
> that look weird or require follow-up. Don't trust p-values other than
> `enrichment_q` and `fisher_p_zeros` — the per-reaction MW_q column
> doesn't exist in these tables and shouldn't be re-computed on hopsy
> samples because they're not biological replicates.
