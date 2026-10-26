# Briefing v2 — Zebrafish embryo metabolic model comparisons

Self-contained description of the dataset and analysis for a downstream read.
Primary question: **how does light regime (BL / D / LD) affect embryonic
metabolism?** Secondary question: **how does metabolism change across
development (24 → 48 → 72 → 96 → 120 hpf)?**

This replaces the previous briefing. Notable changes: per-cell
magnitude-normalized tables are now the primary analytical layer, and
several sanity checks (BL_96 outlier mechanism, peroxisomal switch,
melatonin-chain measurements) have been run and are summarized below.

---

## 1. Experimental setup

Zebrafish embryos reared under three light regimes from fertilization:

- **BL** — constant blue light
- **D** — constant darkness
- **LD** — standard 14 h light / 10 h dark cycle (the "normal" condition)

Samples taken at **24, 48, 72, 96, 120 hpf** → **15 cells** total
(`BL_24`, `BL_48`, …, `LD_120`). Each cell has:
- bulk transcriptome (RNA-seq → zFPKM)
- metabolomic exchange envelope (105 measured intake rates with fc values)
- measured whole-embryo oxygen uptake rate for its hpf

No chemical exposure, no perturbation beyond lighting.

## 2. Pipeline (short)

1. **Base model.** Zebrafish1 (Wang et al.) — 12,909 reactions, 8,344
   metabolites, 9 compartments.
2. **Context extraction.** rFASTCORMICS per cell, using zFPKM
   discretization. Produces 15 extracted models, each with ~4,800–5,600
   reactions.
3. **Flux envelope.** For each cell, `apply_context_bounds()` sets exchange
   bounds from measured metabolomics (lb = −10 × fc on 105 exchanges),
   OCR cap on `MAR09048` (97/173/227/273/313 at 24/48/72/96/120 hpf),
   §11.19-corrected vitamin bounds.
4. **Flux sampling.** hopsy billiard walk, 4 chains × 500 draws, thinning
   300, CycleFreeFlux-corrected to remove thermodynamic loops. Each cell
   samples biomass pinned at 0.9 × its own biomass_max (per-cell floor
   convention).
5. **Per-reaction comparison.** For every reaction present in both cells
   of a contrast, compute non-zero flux statistics, log2 fold change,
   direction, and distribution-shape metrics.
6. **Normalization (new, see §3).** Centered log2FC values added as a
   parallel analytical layer.
7. **Subsystem rollup + enrichment.** Hypergeometric test on counts of
   "notable" reactions per subsystem, BH-corrected per contrast.

## 3. What normalization did and why

An exploratory PCA showed PC1 correlated r = 0.95 with per-cell median
log flux. Since OCR rises 97 → 313 across hpf, per-cell flux magnitude
tracks hpf by construction. Separately, BL_96 has a notably lower
biomass_max (~8.5) than D_96 or LD_96 (~12.5–13), traced to the
isoleucine (MAR09039) uptake constraint — see §7.

**To separate magnitude from pathway-composition signals**, we centered
each cell's log flux vector by subtracting its own median across the
3,905 reactions present in all 15 cells. For pairwise comparisons this
is mathematically equivalent to shifting each contrast's log2FC
distribution by `(median_log_flux_A − median_log_flux_B)`.

Per-cell median shifts (centered − raw log2FC offsets):

```
24 hpf:  BL=-2.74  D=-2.59  LD=-2.59   small spread (0.15)
48 hpf:  BL=-2.39  D=-2.35  LD=-2.48   small spread (0.13)
72 hpf:  BL=-1.97  D=-1.98  LD=-1.90   small spread (0.08)
96 hpf:  BL=-1.11  D=-1.44  LD=-1.75   BIG spread (0.64) — biomass asymmetry
120 hpf: BL=-1.30  D=-1.41  LD=-1.50   small spread (0.21)
```

At 24, 48, 72, 120 hpf centering is a small adjustment. At 96 hpf and
across F2 (within-regime adjacent-hpf developmental contrasts) centering
is substantial because it removes the OCR-driven magnitude baseline.

**Both raw and centered tables are provided.** Centered is the recommended
primary read for pathway-composition questions; raw is the backup when
total metabolic activity (magnitude + composition) is the biology of
interest.

## 4. Direction conventions (do NOT mix these up)

- **F1** (same hpf, cross-condition): `A` is the alphabetically earlier
  condition (BL < D < LD). `direction = +1` means A higher than B.
  `log2FC_centered = +X` means A higher by factor 2^X.
- **F2** (same regime, adjacent hpf): `A` is the earlier hpf.
  `direction = +1` means A higher (= flux falls across development).
  `direction = −1` means B higher (= flux rises across development).

Positive log2FC = A > B. Negative = B > A. Same for `log2FC_mean` (raw)
and `log2FC_centered` (normalized).

## 5. Files to use

Primary (centered): `grace_v4/results/chapter_assets/comparison_tables/`

```
subsystem_summary_F1_centered.csv       light-regime subsystems with enrichment_q
subsystem_summary_F2_centered.csv       developmental subsystems with enrichment_q
per_reaction_flux_F1_centered.csv       light-regime per-reaction
per_reaction_flux_F2_centered.csv       developmental per-reaction

reaction_inventory_F1.csv               presence/absence (unchanged by normalization)
reaction_inventory_F2.csv               same

pca_cells.csv + pca_cells.pdf           raw PCA (magnitude + composition)
pca_cells_centered.csv + pca_cells_centered.pdf  centered PCA (composition only)
pca_cells_comparison.pdf                side-by-side raw vs centered
```

Raw backups (same schema, uncorrected):

```
subsystem_summary_F1.csv
subsystem_summary_F2.csv
per_reaction_flux_F1.csv
per_reaction_flux_F2.csv
```

Column schema for the centered subsystem summary:

```
contrast, cell_A, cell_B, subsystem
n_reps_shared, n_reps_tested, n_reps_up, n_reps_down
n_reps_notable_up, n_reps_notable_down  (|log2FC_centered| >= 1)
frac_notable, coherence                 (|up-down| / (up+down))
mean_log2FC_centered                    PRIMARY effect size
median_log2FC_centered, max_abs_log2FC_centered, max_abs_rxn
verdict                                 flux_A_higher / flux_B_higher / mixed_flux / flat
odds_ratio, enrichment_p, enrichment_q  Fisher exact + BH-corrected q, per contrast
```

Per-reaction table has the same identifiers plus full statistics per
reaction: raw and centered log2FC, direction, log_var_ratio,
Wasserstein distance, Fisher p on zero fractions, carry_both and
loop_flagged flags.

## 6. Statistical framing (important)

Hopsy's 2,000 samples per reaction per cell describe the extracted
model's **feasible flux region**, not biological replicates. Samples
are correlated (uniform draws from a convex polytope), so classical
Mann-Whitney / t / DESeq-style p-values on per-reaction samples are
massively overpowered and are not reported.

The two defensible significance measures:

- **`enrichment_q`** at the subsystem level — Fisher exact test on
  counts of notable reactions per subsystem (classical, well-calibrated,
  tests proportions not samples).
- **`fisher_p_zeros`** at the per-reaction level — Fisher exact on
  zero-fraction asymmetry between cells.

Everything else is magnitude-based: `|log2FC_centered| >= 1` for notable
reactions, `coherence >= 0.7` for coherent subsystem direction,
`n_reps_tested >= 5` to avoid tiny-subsystem noise.

## 7. Sanity checks already run

Previous analysis of these tables identified four candidate light-regime
findings and flagged BL_96 as a possible artefact. The following direct
checks were then run:

1. **96 hpf measured envelope sum.** BL_96 envelope `Σ|lb| = 1236`;
   D_96 = 1308; LD_96 = 1193. Spread 9.7 %, in line with 24/48/72/120
   hpf (4.8 – 7.8 % spread). **BL_96 is NOT in a loose envelope**; the
   earlier "polytope-size artefact" hypothesis is wrong.

2. **BL_96 biomass_max binding analysis.** BL_96 reaches biomass_max =
   8.68 bound by MAR09039 (Exchange of isoleucine, fc 0.818, lb −8.18).
   D_96 and LD_96 reach biomax 13.83 / 14.22 bound by MAR09152
   (α-tocotrienol at the §11.19 numerical floor −0.001). LD_96 has
   nearly identical isoleucine fc to BL_96 (0.807 vs 0.818) but reaches
   biomax 14.22 — the extracted model retains internal isoleucine-
   supplying reactions that BL_96's does not. **BL_96's low biomax is
   a topology-level finding (retained reaction set difference), not a
   measurement artefact.**

3. **Peroxisomal β-oxidation switch at 72 hpf.** The 0 → 47 reaction
   switch LD leads over BL and D is not driven by a single gene
   crossing a threshold. `acox1`, `fasn`, `acaca` show modest LD-higher
   zFPKM (deltas 0.19 – 0.34) but none crosses an rFASTCORMICS
   inclusion threshold on its own. Likely multi-gene + FASTCORE
   connectivity-driven.

4. **Melatonin chain exchanges.** None of the five candidate melatonin-
   chain exchanges (5-HIAA, Indole-3-acetate, Indole, Melatonin,
   6-hydroxymelatonin) is in the 105-metabolite measured envelope.
   Neither is tryptophan. **The 72-hpf BL melatonin depression is a
   pure model inference**, not a metabolomics-input artefact.

## 8. Framing for interpretation

- `enrichment_q < 0.05` is the primary significance gate at the
  subsystem level.
- `|mean_log2FC_centered| >= 1` and `coherence >= 0.7` are the magnitude
  and directional-consistency filters for calling a subsystem a
  pathway-level finding.
- `n_reps_tested >= 5` ensures the subsystem has enough carrying
  representatives for the comparison to be meaningful.
- Verdict categories (`flux_A_higher`, `flux_B_higher`, `mixed_flux`,
  `flat`) encode these filters. `mixed_flux` and `flat` subsystems with
  low `enrichment_q` still deserve attention — they can indicate
  bidirectional or shape-only shifts.
- Reaction-level findings with `loop_flagged = True` should be
  interpreted with caution (possible thermodynamic cycle artefact).

## 9. What to look for

### Primary (light regime, F1 tables)

- Subsystems with `enrichment_q < 0.05` AND
  `verdict in ("flux_A_higher","flux_B_higher")` at specific hpf.
- Cross-hpf patterns: does one condition consistently differ at one hpf?
- Known robust findings from previous analysis of these data:
  - 48 hpf BL fatty-acid activation (BL retains lipid activation that
    D and LD transiently pause between 24 and 48 hpf)
  - 72 hpf constant-light serotonin/melatonin depression (BL < D, LD)
  - 72 hpf LD peroxisomal β-oxidation lead
  - 120 hpf LD-specific sterol/bile acid step (LD > BL, D)

### Secondary (development, F2 tables)

- Subsystems rising in all three regimes (shared developmental program):
  check for `flux_B_higher` verdicts at the same hpf transition across
  BL, D, LD. The 48 → 72 hpf lipid activation program is the strongest
  example (Leukotriene, Arachidonic acid, Fatty acid activation,
  Fatty acid oxidation all rise in all three).
- Lead/lag effects: does one regime complete a transition one stage
  earlier?

### What NOT to over-interpret

- `MW_q` columns are not in these tables (deliberately removed —
  overpowered by sample autocorrelation).
- Subsystems with `n_reps_tested < 5` can show large log2FC but are
  noise-prone.
- `Transport reactions` as a subsystem dominates many enrichments by
  virtue of size; real but rarely biologically specific.
- `Retinol metabolism` and `Drug metabolism` findings can be driven
  by unknown-zone amplification (near-floor gene expression where
  rFASTCORMICS inclusion is connectivity-driven, not transcriptional);
  see §11.21.5 of the parent chapter.

## 10. Suggested analysis prompt for a fresh pass

> You have per-reaction and per-subsystem flux comparison tables for 15
> zebrafish embryo samples (3 light regimes × 5 developmental stages)
> from genome-scale metabolic models sampled via hopsy. Primary and
> centered log2FC values are provided; use the `_centered.csv` files
> as the primary analytical layer. Primary question: which metabolic
> subsystems differ between light regimes, and are the differences
> stage-specific? Secondary question: what metabolic programs unfold
> across development, and do the three regimes follow the same schedule?
>
> The primary significance column is `enrichment_q` (BH-corrected
> hypergeometric on notable reactions). The primary effect size is
> `mean_log2FC_centered` with `coherence` as directional consistency.
> `verdict` gives a categorical read.
>
> Convention: F1 `A` = alphabetically earlier condition (BL<D<LD);
> F2 `A` = earlier hpf; `+log2FC` means A higher.
>
> Four light-regime findings were previously identified (see §9).
> Confirm or refine them with the centered tables, and flag any
> additional patterns the centered analysis surfaces. Do not re-compute
> per-reaction p-values — hopsy samples are not biological replicates
> and classical tests are overpowered.
