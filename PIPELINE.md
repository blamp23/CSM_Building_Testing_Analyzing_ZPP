# v7/v8 pipeline — linear flow

End-to-end sequence from measured data to frozen findings. Each stage lists
inputs, outputs, and the amendment it's governed by. Nothing in this doc
re-defines prior stages — read `FAMILY_PREREG.md` for the amendment trail.

**2026-10-26 reframe:** Stages 6–10 (W1 contrasts → tier export → figures) are
**superseded** for the chapter. The chapter uses a per-reaction DESeq-style
framework instead; see Stage 12 and `v8_analysis/`. Stages 0–5 (model + sample
+ loop + pool + subset) are unchanged and feed both pipelines. Stage 11 (F2
matched-floor resample) also feeds the v8 pipeline. The `v7_final/` deliverable
is retained as a frozen artifact of the earlier Tier A/B/C analysis.

## Snapshot of outputs (frozen 2026-10-13)

`v7_final/` contains the frozen deliverables; `v7_final/addenda/` holds the
Tier 0 add-ons (MAR10693 vent check, robustness LPs). See
`v7_final/FREEZE.md` for md5s and amendment list.

Headline funnel: **1,687 rows in family → 146 at q ≤ 0.10 → Tier A 50
(F1=26, F2=24), Tier B 63, Tier C 33.** Positive controls validated
(Steroid KO Fisher p = 8e-11; KD2 W1 = 2.29, snr ≈ 38).

---

## Stage 0 — base model

`base_model/` + `extraction/`. Wang zebrafish GEM → baked_ocr_{hpf}_v7.mat
(5 per-hpf parents, flux bounds + ε = 10⁻³ cofactor pool).
→ rFASTCORMICS extraction → `v7_models/trans_rfastcormics_{cell}.mat`
(15 context-specific child models). Frozen at tag `v7-frozen`.

## Stage 1 — envelope + OCR + ledger

`functions/apply_context_bounds.py`. Layers onto any model's lb/ub:
- 105 measured exchange fc multipliers per (cond, hpf) from
  `data/phase2_exchange_bounds_long.csv` (§1)
- OCR cap on MAR09048 per hpf (§1)
- Medium opens, ledger closes, Part 1 secretion ub=1000, Part 2 lb=0 (§1)
- medium_uptake demand bounds (§11.19, added 2026-10-12)

## Stage 2 — polytope + rounding

PolyRound: simplify → affine transform → MVE round. Gurobi-required (local
only; not available on Grace venv). Writes
`{cell}_polytope_cache.pkl` + `{cell}_polytope_meta.pkl`.

## Stage 3 — sampling

`v7_analysis/matched_floor/hopsy_sample_direct.py`. hopsy billiard-walk,
4 chains × 500 samples × thinning 300 per cell. Explicit starting-point
Chebyshev + feasibility shrink (patched 2026-10-12 for D_72). Writes
`flux_samples_{cell}_ch{N}.csv`.

(Original 15-cell sample set runs under hopsy_v7_local/; §11.20 matched-
floor resample for 9 cells writes to hopsy_v7_matched_f2/ via Grace.)

## Stage 4 — loop removal

`v7_analysis/loop_removal/cff_gurobi.py`. CycleFreeFlux (Desouki 2015) via
Gurobi; warm-starts between consecutive samples. Writes loopless_{cell}_ch{N}.csv.

## Stage 5 — pooling + enzyme subsets

Pool 4 chains per cell → `pooled_v7/{cell}_pooled.csv` (rxn × 2000).
Enzyme-subset collapse → `subsets_v7/{cell}_subset_membership.csv`
(reaction → representative + subsystem + GPR flag).

## Stage 6 — two-tier W1 contrasts

`v7_analysis/statistics/two_tier_fast_mean_dropzero.py`. 42 contrasts
(27 primary + 15 split-chain neg-control) run via
`v7_analysis/launchers/run_contrasts_mean_dropzero.sh`. Writes
`stats_v7_mean_dropzero/{A}_vs_{B}.csv`.

Method of record = mean normaliser + drop zeros (§11.14). Median variant
runs in parallel for concordance (§11.12).

## Stage 7 — family findings

`v7_analysis/statistics/findings_v7_mean_dropzero.py`. BH-FDR within
family, SNR vs split-chain noise floor (§11.13), envelope-tracking flag
(§11.13), floor-suspect flag (§11.10). Writes `findings_v7_mean_dropzero/`
with `findings_family1.csv`, `findings_family2.csv`, `negcontrol_per_subsystem.csv`.

## Stage 8 — tier export (headline)

`v7_analysis/statistics/export_tier_tables.py` (= s85). Applies:
- median-concordance gate (§11.14)
- matched-floor resolution — F1 via BL_96 (§11.10), F2 via six-contrast
  matched-floor resample (§11.20)
- Tier 0 floor-exposed flag (§11.20)
- replication demotion
- F2 matched-floor demotion to Tier C `floor_artifact` when p > 0.05 (§11.20)
- `pre_correction_retinoid` provenance flag (§11.19 / §11.17)

Writes `tierA_final.csv`, `tierB_final.csv`, `tierC_final.csv`,
`all_rows_in_family.csv`, `tierA_by_subsystem.csv`.

## Stage 9 — Tier 0 structural analysis

`v7_analysis/tier0/`:
- **tier0_2x2.py (§11.17)** — parent vs child FBA per cell + 2×2 topology ×
  envelope swap per F1 pair. Writes `tier0_parent_child_v7.csv` (15 cells),
  `tier0_2x2_v7.csv` (15 pairs × 4 corners + decomposition).
- **bl96_drilldown.py (§11.17)** — greedy single-reaction add-back from
  BL_96 baseline toward LD_96 ceiling. Writes `bl96_lacking_reactions.csv`,
  `bl96_greedy_addback.csv`, `bl96_greedy_gene_zfpkm.csv`.
- **mar10693_vent_check.py (§11.17 footnote, 2026-10-13)** — 15 LPs
  closing MAR10693 ub. Confirmed dual degeneracy (0 % Δ in every cell).
- **robustness.py (§11.17, 2026-10-13)** — 60 LPs: valine / histidine /
  isoleucine / glucose uptake relaxed by 10 % per cell. Writes
  `tier0_robustness.csv`. Nitrogen-bound, histidine-dominated at hpf
  48–72, isoleucine-bound at BL_96, glucose carbon control passes.

## Stage 10 — figures

`v7_analysis/figures/`:
- **figs_easy.py** — fig1 (design + funnel), fig2 (F2 developmental
  heatmap), fig3 (F1 photoperiod heatmap with ◆/▲/▼ glyphs), fig5
  (growth capacity attribution), figS4 (R̂/ESS violins), figS5 (mean vs
  median concordance). Each with `_data.csv`.
- **condition_axis_curves.py** — per-(subsystem × hpf) CDF document, 3
  sub-panels per page (BL-vs-D, BL-vs-LD, D-vs-LD); emits `tierA_table.pdf`
  + markdown reader.

## Stage 11 — matched-floor F2 resample (§11.20)

`v7_analysis/matched_floor/run_pipeline.sh`. 9 cells × 4 chains sampled on
Grace (polytope caches built locally; sampler runs CPU-only billiard walk).
CFF + pool + 6 F2 contrasts locally. Writes
`stats_v7_matched_f2/{A}_vs_{B}.csv` which `export_tier_tables.py` reads.

Verdict: 10 of 22 flagged F2 Tier A rows survive (lipid/eicosanoid family),
12 demote to Tier C (`floor_artifact`; amino-acid 24→48 rows).

## Stage 12 — chapter analytical reframe (2026-10-26, `v8_analysis/`)

**Rationale.** Hopsy samples are correlated draws from one polytope, not
biological replicates. Per-sample W1 contrasts (Stage 6) and per-reaction
Mann-Whitney p-values are both overpowered on 2,000-sample pools. The
reframe drops both and reports exactly one significance measure: a
hypergeometric enrichment q-value on the count of notable reactions per
subsystem, BH-corrected within contrast. Combined with per-cell median
centering of log flux (removing the OCR cap × hpf and biomass_max scaling
baselines), the framework asks whether a subsystem's reactions shift
coherently in direction + magnitude relative to the rest of the model.

**Inputs:** `pooled_v7/{cell}_pooled.csv` (Stage 5) and
`subsets_v7/{cell}_subset_membership.csv` (Stage 5), plus the v7_models
.mat files (Stage 0) for the capacity LPs. The matched-floor F2 pool
(Stage 11) is implicit in the F2 contrasts.

**Scripts (linear order in `v8_analysis/scripts/`):**

| # | Script | Produces |
|---|---|---|
| 01 | `01_substrate_yields.py` | ATP yield per substrate × cell (Table S7) |
| 02 | `02_atp_ceiling_physio.py` | ATP ceiling under physio OCR + envelope (§3.2) |
| 03 | `03_atp_ceiling_glucose.py` | ATP ceiling glucose-only control (§3.2) |
| 04 | `04_comparison_tables.py` | Three-layer F1/F2 tables (raw log2FC) |
| 05 | `05_pca_raw.py` | PCA + initial hypergeometric enrichment |
| 06 | `06_pca_centered.py` | Per-cell-centered PCA + comparison figure |
| 07 | `07_normalize_centered.py` | Centered F1/F2 tables + enrichment (§3.1 primary) |
| 08 | `08_capacity_figure.py` | Xenobiotic clearance LPs + 2×2 figure (§3.3) |
| 09 | `09_density_plots.py` | Developmental + light-regime density ridges (§3.1) |

Scripts 01, 02, 03, 08 need Gurobi. 04–07, 09 are pandas-only.

**Path configuration.** All 9 scripts import from `v8_analysis/scripts/config.py`.
Repo-internal paths resolve from `__file__`; large external data (v7_models,
pooled_v7, subsets_v7) resolve from `$CSM_DATA_ROOT`, defaulting to the
author's OneDrive working location. Set the env var to point at a local
mirror.

**Outputs (committed to git under `v8_analysis/chapter_assets/`):**
`comparison_tables/subsystem_summary_F{1,2}_centered.csv` (§3.1 primary),
`per_reaction_flux_F{1,2}_centered.csv`, `reaction_inventory_F{1,2}.csv`,
`pca_cells*.pdf/csv`, `tables/table_{1..6}.csv`, `table_S{1..15}.csv`,
`xenobiotic_capacity.csv`, `figures/metabolic_capacity_developmental.pdf`.
Density figures (`figures/subsystem_density/`) are regenerable by script 09
and git-ignored.

**Deprecated by Stage 12:**
- Stage 6 outputs (`stats_v7_mean_dropzero/`) — W1 contrasts are no longer
  reported.
- Stage 7 outputs (`findings_v7_mean_dropzero/`) — family findings are not
  cited in the chapter.
- Stage 8 outputs (`tierA_final.csv`, `tierB_final.csv`, `tierC_final.csv`)
  — the tier framework is retained in `v7_final/` as a frozen artifact of
  the earlier analysis, but the chapter does not cite the three-tier split.
- Stage 10 figures — superseded by density ridges (script 09) and the 2×2
  capacity figure (script 08).

---

## Reproducing the whole thing from a clean checkout

```
# Prerequisites: Gurobi, Python 3.13, hopsy, PolyRound,
# scipy, pandas, matplotlib, gurobipy, h5py.

# 0-5: build models + sample + CFF + pool (see v7_analysis/README.md stages)
# 6-7: contrasts + findings
bash v7_analysis/launchers/run_contrasts_mean_dropzero.sh
python v7_analysis/statistics/findings_v7_mean_dropzero.py

# 8: tier export
python v7_analysis/statistics/export_tier_tables.py

# 9: Tier 0
python v7_analysis/tier0/tier0_2x2.py
python v7_analysis/tier0/bl96_drilldown.py
python v7_analysis/tier0/mar10693_vent_check.py
python v7_analysis/tier0/robustness.py

# 10: figures
python v7_analysis/figures/figs_easy.py
python v7_analysis/figures/condition_axis_curves.py

# 11: matched-floor F2 resample (requires Grace sampling first; see §11.20)
bash v7_analysis/matched_floor/run_pipeline.sh
python v7_analysis/statistics/export_tier_tables.py   # re-run to pick up matched-floor
```

Outputs land in `grace_v4/results/` (sampling/stats/findings) and
`v7_final/` (frozen snapshot). The `v7_final/` folder is immutable;
`FREEZE.md` carries md5s.
