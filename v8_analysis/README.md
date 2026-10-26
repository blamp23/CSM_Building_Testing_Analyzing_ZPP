# v8_analysis — chapter analytical pipeline

Stage 12 of PIPELINE.md. Nine linear scripts that take frozen Stage 0–5 outputs
(context-specific models + hopsy flux samples) and produce the chapter's
findings tables, PCA, metabolic capacity LPs, and density figures.

## Why "v8" when the repo is v7-frozen

The v7 pipeline (Stages 6–10 in PIPELINE.md) summarized sampled flux via W1
contrasts and a Tier A/B/C exportable table. The chapter reframe (2026-10-26)
abandons that in favor of a per-reaction DESeq-style framework:

- Per-cell median centering of log flux (removes OCR + biomass_max baseline)
- log2FC on reactions carried by both cells
- Direction + coherence + magnitude thresholds (|log2FC| ≥ 1)
- Hypergeometric enrichment on count of notable reactions, BH-corrected per
  contrast — the ONLY q-value reported, because hopsy draws are not biological
  replicates

Stages 0–5 (model + sampling + loop removal + pooling + subsets) are unchanged
and still feed this pipeline. Stages 6–10 are superseded for the chapter; the
tier tables in `v7_final/` are retained as frozen artifacts of the earlier
analysis.

## Linear script order

| # | Script | Produces | Reads |
|---|---|---|---|
| 01 | substrate_yields.py | ATP yield per substrate × cell (Table S7) | v7_models, LEDGER |
| 02 | atp_ceiling_physio.py | ATP ceiling under physio OCR + envelope (§3.2) | v7_models, bounds, LEDGER |
| 03 | atp_ceiling_glucose.py | ATP ceiling glucose-only control (§3.2) | v7_models, 01 output |
| 04 | comparison_tables.py | Three-layer F1/F2 tables (raw) | POOLED, SUBSETS |
| 05 | pca_raw.py | PCA + hypergeometric enrichment | POOLED, 04 output |
| 06 | pca_centered.py | Per-cell-centered PCA + comparison fig | POOLED, 05 output |
| 07 | normalize_centered.py | Centered F1/F2 tables + re-run enrichment (§3.1 primary) | POOLED, 04 output |
| 08 | capacity_figure.py | Xenobiotic clearance LPs + 2×2 figure (§3.3) | v7_models, 02 + 03 output |
| 09 | density_plots.py | Developmental + light-regime ridges (§3.1 figures) | POOLED, SUBSETS |

## Running

Stage 0–5 inputs are large and stay outside git. Point `CSM_DATA_ROOT` at
the laptop's local Discrete_Models tree, or omit the env var to use the
default path in `scripts/config.py`.

```bash
export CSM_DATA_ROOT=/Users/.../Discrete_Models
python v8_analysis/scripts/01_substrate_yields.py
python v8_analysis/scripts/02_atp_ceiling_physio.py
python v8_analysis/scripts/03_atp_ceiling_glucose.py
python v8_analysis/scripts/04_comparison_tables.py
python v8_analysis/scripts/05_pca_raw.py
python v8_analysis/scripts/06_pca_centered.py
python v8_analysis/scripts/07_normalize_centered.py
python v8_analysis/scripts/08_capacity_figure.py
python v8_analysis/scripts/09_density_plots.py
```

Scripts 01, 02, 03, 08 need Gurobi. Scripts 04–07, 09 are pandas-only.

## Outputs

Committed to git under `chapter_assets/`:
- `comparison_tables/subsystem_summary_F{1,2}_centered.csv` — §3.1 primary
- `comparison_tables/per_reaction_flux_F{1,2}_centered.csv` — atomic layer
- `comparison_tables/reaction_inventory_F{1,2}.csv` — presence/absence
- `comparison_tables/pca_cells*.{csv,pdf}` — PCA (raw + centered)
- `tables/table_{1..6}.csv`, `table_S{1..15}.csv` — chapter tables
- `tables/xenobiotic_capacity.csv` — §3.3
- `figures/metabolic_capacity_developmental.pdf` — §3.2+§3.3 2×2
- `figures/subsystem_density/{developmental,light_regime}/` — §3.1 ridges

## Archive

`archive/` holds ~100 exploratory scripts from the pre-reframe work:
Tier 0 drilldown, lipid-relay, Grace resample remediation, density-plot
iterations v1–v4, etc. Preserved for provenance; not needed to reproduce
the chapter.
