# pipeline/05_chapter — chapter analytical pipeline (Stage 12)

Nine linear scripts that take frozen Stage 00–04 outputs (context-specific
models + hopsy flux samples + pooled samples) and produce the chapter's
findings tables, PCA, metabolic capacity LPs, and density figures.

## Reframe background

The earlier pipeline (now in `legacy/v7_pipeline_deprecated/`) summarized
sampled flux via Wasserstein-1 (W1) contrasts and a Tier A/B/C exportable
table. The 2026-10-26 reframe dropped that because hopsy samples are not
biological replicates — per-sample tests are overpowered. The current
framework instead:

- Per-cell median centering of log flux (removes OCR + biomass_max baseline)
- log2FC on reactions carried by both cells
- Direction + coherence + magnitude thresholds (|log2FC| ≥ 1)
- Hypergeometric enrichment on count of notable reactions, BH-corrected
  per contrast — the ONLY q-value reported

Stages 00–04 and the Stage 11 matched-floor resample feed this pipeline
unchanged. The frozen Tier A/B/C deliverable at `v7_final/` is retained
as a reviewer-citable artifact of the earlier analysis.

## Linear script order

| # | Script | Produces | Reads |
|---|---|---|---|
| 01 | `01_substrate_yields.py` | ATP yield per substrate × cell (Table S7) | v7_models, LEDGER |
| 02 | `02_atp_ceiling_physio.py` | ATP ceiling under physio OCR + envelope (§3.2) | v7_models, bounds, LEDGER |
| 03 | `03_atp_ceiling_glucose.py` | ATP ceiling glucose-only control (§3.2) | v7_models, 01 output |
| 04 | `04_comparison_tables.py` | Three-layer F1/F2 tables (raw) | POOLED, SUBSETS |
| 05 | `05_pca_raw.py` | PCA + hypergeometric enrichment | POOLED, 04 output |
| 06 | `06_pca_centered.py` | Per-cell-centered PCA + comparison fig | POOLED, 05 output |
| 07 | `07_normalize_centered.py` | Centered F1/F2 tables + re-run enrichment (§3.1 primary) | POOLED, 04 output |
| 08 | `08_capacity_figure.py` | Xenobiotic clearance LPs + 2×2 figure (§3.3) | v7_models, 02 + 03 output |
| 09 | `09_density_plots.py` | Developmental + light-regime ridges (§3.1 figures) | POOLED, SUBSETS |

## Running

Stage 00–04 inputs (model files, hopsy samples) are large and stay outside
git. Point `CSM_DATA_ROOT` at the laptop's local Discrete_Models tree, or
omit the env var to use the default path in `scripts/config.py`.

```bash
export CSM_DATA_ROOT=/Users/.../Discrete_Models
cd scripts
python 01_substrate_yields.py
python 02_atp_ceiling_physio.py
python 03_atp_ceiling_glucose.py
python 04_comparison_tables.py
python 05_pca_raw.py
python 06_pca_centered.py
python 07_normalize_centered.py
python 08_capacity_figure.py
python 09_density_plots.py
```

Scripts 01, 02, 03, 08 need Gurobi. 04–07, 09 are pandas-only.

## Outputs

Committed under `chapter_assets/`:
- `comparison_tables/subsystem_summary_F{1,2}_centered.csv` — §3.1 primary
- `comparison_tables/per_reaction_flux_F{1,2}_centered.csv` — atomic layer
- `comparison_tables/reaction_inventory_F{1,2}.csv` — presence/absence
- `comparison_tables/pca_cells*.{csv,pdf}` — PCA (raw + centered)
- `tables/table_{1..6}.csv`, `table_S{1..15}.csv` — chapter tables
- `tables/xenobiotic_capacity.csv` — §3.3
- `figures/metabolic_capacity_developmental.pdf` — §3.2+§3.3 2×2

Density figures (`figures/subsystem_density/`) are regenerable by script 09
and git-ignored.

## Archive

`archive/` holds ~100 exploratory scripts from the pre-reframe work: Tier 0
drilldown, lipid-relay, Grace resample remediation, density-plot iterations
v1–v4, etc. Preserved for provenance; not needed to reproduce the chapter.
See `archive/README.md`.
