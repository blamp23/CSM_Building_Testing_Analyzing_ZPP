# Zebrafish-embryo context-specific metabolic modeling under three light regimes

Fifteen context-specific genome-scale metabolic models (3 photoperiods × 5
developmental stages, 24–120 hpf), extracted with rFASTCORMICS from measured
transcriptomes + metabolomics, sampled with hopsy, and analyzed for shared
developmental programs and light-regime-specific metabolic signatures.

**Chapter output lives in `pipeline/05_chapter/chapter_assets/`.** The nine
scripts that produce it are in `pipeline/05_chapter/scripts/01..09`, run in
that order.

## 60-second reproduction

```bash
# 1. Install Python deps (requires Gurobi + MATLAB R2025b for Stages 0-4)
pip install -r pipeline/05_chapter/requirements.txt

# 2. Point at the local mirror of the Stage 0-5 outputs — model .mat files
#    and hopsy sample pools — too big to ship in git.
export CSM_DATA_ROOT=/path/to/Discrete_Models

# 3. Run the chapter pipeline in order.
cd pipeline/05_chapter/scripts
python 01_substrate_yields.py     # Table S7 — ATP yield sanity check
python 02_atp_ceiling_physio.py   # §3.2 — ATP ceiling under physio OCR
python 03_atp_ceiling_glucose.py  # §3.2 — glucose-only control
python 04_comparison_tables.py    # §3.1 — raw per-reaction F1/F2 tables
python 05_pca_raw.py              # PCA + initial enrichment
python 06_pca_centered.py         # Centered PCA comparison
python 07_normalize_centered.py   # §3.1 primary — centered F1/F2 tables
python 08_capacity_figure.py      # §3.3 — xenobiotic clearance + 2×2 figure
python 09_density_plots.py        # §3.1 — subsystem density ridges
```

Scripts 01, 02, 03, 08 need Gurobi. 04–07, 09 are pandas-only.

## Repo layout

```
pipeline/                All executable stages, numbered in order
  00_bake/               Wang zebrafish GEM + OCR cap per hpf → baked_ocr_{hpf}_v7.mat
  01_extract/            rFASTCORMICS → 15 context-specific .mat models
  02_sampling/           PolyRound polytope rounding + hopsy MCMC + polytope cond.
  03_loop_removal/       CycleFreeFlux (Gurobi) → loopless flux samples
  04_matched_floor/      F2 matched-floor resample for 9 cells (§11.20)
  05_chapter/            THE CHAPTER PIPELINE
    scripts/01..09.py    Linear chapter pipeline (reads pooled samples)
    scripts/config.py    Single source of truth for all file paths
    chapter_assets/      §3.1-3.4 primary outputs (tables + PCA + figures)
    archive/             104 exploratory scripts preserved for provenance

functions/               Shared Python/MATLAB lib — apply_context_bounds
data/                    Shared input data — exchange bounds, gene maps
rules/                   Shared input data — ledger + pool table
docs/                    Small validation CSVs + version notes
v7_final/                Frozen Tier A/B/C deliverable (reviewer-citable artifact)
legacy/                  Everything the chapter does NOT use (see legacy/README.md)
PIPELINE.md              Full stage-by-stage description
```

## What each stage produces

| Stage | Produces | Needs |
|---|---|---|
| 00_bake | 5 base `.mat` files (`baked_ocr_{hpf}_v7.mat`) | Wang GEM, measured OCR per hpf |
| 01_extract | 15 context-specific `.mat` files (`trans_rfastcormics_{cell}.mat`) | 00_bake output, 15 transcriptomes |
| 02_sampling | 15 × 4 hopsy chain CSVs (billiard-walk samples) | 01_extract output, measured envelope |
| 03_loop_removal | 15 × 4 CycleFreeFlux-corrected chain CSVs | 02_sampling output, Gurobi |
| 04_matched_floor | 9-cell F2 resample pool at matched biomass floors | 01_extract output, Grace cluster |
| **05_chapter** | **§3.1–3.4 tables, PCAs, 2×2 figure, density ridges** | **03 + 04 pooled samples** |

Full details, inputs/outputs, and the amendment governing each stage are in
`PIPELINE.md`.

## Dependencies

- Python 3.13+ — pandas, numpy, scipy, matplotlib, scikit-learn, h5py, hopsy, PolyRound, gurobipy
- MATLAB R2025b + COBRA Toolbox (stages 00 and 01 only)
- Gurobi (academic license; used by 02, 03, 05/scripts/01/02/03/08)
- R 4.4+ with tidyverse (optional — only for legacy `legacy/older_investigations/fba_fva/` R scripts)

## Reproducing from scratch

- **With the model files on hand:** export `CSM_DATA_ROOT` and run the chapter
  scripts (01..09 in `pipeline/05_chapter/scripts/`). That reproduces every
  chapter output.
- **From Stage 00:** run `pipeline/00_bake/01_bake_ocr.m` →
  `pipeline/00_bake/02_build_ledger.py` →
  `pipeline/00_bake/03_apply_v7_boundary.m` →
  `pipeline/01_extract/04_extract_v7.m` to regenerate the 15 context-specific
  models from the Wang GEM. Then sample (`pipeline/02_sampling/`), apply CFF
  (`pipeline/03_loop_removal/`), pool, and run `pipeline/05_chapter/scripts/`.

## Tags

- `v7-frozen` (2026-10-01) — pre-reframe pipeline + initial v7 model set
- `chapter-v1` (2026-10-26) — chapter outputs + flat `pipeline/` restructure
