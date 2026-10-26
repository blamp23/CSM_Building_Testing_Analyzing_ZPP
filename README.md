# Zebrafish-embryo context-specific metabolic modeling under three light regimes

Fifteen context-specific genome-scale metabolic models (3 photoperiods × 5
developmental stages, 24–120 hpf), extracted with rFASTCORMICS from measured
transcriptomes + metabolomics, sampled with hopsy, and analyzed for shared
developmental programs and light-regime-specific metabolic signatures.

**Chapter output lives in `v8_analysis/chapter_assets/`.** The nine scripts
that produce it are in `v8_analysis/scripts/01..09`, run in that order.

## 60-second reproduction

```bash
# 1. Install Python deps (requires Gurobi + MATLAB R2025b for Stages 0-5)
pip install -r v8_analysis/requirements.txt

# 2. Point at the local mirror of the large binary artifacts (model .mat files
#    and hopsy sample pools) that are too big to ship in git.
export CSM_DATA_ROOT=/path/to/Discrete_Models

# 3. Run the chapter pipeline, in order.
cd v8_analysis/scripts
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

## What's where

```
base_model/            Stage 0 — bake Wang zebrafish GEM with OCR cap per hpf
extraction/            Stage 0 — rFASTCORMICS → 15 context-specific models
functions/             Stage 1 — apply_context_bounds (envelope + OCR + ledger)
data/                  Stage 1 inputs — exchange bounds, gene maps
rules/                 Stage 1 inputs — ledger + pool table
sampling/              Stages 2-3 — PolyRound + hopsy MCMC
v7_analysis/
  preprocessing/       Stage 2 — LD_72 polytope conditioning
  loop_removal/        Stage 4 — CycleFreeFlux via Gurobi
  matched_floor/       Stage 11 — matched-floor F2 resample pipeline
v8_analysis/           Stage 12 — CHAPTER PIPELINE
  scripts/01..09       the 9-script linear pipeline
  chapter_assets/      §3.1-3.4 primary outputs (tables + figures)
  archive/             104 exploratory scripts (not needed to reproduce)
v7_final/              Frozen Tier A/B/C deliverable from the earlier analysis
                         (retained as reviewer-citable artifact; not cited by chapter)
docs/                  Small validation CSVs + version notes
legacy/                Everything the chapter does NOT use. See legacy/README.md.
PIPELINE.md            Full stage-by-stage description of what runs when
```

## Which stage does what

| Stage | Produces | Where |
|---|---|---|
| 0 | 15 context-specific `.mat` models | `base_model/`, `extraction/` |
| 1 | Envelope + OCR + ledger applied to each model | `functions/apply_context_bounds.py` |
| 2 | 15 rounded flux polytopes | `sampling/`, `v7_analysis/preprocessing/` |
| 3 | 15 hopsy sample files (4 chains × 500 draws each) | `sampling/scripts/`, `sampling/slurm/` |
| 4 | CycleFreeFlux-corrected (loopless) samples | `v7_analysis/loop_removal/` |
| 5 | Pooled samples + enzyme-subset membership | sampling-pipeline outputs |
| 11 | Matched-floor F2 resample for 9 cells | `v7_analysis/matched_floor/` |
| **12** | **Chapter outputs — tables, PCAs, figures, density ridges** | **`v8_analysis/scripts/01..09`** |

Full details, inputs/outputs, and the amendment governing each stage are in
`PIPELINE.md`. Stages 6–10 (W1 contrasts, Tier A/B/C export, pre-reframe
figures) are no longer part of the chapter and live in `legacy/`.

## Dependencies

- Python 3.13+ with: pandas, numpy, scipy, matplotlib, scikit-learn, h5py, hopsy, PolyRound, gurobipy
- MATLAB R2025b + COBRA Toolbox (Stage 0 only)
- Gurobi (academic license; used by Stages 2, 4, and v8 scripts 01/02/03/08)
- R 4.4+ with tidyverse (optional — only for the legacy `fba_fva/` R scripts)

## Reproducing from scratch

- **With the model files on hand:** export `CSM_DATA_ROOT` and run
  `v8_analysis/scripts/01..09`. That reproduces every chapter output.
- **From Stage 0:** run `base_model/01_bake_ocr.m` → `02_build_ledger.py` →
  `03_apply_v7_boundary.m` → `extraction/04_extract_v7.m` to regenerate the
  15 context-specific models from the Wang GEM + the measured exchange bounds.
  Then sample (Stages 2–3), apply CFF (Stage 4), pool (Stage 5), and run
  `v8_analysis/scripts/01..09`.

The `v7-frozen` tag marks the pre-reframe pipeline state. The `chapter-v1` tag
marks the state used to produce the chapter.

## Tags

- `v7-frozen` (2026-10-01) — pre-reframe pipeline + initial v7 model set
- `chapter-v1` (2026-10-26) — current chapter output + v8 pipeline added
