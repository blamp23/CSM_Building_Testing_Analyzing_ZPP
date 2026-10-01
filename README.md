# CSM_Building_Testing_Analyzing_ZPP

Context-specific metabolic model pipeline for zebrafish embryos under three photoperiods (BL / D / LD × 5 developmental stages = 15 cells).

**Model set**: v7 (frozen 2026-10-01). See `FAMILY_PREREG.md` §10 for the pre-registration and amendment log.

## Repo layout

```
data/           5 input files (metabolomics bounds, gene labels, mappings)
rules/          v7 curation ledger + pool table (immutable in v7)
functions/     apply_context_bounds (shared Python + MATLAB twin) + COBRA helpers
base_model/     Phase 1 — bake OCR cap, build ledger, apply boundary rule
extraction/     Phase 2 — rFASTCORMICS extraction × 15 cells
validation/     Phase 3 — closure, ATP, parent≥child, s128c gate, Tier 0, sizes
fba_fva/        Phase 4 — pFBA + ll-FVA (Frame A + matched floor), augmentation, analysis
sampling/       Phase 5 — hopsy MCMC + CFF + two-tier (local prep + Grace SLURM)
statistics/     Phase 6 — gamma-tail BH, neg-control, findings annotation
docs/           version notes, pipeline doc
FAMILY_PREREG.md   pre-registration + §10 v7 amendment
```

## Setup

```bash
export V7_REPO=/path/to/compute/workspace   # dir containing baked_ocr_*_v7.mat
pip install -r requirements.txt             # (TODO)
```

Scripts read `$V7_REPO` to locate the large binary artefacts (base model `.mat` files, extracted cells, PolyRound caches) that are not tracked in git.

## Pipeline

Numbered scripts run in order. Each stage writes outputs that the next stage consumes.

1. `base_model/01_bake_ocr.m` — Wang zebrafish GEM + OCR cap per hpf → 5 `baked_ocr_{hpf}.mat`
2. `base_model/02_build_ledger.py` — automatic class-1/2/3 classification → `rules/v7_ledger.csv`
3. `base_model/03_apply_v7_boundary.m` — ε=1e-3 pool coef + ledger closes + Part 1/Part 2 + medium → 5 `baked_ocr_{hpf}_v7.mat`
4. `extraction/04_extract_v7.m` — rFASTCORMICS × 15 cells (parfor) → 15 `trans_rfastcormics_{cell}.mat`
5. **Validation gates** (`validation/05_*` through `12_*`): closure, ATP tests, parent ≥ child, s128c gate, Tier 0 constraint decomposition + ε/BOUND_EX scans, sizes, phosphocholine note.
6. `fba_fva/13_*` through `18_*` — deterministic analysis (pFBA, ll-FVA Frame A + matched floor, R analysis + figures).
7. `sampling/preprocess_local.sh` — PolyRound locally with Gurobi, caches 15 polytopes.
8. `sampling/slurm/submit.sh` — Grace SLURM chain: sample → aggregate → CFF → subsets → verify → two-tier.
9. `statistics/26_*` through `28_*` — gamma-tail p + BH + split-chain neg control + Tier A/B/C annotation on the samples.

See `FAMILY_PREREG.md` §10 for the frozen configuration (BOUND_EX=10, ε=1e-3, biomass floor 0.9×, OCR per hpf, pre-registered family filter, etc.).

## Dependencies

- Python 3.13+ (scipy, numpy, pandas, h5py, hopsy, PolyRound, gurobipy)
- MATLAB R2025b + COBRA Toolbox
- Gurobi (academic license)
- R 4.4+ (tidyverse, ggplot2, data.table) for `fba_fva/16_*–18_*`
