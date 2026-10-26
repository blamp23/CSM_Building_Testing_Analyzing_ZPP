# legacy/

Preserved for provenance. **Nothing in here is needed to reproduce the chapter.**

If you cloned this repo to re-run the analysis, skip this directory. If you're
reviewing the project's history or asking "why did they do X instead of Y,"
this is where the earlier work lives.

## What's here

### `v7_pipeline_deprecated/` — Stages 6–10 of the earlier W1/Tier A pipeline

The repo was originally built around a two-tier Wasserstein (W1) contrast
framework (Stage 6 → statistics/ → findings → Tier A/B/C export → figures).
The 2026-10-26 chapter reframe dropped that framework after it became clear
that hopsy samples are not biological replicates and per-sample W1 tests are
overpowered. The scripts still work; they just don't appear in the chapter.

| Subdir | What it did | Superseded by |
|---|---|---|
| `statistics/` | W1 contrasts (3 normaliser variants) + family findings | `pipeline/05_chapter/scripts/04`, `07` |
| `figures/` | fig1–5 + figS4/S5 heatmaps / developmental axis | `pipeline/05_chapter/scripts/08`, `09` |
| `tier0/` | Parent/child FBA + BL_96 greedy addback + MAR10693 vent | Finding is now carried by §3.1/§3.2 isoleucine topology thread |
| `positive_control/` | Steroid KO + KD2 positive controls on LD_72 | Chapter does not report positive-control validation |
| `launchers/` | run_contrasts_*.sh batch drivers for 42 contrasts | Not needed; no W1 contrasts in chapter |
| `plots/` | KDE / CDF / trajectory alternatives | `pipeline/05_chapter/scripts/09` density ridges |

The frozen Tier A/B/C deliverable from this pipeline is retained at the repo
root in `v7_final/` — see its `FREEZE.md` for md5s.

### `older_investigations/` — pre-W1 investigation code

Early mixed-language scripts (MATLAB + R + Python) that preceded the v7
pipeline. Numbered 05–28 at root level originally.

| Subdir | What it did |
|---|---|
| `validation/` | 05–12: closure, ATP base/per-cell, parent≥child, gate essentiality, Tier 0, sizes, phosphocholine |
| `fba_fva/` | 13–17: pFBA, ll-FVA (Frame A + matched floor), augmentation, R pairwise deltas |
| `statistics/` | 26–28: gamma-tail BH, split-chain neg control, findings annotation |

These produced the small CSVs in `docs/` (`validation_atp_per_cell.csv` etc.)
that the chapter does not cite.

### `docs/`

- `FAMILY_PREREG.md` — 51KB pre-registration + §10 v7 amendment trail +
  §11.1–§11.22 amendments. Reviewer-facing, cites the W1 analysis.
- `v7_methods_for_reviewer.md` — 20KB methods walkthrough for the W1 pipeline.

Both are retained for anyone asking what was pre-registered and how the earlier
analysis was done. The chapter's methods section covers the v8 reframe and does
not depend on these.
