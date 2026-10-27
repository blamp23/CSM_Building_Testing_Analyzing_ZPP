# Finding 1 — shared lipid activation program, 48 to 72 hpf

Four PPT-ready figures (16:9, 300 DPI, PDF + PNG) produced by
`pipeline/05_chapter/scripts/10_ppt_finding1.py`.

Finding 1 asks: do leukotriene, arachidonic acid, eicosanoid, fatty-acid
oxidation and fatty-acid activation subsystems all shift the same way
across the 48→72 hpf transition in all three regimes? Yes for the first
three (q passes in all regimes); yes with lower coherence for FA oxidation
(q passes only in BL); the FA activation signal is covered by Finding 2.

| File | What it shows |
|---|---|
| `1a_table` | Per-subsystem × per-regime effect size, direction, coherence and q at 48→72 |
| `1b_bars` | Grouped bar chart of mean centered log₂FC (positive = rose from 48 to 72 hpf); non-significant bars are hatched |
| `1c_trajectory` | Mean centered log₂FC across all four hpf transitions (24→48, 48→72, 72→96, 96→120), one subplot per subsystem, three regimes overlaid; 48→72 band highlighted |
| `1d_density_leukotriene` | Density ridges of sampled flux across the five hpf stages, three panels (one per regime) — visualises the magnitude rise at 72 vs 48 |
| `1d_density_arachidonic` | Same as above for arachidonic acid metabolism |

Sign convention: y-axis is **log₂FC (later vs earlier hpf)**, so positive =
rose across the transition. This flips the sign of the raw F2 CSV, which
stores earlier-hpf-as-A.

Rebuild: `python pipeline/05_chapter/scripts/10_ppt_finding1.py`.
