# Figure B — Findings heatmap

Two outputs, built by `pipeline/05_chapter/scripts/11_master_figures.py`.

## Chapter figure

**`B_findings_heatmap.{pdf,png}`** — the 14 subsystems the chapter
findings live in, grouped by thread (bracketed labels on the right).

Rows (in chapter order):
1. Leukotriene metabolism
2. Arachidonic acid metabolism
3. Eicosanoid metabolism
4. Fatty acid oxidation
5. Fatty acid activation (cytosolic)
6. Beta oxidation of even-chain fatty acids (peroxisomal)
7. Beta oxidation of unsaturated fatty acids (n-9) (peroxisomal)
8. Fatty acid desaturation (even-chain)
9. Omega-3 fatty acid metabolism
10. Bile acid biosynthesis
11. Steroid metabolism
12. Serotonin and melatonin biosynthesis  ← chapter exception (n=4)
13. Phenylalanine, tyrosine and tryptophan biosynthesis
14. Tyrosine metabolism

Columns: `BL D LD` within each stage block (24, 48, 72, 96, 120 hpf).

**Cell value.** For each cell and each representative, centered log₂|flux|
= log₂(mean non-zero flux) − cell median log₂|flux|. Subsystem value =
mean over representatives present in that cell. Row-z-scored across the
15 cells so every row has the same scale. Values are from the F1
(per-cell-floor) samples, not matched-floor F2 samples.

**Grey cells.** Subsystem absent from that cell per
`reaction_inventory_F1.csv` (peroxisomal β-oxidation at 24 and 48 hpf;
BL_72 and D_72 for even-chain peroxisomal; LD_72 for same subsystem is
filled, which is Finding 4).

**Markers.** Dots between regime columns at a stage = F1 enrichment q <
0.05 with n ≥ 5. Arrows between stage blocks within a regime row = F2
enrichment q < 0.05 with n ≥ 5. F2 markers come from the matched-floor
F2 samples, so they are annotations from a different sample set than
the heatmap colours.

**Caveat.** Row-z-scoring makes a 4-reaction subsystem (melatonin) look
as confident as a 48-reaction one. Each row label carries its max cell
n so the reader sees it.

## Supplementary master

**`SB_all_subsystems_heatmap.{pdf,png}`** — same column structure, same
cell value computation, every subsystem with ≥5 representatives in ≥1
cell. 83 rows. Chapter rows highlighted in bold. This is the figure the
chapter 14 are selected from; Figure B is literally a row-selection of
SB plus the chapter row ordering + thread brackets.
