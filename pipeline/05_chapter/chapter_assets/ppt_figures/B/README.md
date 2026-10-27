# Figure B family — findings heatmaps

Four outputs, built by `pipeline/05_chapter/scripts/11_master_figures.py`.
Two chapter figures (14 subsystems), two supplementary masters (every
qualifying subsystem). The chapter figures are literal row-selections of
the masters.

Replaces the earlier per-cell heatmap (removed): that version mixed F1
and F2 markers in one display and used confusing arrow glyphs. The split
below mirrors how the chapter itself talks about the findings — one
section per contrast family.

## Chapter figures

**`B_F1_regime_pair_heatmap.{pdf,png}`** — F1 (light-regime) findings.
Three facets side-by-side: BL vs D | BL vs LD | D vs LD. Rows =
chapter-14 subsystems. Columns = 5 hpf stages. Cell value = mean
centered log₂FC between the two cells for that (subsystem, pair, hpf)
contrast (positive = alphabetically earlier regime higher). Stars inside
cells: `*` for q<0.05, `**` for q<0.001. Grey = subsystem had <5
representatives in that contrast (test not powered).

**`B_F2_transition_heatmap.{pdf,png}`** — F2 (developmental) findings.
Four facets side-by-side: 24→48 | 48→72 | 72→96 | 96→120 hpf. Rows =
chapter-14 subsystems. Columns = 3 regimes (BL/D/LD). Cell value =
mean centered log₂FC with sign flipped so **red = rose from earlier
to later hpf**. Stars as above.

## Supplementary masters

**`SB_F1_regime_pair_heatmap.{pdf,png}`** — same facet structure as
Figure B, but every subsystem with ≥5 reps tested in at least one F1
contrast (~80 subsystems).

**`SB_F2_transition_heatmap.{pdf,png}`** — same facet structure as
Figure B2, but every subsystem with ≥5 reps tested in at least one F2
contrast.

## How each finding reads

| Finding | Figure | Facet | Where |
|---|---|---|---|
| 1 — shared lipid rise 48→72 | B2 | 48→72 | Leukotriene / Arachidonic / Eicosanoid rows, red ** in BL+D+LD |
| 2 — BL_48 FA activation | B / B2 | B: BL vs D (48 hpf) red **  ·  B2: 24→48 (D+LD blue **) | |
| 3 — BL_72 melatonin dip | B | BL vs D (72 hpf) blue * and BL vs LD (72 hpf) blue * | |
| 4 — LD peroxisomal lead at 72 | B2 | 48→72 and 96→120 (D red ** at 96→120 = D catching up) | |
| 5 — BL_96 eicosanoid stall | B | BL vs D (72 hpf) FA desaturation blue * | |
| 6 — LD_120 bile+steroid | B | BL vs LD and D vs LD at 120 hpf, bottom-right ** | |
| 7 — BL_120 Phe/Tyr/Trp | B | BL vs D (120 hpf) ** on Phe/Tyr/Trp + Tyrosine rows | |

## Caveat

Row n (max tested across the contrasts shown in that panel) is in the
row label. A 4-rep subsystem (Serotonin and melatonin biosynthesis)
sits next to a 165-rep subsystem (Fatty acid oxidation) and the
colour scale treats them identically. The row n tells the reader the
power behind each colour.
