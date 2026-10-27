# Figure B′ — Findings trajectories

Two outputs, built by `pipeline/05_chapter/scripts/11_master_figures.py`.

Line-plot twin of Figure B. Same rows, each as its own small panel; x =
hpf, y = mean centered log₂|flux| for the subsystem in that cell (NOT
z-scored — raw values are preserved so panels are comparable to the
density PDFs in `../../figures/subsystem_density/`). Three lines per
panel, one per regime (BL blue, D grey, LD orange). Open white circle
at y=0 means the subsystem is absent from that cell.

## Chapter figure

**`B_prime_findings_lines.{pdf,png}`** — the 14 subsystems from Figure B,
3 columns × 5 rows. Each finding is the point where one regime leaves
the other two (Finding 2: BL_48 FA activation; Finding 3: BL_72
melatonin dip; Finding 4: LD_72 peroxisomal spike; Finding 5: BL_96
eicosanoid/desaturation stall; Finding 6: LD_120 bile and steroid
climb; Finding 7: BL_120 Phe/Tyr/Trp rise; Finding 1 is the shared
48→72 lipid peak across all three regimes).

## Supplementary master

**`SB_prime_all_subsystems_lines.{pdf,png}`** — the same panel for every
one of the 83 qualifying subsystems in SB, sorted by metabolic class.
Several pages tall. Figure B′ is a row selection of SB′.
