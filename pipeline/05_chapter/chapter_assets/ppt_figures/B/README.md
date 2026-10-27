# Figure B family — findings heatmaps per thread

Seven outputs, built by `pipeline/05_chapter/scripts/11_master_figures.py`.
Five per-thread chapter documents (one PDF per thread) + two giant
supplementary masters.

Significance: single threshold at **q < 0.05**, marked by a **maroon
box border** around the cell. No stars, no second threshold. Grey = the
subsystem had fewer than 5 representatives tested in that contrast
(test not powered).

## Chapter figures — one per thread

Each document contains the subsystems of one thread with both panels:
- Top: F1 (light regime) — 3 facets side-by-side, cells = regime-pair log₂FC at each hpf
- Bottom: F2 (hpf transition) — 4 facets side-by-side, cells = transition log₂FC (sign flipped: red = rose)

| Document | Thread | Subsystems |
|---|---|---|
| `B_group_lipid_backbone.{pdf,png}` | Lipid backbone | Leukotriene; Arachidonic acid; Eicosanoid; Fatty acid oxidation; Fatty acid activation (cytosolic) |
| `B_group_peroxisomal.{pdf,png}` | Peroxisomal | β-oxidation even-chain (peroxisomal); β-oxidation n-9 (peroxisomal) |
| `B_group_eicosanoid_stall.{pdf,png}` | Eicosanoid stall | Fatty acid desaturation (even-chain); Omega-3 fatty acid metabolism |
| `B_group_sterol_bile.{pdf,png}` | Sterol & bile | Bile acid biosynthesis; Steroid metabolism |
| `B_group_signaling_aa.{pdf,png}` | Signaling / AA | Serotonin and melatonin biosynthesis; Phe, Tyr and Trp biosynthesis; Tyrosine metabolism |

## Supplementary masters — every qualifying subsystem

| Document | What |
|---|---|
| `SB_F1_regime_pair_heatmap.{pdf,png}` | F1 facets: BL vs D | BL vs LD | D vs LD. Rows = ~80 subsystems with ≥5 reps in ≥1 F1 contrast |
| `SB_F2_transition_heatmap.{pdf,png}` | F2 facets: 24→48 | 48→72 | 72→96 | 96→120. Rows = same ~80 subsystems |

## How each finding reads

| Finding | Document | Where |
|---|---|---|
| 1 — shared lipid rise 48→72 | `B_group_lipid_backbone` | F2 48→72 panel: maroon boxes on Leukotriene / Arachidonic / Eicosanoid rows across all 3 regimes |
| 2 — BL_48 FA activation | `B_group_lipid_backbone` | F1 BL vs D & BL vs LD panels (48 hpf) red + maroon on FA activation (cyto) row |
| 3 — BL_72 melatonin dip | `B_group_signaling_aa` | F1 BL vs D & BL vs LD panels (72 hpf) blue + maroon on Serotonin/melatonin row |
| 4 — LD peroxisomal lead at 72 | `B_group_peroxisomal` | F1 D vs LD panel (72 hpf) dark blue + maroon on both peroxisomal rows; F2 96→120 D panel maroon (D catches up) |
| 5 — BL_96 eicosanoid stall | `B_group_eicosanoid_stall` | F1 BL vs D (72 hpf) FA desaturation blue + maroon |
| 6 — LD_120 bile+steroid | `B_group_sterol_bile` | F1 BL vs LD & D vs LD panels (120 hpf) dark blue + maroon on Bile acid; F2 96→120 LD panel maroon on Steroid |
| 7 — BL_120 Phe/Tyr/Trp | `B_group_signaling_aa` | F1 BL vs D panel (120 hpf) red + maroon on Phe/Tyr/Trp and Tyrosine rows |

## Caveat

Row n (max reps tested across the contrasts shown in that panel) is in
the row label. A 4-rep subsystem (Serotonin and melatonin biosynthesis)
sits next to a 165-rep subsystem (Fatty acid oxidation) with the same
colour scale and the same maroon-box rule. The row n tells the reader
the power behind each colour.
