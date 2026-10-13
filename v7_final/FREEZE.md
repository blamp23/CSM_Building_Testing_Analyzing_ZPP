# v7_final — FREEZE 2026-10-13

**Nothing in this folder is edited after this freeze.**

## Funnel (method of record: mean + drop zeros)

- Rows entered the §1 family: **1687** (across 27 contrasts)
- Rows with q_bh ≤ 0.10: **146**
- Tier A: **50** (F1=26, F2=24)
- Tier B: **63** (F1=26, F2=37)
- Tier C: **33** (F1=12, F2=21)  — includes 12 `flag_reason = floor_artifact`

## Tier A by subsystem

| subsystem | n | F1 | F2 |
|---|---|---|---|
| Arginine and proline metabolism | 7 | 5 | 2 |
| Leukotriene metabolism | 6 | 2 | 4 |
| Phenylalanine, tyrosine and tryptophan biosynthesis | 6 | 4 | 2 |
| Pyrimidine metabolism | 6 | 4 | 2 |
| Arachidonic acid metabolism | 4 | 0 | 4 |
| Fatty acid activation (cytosolic) | 4 | 0 | 4 |
| Nucleotide metabolism | 4 | 3 | 1 |
| Valine, leucine, and isoleucine metabolism | 3 | 2 | 1 |
| Acylglycerides metabolism | 2 | 2 | 0 |
| Bile acid biosynthesis | 2 | 2 | 0 |
| Fatty acid oxidation | 2 | 0 | 2 |
| Glycolysis / Gluconeogenesis | 2 | 0 | 2 |
| Steroid metabolism | 2 | 2 | 0 |

## pre_correction_retinoid flag (§11.17 / §11.19 provenance)

- Tier A: 0  (0 expected)
- Tier B: 3  (3 expected)
- Tier C: 3  (3 expected)

## Amendments included (§11.1–§11.20)

- **§11.1** — F1 family — between-condition contrasts at the same hpf
- **§11.2** — F2 family — within-condition adjacent-hpf contrasts
- **§11.3** — Negative control — split-chain within each cell
- **§11.4** — Per-row inclusion rules (carry_both, loop_flagged, etc.)
- **§11.5** — Robustness tiers v1 (superseded by §11.13)
- **§11.9** — Pool reactions excluded from stats
- **§11.10** — Floor-suspect flag + BL_96 matched-floor resample (F1)
- **§11.11** — Frame A retinol 48→72 limitation
- **§11.12** — Median normaliser sensitivity
- **§11.13** — SNR tiering + envelope-tracking flag (supersedes §11.5)
- **§11.14** — Zero handling: mean + drop zeros is method of record; EPS retracted
- **§11.15** — Positive control: Steroid KO Fisher p = 8e-11
- **§11.16** — Positive control: Steroid KD2 W1 = 2.29, snr ≈ 38
- **§11.17** — Tier 0 topology × envelope decomposition + BL_96 drill-down (+ MAR10693 degenerate-vent footnote + Tier 0 robustness)
- **§11.19** — medium_uptake demand bounds (ATP session)
- **§11.20** — F2 matched-floor resample + demotion rule

## MAR10693 Tier 0 vent check (§11.17 footnote)

- 15 LPs with MAR10693.ub = 0: biomass_max change **0.00% in every cell**
- Verdict: dual degeneracy — one of a family of interchangeable gene-less secretion stubs
- Transporter MAR10935 is 1:1 H⁺ symport (not Na⁺)
- Data: `grace_v4/results/v7_final_addenda/mar10693_biomax_table.csv`

## Tier 0 robustness — 60 LPs with 10% uptake relaxation (§11.17)

- **Glucose (carbon control):** max |Δ| = 0.07 % — passes
- **Valine:** max |Δ| = 0.007 % — not a bottleneck anywhere
- **Histidine:** +4.4 to +5.5 % at BL_48, D_48, D_72, LD_48, LD_72; +2.4 % at LD_96 — mid-development ceiling
- **Isoleucine:** +8.2 % at BL_96 specifically; +0.7 % at BL_120
- Verdict: nitrogen-bound, **histidine-dominated in the 48–72 window, isoleucine-bound at BL_96**
- Data: `grace_v4/results/v7_final_addenda/tier0_robustness.csv`

## md5 checksums

| file | md5 | bytes |
|---|---|---|
| `FAMILY_PREREG.md` | `a987fca8bb711952467802ecaa79e8cf` | 51843 |
| `all_rows_in_family.csv` | `f00378ad1ca9927a68b0a022befc094b` | 559584 |
| `heat_F1.pdf` | `18cc2e093155c96f2d93f19f38ab157d` | 28727 |
| `heat_F1.png` | `6f44e825094db6560a07a7867ff4a403` | 273315 |
| `heat_F2.pdf` | `37d2f1354973f462d236baf22ee04da9` | 35137 |
| `heat_F2.png` | `674aed7f82768ff1308978653f2aad7e` | 519430 |
| `negcontrol_all_rows.csv` | `957f832b6ad3a3406db11d77870580d8` | 167331 |
| `tierA_by_subsystem.csv` | `ff24e934c9795f999e06a0afb808477f` | 1844 |
| `tierA_final.csv` | `b81056a83607d703354f38b8c7261133` | 19973 |
| `tierB_final.csv` | `bc625ae7c375da36b38451eca29723cd` | 21750 |
| `tierC_final.csv` | `14906bf0c96f17c0e200faa071093d94` | 15247 |