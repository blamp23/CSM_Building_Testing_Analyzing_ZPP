# v7 frozen — 2026-10-01

The v7 base, extractions, ledger, and bounds function are frozen as of
2026-10-01. No change permitted without a failed gate.

## Base
- `ocr_anchored_extraction/models/baked_ocr_{24,48,72,96,120}_v7.mat`
  - MAR00022 pool: 4 reduced-redox substrates coef=0, others unchanged
  - MAR00021 pool coef (MAM01602c) = −1e-3 (ε)
  - 20 ledger close_exchange rows: lb=ub=0
  - Part 1 (secretion uncap): ub=1000 for unmeasured non-ledger non-medium
  - Part 2 (uptake close): lb=0 for same set
  - 17 medium rows: lb=-1000, ub=1000

## Extracted
- `reviewer_packet_v4_slim/v7_models/trans_rfastcormics_*.mat` (15 cells)
- Sizes: see `s147_sizes_v6_v7.csv`

## Shared bounds function
- `functions/apply_context_bounds.{py,m}`
- lb = -BOUND_EX × fc, ub = +BOUND_EX × fc on 105 (BOUND_EX=10, fallback fc=1.0)
- OCR cap on MAR09048 at -OCR(hpf)
- medium_uptake not in 105: lb=-1000
- Assertions: closes closed; Part 1+2 rows at [0,1000]; medium lb<0; exactly N measured scaled

## Ledger
- `reviewer_packet_v4_slim/s130_v7_ledger.csv` (75 rows)

## Validation
- `s138_atp_v7_base.csv` — base ATP tests 1–5 (all pass)
- `s141_v7_atp_per_cell.csv` — per-cell ATP tests 1 & 3 (all pass)
- `s146_tier0_*.csv` — Tier 0 constraint decomposition + ε/BOUND_EX scans
- `s128c_uptake_essential_v7.csv` — gate pass (class-1 only + biomass pseudo)
- `s128c_secretion_essential_v7.csv` — 5-dA LIAS byproduct sink separately
- `s147_sizes_v6_v7.csv` — model sizes
- `s148.log` — phosphocholine note (limitations)

## Amendment
`reviewer_packet_v4_slim/FAMILY_PREREG.md` §10 (full text).
