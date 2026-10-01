# Pre-registration of the Tier 2 test family and multiplicity correction

**Date locked:** 2026-09-28
**Status:** written *before* the matched-floor LD_72 job and the replicate ×
threshold biomass_max sweep complete. Both are running on 2026-09-28 but no
result from either has been read.
**Purpose:** the two negative-control exclusions from the split-chain check
(Leukotriene metabolism, Glycosphingolipid ganglio series) are only legitimate
if the exclusion rule was fixed before the positive results were selected.
Same reasoning applies to BH: it decides which hits exist, so it belongs before
the results draft, not in it.

## 1. Test family definition

A subsystem-level contrast enters the Tier 2 family iff **all** of:

1. `carry_both >= 5` — at least 5 representatives carry non-zero flux in both
   cells at the median-|v| > 1e-4 cutoff.
2. `loop_flagged / N_tier2 < 0.30` — fewer than 30% of the Tier-2-eligible
   representatives had max-|v| > 100 (the loop-flagged threshold).
3. Not the Exchange, Artificial, or Transport subsystem (structural, not
   biological — per reviewer 2026-09-28 note on Exchange/demand heatmap row).

Everything else — including any Fisher/Tier 1 tests, per-representative
stats, and diagnostic tables — is reported but does not receive multiplicity
correction and is not counted as a "finding."

## 2. Reported p-value per subsystem

Depends on `N_tier2`:

- **`N_tier2 >= 15`**: empirical p from the 2000-draw null (current pipeline).
  Reported as `p_emp`.
- **`N_tier2 < 15`**: gamma fit to the 2000 W1 null draws. The gamma is fit
  by MLE on the strictly positive W1 draws (W1 is nonnegative, right-skewed;
  gamma is the standard tail model). The reported p is 1 - Gamma_CDF(W1_obs;
  shape, scale). Reported as `p_gamma`. `p_emp` is also reported alongside
  as a reference but does not drive the family test.

Rationale: neg control on 2026-09-28 showed Glycosphingolipid ganglio series
at z = 5.22 with N_tier2 = 11 and 0 loops — a small-N tail failure of the
empirical null, not a biology signal. The gamma tail is what the reviewer
called for on 2026-09-28.

## 3. Multiplicity correction

Benjamini-Hochberg over the family defined in §1, using the p-value defined
in §2.

- Family size: reported alongside every corrected p (`family_N`).
- Corrected values reported as `q_bh`.
- Discovery threshold: `q_bh <= 0.10` (standard for exploratory omics).

## 4. What counts as a "finding" in v5

A subsystem is reported as a Tier 2 finding iff:

1. It is in the family (§1).
2. `q_bh <= 0.10` (§3).
3. Shift direction and W1_obs are consistent with the effect claimed
   (per reviewer 2026-09-28: not just a z threshold, but the paired W1
   and direction).

The v4 packet's ranking of hits by raw z is a pre-multiplicity ranking and
is not treated as a claim of significance until §1–§4 are applied.

## 5. What this pre-registration means for the v4 packet

Under the rules above, the v4 packet's top hits pre-audit as follows
(family membership only; final q-values await matched-floor + sweep):

| finding | N_tier2 | loops | in family? |
|---|---:|---:|:---:|
| D_48 Phe/Tyr/Trp                 | 25 |  1 | yes |
| D_72 Phe/Tyr/Trp                 | 38 |  1 | yes |
| BL_72 Glycerolipid               | 26 |  0 | yes |
| D_96 cholesterol esters          | 16 |  0 | yes |
| BL_96 Glycosphingolipid          |  8 |  0 | yes (N≥5 per §1) |
| BL_48 FA activation              | 66 |  6 | yes |
| BL_24 Starch/sucrose             |  8 |  0 | yes (N≥5 per §1) |

**Amendment 2026-09-28b (per reviewer):** the earlier draft of this section
erroneously marked BL_96 Glycosphingolipid (N=8) and BL_24 Starch/sucrose
(N=8) as "borderline" and outside the family. That contradicted §1
(`carry_both >= 5`). Both are in-family. Table above corrected; final
findings in `two_tier_v5_nonzero_gamma_bh.csv` includes BL_96
Glycosphingolipid metabolism at N=8, z=5.56, q=0.013 as a confirmed finding.

## 6. Amendments

Any modification of §1–§4 after the matched-floor and sweep results are read
must be logged with reason and date in this file, dated, and disclosed in
the response letter. No stealth changes.

## 7. Tier 0 finding: envelope-dominated (result addendum, 2026-09-28)

The replicate × threshold sweep and topology × envelope 2×2 decomposition
completed 2026-09-28 12:00. Sweep results are recorded here for provenance
because they inform how Tier 0 is framed in v5.

- **Replicate × threshold sweep** (per-replicate zFPKM at inex caps −4/−3/−2,
  4 cells × 3 reps × 3 thresholds = 36 extractions, extraction bounds only):
  every extraction reached biomass_max = 9.524. D_72 vs LD_72 gap = 0.000
  in 0/27 pairings. Confirms extraction topology sensitivity to labels is
  small.
- **Topology × envelope 2×2** (packet extracted models, sampling-style
  bounds = BOUND_EX × metabolomics fc + OCR cap):
  - D_72 vs LD_72 at 72 hpf: total gap = 1.515. Envelope effect = 1.515
    (100% of gap). Topology effect = 0.000.
  - BL_72 vs LD_72: gap = 0.812. Envelope 100%, topology 0.
  - D_48 vs LD_48: gap = 0.094. Envelope ~90%, topology ~10%.
  - BL_48 vs LD_48: gap = 0.101. Envelope ~95%, topology ~5%.

**Conclusion:** the packet's Tier 0 biomass_max differences are driven
almost entirely by the metabolomics-derived exchange envelope, not by the
extracted topology. In v5, Tier 0 is framed as "the model transmits a
measured metabolomics difference" and the next step is identifying which
specific exchange is binding via the LP dual (task #114). This does not
change §1–§4 (Tier 2 family, gamma tail, BH), which apply to sampled flux
distributions under the same envelope for both cells within a contrast.

## 8. Amendment 2026-09-29a: v6 extraction change (RFK/FLAD1/riboflavin_ex forced)

**What changed.** Three MAR IDs were added to the fastcore L1 "hardcoded
utilities" list before re-running extraction on all 15 cells:

- `MAR06506` — riboflavin kinase (RFK): riboflavin + ATP → FMN
- `MAR06508` — FAD synthase (FLAD1): FMN + ATP → FAD
- `MAR09143` — riboflavin exchange (extracellular)

New extractions live in `reviewer_packet_v4_slim/ribo_forced_models/`. v5
extractions in `ocr_anchored_extraction/models/` remain untouched.

**Why.** The v5 LP-dual analysis (task #114) identified `MAR01939` (FAD
exchange) as the sole binding constraint on biomass_max at 72 hpf. The
MAR-ID audit (`s121_mar_audit.csv`) showed that RFK and FLAD1 were both
dropped from every extracted model at 72, 96, and 120 hpf (partial keep
at 48; both present at 24). Without RFK/FLAD1 the cells cannot synthesise
FAD from riboflavin, forcing dependence on `MAR01939` and making
biomass_max a linear function of that one metabolomics fc. This is a
Tier 0 modelling artefact of the transcriptome-driven pruning, not
biology, and the reviewer (Q2 of 2026-09-28) directed us to force these
enzymes into the core and rerun.

**Motivated by Tier 0 diagnostic ONLY.** The change was proposed before
any Tier 2 v6 results existed. It is not motivated by any Tier 2
outcome. Same three reactions forced symmetrically across all 15 cells
— not selectively into any condition or hpf.

**Scope of impact.**
- §1 (family filter definition): unchanged.
- §2 (gamma-fit primary p for N<15): unchanged.
- §3 (BH-Hochberg over family, q≤0.10): unchanged.
- §4 (`is_finding` definition): unchanged.
- §5 (v4 audit table): historical, unchanged.
- §7 (Tier 0 envelope-dominated conclusion): now qualified — the envelope
  DOES dominate the biomass_max gap under v5 extractions, but part of
  that envelope effect was itself a modelling artefact of RFK/FLAD1
  dropout. Corrected biomass_max values under v6 sampling bounds are:
  D_72 = 8.91 (vs v5 7.08), LD_72 = 9.52 (vs v5 8.60), BL_72 = 8.12 (vs
  v5 7.78). Residual D_72 vs LD_72 gap = 6.4% (vs v5 18%). Under v6 the
  binding constraint on biomass_max is a different metabolite (identity
  pending v6 LP dual).

**Companion audit.** The v5-vs-v6 concordance column in
`two_tier_v6_gamma_bh.csv` (pending Grace) makes the effect of the
extraction change auditable per subsystem × contrast. Findings that
appear in both v5 and v6 are extraction-robust; findings that appear only
in v5 or only in v6 are extraction-dependent.

**Systematic follow-up.** Two-step audit. First pass —
`s128b_biomass_component_audit_fba.py` — enumerated MAR00021 substrates
(32) + MAR00022 substrates (16) and knocked out their exchanges. On v5 it
flagged MAM02750 (PI pool) essential in all 15 cells but missed FAD.
Post-mortem: the pool substrate is FADH2 (MAM01803), not FAD (MAM01802),
so a substrate-driven enumeration doesn't hit FAD's exchange. Replacement
—`s128c_exchange_essentiality.py` — knocks out **every exchange** in each
extracted cell (no precursor-list gating) and annotates each exchange as
(a) direct biomass substrate, (b) MAR00022 substrate, (c) member of a
curated chemical family (FAD/FADH2, FMN, riboflavin, NAD(P)/H). Any
exchange with ratio < 0.01 in ≥1 cell is an artefact candidate.
Cross-referenced with the v6 LP dual: if the new binding constraint at
D_72 hits this list, that's another instance of the same pattern and
warrants force-into-core treatment; if not, the residual 6.4% is a
genuine envelope difference and gets the 2×2 topology-envelope
decomposition plus per-replicate metabolomics check on the responsible
metabolite.

## 9. Sampler and pipeline notes (informational, no rule change)

- **PolyRound LP backend** was mixed across v6 preprocesses due to
  numerical failures of glpk on 3 tighter-floor polytopes (D_72, D_96,
  LD_120). Those 3 were preprocessed locally on Gurobi; the other 9 used
  glpk on Grace. Rounding affects sampler efficiency, not target
  distribution, but this is noted here for methodological transparency.
- **v6 biomass_lb rose substantially** in cells where biomass_max moved
  (e.g., D_72 from 6.37 → 8.02). Thinner polytopes → mixing may be
  harder. `s41d_aggregate` Rhat/ESS on v6 samples must be read before
  the two-tier output. A v6 split-chain neg control (analog of s120)
  will confirm the noise floor holds under the corrected pipeline.

## 10. Amendment 2026-10-01 — v7 freeze

The v6 flavin-forced extraction and its downstream analyses are retired.
This amendment freezes a new base (`baked_ocr_{24,48,72,96,120}_v7.mat`)
and new extractions (`v7_models/trans_rfastcormics_{cond}_{hpf}.mat`) as
the final model set for the chapter. Nothing downstream of the ledger or
the base is permitted to change without a failed gate — any new finding
goes in the limitations section, not the pipeline.

### 10.1 Curation ledger

`reviewer_packet_v4_slim/s130_v7_ledger.csv` is the authoritative rule.
75 rows across 3 classes × 4 actions:

| action | class 1 | class 2 | class 3 | total |
|---|---:|---:|---:|---:|
| keep_exchange            | 32 |  0 |  0 | 32 |
| medium_uptake            |  4 |  0 |  0 |  4 |
| ysl_lipid_delivery       |  6 |  0 |  0 |  6 |
| close_exchange           |  0 |  7 | 13 | 20 |
| none (no exchange in base) |  4 |  5 |  4 | 13 |

Class definitions:
- **Class 1** — true vitamins / yolk-derived essentials the embryo cannot
  synthesise; essential AAs; bulk inorganics; apolipoproteins and
  lipoprotein assemblies for yolk-syncytial lipid delivery.
- **Class 2** — the embryo can synthesise this from a class-1 precursor
  (pathway exists in Zebrafish1 with a GPR); close the exchange.
- **Class 3** — boundary artefact (pool metabolite, byproduct, xenobiotic);
  close the exchange.

Measurement override (per Q2 2026-09-29): any exchange in the 105 metabolomics
measured set is class 1 by convention, regardless of synthesis potential.

### 10.2 Pool coefficient ε = 1e-3

The pool pseudo-metabolite (MAM01602c) coefficient in biomass (MAR00021) was
changed from Wang's placeholder −1 to −ε = −1e-3. MAR00022 pool
stoichiometry otherwise unchanged (4 reduced-redox substrates — FADH2,
NADH, NADPH, ubiquinol — set to 0 per v7; the other 12 pool substrates
keep their Wang coefficients).

ε scan {1e-2, 1e-3, 1e-4} on all 15 extracted v7 models under
apply_context_bounds confirms biomass_max is flat across ε: ratios
1.0000–1.0001 for 12/15 cells; three cells show slight variation
(BL_48: 1.0144, D_48: 1.0145, BL_96: 1.0066). Pool coefficient does not
drive biomass_max. ε=1e-3 adopted.

### 10.3 Boundary rule (Part 1 + Part 2)

v7 base, per `s135_build_v7_base_full.m`, applied to all 5 hpf baselines:
- **Part 1 (secretion uncapped)**: for every exchange not in the 105, not
  in the ledger, not in the medium list, not O2: ub = 1000. Fixes the
  5-deoxyadenosine artefact and every analogous byproduct sink.
- **Part 2 (unmeasured organic uptake closed)**: same set: lb = 0.
- **Ledger closes**: lb = ub = 0 for the 20 close_exchange rows.
- **Medium (17 bulk inorganics / water / H+ / gases)**: lb = -1000,
  ub = 1000.

### 10.4 Free supply for unmeasured class-1 medium_uptake

Amendment 2026-10-01: in `apply_context_bounds`, lb = -1000 for every
`medium_uptake` ledger row not in the 105 (α-tocotrienol MAR09152,
γ-tocotrienol MAR09154, β-carotene MAR09276, inositol MAR09361).
Rationale: unmeasured class-1 essentials are yolk-derived and should not
default to the ±10 relative-envelope cap; otherwise an un-measured
vitamin becomes the LP's growth ceiling as an artefact of our choice
not to measure it.

### 10.5 Envelope convention BOUND_EX × fc

`apply_context_bounds` (Python + MATLAB twins, in `functions/`): for the
105 measured exchanges, lb = -BOUND_EX × fc, ub = +BOUND_EX × fc.
**BOUND_EX = 10 is a relative-units convention**; the envelope is a scale
centred on fc = 1.0. Only cross-cell ratios are physically interpretable
— the absolute "mmol/gDW/hr" magnitude is not calibrated against a
measured growth rate for v7. State this once in methods.

BOUND_EX scan {3, 10, 30} on all 15 extracted models: biomass_max
responds to the scale as expected (bex=3 → 2.6–6.1; bex=10 → 8.6–10.0;
bex=30 saturates at 8.6–10.0). 30/3 ratios 1.64–3.82. Cross-cell biomass
ratios at bex=10 (vs LD_72 reference) 0.86–1.00. Interpretable spread.

### 10.6 OCR unchanged

MAR09048 cap at −OCR(hpf), unchanged from v6:
24/48/72/96/120 hpf → 97/173/227/273/313 µmol O₂/g DW/hr.

### 10.7 Extraction convention

rFASTCORMICS ran under default ±10 on the 105 (no fc applied at
extraction time); fc envelope is applied post-extraction via
`apply_context_bounds` for every downstream analysis. This is a
documented deviation from applying fc at extraction time; it means the
extracted topology is calibrated to the maximum envelope the model will
see, so no fc-feasibility conflict arises during extraction.

### 10.8 v7 base validation (ATP sanity on D_72 base, MATLAB + Gurobi)

| test | expected | observed | gate |
|---|---|---:|---|
| 1. ATP from nothing, O₂ closed      | 0        | 0.000 | pass |
| 2. ATP from nothing, O₂ open        | 0        | 0.000 | pass |
| 3. aerobic yield, glucose = 1       | ≈30–32   | 31.500 | pass (≤32) |
| 4. anaerobic yield, glucose = 1     | ≤2       | 2.000 | pass |
| 5. full v7 medium, biomass unfixed  | report   | ATPD = 1000 (ub), O₂ = −226.63 | report-only |

### 10.8b ATP tests 1 & 3 per extracted v7 model (s141)

All 15 cells pass both gates (test 1 = 0, test 3 ≤ 31.5). Three 24 hpf
cells show test 3 = 30.75 vs 31.50 for later hpf. **This is an
extraction-topology difference, not an OCR effect** — test 3 is run with
O₂ unlimited (`lb = -1000` on MAR09048), so the per-hpf OCR cap does
not apply. The 24 hpf extractions have pruned a step whose later
counterparts recover the final 0.75 ATP/glucose.

| cell | test 1 | test 3 |
|---|---:|---:|
| BL_24 / D_24 / LD_24  | 0 | 30.75 |
| BL_48–120 | 0 | 31.50 |
| D_48–120  | 0 | 31.50 |
| LD_48–120 | 0 | 31.50 |

### 10.9 Model sizes v6 → v7

See `s147_sizes_v6_v7.csv`. v7 extractions have ~700–900 fewer reactions,
~450–630 fewer metabolites, 20–63 fewer genes per cell than v6. The
curation prunes close_exchange transports and their downstream orphans
plus the 4 zeroed pool substrates and their satellites, while retaining
nearly all gene-supported network.

### 10.10 Frozen configuration (no re-opens without a failed gate)

- Internal reaction default bounds: ±1000 (unchanged from Wang).
- BOUND_EX = 10 (relative-units convention).
- fc from `phase2_exchange_bounds_long.csv` applied post-extraction.
- Unmeasured class-1 medium_uptake rows: lb = -1000 (free supply).
- Biomass fraction floor: 0.9 × biomass_max (Frame A for within-condition;
  per-hpf matched floor for between-condition).
- zFPKM threshold: unchanged (per phase3_transcriptomics).
- Sampler settings: hopsy billiard walk, 4 chains per cell, Rhat/ESS gate.
- OCR mapping: 97/173/227/273/313 for 24/48/72/96/120 hpf.
- Pool coefficient ε = 1e-3 in MAR00021.
- Ledger (s130_v7_ledger.csv) is immutable.

### 10.11 Limitations

- **Phosphocholine (MAM02738, MAR09845)** — unmeasured; v7 base route is
  present but gapped at certain hpf. 7 cells show it as uptake-essential
  (see s148 report). Not re-opened; treated as a known gap. Growth
  shortfall in those cells is attributed to this gap in reporting.
- **5-deoxyadenosine (MAM01098, MAR09298)** — LIAS byproduct with no
  internal consumer in Zebrafish1; Part 1 secretion uncap (ub=1000) lets
  it exit. Appears as "secretion-essential" in the s128c gate split; not
  a curation failure.
- **Wang biomass coefficients** — the amino-acid coefficients and the
  11 remaining pool substrate coefficients are inherited from Wang. Not
  re-estimated against measured embryo composition for v7.
- **Unmeasured vitamins at free supply** — tocotrienols, β-carotene,
  inositol default to lb=-1000 as class-1 essentials. Any between-cell
  growth difference attributable to one of these is a modelling choice,
  not a measurement.
- **Extraction fc convention** — fc applied post-extraction. The
  extracted topology is calibrated to ±10 default, not to measured
  envelope.
- **Relative-units envelope** — only cross-cell biomass_max ratios at a
  fixed BOUND_EX are physically interpretable.

### 10.12 Launch gate (closed)

The v6 Grace hopsy chain and its afterok dependencies are to be
cancelled. v7 models are shipped. Enzyme subsets recomputed on v7.
PolyRound preprocessing (local, Gurobi) over all 15 cells. Sampling
chain with unchanged settings. v5→v7 concordance column added to the
findings table.
