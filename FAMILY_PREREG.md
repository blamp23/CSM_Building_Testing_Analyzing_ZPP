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

## 11. Amendment 2026-10-03 — v7 contrast families and findings pipeline

The v4 packet (§5) sampled 10 cell-pairs under the matched-floor regime
and had a single family. v7 sampled all 15 cells (local hopsy, Fix A for
LD_72) and runs a wider contrast panel that must be declared **before**
any q-value is read.

### 11.1 Family 1 — between-condition at each hpf (15 contrasts)

For each hpf in {24, 48, 72, 96, 120}, three contrasts:

- BL vs D
- BL vs LD
- D vs LD

Tests whether lighting condition changes subsystem flux distribution
at a given developmental stage.

### 11.2 Family 2 — within-condition developmental (12 contrasts)

For each condition in {BL, D, LD}, four adjacent-hpf contrasts:

- 24 vs 48
- 48 vs 72
- 72 vs 96
- 96 vs 120

Tests whether a subsystem's flux distribution shifts between adjacent
developmental stages within a lighting condition.

### 11.3 Negative control — split-chain within each cell (15 contrasts)

For each of the 15 cells, half-A (chains 1+2, 1000 samples) vs half-B
(chains 3+4, 1000 samples) through the identical two-tier pipeline.

- Reported **separately**, never pooled with Family 1 or Family 2 for
  multiplicity correction.
- Expected under the null: median |z| ≈ 0, 95th-pct |z| ≈ 2, zero
  subsystems reaching q ≤ 0.10 under Family-1/2 BH if the test is
  calibrated.
- Used to compute an empirical FDR floor and a per-subsystem
  neg-control-suspect flag (any split-chain |z| > 3 ⇒ flagged for Tier
  C in Families 1/2).

### 11.4 Per-row rules (unchanged from §1–§4, restated for v7)

1. Row enters the family iff `carry_both >= 5`, `loop_flagged /
   N_tier2 < 0.30`, subsystem ∉ {Exchange, Artificial, Transport}.
2. `p_primary = p_gamma` if `N_tier2 < 15`, else `p_emp`. Gamma MLE on
   the strictly positive W1 null draws; `p = 1 - Gamma_CDF(W1_obs;
   shape, scale)`. Rationale: §2, re-affirmed.
3. BH-FDR applied **globally within each family** (not per contrast).
   `q_bh <= 0.10` is a finding.
4. Reported alongside q_bh: `W1_obs`, `z`, `direction` (sign of mean
   log|v| difference), `asymmetric_zeros` flag (`|zero_frac_A −
   zero_frac_B| > 0.02`), `p_fisher_zeros`, `neg_control_suspect`.

### 11.5 Robustness tiers for Family 1 and Family 2 findings

Headlines come **only** from Tier A.

- **Tier A** — mid-N (`N_tier2 >= 15`), zero-symmetric (no asymmetric-
  zeros flag), not neg-control-suspect, appears in `>= 2` contrasts of
  the same type (adjacent hpf within a condition, or two conditions at
  one hpf) at q ≤ 0.10.
- **Tier B** — passes family + q ≤ 0.10 but single-contrast, or small
  N (`N_tier2 < 15` driven by p_gamma).
- **Tier C** — passes family + q ≤ 0.10 but flagged: asymmetric zeros,
  neg-control-suspect subsystem, or loop-dense.

### 11.6 Outputs

- `findings_family1.csv` — per-row Family 1 results (all 15 contrasts
  concatenated), with q_bh, tier, concordance vs v5.
- `findings_family2.csv` — same for Family 2 (12 contrasts).
- `negcontrol_summary.csv` — median |z|, 95th pct |z|, count at q ≤
  0.10 (empirical FDR estimate), per-subsystem neg-control-suspect
  table.

### 11.7 Concordance against v5

Each row in `findings_family1.csv` and `findings_family2.csv` carries
a concordance column against `reviewer_packet_v4_slim/
s124_findings_annotated.csv` (v5): `found_in = {v5_only, v7_only,
both}` plus `v5_z` and `v5_q` where applicable. Match key is
`(subsystem, contrast)`.

### 11.8 No stealth changes

Per §6, any further modification of §11 after `findings_family1.csv`
or `findings_family2.csv` is read must be logged here with reason and
date and disclosed in the response letter.

### 11.9 Amendment 2026-10-03 — "Pool reactions" excluded as artificial

"Pool reactions" (MAR00021 biomass precursor pool, MAR00022 cofactor
pool, and their derived pseudo-reactions in the Human-GEM / Zebrafish-
GEM lineage) are **pseudo-reactions used to assemble lumped substrate
pools**, not biochemical pathway members. They were not an explicit
exclusion in §1 but belong with Exchange / Artificial / Transport
on the same grounds: structural, not biological.

Applied retroactively. First-run findings had "Pool reactions" at
rank 10 in Family 1 (BL_120 vs D_120, N_tier2 = 5, z = 3.83,
q = 0.045). That row is dropped from Family 1; q-values are
re-computed by BH over the remaining 956 rows. Documented in
`findings_family1.csv` as `subsystem_excluded = True` for provenance.

### 11.10 Amendment 2026-10-03 — matched-floor resample at 96 hpf

Any q <= 0.10 hit at 96 hpf involving BL_96 is flagged
`floor_suspect = True` pending the matched-floor resample.

Rationale: BL_96 is isoleucine-bound at biomass_max = 8.61, while
D_96 = 12.3 and LD_96 = 12.7. Under Frame A, each cell's biomass
floor is 0.9 x bmax, so BL_96's floor is ~40 % below the other two.
The polytope shape differs by construction and this is the same
confound that produced the fake 72 hpf coordination signal in the
FVA diagnostic report.

Matched-floor resample launched 2026-10-03 for BL_96, D_96, LD_96
at biomass_lb = 7.749 (= 0.9 x 8.61) using s41d
`--biomass_lb 7.749` (Frame B). Three contrasts rerun on the
matched-floor samples: BL vs D, BL vs LD, D vs LD at 96 hpf.
`findings_family1_matched96.csv` written next to
`findings_family1.csv`. Any Family 1 or Family 2 claim that cites
BL_96 is held until the matched-floor q is read.

### 11.12 Amendment 2026-10-03 — median normaliser (sensitivity + rationale)

§2 as originally pre-registered uses the per-representative mean across
A and B as the normaliser for the log-ratio statistic:

    mu_r = (AA_r.mean() + AB_r.mean()) / 2
    log_ratio = log2(|v| / mu_r)   (with EPS=1e-30 added to |v|)

Sampled |v| distributions are right-skewed, so the mean sits above the
bulk. Every rep's log-ratio distribution is therefore centred somewhere
below zero by an amount proportional to its skew, and "0" on the x
axis no longer reads as "typical flux." The pooled subsystem
distribution is harder to interpret, and KDEs across many reps look
muddied.

Switching the normaliser to the **median of pooled non-zero values**:

    mu_r = median(|v| across A and B, where |v| >= 1e-10)
    log_ratio = log2(|v| / mu_r)   for samples with |v| >= 1e-10;
                                    zero samples dropped from the W1 input.

Rationale: median centres every rep at zero on the pooled scale, so
the subsystem KDE becomes an interpretable mixture and "0" reads as
"typical flux." Both A and B share the same per-rep normaliser, so
the W1 moves only a little; the ranking of hits should move less.

**Not a free swap:** this is a §2 amendment. The findings produced
under the mean normaliser (`findings_family1.csv`, `findings_family2.csv`)
remain the record of what the prereg actually prescribed. The median
rerun is reported as a **sensitivity analysis** alongside, with a
concordance table:

- Rows that are findings under **both** normalisers remain Tier A/B as
  classified.
- Rows that are **mean-only** or **median-only** findings are
  downgraded to **Tier C** regardless of z or N_tier2, with a
  `normaliser_disagree` flag.

Scripts:
- `grace_v4/scripts/two_tier_subsystem_test_fast_median.py`
- `grace_v4/scripts/s80_findings_v7_median.py`
- outputs → `results/stats_v7_median/`, `results/findings_v7_median/`
- concordance → `results/findings_v7/concordance_mean_vs_median.csv`

### 11.13 Amendment 2026-10-03 — SNR tiering + envelope-tracking flag

Two replacements of previous §11 tiering devices.

**(A) noise_floor_W1 + SNR replaces the |z| > 3 neg-control-suspect flag.**

For every subsystem we now compute:

    noise_floor_W1[sub] = 95th percentile of W1_obs across the 15
                          split-chain contrasts for that subsystem.

For every Family 1 / Family 2 row we add:

    snr = W1_obs / noise_floor_W1[sub]

This is a direct effect-size-above-method-noise ratio, instead of the
earlier |z| > 3 proxy. SNR replaces the `neg_control_suspect` field in
the tiering decision. Retained in the CSV as `snr` and
`noise_floor_W1` columns. New tier rule:

- **Tier A**: q ≤ 0.10 AND `snr ≥ 3` AND N_tier2 ≥ 15 AND NOT
  asymmetric_zeros AND NOT envelope_tracking AND `n_sig_contrasts ≥ 2`
- **Tier B**: q ≤ 0.10 AND `snr ≥ 2` AND NOT asymmetric_zeros
  AND NOT envelope_tracking (and doesn't satisfy Tier A)
- **Tier C**: q ≤ 0.10 AND (`snr < 2` OR asymmetric_zeros OR
  envelope_tracking)

Rationale: a hit that's only 1-1.5x above the test's own split-chain
noise floor is not a biology claim; it's the floor itself. 3× is a
conservative SNR threshold that holds across the n_tier2 spectrum
(small-N and large-N subsystems both get normalised to their own
split-chain scale).

**(B) envelope_tracking flag.**

For every Tier A/B row we now report:

    ex_fc_summary = list of {exchange : fc_A : fc_B : fc_ratio}
                    where exchange is any ex_rxn whose single
                    metabolite appears in any reaction of the
                    subsystem (either cell's extracted model).
    envelope_tracking = True iff
        (direction = +1 AND any fc_A / fc_B >= 2) OR
        (direction = -1 AND any fc_A / fc_B <= 0.5)

Rationale: v7 scales each measured exchange's envelope by the
cell-specific fc multiplier. If a subsystem's reactions are tightly
coupled to a measured exchange whose fc ratio across the two cells
exceeds 2, the model's directional "finding" may just be the measured
envelope difference reasserting itself through the FBA LP, not an
emergent topology effect. Flagging as Tier C.

### 11.14 Amendment 2026-10-05 — pre-reg §2 zero handling corrected

§2 as written pre-registers the W1 statistic as `log₂(|v| / mu + ε)`
with `ε = 1e-30`, so samples with `|v| = 0` enter the test at
`log₂(ε) = −99.66`. This was carried over from an early draft and
contradicts the v5 convention, where non-zero W1 is the method of
record and a separate per-row Fisher-exact on zero counts (§4
`p_fisher_zeros`) handles the zero mass as a distinct test. The
reviewer's 2026-10-05 read correctly identified that the 18 → 63 gap
between the §2 pre-reg run (mean + EPS-shift) and the §11.12
sensitivity run (median + drop zeros) is almost entirely explained by
zero handling, **not** by the mean-vs-median normaliser.

Correction: non-zero W1 is the method of record. Samples with `|v|
< 1e-10` are dropped from the W1 input and the zero fraction per
condition is reported alongside with Fisher-exact on the 2×2 table
(zeros vs non-zeros × A vs B). This matches v5.

To separate the two effects empirically, a third pipeline is run:
mean normaliser **with zeros dropped**. Scripts:

- `grace_v4/scripts/two_tier_subsystem_test_fast_mean_dropzero.py`
- `grace_v4/scripts/s80_findings_v7_mean_dropzero.py`
- outputs → `results/stats_v7_mean_dropzero/`,
  `results/findings_v7_mean_dropzero/`

The three variants are:

| variant | mu_r | zeros | role |
|---|---|---|---|
| `findings_v7/` | pooled mean of A, B | EPS-shift to −100 | pre-reg §2 as written (artifact) |
| `findings_v7_median/` | pooled median non-zero | dropped | §11.12 sensitivity |
| `findings_v7_mean_dropzero/` | pooled mean non-zero | dropped | **method of record (v5)** |

Three-way concordance written to
`results/findings_v7/concordance_three_way.csv`. Rows surviving as
Tier A/B in the method-of-record pipeline are the ones carried into
the response. The other two variants are reported as sensitivities
with the amendment trail.

### 11.15 Amendment 2026-10-05 — positive control (corrected 2026-10-05)

To demonstrate pipeline sensitivity (not only specificity, which the
split-chain negative control already shows), two in-silico
perturbations are applied to LD_72 (clean under all §11 flags, Fix A
applied for MVE rounding).

**Target subsystem**: Steroid metabolism (not on any current hit
list; **28 reactions in the LD_72 extracted model**, mid-sized, not
biomass-essential — the biomass_max under raw extraction bounds is
10.02 before and 10.02 after a full KO, so the perturbation is
non-lethal; the fc-envelope biomass_max for LD_72 is 12.85, which is
the number the sampler sees. The 10.02 figure is logged here only to
document that the KO itself does not collapse growth.).

**(A) Full knockout — Tier 1 (zero-fraction Fisher) validator.**

- lb = ub = 0 for every reaction in Steroid metabolism in the LD_72
  Fix A tightened model. Re-run FVA tightening on the resulting
  model (steroid-zero may make neighbouring reactions' prior tight
  bounds infeasible) to recover a sampling-ready polytope.
- Resample 4 chains, 500 effective samples each, thinning 300.
- Contrast LD_72 vs LD_72_KO_steroid through the method-of-record
  pipeline.
- **Correct expected signature** (the earlier draft of this
  amendment pre-registered a Tier 2 W1 hit here, which cannot
  happen by construction: with lb = ub = 0 every steroid rep has
  zero flux in the KO cell, carryA = 0, so the subsystem fails the
  §1 `carry_both ≥ 5` filter and never enters Tier 2):
  1. **Tier 1 (p_fisher_zeros)** hit on Steroid metabolism — the
     zero fraction in the KO cell is 1.0 vs ~0.0 in the original.
  2. **Tier 2 (W1)** hits in metabolically connected subsystems
     (sterol transport, cholesterol biosynthesis, bile acid
     biosynthesis — whatever has steroid metabolism's substrates
     or products in its reactions) via the LP redistributing flux
     away from the newly-blocked pathway.
  3. Steroid metabolism itself appears in the Tier 2 output with
     `carry_both < 5 → excluded` and `p_fisher_zeros` very small —
     documented as filtered-out, not interpreted.

**(B) Knockdown — Tier 2 (W1) validator.**

- For each reaction in Steroid metabolism, under LD_72 Frame A
  bounds: run FVA to get `[FVA_min, FVA_max]`. Cap the reaction's
  envelope to **10 % of the Frame A FVA range, preserving sign**:

      new_ub = 0.1 × max(FVA_max, 0)
      new_lb = 0.1 × min(FVA_min, 0)

  (so a forward-only reaction keeps its forward direction but is
  capped at 10 %; a reversible reaction keeps its sign but has both
  ends tightened to 10 %). Then re-run full-model FVA tightening on
  the knockdown model.
- Resample 4 chains, 500 samples, thinning 300.
- Contrast LD_72 vs LD_72_KD_steroid.
- **Expected signature**: Steroid metabolism itself fires as a Tier
  2 W1 hit with `snr ≫ 3`, direction = ↓, and the subsystem passes
  the §1 filter (reactions still carry non-zero flux at 10 %
  magnitude). Downstream subsystems may also fire but less
  dramatically than under full KO.

**(B) is the Tier 2 validator — the finding machinery most of the
response relies on.** It runs first.

Outputs:
- `results/stats_v7_positivecontrol/LD_72_vs_LD_72_KO_steroid.csv`
- `results/stats_v7_positivecontrol/LD_72_vs_LD_72_KD_steroid.csv`
- `results/findings_v7_positivecontrol/report.txt` (both)

If the knockdown does not fire a Tier 2 W1 hit on Steroid metabolism
itself with snr ≫ 3, the method's Tier 2 sensitivity floor is above
a 10× flux magnitude shift and the "few findings" reading of the
real-data run is uninformative; both results — the method's
sensitivity ceiling and the finding set — are reported in that case
rather than suppressed.

### 11.16 Amendment 2026-10-05 — KD positive control redesign

The §11.15(B) knockdown as originally spec'd capped each Steroid
metabolism reaction at 10% of its Frame A FVA max. Ran, test fired
`W1_obs = 0.015, z = -0.26, p = 0.58` on Steroid metabolism itself —
no signal. Diagnosis after inspecting KD vs original sampled fluxes:

    rxn        KD_cap_ub   mean |v| in sampled   headroom
    MAR01644   2.80        0.049                 57× above cap
    MAR01924   100.0       0.294                340×
    MAR04766   69.25       8.2e-05               ~800,000×

The caps are correctly applied to the mat file, but FVA_max is the
absolute maximum each reaction *could* achieve in principle, while
actual sampled fluxes live in the polytope *interior* far below those
extremes. 10% of FVA_max is still enormously larger than the typical
sampled magnitude, so the caps never bind and the KD polytope is
functionally the same as the original.

Methodological finding: for flux-sampling pipelines, FVA-max-based
knockdown is the wrong operationalisation. The correct cap must be
tied to the *sampled* distribution, not the polytope's absolute bounds.

**(B2) Revised Tier 2 knockdown.** For each Steroid metabolism
reaction, read the original LD_72 pooled samples and compute
`max_abs_v = max(|v_r|)` across all 2000 samples. Set:

    new_ub = +0.1 * max_abs_v   if any sampled v_r > 0 else 0
    new_lb = -0.1 * max_abs_v   if any sampled v_r < 0 else 0

This caps each steroid reaction at 10% of the maximum magnitude it
actually achieved in sampling of the original polytope, preserving
sign direction. Re-FVA tighten; sample 4 chains. Scripts:

- `grace_v4/scripts/s82b_build_kd_percentile.py`
- outputs → `reviewer_packet_v4_slim/v7_models_kd2_steroid/`,
  `results/hopsy_v7_kd2_steroid/`

Pre-registered expected signature (unchanged from §11.15(B)): Steroid
metabolism itself fires at Tier 2 W1 with snr ≫ 3, direction = ↓.
If (B2) also fails to fire, the pipeline's Tier 2 sensitivity ceiling
is above a per-reaction 10× flux magnitude reduction — reported in
the response alongside the finding set.

### 11.11 Known Frame A limitation — retinol at 48->72 across all conditions

Retinol metabolism fires at q <= 0.10 at the 48->72 transition in
BL, D, AND LD (Family 2, all Tier C, all neg-control-suspect, all
direction = decreasing). This is the biomass-coupled retinoid subset
(bco1-dependent in D; rbp4/stra6-dependent in BL and LD) reacting
to Frame A's per-cell biomass floor shifting between hpf — the
floor goes up with developmental stage because bmax goes up, and
the retinoid subset is tightly coupled to biomass_lb through the
vitamin A derivatives pool (MAM03139c, previously identified as
D's biomass bottleneck in `project_vitamin_a_finding.md`). This is
reported as a Frame A limitation in the response letter, not as a
developmental finding.


### 11.17 Amendment 2026-10-07 — Tier 0 topology × envelope decomposition under v7

The v6-era "envelope ≈ 100% of parent→child gap" statement is superseded.
Rerun under the final v7 bounds (apply_context_bounds; Frame A off; no biomass
floor; envelope + OCR + medium + ledger only) now gives a pair-level verdict.

**Scope.** `grace_v4/scripts/s150_tier0_2x2_v7.py` + `s151_tier0_fixups_bl96.py`.
Plain FBA via Gurobi, no sampling.

**Part 1 — parent vs child, per cell (15 + 15 LPs).**
`tier0_parent_child_v7.csv` with top-3 exchange constraints by |dual × bound|
per cell (lb = 0 irreversibility rows excluded). Parent biomass_max values
match the previously published `tier0_biomass_max.csv` (26.6, 30.6, 31.2,
31.2, 32.1, …) — the parent and child columns were computed under the same
bound set. Across cells the parent is pinned on the three plasma-protein
exchanges MAR04117 (fibrinogen), MAR04129 (plasminogen), MAR04138
(prothrombin) with duals ≈ 1.4×10⁵ each — a consequence of the small biomass
stoichiometric coefficient on these proteins, not biology. Child models are
pinned on much smaller contributions (10²–10⁴), **except BL_96** whose child
is pinned on MAR09052 (chylomicron-remnant exchange) with a dual contribution
of 2.70×10⁵ — ≈ 50× the next-highest cell. This is the Tier 0 signature of
BL_96's bottleneck.

*Footnote on MAR10693 (added 2026-10-13, §11.20 addendum).* The gene-less
Trp-Gly-Leu secretion exchange MAR10693 appears as the top-1 or top-2
ub-side binder (ub = 1000, Part 1 secretion default) in 9 of 15 cells with
|dual × bound| contributions 100–840. Setting MAR10693.ub = 0 and
re-solving FBA under `apply_context_bounds` changes biomass_max by **0.00%
in all 15 cells** (`grace_v4/results/v7_final_addenda/mar10693_biomax_table.csv`).
A constraint that carries the largest dual contribution and whose removal
changes the optimum by exactly zero is the signature of **dual degeneracy**:
there is a family of gene-less secretion stubs (the chylomicron / VLDL
remnant exchanges and the small-peptide exchanges that the ATP session
enumerated as adjacent top-duals, e.g. MAR09052, MAR10709, MAR10188,
MAR10543) that are interchangeable vents in this LP, and the solver reports
whichever one the basis happens to land on. The underlying topology in e.g.
LD_24 is: intracellular condensation MAR11177 (`tryptophan[c] +
Glycylleucine[c] <=> H2O[c] + TrpGlyLeu[c]`) → H⁺-symport transporter
MAR10935 (`H+[e] + TrpGlyLeu[e] <=> H+[c] + TrpGlyLeu[c]`) → secretion
exchange MAR10693 (`TrpGlyLeu[e] -> ∅`). MAR10693 is one of a degenerate
set of gene-less secretion stubs: closing any single member is absorbed by
the others; **none affects biomass_max**. Tier 0 duals that cite these
stubs should be read as pointers that found the degenerate set, not as
operating costs. See the Tier 0 robustness perturbation table below for
what biomass_max actually responds to.

*Tier 0 robustness — 60 LPs (added 2026-10-13).* Because the top-ub duals
are degenerate (above), the "nitrogen-bound" characterisation of Tier 0
rests on a non-degenerate test: relaxing an uptake bound by 10% and asking
whether biomass_max actually moves. For each of 15 cells × 4 measured
exchanges (valine MAR09046, histidine MAR09038, isoleucine MAR09039,
glucose MAR09034 as a carbon-side control), FBA under
`apply_context_bounds`, then `lb(exchange) := 1.10 × lb(exchange)` and
re-solve (`grace_v4/results/v7_final_addenda/tier0_robustness.csv`):

- **Glucose (carbon control):** max |Δ| = 0.07 % across 15 cells; nowhere
  near biomass-binding. Confirms the Tier 0 capacity is not carbon-limited.
- **Valine:** max |Δ| = 0.007 %, zero in every cell. Does not bind.
- **Histidine:** moves biomass by **4.4–5.5 %** at BL_48, D_48, D_72, LD_48,
  LD_72 and **+2.4 %** at LD_96 — a mid-development histidine ceiling.
- **Isoleucine:** moves biomass by **+8.2 %** at BL_96 and **+0.7 %** at
  BL_120. BL_96 is isoleucine-bound specifically; nowhere else does
  isoleucine bind materially.

Honest Tier 0 characterisation: the capacity ceiling is nitrogen-limited,
driven by **histidine in the hpf 48–72 window across all three conditions**
(LD_96 included) and by **isoleucine specifically at BL_96** (consistent
with §11.17's finding that BL_96 is a topology-dominated outlier). Valine
is not a Tier 0 bottleneck under v7 bounds. The duals pointed to a
degenerate secretion family and were the diagnostic entry point; the
perturbation table is what the operating ceiling claim rests on.

**Part 2 — 2×2 topology × envelope swap, per pair (15 pairs × 4 corners = 60
LPs).** `tier0_2x2_v7.csv` reports the four corners in absolute biomass units
plus the topology effect, envelope effect, interaction, and total diagonal
gap. Shares are reported **only** where `|total_gap| > 5% of the smaller
corner`; otherwise the row is marked `no material gap`. No per-hpf
averaging. Verdict per pair:

| pair | gap (hr⁻¹) | top effect | env effect | dominance |
|---|---|---|---|---|
| BL_24 vs {D,LD}_24, D_24 vs LD_24 | ±0.04 to 0.09 | — | — | no material gap (gap < 1%) |
| BL_48 vs {D,LD}_48, D_48 vs LD_48 | ±0.00 to 0.09 | — | — | no material gap (gap < 1%) |
| BL_72 vs D_72 | +1.57 | 0 | +1.59 | envelope (101%) |
| BL_72 vs LD_72 | +0.23 | — | — | no material gap (1.8%) |
| D_72 vs LD_72 | –1.35 | 0 | –1.37 | envelope (102%) |
| **BL_96 vs D_96** | **–4.04** | **–3.43** | –0.61 | **topology (85%)** |
| **BL_96 vs LD_96** | **–4.39** | **–4.06** | –0.33 | **topology (92%)** |
| D_96 vs LD_96 | –0.36 | — | — | no material gap (2.8%) |
| BL_120 vs D_120 | –0.11 | — | — | no material gap (0.9%) |
| BL_120 vs LD_120 | –0.84 | –0.06 | –0.78 | envelope (93%) |
| D_120 vs LD_120 | –0.73 | +0.04 | –0.77 | envelope (106%) |

**9 of 15 pairs** carry no material gap. Of the **6 material pairs**, 4 are
envelope-dominated (72_B vs D, 72_D vs LD, 120_BL vs LD, 120_D vs LD) and
2 are topology-dominated (both BL_96 pairs). `n_unswappable_exchanges = 0`
for every pair — every measured exchange appears in both members'
extractions, so the envelope swap is clean.

**Part 3 — BL_96 drill-down.** `bl96_lacking_reactions.csv` lists the
54 transport reactions and 43 exchange/demand reactions in LD_96 but not
BL_96, with substrate/product metabolite names and compartments. A greedy
single-reaction add-back from BL_96's 8.535 baseline toward the 12.255
all-173-extras ceiling (`bl96_greedy_addback.csv`) hits a plateau at **88.6%
recovery after four meaningful picks** (steps 5–6 add < 10⁻⁵, i.e. degenerate
alternates):

| step | rxn | subsystem | biomass_max | fraction of recovery |
|---|---|---|---|---|
| 1 | **MAR00466** | Transport reactions | 11.737 | **86.1%** |
| 2 | MAR11261 | Transport reactions | 11.792 | 87.5% |
| 3 | MAR00796 | Transport reactions | 11.805 | 87.9% |
| 4 | MAR13032 | Exchange/demand reactions | 11.832 | 88.6% |

The single most-restoring reaction is **MAR00466: PS pool (phosphatidylserine)
transport from Cytosol to Endoplasmic reticulum** — one reaction alone takes
BL_96 from 8.5 to 11.7, 86% of the available recovery. The next picks are
hexanoyl-carnitine import (MAR11261), sphingomyelin pool
Cytosol → Lysosome transport (MAR00796), and TAG-chylomicron exchange
(MAR13032). All four picks are **phospholipid / sphingolipid / fatty-acid
trafficking reactions**. The remaining 11.4% of recovery requires a coupled
multi-reaction add that greedy cannot find in a single step.

**zFPKM check — biology vs extraction.** `bl96_greedy_gene_zfpkm.csv`
confirms that **all six picked reactions have an empty GPR** (no gene
assignment in Human-GEM / Wang). The direct question "are these genes low in
BL_96 RNA-seq?" collapses: there are no genes to check. FASTCORE's decision
to prune these transporters is therefore a **topology cascade from adjacent
GPR-mapped reactions** that were not activated by BL_96's expression. The
clean verdict: **the BL_96 biomass gap is extraction, not direct expression
of these transporters**. Identifying which adjacent lipid-handling GPR-mapped
reactions FASTCORE dropped in BL_96 (vs retained in LD_96) is the next drill
level and belongs to a separate analysis — the Tier 0 scope stops at
"topology-dominated, driven by PS/SM/FA trafficking reactions that FASTCORE
pruned."

**Supersession.** The v6-era "envelope ≈ 100%" claim is retracted. v7 shows
envelope dominates at hpf 24/48/72/120 (where there is material gap), but at
hpf 96 the gap is topology-dominated for both BL-involving pairs. This is a
specific, reproducible finding and consistent with the subsystem-level v7
sampling story that singled out BL_96 as the outlier cell.

**Partial supersession (2026-10-12, §11.19).** The β-carotene / inositol
free-supplied-at-`lb=-1000` bug uncovered in §11.19 inflates the apparent
topology effect in both BL_96 pairs. Under the §11.19-corrected bounds
(`medium_uptake_demand_bounds.csv`), the BL_96 pair-level numbers become:

| pair | topology (OLD) | topology (NEW) | envelope (unchanged) |
|---|---:|---:|---:|
| **BL_96 vs D_96** | −3.43 | **−0.887** | −0.717 |
| **BL_96 vs LD_96** | −4.06 | **−0.888** | +0.084 |

The topology component shrinks by ≈ 75% in both pairs. The envelope component
is unchanged. Under corrected bounds the pair-level verdict changes:

- **BL_96 vs D_96** total gap ≈ −1.6, topology ≈ −0.89, envelope ≈ −0.72:
  roughly 55/45 topology/envelope. **No longer a cleanly topology-dominated
  pair** — the correction converts it to a mixed case.
- **BL_96 vs LD_96** total gap ≈ −0.80, topology ≈ −0.89, envelope ≈ +0.08:
  topology still dominates (≈ 100% of a now-smaller gap), but the gap itself
  is five times smaller than previously reported.

The PS-transporter add-back finding stands qualitatively but is restated
against the new gap. BL_96 biomass_max is unchanged under correction
(8.535 → 8.535); the ceiling (BL + all 173 missing reactions, 12.255) is
unchanged; the four-pick 88.6%-recovery curve is unchanged
(`bl96_greedy_addback.csv` reruns identically because the add-back model
does not involve the 4 free-supplied class-1 exchanges on its restoring
reactions). What **does** change is the interpretation of what BL_96's
3.72-unit gap *costs on the pair contrast*: previously the pair-level gap
to D_96 and LD_96 was ≈ −4 units of apparent topology, now it is ≈ −0.9.
Most of the gap-of-interest in the pair comparison was β-carotene
inflation of the non-BL_96 member of the pair, which under the old bounds
was running a retinoid-fed biomass boost that BL_96 did not match.

The pre-correction tables remain in `grace_v4/results/tier0_v7/tier0_2x2_v7.csv`;
the corrected tables are in `grace_v4/results/atp_sweep_v7/tier0_C3_2x2_before_after.csv`.
Any Tier A/B/C row in Vitamin A metabolism, Retinol metabolism, or
Carotenoid-touching subsystems carries a `pre_correction_retinoid` flag
pending the LD_72 spot-resample contrast (§11.19).

Outputs:
- `grace_v4/results/tier0_v7/tier0_parent_child_v7.csv` (15 cells, top-3 duals)
- `grace_v4/results/tier0_v7/tier0_2x2_v7.csv` (15 pairs, 4 corners, decomposition, dominance; pre-§11.19)
- `grace_v4/results/atp_sweep_v7/tier0_C3_2x2_before_after.csv` (post-§11.19, 15 pairs, old vs new)
- `grace_v4/results/tier0_v7/bl96_lacking_reactions.csv` (97 missing T + E reactions, mets, compartments)
- `grace_v4/results/tier0_v7/bl96_greedy_addback.csv` (6-step greedy, 88.6% recovery — unchanged by §11.19)
- `grace_v4/results/tier0_v7/bl96_greedy_gene_zfpkm.csv` (no GPRs → extraction, not expression)

### 11.20 Amendment 2026-10-13 — F2 matched-floor rule + resample

**Rule (prospective, applied in `s85_export_tier_tables.py`):**
An F2 Tier A row requires `p_matched_floor ≤ 0.05` when
`max(child_biomax_A, child_biomax_B) / min(…) > 1.15`. Otherwise it is
demoted to Tier C with `flag_reason = "floor_artifact"`.

**Why.** Within-condition F2 contrasts frequently compare adjacent-hpf cells
whose biomass_max differs enough to change polytope SHAPE independently of
biology. Across all 15 cells, biomass_max drops from ~11 at 24 hpf to ~9.4 at
48 hpf (the Tier 0 BL/D/LD "u-shape" minimum), then recovers to 11–13 by
72 hpf. Any 24↔48 or 48↔72 contrast therefore crosses that 11 → 9.4 drop and
can carry a 1.17–1.38× polytope-size ratio. In §11.10 we flagged this but
only resampled BL_96 (F1). §11.20 extends matched-floor resampling to the F2
axis.

**Resample.** Grace job 20053143 (launched 2026-10-12 10:41 CDT, completed
2026-10-12 ~14:50 CDT after D_72 fix via `starting_point=` patch in
`s41d_hopsy_sample_direct.py`). 9 cells × 4 chains × 500 samples at matched
biomass_lb per condition (BL=8.398, D=8.480, LD=8.479). CFF locally; pooled;
six F2 contrasts via `two_tier_subsystem_test_fast_mean_dropzero.py` →
`results/stats_v7_matched_f2/`.

**Note on bounds used.** The 20053143 resample sampled from polytope caches
built by the 2026-10-06 local preprocess, which called
`apply_context_bounds` BEFORE §11.19 added demand bounds (2026-10-12 17:29).
The resample therefore uses pre-§11.19 bounds — the same bounds as the main
v7 sample set, which is the correct design: like-for-like matched-floor
comparison, not a defect. The LD_72 spot check in §11.19 is what licenses
the whole pre-correction sample set. The §11.20 verdict applies to
polytope-shape confounds under the pre-§11.19 ground truth; the retinoid
family is flagged via the new `pre_correction_retinoid` column (true iff
subsystem ∈ {Retinol metabolism, Vitamin A metabolism, Carotenoid
biosynthesis}).

**Verdict (22 flagged F2 Tier A rows, split 10 / 12):**

Surviving (10): Arachidonic acid (BL 48→72, LD 48→72); FA activation
(cytosolic) at D 24→48, D 48→72, LD 24→48, LD 48→72 — all four FA-activation
rows survive; Fatty acid oxidation (BL 48→72, LD 48→72); Leukotriene
(D 48→72, LD 48→72).

Demoted to Tier C `floor_artifact` (12): Arachidonic acid at D 48→72;
Aromatic AA (phenylalanine/tyrosine/tryptophan biosynthesis) across
BL/D/LD 24→48 — all three; BCAA (valine/leucine/isoleucine) across BL/D/LD
24→48 — all three; Fatty acid oxidation at D 48→72; Leukotriene at BL 24→48
and BL 48→72; Nucleotide at LD 24→48; Pyrimidine at D 48→72.

**Headline impact.** Tier A 62 → 50; Tier C 21 → 33 (12 added with
`flag_reason = "floor_artifact"`, 0 of which are retinoid). The amino-acid
24→48 story (Aromatic AA + BCAA, 6 of the 12 demotions) was a polytope-shape
artefact of the 24→48 biomass_max drop. The FA-handling / eicosanoid
developmental signal is robust under matched floors.

**Columns added by s85.** `W1_matched_floor`, `p_matched_floor` (NaN unless
F2 row with a matched-floor result); `pre_correction_retinoid` (bool, set
from the subsystem list — survives any in-place file edit by downstream
scripts); `flag_reason` (Tier C only).

**Outputs.** `grace_v4/results/stats_v7_matched_f2/*.csv` (6 contrasts),
`grace_v4/results/findings_v7_mean_dropzero/{tierA,tierB,tierC}_final.csv`
+ `all_rows_in_family.csv` (rebuilt), figures in
`grace_v4/results/figures_v7/`.
