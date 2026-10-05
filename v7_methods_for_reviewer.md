# v7 Analysis Methods and Current Findings

A methods document accompanying the v7 revision. Written to be readable
on its own; references `FAMILY_PREREG.md` for pre-registered statistical
rules and `base_model/README.md` for v7 model construction.

Primary question: do subsystem-level flux distributions differ (a)
between lighting conditions BL / D / LD at a given developmental stage,
and (b) between adjacent developmental stages within a lighting
condition. We approach this with context-specific genome-scale metabolic
models, flux sampling of the feasible polytope, loopless post-processing,
and a two-tier Wasserstein statistic that is pre-registered and audited
by multiple sensitivity analyses before any claim is reported.

Nothing in this document was decided after reading q-values; every
analytical choice was logged in `FAMILY_PREREG.md` before the stats were
read, and dated amendments were added when choices had to change
mid-analysis.

---

## 1. The dataset

Zebrafish embryos sampled at 24, 48, 72, 96, 120 hpf under three
lighting regimes: BL (constant blue), D (constant dark), LD (12 h / 12 h
light-dark cycle; **LD is the physiological reference**, D and BL are
experimental deviations). This gives 15 cell contexts.

Per cell we extract a context-specific GEM from the Wang Zebrafish-GEM
using rFASTCORMICS gated on per-cell RNA-seq (zFPKM). All 15 cells share
the same base reaction/metabolite namespace; what differs across cells
is which reactions are active, which exchanges carry flux at what
magnitude, and the biomass_max achievable.

---

## 2. v7 model construction (summary; see `base_model/README.md`)

The v7 "frozen" base is defined by a single tunable scalar ε = 1e-3
applied to the pool-assembly pseudo-reaction (MAR00021 → MAM01602c
coefficient −ε rather than Wang's −1). Every other model-level choice —
the curation ledger (`v7_ledger.csv`), the Part 1 / Part 2 boundary
rule, the medium whitelist, the OCR cap, and the BOUND_EX × fc envelope
convention — is predetermined by the ledger, with a one-time §10
amendment logged for the whole v7 freeze.

Each of the 15 cells is extracted independently with rFASTCORMICS and
"baked" to the v7 ledger. Per-cell biomass_max under Frame A ranges from
≈ 8.5 (BL_96, isoleucine-bound) to ≈ 12.9 (LD_96). The envelope on each
measured exchange is `[−BOUND_EX·fc, +BOUND_EX·fc]`; unmeasured
exchanges inherit the ledger's class-1/2/3 defaults.

All 15 extractions pass four ATP sanity checks and a parent-model
sanity gate.

---

## 3. Sampling

For every cell we build the Frame A polytope (apply_context_bounds with
biomass_lb = 0.9 × bmax_cell, OCR cap on MAR09048, measured-exchange
envelope, medium whitelist). We preprocess with PolyRound (simplify →
transform → maximum-volume ellipsoid rounding) using Gurobi as the LP
back-end. We then sample with hopsy's billiard walk, four independent
chains of 500 effective samples at thinning 300, seeded identically
across cells. Total: 60 chains, 15 cells × 4 chains × 500 samples each.

**LD_72 polytope conditioning.** The LD_72 Frame A polytope has 16
reactions in the retinoid / vitamin D / vitamin E network with ranges
of 10⁻⁶–10⁻⁴ around a biomass-coupled solution; this created slit
directions MVE rounding could not resolve. "Fix A" (`s150_fix_A_ld72.py`)
tightens every reaction's bound to [FVA_min − 1e-4, FVA_max + 1e-4]
under biomass_lb = 11.21 before PolyRound; the samples are still from
the Frame A feasible region but on a numerically-tractable polytope.
Rounding converged in 79 minutes. Zero LP failures across all four
LD_72 chains.

---

## 4. Loop removal

Flux-balance sampling can accept thermodynamically infeasible internal
loops. We post-process every sample with CycleFreeFlux (Desouki et al.
2015): per sample, pin exchanges and biomass to the sampler's value,
minimise the L1 norm of internal reactions subject to steady state and
sign-preserving internal bounds. We implemented this with Gurobi
(`s51b_cyclefreeflux_gurobi.py`); each sample is one LP, warm-started
across adjacent samples via Gurobi dual simplex.

Loop removal rate: 87.75 – 90.30 % of internal L1 removed, consistent
across all 15 cells. Zero LP failures across all 30 000 samples. We
verified bit-equivalent objective value against the scipy-HiGHS
reference implementation on BL_24 chain 1 (relative difference 3 × 10⁻⁹
on total L1); where the two solvers disagree they do so on which vertex
of the optimal L1 face they land on, not on the L1 value — expected for
L1-min with flat optimal faces.

---

## 5. The two-tier Wasserstein statistic

For a contrast between cell A and cell B, for each enzyme-subset
representative r, we compute per-sample log-ratio

    log₂( |v_r| / mu_r )

where mu_r is a per-representative normaliser. Under the pre-registered
spec (§2) this is `mu_r = (mean|v_r|_A + mean|v_r|_B) / 2` and samples
with `|v_r| < 1e-10` are shifted to `|v| + ε` (ε = 1e-30) rather than
dropped. Under the §11.12 sensitivity amendment it is `mu_r = median`
over non-zero values pooled across A and B, with zeros dropped.

For a subsystem S passing the §1 filter, we pool the subsystem's reps'
log-ratios in cell A and cell B separately and compute the
Wasserstein-1 distance between the two empirical distributions over 201
quantiles:

    W1_obs = mean | q_A(p) − q_B(p) |  for p ∈ [0, 1] with 201 grid points.

The null distribution is 2000 random subsample draws of the same size
from the subsystem-neutral background (`gpr & ¬exchange & ¬loop &
carry_both`). For N_tier2 ≥ 15 we take `p_emp = (1 + # draws ≥ W1_obs)
/ (ndraw + 1)`. For N_tier2 < 15 we fit a gamma MLE to the strictly
positive null draws and use `p = 1 − Γ_cdf(W1_obs; shape, scale)` — the
reviewer's 2026-09-28 small-N tail correction, pre-registered at §2.

---

## 6. §1 filter (what enters the family)

A subsystem contrast enters the Tier-2 family **before** any q-value is
computed iff:

- `carry_both ≥ 5`: at least 5 reps carry non-zero flux (median > 1e-4)
  in both cells.
- `loop_flagged / N_tier2 < 0.30`: fewer than 30 % of carrying reps had
  `max|v| > 100` at any point in their 2000 samples.
- Subsystem is not in {Exchange reactions, Artificial reactions,
  Transport reactions, Isolated, Pool reactions, UNMAPPED}. The first
  three were §1 originally; Isolated and Pool reactions were added by
  §11.9 as structurally identical pseudo-reaction classes.

---

## 7. Pre-registered families

Added as §11.1 and §11.2 before the stats were read:

**Family 1 (between-condition at each hpf, 15 contrasts.)**
For each hpf in {24, 48, 72, 96, 120}: BL vs D, BL vs LD, D vs LD.

**Family 2 (within-condition adjacent hpf, 12 contrasts.)**
For each condition in {BL, D, LD}: 24 vs 48, 48 vs 72, 72 vs 96, 96 vs
120.

**Negative control (15 contrasts.)** For each cell, split-chain A
(chains 1 + 2, 1000 samples) vs B (chains 3 + 4, 1000 samples). Run
through the identical pipeline. Reported separately, never pooled with
Family 1 or Family 2 for multiplicity correction.

Multiplicity correction: Benjamini-Hochberg applied **globally within
each family** (not per contrast). Discovery threshold q_BH ≤ 0.10. The
negative control establishes the empirical FDR floor against which the
family's hits are read.

---

## 8. §11.13 SNR + envelope tracking

Replaces the earlier `|z| > 3 → neg_control_suspect` flag with two
quantitative devices:

**noise_floor_W1 per subsystem** = 95th percentile of W1_obs across the
15 split-chain contrasts for that subsystem.

**snr = W1_obs / noise_floor_W1** per family row. This is a direct
effect-size-above-method-noise ratio, normalised to each subsystem's
own split-chain scale.

**envelope_tracking.** For every Tier A/B row we enumerate the
exchange reactions whose single metabolite appears in any reaction of
the subsystem (either cell's extracted model, union). For each such
exchange we compare fc_A / fc_B. If the direction of the W1 shift
matches the fc envelope direction with a ratio ≥ 2 (or ≤ 0.5), we flag
the row — because the model's "finding" may just be the measured
envelope reasserting itself through the LP.

---

## 9. Tier definitions

Only rows with q_BH ≤ 0.10 enter the tier system.

- **Tier A — headline.** snr ≥ 3, N_tier2 ≥ 15, zero-symmetric
  (|Δ zero_frac| ≤ 0.02), not envelope-tracked, and the same subsystem
  × direction fired in ≥ 2 contrasts of the same family.
- **Tier B — reportable, constrained claim.** snr ≥ 2, not asymmetric,
  not envelope-tracked, fails one Tier A requirement (snr 2–3, or
  N_tier2 < 15, or single-contrast).
- **Tier C — flagged, do not interpret in isolation.** snr < 2 OR
  asymmetric_zeros OR envelope_tracking, including cases where
  neg-control variance in that subsystem is roughly the magnitude of
  the "signal."

Tier C is not evidence against biology; it is evidence that this
pipeline cannot discriminate signal from one of the specified
artefacts for that row.

---

## 10. Floor-suspect flag (§11.10)

Each cell's Frame A biomass floor is 0.9 × its own bmax. For the three
96-hpf cells this gives BL_96 at 7.75 vs D_96 at 11.07 vs LD_96 at
11.63 — a 40 % envelope-size difference by construction. Between-cell
contrasts at 96 hpf involving BL_96 therefore share a polytope-size
confound.

Every q ≤ 0.10 row where the contrast involves BL_96 is flagged
`floor_suspect = True`. Before any such row is reported as a biology
finding we resample BL_96, D_96, LD_96 at a shared floor
(biomass_lb = 7.749 for all three; cells have their native envelopes
otherwise), run those three contrasts, and compare. Only rows that
survive the matched-floor rerun with q ≤ 0.10 are kept.

---

## 11. Sensitivity analyses (full re-runs)

**Median normaliser (§11.12).** Rerun the 42 contrasts with mu_r =
median of pooled non-zero values and samples with `|v| < 1e-10`
dropped (not EPS-shifted). This is logged as a §2 amendment because the
mean-normaliser pre-registration remains the method of record; the
median rerun is a sensitivity report.

**Concordance** of findings_mean and findings_median is written as
`concordance_mean_vs_median.csv`. Rows that are findings under both
normalisers keep their classified tier. Rows that are mean-only or
median-only findings are demoted to Tier C with a `normaliser_disagree`
flag.

**Matched-floor 96 hpf (§11.10).** New sampling of BL_96, D_96, LD_96
at biomass_lb = 7.749 (12 chains). Run through CFF Gurobi. Three
contrasts rerun. Any row previously `floor_suspect` is updated to its
matched-floor q value, with a `matched_floor_concordance` column
indicating whether the row holds up.

---

## 12. Visualisation

Four companion PDFs. Each has the same organisation — rows are
subsystems, columns are the five hpf, red border on any (subsystem,
hpf) that was significant in at least one contrast involving that hpf
— so the three visualisations can be read side-by-side.

### 12.1 Density KDEs (mean + median)

Files: `density_sig_subsystems_v2.pdf`, `density_sig_subsystems_median.pdf`
(and the paginated `..._all_subsystems_*` siblings).

Three overlaid KDEs per panel: BL (blue), D (red), LD (green). The
quantity on the x-axis is the W1 test's input — `log₂(|v| / mu_r)`
over all `(reaction, sample)` pairs in the subsystem — with mu_r =
pooled mean of per-reaction values across all three conditions. Zero
fractions per condition are annotated in each panel.

**Reads clearly for subsystems with coordinated shift** (all reactions
move the same direction). It also *under-reads* subsystems where the
shift lives in the tails of the per-reaction distribution: KDE
smoothing concentrates visible density at the mode, so if two
distributions have identical modes but very different 95th-percentile
behaviour, the KDE curves overlap even though the W1 is large. This
is what we saw with Fatty acid activation (cytosolic) BL vs LD at 48
hpf — see 12.2.

### 12.2 Empirical CDFs

File: `cdf_sig_subsystems_median.pdf`.

Same statistic as the KDE but plotted as empirical cumulative
distribution — sorted log-ratio on x, cumulative fraction on y. **Tail
differences are visible here.** Where the KDEs overlap at the mode
but the test fires, the CDFs show the three stepped curves spreading
apart in the extremes. This is particularly useful for diagnosing why
a q ≤ 0.10 panel looks "the same" under KDE — if its CDFs separate
in the upper 10 % quantile but not around the median, the W1 is being
driven by a right-tail subset of reactions.

### 12.3 Per-reaction trajectory plot — the "which reactions moved" view

File: `trajectory_sig_subsystems.pdf`.

Each panel plots **one point per reaction per condition** on a shared
y-axis of `log₂(median |v|)` (median taken across the 2000 samples
for that reaction in that cell). The three conditions BL, D, LD are
placed on the x-axis in that order, and for every reaction the three
points for that reaction are **connected by a thin line**. The thick
black line is the subsystem-level mean of per-reaction
`log₂(median |v|)`.

How to read it:

- **A reaction whose thin line slopes upward from left to right** has
  higher log₂(median flux) in LD than in BL at that hpf.
- **A reaction whose line slopes downward** has higher flux in BL.
- **A reaction whose line is flat** has similar flux across all three
  conditions at that hpf.
- The thick black overlay is the subsystem's mean trajectory. If it
  is tilted, the subsystem is shifting coordinately in that direction
  at that hpf. If it is flat but the thin lines fan out widely, the
  subsystem has heterogeneous per-reaction responses — some reactions
  go up and some go down, cancelling in the mean — which is exactly
  the pattern the W1 statistic is sensitive to and which KDE and even
  CDF overlays can hide, because W1 does not require the two
  distributions to have different means to differ.
- The "n = …" annotation in each panel is the number of reactions
  surviving the loop filter in that panel (one reaction per enzyme-
  subset representative, loop-flagged ones excluded).
- The red border on a panel means at least one contrast involving
  that hpf was at q ≤ 0.10 under the mean-normaliser pipeline. This
  does **not** say which pairwise contrast fired — see the findings
  CSVs for that.

This figure is the clearest answer to the question "if the test says
this subsystem is different, which reactions is it actually
disagreeing about?" If a panel has a red border and the thick line is
visibly tilted, the subsystem shift is coordinated. If a red-bordered
panel has a flat thick line but visibly fanned-out thin lines,
different reactions in the subsystem are shifting in different
directions and the W1 is picking up the heterogeneity itself as a
distributional difference.

### 12.4 Reading all four plots together

For any Tier A row we recommend reading in this order:

1. **CDF** of that subsystem at that hpf — does the right or left tail
   diverge?
2. **Trajectory** of that subsystem at that hpf — do the thin lines
   slope coherently, or fan out?
3. The row's **`ex_fc_summary`** column in the findings CSV — are
   there measured exchanges with fc ratios ≥ 2 (or ≤ 0.5) in the
   same direction as the shift?

A Tier A row where (1) and (2) show a visible, coordinated shift and
(3) is empty is the strongest claim this pipeline can make.

---

## 13. What the pipeline will and will not do

**It will:** fire on reproducible per-subsystem flux-distribution
differences at defined q level, calibrated against a per-subsystem
split-chain noise floor and a per-contrast envelope sensitivity.

**It will not:** distinguish a biological subsystem shift from the
measured-envelope shift the LP recomputes to be consistent with it.
Hence `envelope_tracking → Tier C`.

**It will not:** reliably distinguish small per-reaction shifts in a
large subsystem (N_tier2 > 50) from the test's statistical floor when
the W1 is < 2× the subsystem's split-chain 95th percentile. Hence
`snr < 2 → Tier C`.

**It will not:** be interpreted as a per-reaction mechanistic finding
from subsystem-level significance alone. Reaction-level checks are
done in the trajectory plots and the per-reaction stats in the
findings CSVs' `ex_fc_summary` column.

---

## 14. Current status (reading date: 2026-10-05)

With the mean-normaliser pre-registered pipeline:

- 976 rows pass §1 filter in the negative control.
- 0 rows at q ≤ 0.10 in the negative control; median |z| = 0.40;
  95th-pct 2.19. Empirical FDR = 0.0 at q = 0.10. ✓
- Family 1: 948 rows pass §1; 18 at q ≤ 0.10 (4 Tier A).
- Family 2: 739 rows pass §1; 15 at q ≤ 0.10 (4 Tier A).

Under the §11.12 median-normaliser sensitivity:

- Negative control empirical FDR climbs from 0.0 to 0.005
  (5 of 976 rows at q ≤ 0.10). Still well below α = 0.10.
- Family 1: 63 at q ≤ 0.10; Family 2: 78 at q ≤ 0.10.
- The median normaliser is substantially more liberal; this is
  reported in the response, not glossed over.

Concordance of 27 rows **in both methods**, 6 mean-only, 114
median-only, zero direction disagreements. Under the §11.12 rule
(mean-only and median-only → Tier C), six subsystem × contrast rows
are both-methods Tier A, with floor_suspect status noted:

| contrast | subsystem | floor_suspect |
|---|---|:---:|
| BL_120 vs D_120 | Phe/Tyr/Trp biosynthesis | — |
| BL_96 vs D_96 | Phe/Tyr/Trp biosynthesis | ⚠ |
| BL_96 vs LD_96 | Arginine / proline | ⚠ |
| D_96 vs LD_96 | Arginine / proline | — |
| LD_24 vs LD_48 | Fatty acid activation (cytosolic) | — |
| D_24 vs D_48 | Fatty acid activation (cytosolic) | — |

Four rows are **both-methods Tier A AND not floor-suspect** — the
truly rock-solid set: Phe/Tyr/Trp biosynthesis at BL_120 vs D_120,
Arginine/proline at D_96 vs LD_96 (matches v5 at v5_z = 2.71), and
Fatty acid activation (cytosolic) at LD_24 → LD_48 and D_24 → D_48.

The matched-floor 96-hpf resample is complete (12/12 samples on disk);
CFF and the three matched-floor contrasts are the final step before the
floor-suspect rows are converted to Tier A/B/C.

---

## 15. Known limitations

- **Retinol metabolism** fires at q ≤ 0.10 at the 48→72 transition in
  all three conditions (Family 2, all Tier C, all neg-suspect, all
  direction = ↓). This is the biomass-coupled retinoid subset reacting
  to Frame A's per-cell biomass floor shifting as bmax rises with hpf;
  not a developmental finding. Reported as a Frame A limitation.
- **ISYNA1 absent** from the Zebrafish-GEM. Inositol is in the medium
  whitelist with NRC 2011 justification.
- **Phosphocholine not producible internally** in any extracted model.
  Documented as a known model limitation (§10.11).
- **Fatty acid oxidation at BL_24 vs LD_24** (Tier C, snr barely above
  neg-control floor). The method's own split-chain noise for this
  subsystem is ≈ 0.48, observed W1 = 0.62 — not a biology claim.

---

## 16. What claims we will make

The response will headline only Tier A rows that:

1. are at q ≤ 0.10 in **both** the mean-normaliser (§2) and the
   median-normaliser (§11.12) pipelines;
2. have snr ≥ 3, N_tier2 ≥ 15, zero-symmetric, not envelope-tracked,
   not floor-suspect (or floor-suspect but validated under the matched-
   floor rerun);
3. replicated in ≥ 2 same-family contrasts under both normalisers.

Tier B rows (strong single-contrast or small-N) will be discussed in
results, not abstract. Tier C rows will appear only in supplementary
with their flag reasons annotated.

Every amendment is dated and logged in `FAMILY_PREREG.md`; every
sensitivity rerun is in a separate output directory alongside the
pre-registered one for direct comparison.

— end of document —
