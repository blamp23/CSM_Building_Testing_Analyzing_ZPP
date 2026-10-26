# Handoff packet — Results + Discussion writing

For: a downstream AI session that will write the Results and Discussion
sections of a thesis chapter on light-regime effects in zebrafish embryo
metabolic modeling.

Hand over this document, the BRIEFING file, the data files listed in
§3, and the density-plot figure folders.

---

## 1. Overall task

Write a complete **Results** section (3 analyses, up to 3 pages each —
~8–10 pages total) and a complete **Discussion** section (1–1.5 pages).

The chapter's thesis: context-specific metabolic models (CSMs) extracted
from transcriptomics + metabolomics predict stage-specific, light-regime-
specific metabolic differences in zebrafish embryos reared under three
light conditions (constant blue, constant dark, normal 14/10 light-dark
cycle) at five developmental stages (24, 48, 72, 96, 120 hpf).

Three analyses contribute to that claim:

1. **Statistical story** — per-reaction and subsystem flux comparisons
   identify five stage-specific light-regime findings and a shared
   developmental program through 72 hpf. This is the chapter's main
   contribution. ~3 pages.
2. **Metabolic capacity across development: ATP** — maximum ATP ceilings
   under physiological OCR scale ~2× across development; three regimes
   track together. ~1.5 pages.
3. **Metabolic capacity across development: xenobiotic clearance** —
   maximum detox clearance rates (glucuronide and glutathione routes)
   are stage-specific and modestly condition-dependent. ~1.5 pages.

Plus a **short section** (~0.5 page) on clock gene expression at 120 hpf
linking transcriptional light output to the §1 metabolic findings.

Plus the **Discussion** (~1.5 pages).

Total target: ~10–13 pages of prose + table/figure references.

---

## 2. Required Results structure

### 3.1  Context-specific flux comparison — five light-regime findings on a shared developmental program

Up to 3 pages. One paragraph per finding below, prefaced by one orienting
paragraph on method.

**Opening paragraph (method orientation).** State what the comparison
layer measures: per-reaction log2 fold change of non-zero flux between
cells, aggregated to subsystem level via counts of notable reactions
(|log2FC| ≥ 1) and coherence (fraction moving in a consistent direction),
with per-cell median centering applied uniformly to separate pathway
composition from magnitude scaling. Significance is reported as
hypergeometric q-values (BH-corrected) on reaction counts, not on
per-reaction tests. Reference the comparison tables.

**Finding 1 — Shared developmental lipid activation program (48 → 72 hpf,
all three regimes).** Leukotriene metabolism, arachidonic acid metabolism,
fatty acid oxidation, and fatty acid activation all rise coherently
during the 48 → 72 hpf transition in BL, D, and LD, with mean centered
log2FC ~1.1–1.8 and coherence ≥ 0.9. Reference developmental density
plots in `figures/subsystem_density/developmental/`. Interpret as yolk
mobilization + onset of β-oxidation.

**Finding 2 — 48 hpf BL retains a lipid-activation signature the other
regimes drop (strongest hit, q ≈ 1e-25).** Fatty acid activation
(cytosolic) is dramatically higher in BL_48 than in D_48 or LD_48;
F2 trajectories show D and LD drop this subsystem 2–3× between 24 and
48 hpf and recover at 48 → 72, while BL holds its 24 hpf level flat.
Reference `Fatty_acid_activation_cytosolic.pdf` in both plot folders.

**Finding 3 — 72 hpf constant-light melatonin depression.** Serotonin
and melatonin biosynthesis flux is lower in BL_72 than in either D_72
(mean log2FC −2.55, q = 0.001) or LD_72 (−1.83, q = 0.002), despite the
pathway being switched on in all three regimes between 24 and 48 hpf.
Note the small subsystem size (4 reps) and the fact that no melatonin-
chain exchange is in the measured metabolomic envelope — this is a model
inference, not a metabolomics-input artefact.

**Finding 4 — 72 hpf LD peroxisomal β-oxidation lead.** Peroxisomal β-
oxidation reactions switch on in LD_72 but not in BL_72 or D_72; by
96 hpf all three regimes have switched them on. LD leads by one stage.
Report as both an inventory finding (reaction presence/absence) and a
flux finding at 96 hpf where LD_96 remains higher than D_96 on 5/5 reps.
Note that no single transcriptomic gene drives the switch; it's
multi-gene plus rFASTCORMICS connectivity.

**Finding 5 — 96 hpf BL stalls and reverses the eicosanoid program.**
After centering, BL_96 is lower than both D_96 and LD_96 on Leukotriene
(24/24), Arachidonic acid (6/6), fatty acid desaturation (11/11), and
omega-3 fatty acid metabolism (13/15), with mean centered log2FC −0.9
to −1.1 and coherence ≥ 0.95. F2 trajectories show BL fluxes drop at
72 → 96 while D and LD rise; BL then catches up at 96 → 120 while LD
comes back down. Interpret as a BL-specific stall-and-reversal that
resolves one stage later. Reference `Arachidonic_acid_metabolism.pdf`
and `Leukotriene_metabolism.pdf` in `light_regime/` folder.

**Finding 6 — 120 hpf LD-only sterol and bile acid step.** Bile acid
biosynthesis is higher in LD_120 than in both D_120 (19/21 reps,
q < 0.001) and BL_120 (20/24, q < 0.001); steroid metabolism is higher
in LD_120 than in both D (17/18, q = 0.002) and BL (17/18). F2 shows
these subsystems rise at 96 → 120 only in LD; BL and D are flat.
Interpret as the metabolic signature of entrained clock output at
120 hpf, when the circadian clock is fully expressed — the two
constant regimes cluster together against the cycling one.

**Short caveat paragraph.** State that two subsystems (Retinol metabolism,
Drug metabolism) show light-regime findings in the comparison tables but
are flagged for unknown-zone amplification — a specific gene sitting at
rFASTCORMICS's expression detection floor in all cells, where inclusion
is driven by connectivity rather than transcription. Those findings are
reportable at the model level but not interpretable as transcriptional
regulation.

### 3.2  Metabolic capacity across development — ATP

Up to 1.5 pages. Three paragraphs.

**Opening.** Metabolic flux magnitude is bounded by the physiological
oxygen uptake rate (OCR), which rises with developmental stage
(97 → 173 → 227 → 273 → 313 µmol O₂/g DW/hr at 24/48/72/96/120 hpf).
State what the ATP capacity test measures: maximum ATP maintenance flux
achievable under the measured metabolomic envelope at the hpf-specific
OCR cap, with no biomass constraint applied. One LP per cell.

**Result.** ATP ceiling scales linearly with developmental OCR, rising
from ~1,000 ATP-equivalents at 24 hpf to ~2,100 at 120 hpf — a 2.1× rise
on a 3.2× OCR rise, with the ratio reflecting a mix of oxidative
phosphorylation (~5 ATP per O₂), substrate-level phosphorylation via
glycolysis, and anaerobic contributions. All 15 cells are OCR-cap-limited
at the optimum (confirmed by LP dual values). Cross-regime spread at
matched hpf is 1–8 %, reflecting that OCR is defined per developmental
stage and is identical across BL, D, and LD at matched hpf by
construction. Any residual spread reflects retained wiring (e.g., how
efficiently each cell's extracted model routes carbon through the
TCA cycle and electron transport chain).

**Interpretation.** The models predict that bulk energetic capacity is
set by developmental stage, not by light regime. The flat cross-regime
baseline is itself a useful result: any susceptibility differences
observed in a future crossed light × exposure experiment would be
attributable to the chemical rather than to light-driven differences
in baseline energetics.

Reference `atp_ceiling_physio_ocr.csv` and the per-substrate yields in
`table_S7_substrate_yields.csv`.

### 3.3  Metabolic capacity across development — xenobiotic clearance

Up to 1.5 pages. Three paragraphs mirroring §3.2 structure.

**Opening.** State what the xenobiotic clearance capacity test measures:
maximum rate at which each cell can take up a generic foreign chemical,
hydroxylate it via a NADPH + O₂ step (CYP-like), conjugate it with
either UDP-glucuronate or glutathione, and export the conjugate, while
maintaining growth at 90 % of maximum. One LP per cell per route.
Note that the chemical has no identity — this is a cofactor-drain
test with real stoichiometry, not a toxicity model.

**Result.** The glucuronide clearance route scales from ~30 at 24 hpf
to ~80 mmol/gDW/hr at 120 hpf, rising smoothly with development and
largely tracking across regimes. One notable exception: BL_96 reaches
~124 mmol/gDW/hr — well above D_96 (62) and LD_96 (58). This is a
consequence of BL_96's lower biomass_max (bound by measured isoleucine
uptake, see Discussion); the lower growth demand frees carbon for
Phase II work. The glutathione clearance route is more variable across
both stage and regime: BL_24 at 38 is nearly 4× D_24 and LD_24 (10 and 9),
and D_120 drops to 7 (the lowest value in the whole dataset).

**Interpretation.** These are baseline capacity measurements on unexposed
embryos; the transcriptome driving the extraction carries no induction
signature. What the models report is the ceiling that measured intake
and retained topology permit, not what the embryo would mobilize under
acute exposure. A crossed light × exposure design would be required to
test the response; the near-flat ATP baseline across regimes means such
a design would attribute susceptibility differences to the chemical,
not to pre-existing baseline differences in detoxification capacity.

Reference `xenobiotic_capacity.csv`.

### 3.4  Clock gene expression at 120 hpf

Short section (0.5 page, 3–4 sentences). Describe the 29-gene panel at
120 hpf and the non-monotonic split (5 genes with BL ≤ D < LD; 3 genes
with D < BL < LD). Connect to the §3.1 Finding 6: the LD-specific
sterol / bile acid step is a candidate metabolic consequence of the
entrained clock output, since the two constant regimes (BL and D)
cluster together against the cycling one for both the transcriptional
and metabolic signatures.

Reference `table_6_clock_gene_zfpkm_120.csv` and `table_S8_clock_panel_full.csv`.

---

## 3. Discussion section — target structure

~1–1.5 pages.

**Paragraph 1 — Summary.** The models reveal that light history produces
stage-specific metabolic signatures while leaving bulk energetic capacity
unchanged. Name the five findings in biological terms.

**Paragraph 2 — Biological context.** For each finding, briefly place in
context of what's known about zebrafish embryonic metabolism and
circadian biology. **Do NOT invent citations.** Mark places where a
reference will be added with `[CITE: topic]`. Specifically consider:
yolk mobilization at 48–72 hpf, peroxisomal β-oxidation ontogeny,
melatonin biosynthesis regulation by light, bile acid / steroid
developmental onset, clock gene entrainment at 120 hpf.

**Paragraph 3 — Methodological strength and limits.** The chapter's main
methodological contribution is a per-reaction + subsystem comparison
framework with per-cell median centering that separates pathway
composition from biomass-driven magnitude scaling. Specifically mention:
(a) that classical per-sample p-values are not meaningful for hopsy-
sampled data because samples are not biological replicates;
(b) that enrichment q-values on reaction counts are the one defensible
significance measure; (c) that per-cell biomass_max asymmetry at 96 hpf
is a real feature of the extraction (traced to BL_96's isoleucine
uptake constraint and the absence of internal isoleucine-supplying
reactions in BL_96's retained topology) rather than a sampling artefact;
(d) that unknown-zone amplification affects two specific subsystems
(Retinol metabolism and Drug metabolism) where findings are not
interpretable as transcriptional regulation.

**Paragraph 4 — Baseline vs response, and the crossed design.** The
chapter's measurements are capacity, not response — the embryos were
never exposed. The flat ATP baseline across regimes is useful precisely
because it means a future crossed light × exposure experiment would
attribute any susceptibility differences to the chemical rather than
to light-driven baseline differences. The natural next experiment is
that crossed design: expose embryos at the five hpf under each of the
three light conditions, extract models on the exposed transcriptomes,
and run the comparison framework against a response rather than a
baseline.

---

## 4. Files to hand over

```
grace_v4/results/chapter_assets/
    BRIEFING_FOR_ANALYSIS_v2.md         ← context, pipeline, conventions, sanity checks
    HANDOFF_PACKET.md                    ← this file (writing instructions)

    comparison_tables/
        subsystem_summary_F1_centered.csv       ← primary F1 subsystem
        subsystem_summary_F2_centered.csv       ← primary F2 subsystem
        per_reaction_flux_F1_centered.csv       ← atomic F1 data
        per_reaction_flux_F2_centered.csv       ← atomic F2 data
        reaction_inventory_F1.csv               ← presence/absence
        reaction_inventory_F2.csv               ← presence/absence
        pca_cells_centered.pdf                  ← PCA figure (centered)
        pca_cells.pdf                           ← PCA figure (raw)

    tables/
        xenobiotic_capacity.csv                 ← §3.3 data
        table_6_clock_gene_zfpkm_120.csv        ← §3.4 data
        table_S7_substrate_yields.csv           ← §3.2 ATP sweep context
        table_S8_clock_panel_full.csv           ← §3.4 supplementary

    figures/
        subsystem_density/
            developmental/                      ← per-regime hpf trajectory plots
                (21 PDFs, one per subsystem)
            light_regime/                       ← cross-regime density plots
                (13 PDFs, one per subsystem with F1 finding)

    (also, elsewhere in the repo:)
grace_v4/results/v7_susceptibility/
    atp_ceiling_physio_ocr.csv              ← §3.2 data
```

---

## 5. Instructions for the writing pass

**Voice:** scientific thesis chapter, confident where data supports it,
honest about limits. Flowing prose paragraphs — no bulleted lists in
the body. First-person plural ("we") is acceptable; passive constructions
are acceptable; mix as feels natural.

**Length:** target ~10–13 pages total for Results + Discussion. Do not
pad. If a finding needs fewer than the allotted paragraphs, use fewer.
If a finding needs more and you can justify the biology, use more
(within the 3-page per analysis limit).

**Numbers:** report specific numbers for every finding (log2FC, q,
coherence, n_reps). Pull them from the comparison tables, not from
memory. Do not round aggressively — one or two significant figures is
fine for log2FC; use scientific notation for q-values.

**Figure references:** cite the density plot folder files by name
(e.g., "see `figures/subsystem_density/light_regime/
Fatty_acid_activation_cytosolic.pdf`"). Cite tables by their CSV name.
Assume the reader will have the figures and tables in front of them.

**Direction conventions (critical — do NOT mix up):**
- F1: A = alphabetically earlier condition (BL < D < LD).
  Positive log2FC = A higher than B.
- F2: A = earlier hpf. Positive log2FC = A higher = flux falls across
  development. Negative = flux rises.

**What NOT to invent:**
- Citations. Mark with `[CITE: topic]` where a reference should go.
- Numbers not in the tables. If you're unsure, say so.
- Causal mechanistic claims that go beyond "the model predicts X"
  unless there is corroborating biological evidence that supports it
  from the briefing or tables.

**What NOT to re-report in the Discussion:** numbers already in the
Results section. The Discussion interprets, it does not re-tabulate.

**Scope statements to include:**
- Baseline capacity, not response (recur in §3.2, §3.3, Discussion).
- Hopsy samples are not biological replicates (methods summary, methods
  limits in Discussion).
- The lipid-relay decomposition (Table S5) is incomplete; do not reference
  it as a chapter finding.

---

## 6. Suggested opening prompt for the writing session

> "You previously analyzed per-reaction and subsystem flux comparison
> tables for a zebrafish embryo metabolic modeling study. You have the
> context from BRIEFING_FOR_ANALYSIS_v2.md and the writing instructions
> from HANDOFF_PACKET.md.
>
> Write the complete Results section (sections 3.1, 3.2, 3.3, 3.4) and
> the complete Discussion section for this thesis chapter. Follow the
> structure specified in HANDOFF_PACKET.md §2 and §3. Report specific
> numbers from the comparison tables. Reference density plots and tables
> by filename. Do not invent citations — use `[CITE: topic]` markers.
> Do not re-report numbers in the Discussion. Target ~10–13 pages total.
>
> Begin with Section 3.1. Pause between major sections so I can comment
> if needed."
