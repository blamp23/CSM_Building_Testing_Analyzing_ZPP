# 2 Methods

The pipeline is implemented in the public repository
`CSM_Building_Testing_Analyzing_ZPP` (archived at
`github.com/blamp23/CSM_Building_Testing_Analyzing_ZPP`, tag `chapter-v1`).
Six executable stages live under `pipeline/` and shared assets under
`functions/`, `data/` and `rules/`. The sequence is:

| Stage | Role |
|---|---|
| **00_bake** | Build five OCR-anchored base models from Zebrafish1 |
| **01_extract** | rFASTCORMICS extraction, 15 context-specific models |
| **02_sampling** | PolyRound rounding + hopsy billiard-walk sampling |
| **03_loop_removal** | CycleFreeFlux correction of sampled fluxes |
| **04_matched_floor** | Resample at matched biomass floors for the 9 cells that participate in the six developmental contrasts |
| **05_chapter** | ATP and xenobiotic capacity LPs + per-reaction comparison framework |

Transcriptomics preprocessing, metabolomics processing and the boundary
helper (`functions/apply_context_bounds.py`) feed every stage. Each
extracted model is probed with an independent single-substrate ATP
capability sweep (§2.9). All code and curation data are tagged frozen at
`chapter-v1`; large binary artefacts (base and extracted `.mat` files,
hopsy pools) live outside git. A supplementary mapping of each methods
step to its implementation script is in `MANIFEST.md` in the same
directory as this file.

## 2.1 Experimental data

*Danio rerio* embryos were raised under three light regimes: constant
blue light (BL), constant dark (D), and 12:12 light/dark (LD). Embryos
were collected at 24, 48, 72, 96 and 120 hpf, with three biological
replicates per regime × stage for RNA-seq and metabolomics.

Oxygen-consumption rate (OCR) is reported per hpf in Table 1. Measured
values anchor the earliest timepoints; later values are stage-wise
interpolations (`Data/Metabolomics/physiology_by_timepoint.csv` on the
analysis host; the five values become the cap on the oxygen exchange
`MAR09048` in §2.4).

## 2.2 Transcriptomics

Reads were quality-trimmed, aligned to GRCz11 and quantified per gene as
fragments per kilobase of transcript per million mapped reads:

$$\mathrm{FPKM}_{g,s} = \frac{10^{9}\,C_{g,s}}{L_{g}\,N_{s}}$$

for gene $g$ in sample $s$, with read count $C$, gene length $L$ and
mapped-read total $N$. FPKM corrects for gene length and depth but not
for distribution shifts between samples, which matters because extraction
needs a per-gene expressed-or-not call. A fixed cutoff calls different
fractions of the genome expressed depending on library yield, a
percentile cutoff forces every sample to the same expressed count, and
raw FPKM as a continuous score (GIMME) conflates abundance with activity.
We used zFPKM (Hart et al. 2013) instead.

Across a sample's genome, $\log_{2}(\mathrm{FPKM})$ is bimodal, with a
noise mode on the left and an approximately Gaussian expressed mode on
the right. For each sample we estimate the density by kernel density
estimation (excluding $\mathrm{FPKM}=0$), locate the expressed-mode peak
$\mu(s)$ and width $\sigma(s)$, and standardise:

$$\mathrm{zFPKM}_{g,s} = \frac{\log_{2}\mathrm{FPKM}_{g,s}-\mu_{s}}{\sigma_{s}}$$

A zFPKM near 0 is a typical expressed gene; below −3 is noise.
rFASTCORMICS uses these scores for a three-level call. Genes at or above
`exp_threshold` (default 0) are expressed and enter the FASTCORE core,
genes at or below `inex_threshold` (default about −3) are removal
candidates, and genes in between are left for flux consistency to decide.
Thresholds are per sample (`zfpkm_thresholds_by_cell.csv`); in our data
`exp_threshold` sits near 0, `inex_threshold` runs from about −1.2 to
−2.5, and expressed-gene counts vary 10–15 % across cells.

Gene calls reach reactions through gene-protein-reaction (GPR) rules,
evaluated as OR (maximum across isoforms) and AND (minimum across
subunits). The symbol-to-Ensembl map
(`data/gene_map_ensdarg_to_model.csv`) is many-to-many through
paralogues and renames; where several Ensembl IDs match one symbol we
take the maximum zFPKM. Zebrafish1 contains 2,697 unique gene symbols
across its GPR rules.

**Unknown-zone caveat.** Reactions tied to a gene sitting in the
uncertain zone of the three-level call enter each extracted model
through the flux-consistency rule, not through transcription. Two
subsystems in our analysis, retinol metabolism and drug metabolism,
contain reactions pulled into every model this way through a single
gene that sits near the detection floor in all 15 cells. Enrichment
hits in these subsystems describe the extraction behaviour, not
transcriptional regulation by light, and are flagged accordingly in
§3.1 (`figure_S_unknown_zone.pdf` plots the gene location on the zFPKM
density).

## 2.3 Metabolomics envelope

Metabolomics enters the models as bounds on exchange reactions. For
metabolite $m$ in sample $s$, the measured abundance is referenced to
the fifteen-cell median:

$$\mathrm{fc}_{m,s} = \frac{\tilde{I}_{m,s}}{\operatorname{median}_{s}\tilde{I}_{m,s}}$$

The median is robust to single-sample outliers. The result is relative
abundance rather than absolute concentration, which is why the model's
flux units are relative too. The processed table
`data/phase2_exchange_bounds_long.csv` holds 1,575 rows: 105 exchanges
× 15 cells. The 105 measured metabolites include the 20 proteogenic
amino acids plus taurine, citrulline, ornithine, GABA and
trimethyl-lysine; hexoses, glycolytic and TCA-cycle intermediates
(lactate, pyruvate, citrate, isocitrate, 2-oxoglutarate, succinate,
fumarate, malate); nucleosides, nucleotides and bases; free fatty
acids (saturated and unsaturated, 10:0 through 24:1); medium- and
long-chain acylcarnitines; water-soluble cofactors (FAD, GSH, GSSG,
NAD, nicotinamide, riboflavin, pantothenate, UDP-glucose,
UDP-glucuronate); choline and its small derivatives (phosphocholine,
glycerophosphocholine, CDP-choline, acetylcholine); carnitine; a few
glycerolipid precursors (CDP-ethanolamine, LysoPC, LysoPE,
glycerol-3-phosphate); and a handful of small acids (picolinic,
argininosuccinic, N-acetyl-aspartate, N-acetyl-neuraminic, creatine).
Each maps to exactly one Zebrafish1 boundary exchange. The remaining
~1,560 boundary exchanges are unmeasured and are handled by the
curation ledger (§2.6). **Bile acids, ketone bodies, retinoids,
sterols and vitamin-D derivatives are not in the measured envelope**,
so findings on those subsystems (notably §3.1 Findings 3 and 6) are
pure model inferences — no envelope input can drive them. The full
list is `data/phase2_exchange_bounds_long.csv`.

Each measured exchange $j$ in condition $(c,h)$ is bounded by its
fold-change:

$$\mathrm{lb}_{j} = -\mathrm{BOUND\_EX}\cdot\mathrm{fc}_{j,c,h},\quad \mathrm{ub}_{j} = +\mathrm{BOUND\_EX}\cdot\mathrm{fc}_{j,c,h}$$

with `BOUND_EX = 10` fixed across every cell. The bounds are symmetric
because a concentration measurement does not reveal flux direction. The
magnitude follows the fold-change, so a cell with doubled extracellular
lactate gets twice the lactate-exchange headroom. The absolute scale is
a convention, not a measurement — no single bound means anything on its
own. Because every cell uses the same `BOUND_EX`, the ratio between two
cells' bounds on the same exchange is the ratio of what was measured,
and those ratios are the only quantities here that should be read as
physical.

The envelope is applied **after** extraction, not during it. Extraction
decides which reactions a cell keeps; the envelope decides how fast
each one can run. Running extraction with cell-specific bounds would
let metabolomics influence the first decision as well as the second: a
cell with a low measured fold-change on some substrate would carry a
tight bound, and a reaction downstream of it could fail the
flux-consistency check in that cell alone, not because any gene was
off, but because the bound left it too little room. The two models
would then differ in topology, and a comparison between them would
report a network difference whose only cause was a concentration.

## 2.4 Reference GEM and OCR-anchored base (`pipeline/00_bake/`)

Zebrafish1 (Wang et al. 2021) provides $S$, bounds, annotations and GPR
rules — 12,909 reactions and 8,344 metabolites across 9 compartments. Of
those reactions, 7,861 (60.9 %) carry a GPR and 5,048 (39.1 %) carry
none; the GPRs reference 2,697 unique gene symbols.

Three scripts in `pipeline/00_bake/` build the per-stage bases:

```
01_bake_ocr.m           Wang GEM + O2/H2O2 leak fixes + OCR cap  →  baked_ocr_{hpf}.mat
02_build_ledger.py      base parent + s128c v5 flags             →  rules/v7_ledger.csv
03_apply_v7_boundary.m  bake + ledger + ε_pool                   →  baked_ocr_{hpf}_v7.mat
```

Step 01 corrects twenty curated oxygen/peroxide leak reactions
(`o2_leak_fixes.csv`), sets every exchange reaction uniformly to
`lb = −10, ub = +10` as a dimensionless envelope radius, and caps the
oxygen exchange `MAR09048` at the stage-wise OCR:

| hpf | OCR (µmol O₂ gDW⁻¹ h⁻¹) |
|---|---|
| 24 | 96.79 |
| 48 | 173.35 |
| 72 | 226.63 |
| 96 | 273.07 |
| 120 | 313.00 |

The 24 and 48 hpf values are direct Seahorse measurements; the 72, 96
and 120 hpf values are stage-wise interpolations from those
measurements. All five are the best available estimate of oxygen
uptake at the corresponding developmental stage and are used as the
hard cap on `MAR09048` throughout the pipeline.

Step 02 builds the curation ledger (§2.6). Step 03 applies the ledger
to the base with four actions:

1. Zero the four reduced-redox substrates (FADH₂, NADH, NADPH,
   ubiquinol) in the cofactor pool `MAR00022`. The right representation
   of a flavin or NAD(P) demand is through riboflavin or niacin uptake
   and internal pool assembly, not through a reduced redox state as a
   growth substrate.
2. Set the pool pseudo-metabolite `MAM01602c` coefficient in the biomass
   reaction `MAR00021` to $\varepsilon_{\mathrm{pool}} = 10^{-3}$. This
   keeps the cofactor pool a trace demand on growth stoichiometry
   rather than a growth ceiling (an $\varepsilon_{\mathrm{pool}}$
   sensitivity scan over $\{10^{-2}, 10^{-3}, 10^{-4}\}$ varies
   biomass_max by a factor of 1.06–1.23 across cells, placing the pool
   as a modest contributor at $10^{-3}$).
3. Close twenty class-2/class-3 exchange reactions (§2.6).
4. Apply a Part 1 / Part 2 default to the remaining ~1,500 unmeasured
   non-ledger exchanges: `ub = +1000` for secretion (prevents byproduct
   sink artefacts) and `lb = 0` for uptake (prevents LP shortcuts
   through unmeasured substrates).

The five `baked_ocr_{hpf}_v7.mat` files are the input to extraction.

## 2.5 Context extraction with rFASTCORMICS (`pipeline/01_extract/`)

Extraction starts from the base network and a core set $C$ of reactions
that must stay active, and asks for the smallest flux-consistent
network that contains them. FASTCORE (Vlassis, Pacheco & Sauter 2014)
solves this via a sequence of linear programs that iteratively activate
core reactions subject to flux consistency; the specific LP sequence
and convergence proof are in their paper. We use the published
rFASTCORMICS implementation (Pacheco et al. 2019) with the standard
flux-is-active threshold $\varepsilon_{\mathrm{fc}} = 10^{-3}$. Note
that $\varepsilon_{\mathrm{fc}}$ (this flux threshold) and
$\varepsilon_{\mathrm{pool}}$ (the biomass-pool coefficient in §2.4)
are the same numeric value by design but govern different things;
sensitivity scans on either are in
`legacy/docs/FAMILY_PREREG.md`.

FASTCC (Vlassis et al. 2014) runs first to remove reactions that
cannot carry flux under any bounds. rFASTCORMICS extends FASTCORE with
the three-level gene call of §2.2: expressed genes put their reactions
in the core, unexpressed genes mark theirs for removal, uncertain genes
leave the decision to flux consistency. A gene is called by the
majority of its three replicates. If the core comes out disconnected,
the method moves borderline reactions from unexpressed to uncertain
one at a time until it connects, and logs each one. We ran it per
sample, giving one model per cell
(`pipeline/01_extract/04_extract_v7.m` → 15 ×
`trans_rfastcormics_{cell}.mat`).

Some boundary reactions are decided by the ledger rather than by
expression. Protected into the core regardless of transcriptomic
status: class-1 essential-nutrient uptakes (vitamins, essential amino
acids, bulk inorganics), yolk-derived lipoproteins not synthesised de
novo, the biomass reaction (`MAR00021`), and ATP maintenance
(`MAR03964`). Closed into the exclusion set: class-2 rows that are
synthesisable from a class-1 precursor and so do not warrant an uptake
shortcut. The ledger (`rules/v7_ledger.csv`, 75 rows) is the audit
record for every one of these decisions.

**Table 1. Extracted model sizes and stage-wise OCR cap.**

| cell | reactions | metabolites | genes | OCR (µmol O₂ gDW⁻¹ h⁻¹) |
|---|---|---|---|---|
| BL_24 | 5,136 | 3,722 | 1,881 | 96.79 |
| BL_48 | 5,001 | 3,737 | 1,863 | 173.35 |
| BL_72 | 5,245 | 3,903 | 1,897 | 226.63 |
| BL_96 | 5,528 | 4,115 | 1,909 | 273.07 |
| BL_120 | 5,504 | 4,103 | 1,913 | 313.00 |
| D_24 | 5,097 | 3,704 | 1,870 | 96.79 |
| D_48 | 4,841 | 3,690 | 1,878 | 173.35 |
| D_72 | 5,213 | 3,876 | 1,888 | 226.63 |
| D_96 | 5,464 | 4,066 | 1,910 | 273.07 |
| D_120 | 5,559 | 4,116 | 1,917 | 313.00 |
| LD_24 | 4,815 | 3,705 | 1,873 | 96.79 |
| LD_48 | 4,862 | 3,710 | 1,877 | 173.35 |
| LD_72 | 5,327 | 3,956 | 1,886 | 226.63 |
| LD_96 | 5,552 | 4,093 | 1,926 | 273.07 |
| LD_120 | 5,592 | 4,163 | 1,919 | 313.00 |

## 2.6 Boundary application and curation ledger (`functions/apply_context_bounds.py` + `rules/`)

Three constraints set the boundary of each cell's polytope. The 105
measured exchanges take the envelope of §2.3, scaled to that cell's
fold-changes. Oxygen uptake is capped at the stage-wise OCR from §2.4.
Four nutrients, α- and γ-tocotrienol, β-carotene and inositol, are
limited to what the biomass reaction stoichiometrically needs rather
than being supplied freely, because unlimited supply of any of them
gave the model an energy source the embryo does not have (the ATP
yield check below rejected those models). The specific demand bound
for each is in `rules/medium_uptake_demand_bounds.csv`. A consistency
check then confirms that every closed exchange is closed, every medium
component is open, and every measured exchange was scaled.

The curation ledger (`rules/v7_ledger.csv`, 75 rows) records every
boundary decision not derivable from the data, with class, action and
justification. Six entries recur in the results:

- **FAD** is allowed in as a nutrient. The two enzymes that make it
  (RFK, FLAD1) drop out of several extractions, so without the uptake
  those cells could not make flavin at all; FAD is measured, so its
  uptake is bounded by the envelope like any other measured exchange.
- **Lipoprotein deliveries** (HDL, LDL, VLDL, chylomicron, apoA1,
  apoB100) are open, since a yolk-fed embryo receives them and does
  not make them.
- **5-deoxyadenosine**, a by-product of radical-SAM reactions that
  nothing in the network consumes, is given a closed sink so it has a
  mass-balanced place to go.
- The small coefficients in the biomass pool are set to $10^{-3}$ (the
  same numeric value as $\varepsilon_{\mathrm{pool}}$), so biomass is
  solvable without those trace components inflating its yield.
- **β-carotene** uptake, identified as a phantom energy source during
  QC, is bounded to biomass demand (per
  `rules/medium_uptake_demand_bounds.csv`).
- **Carotenoid biosynthesis** gets no boundary action, because
  vitamin-A derivatives are internal to the network and have no
  exchange.

**ATP yield sanity check.** With bounds in place,
`apply_context_bounds` closes all uptake, reopens the medium, the four
limited nutrients, one unit of glucose and oxygen capped at the
stoichiometric requirement for complete glucose oxidation
($\mathrm{lb} = -6$), and maximises ATP maintenance. Every extracted
model returns 31.5 ATP (the complete-oxidation yield), against 351 for
the unextracted base under the same test. Limiting the four nutrients
and extracting on expression together remove the reactions the phantom
energy ran through. The per-substrate yield sweep in §2.9 generalises
this check to 44 candidate substrates.

## 2.7 Flux sampling (`pipeline/02_sampling/` + `pipeline/04_matched_floor/`)

For each cell, biomass_max is found by FBA under the §2.6 bounds, and
a floor of $0.9 \cdot \mathrm{biomass}_{\max}$ is imposed on the
biomass reaction, so every sampled state is close to maximum growth.
The 0.9 is a convention; rather than defend it by precedent, its
consequence is controlled directly. For the twelve developmental (F2)
contrasts between adjacent stages within a regime, both cells are
re-sampled at the lower cell's floor so a difference between them
cannot come from one polytope simply being larger than the other. The
resample covers the nine cells that participate in the six F2
contrasts flagged for biomass_max asymmetry under the per-cell floor
(`pipeline/04_matched_floor/run_pipeline.sh`): matched biomass_lb per
regime 8.398 (BL), 8.480 (D), 8.479 (LD), with three cells per regime
spanning the stages where the asymmetry exceeded 15 %.

The polytope is badly shaped for sampling, with some directions
millions of times longer than others, so it is rounded first.
PolyRound (Theorell et al. 2022) removes redundant constraints,
projects the problem into the null space of $S$ (taking it from
~5,300 to ~1,500 dimensions), and applies an affine map that makes the
shape close to isotropic.

We sample the rounded polytope with the billiard walk in hopsy
(Paul et al. 2024, over HOPS; Jadebeck et al. 2021). Hit-and-run
samplers pick a random direction and stop at the first wall; in that
many dimensions, almost every random direction meets a wall almost at
once, so consecutive samples are nearly identical and mixing time
grows exponentially with dimension. Both CHRR and ACHR failed to mix
on these polytopes, returning effective sample sizes in the single
digits from hundreds of draws. The billiard walk picks a random
direction and a random flight distance, travels in a straight line,
and reflects on each face until the distance is spent. A single step
can bounce many times, so a direction that started out pointing across
a thin axis ends up carrying the walker a long way along a long one.
Mixing time grows in proportion to dimension rather than exponentially.

Each cell is sampled with four chains of 150,000 proposed steps
thinned to 500 retained samples per chain (thinning factor 300), for
2,000 samples per cell. Chains start from the Chebyshev centre of the
rounded polytope. A fixed seed is set per chain. Rank-$\hat{R}$
(Vehtari et al. 2021) has a median of 1.013–1.023 across cells, with
98.3–100 % of reactions below 1.10, and median bulk ESS is 155–249.
That clears the usual $\hat{R} < 1.05$ and ESS > 100 bars but not the
stricter $\hat{R} < 1.01$ and ESS > 400. Every claim in §3 is therefore
made at the subsystem level, where a distribution built from many
reactions averages out per-reaction noise, and no per-reaction
confidence intervals are reported.

## 2.8 Loop removal and pooling (`pipeline/03_loop_removal/`)

A sampled flux vector satisfies $Sv = 0$ and the bounds, but it can
still carry flux around internal loops that exchange nothing with the
outside and would need no thermodynamic driving force to run. Those
loops are not biology; they are a property of $S$ that the sampler
cannot see. CycleFreeFlux (Desouki et al. 2015) removes them from each
sample by finding the flux vector that keeps the same exchange and
biomass fluxes and the same sign on every internal reaction, with the
least total internal flux:

$$\begin{aligned}
\min_{v}\;&\sum_{j\in \mathrm{internal}} |v_{j}|\\
\text{s.t.}\;&Sv=0,\quad \mathrm{lb}\leq v\leq \mathrm{ub},\\
&v_{j}=v_{j}^{(0)}\quad (j\in \mathrm{exchange}\cup\{\mathrm{biomass}\}),\\
&\operatorname{sign}(v_{j})=\operatorname{sign}(v_{j}^{(0)})\quad (j\in \mathrm{internal})
\end{aligned}$$

The absolute values are handled by splitting each internal flux into
positive and negative parts, which makes the problem a linear program.
Across all sixty chains the correction removes 89.4 ± 0.1 % of
internal flux magnitude, which is how much of what the sampler
produced was loop flux. Every analysis uses the corrected samples.

**Pooling at the representative level.** The ~5,300 reactions in a
model are not independent. Many run in fixed proportion to one another
because the stoichiometry allows no other way, and treating each as a
separate observation would count the same flux several times. We group
reactions whose fluxes are forced into perfect proportion by $S$,
columns with cosine similarity above 0.9999 on a basis of the null
space, and report one representative per group, chosen as the reaction
with the highest mean absolute flux (ties broken by reaction ID).
Grouping depends on $S$ alone, not on bounds, so a given reaction
belongs to the same subset in every cell that retains it. For LD_72
this takes 5,300 reactions down to 2,649 representatives. The four
chains are pooled per cell into a representatives × 2,000 matrix
(`pooled_v7/{cell}_pooled.csv`) with the subset membership map
(`subsets_v7/{cell}_subset_membership.csv`). **Every count in §2.10
and in the subsystem tables of §3 is over representatives.** When the
Results refer to "reactions" they mean representatives unless noted.

## 2.9 Capability probes (`pipeline/05_chapter/scripts/01..03, 08`)

Four linear programs probe each extracted model's capacity under its
own bounds.

**ATP yield sweep (`01_substrate_yields.py`).** For each model and
each of 44 candidate substrates (hexoses, three-carbon intermediates,
acetate, ketone bodies, long-chain fatty acids, nucleosides, the twenty
proteogenic amino acids, choline and inositol), we start from the §2.6
bounds, close every organic uptake, open the 17 inorganic medium
components at `lb = −1000`, hold the four limited nutrients at their
demand bounds, force the substrate in at exactly one unit
(`lb = ub = −1`), cap oxygen at the stoichiometric requirement for
complete oxidation of one unit of the test substrate, leave biomass
free, and maximise ATP maintenance (`MAR03964`). A substrate counts
as usable if its aerobic yield exceeds 0.5 ATP.

Capping oxygen at the stoichiometric requirement, rather than at the
physiological OCR, is important for this diagnostic: it isolates the
model's intrinsic ability to oxidise the test substrate from any
coupling to the full measured envelope. All 15 extracted models
return 31.5 ATP per glucose (textbook ceiling); palmitate yields
104.5–108.5 across cells against 110.8 in the unextracted base. Nine
of the 44 substrates (acetate, glycerol, the ketone bodies, DHA, EPA)
are not usable by any extracted model because the extraction pruned
their degradation routes; histidine, tyrosine, methionine, tryptophan
and lysine yield a fraction of their base value in some or all cells
for the same reason. Every central-carbon substrate oxidises at the
full base-model yield in all 15 cells. Yields are in
`chapter_assets/tables/table_S7_substrate_yields.csv`.

**ATP ceiling under physiological OCR (`02_atp_ceiling_physio.py`).**
For each cell, maximum flux through ATP maintenance under the full
§2.6 bounds (measured envelope + OCR cap + ledger), under three
biomass regimes: free biomass (primary measurement), per-cell floor
(biomass pinned at $0.9 \cdot \mathrm{biomass}_{\max}$ for that cell),
and matched floor (both cells of a contrast pinned at the lower of
their two floors). Dual values on the first-order KKT conditions
identify which bound sets the ceiling. In the free-biomass regime the
dual on `MAR09048` (oxygen) is largest in every cell; the next three
duals, all on measured exchanges, are an order of magnitude smaller.

**Glucose-only ATP ceiling (`03_atp_ceiling_glucose.py`).** Variant of
the physio-OCR ceiling that closes the 43 non-glucose carbon exchanges
in the yield-sweep panel (`lb = 0`) and opens glucose (`MAR09034`) to
`lb = −1000`. The full measured envelope and OCR cap are preserved;
secretion bounds are unchanged. Separates the question of whether a
model knows how to make ATP from glucose from the question of what its
throughput is under the measured envelope.

**Xenobiotic clearance (`08_capacity_figure.py`).** A generic
xenobiotic $X$ is added to each extracted model with four reactions:
uptake, a Phase I hydroxylation that consumes NADPH and O₂ (a CYP-like
step), a Phase II conjugation of the hydroxylated product with either
UDP-glucuronate or glutathione, and export of the conjugate. For each
cell and each conjugation route we maximise $X$ uptake with biomass
pinned at $0.9 \cdot \mathrm{biomass}_{\max}$ under the §2.6 bounds
and OCR cap. $X$ has no identity and no toxicity. The LP measures the
cofactor drain (NADPH, O₂, UDP-glucuronate or glutathione) each model
can divert to clearance once 90 % of its growth demand is met.

## 2.10 Comparison framework (`pipeline/05_chapter/scripts/04, 07`)

For each contrast between two cells $A$ and $B$ we ask one question of
every representative both models carry: how much does mean non-zero
sampled flux differ, and does the difference cluster in particular
subsystems?

**Per-reaction $\log_{2}$ fold change.** For representative $r$
carried by both cells, with mean absolute non-zero flux
$\bar{v}_{r,A}$ and $\bar{v}_{r,B}$:

$$\log_{2}\mathrm{FC}^{\mathrm{raw}}_{r} = \log_{2}\bar{v}_{r,A} - \log_{2}\bar{v}_{r,B}$$

Representatives with $\max|v_{r,c}| > 100$ in either cell are flagged
as loop-carrying and **excluded from the tested set** (not merely
annotated). Pairs where $\bar{v}_r < 10^{-10}$ in either cell are
excluded as numerically zero.

**Per-cell median centering.** Each cell's log-flux vector is
centered on its own median across the 3,905 representatives shared by
all 15 cells:

$$m_{c} = \mathrm{median}_{r}\,\log_{2}\bar{v}_{r,c}\quad (r \in \mathrm{shared})$$

$$\log_{2}\mathrm{FC}_{r} = \log_{2}\mathrm{FC}^{\mathrm{raw}}_{r} - (m_{A}-m_{B})$$

The shift $m_{A}-m_{B}$ removes the magnitude baseline that scales
with the OCR cap (which rises with hpf by construction) and with each
cell's $\mathrm{biomass}_{\max}$ (which varies by extracted topology).
What remains is pathway composition: a subsystem's reactions shift
relative to the rest of the model, not relative to a common zero.

**Reaction gates.** A representative is called *notable* if
$|\log_{2}\mathrm{FC}_{r}| \geq 1$, carried by both cells, and not
loop-flagged.

**Subsystem rollup (every subsystem).** Every subsystem with at least
one tested representative in a contrast gets a row in
`subsystem_summary_F{1,2}_centered.csv` with:

- $n_{\mathrm{tested}}$, $n_{\mathrm{up}}$, $n_{\mathrm{down}}$:
  direction counts after centering
- $n_{\mathrm{notable\_up}}$, $n_{\mathrm{notable\_down}}$
- coherence $= |n_{\mathrm{up}} - n_{\mathrm{down}}| / (n_{\mathrm{up}} + n_{\mathrm{down}})$
- mean, median, and max-magnitude $\log_{2}\mathrm{FC}$

**Five-representative interpretation floor.** A row is reported as a
finding in §3 only if it has at least five tested representatives.
This is an interpretation filter, not a rollup filter: smaller rows
remain in the CSV but are not cited as findings. One exception is
flagged in §3.1 (the four-reaction melatonin finding at 72 hpf),
reported with that caveat.

**Hypergeometric enrichment.** For each contrast, Fisher's exact
one-sided test asks whether notable representatives are
over-represented in each subsystem relative to the rest of the tested
set. For a subsystem with $a$ notable and $b$ non-notable
representatives, out of $a+c$ notable and $b+d$ non-notable across the
contrast, the test reports the one-sided $p$-value; Benjamini–Hochberg
correction within contrast gives the reported enrichment $q$-value.
This is the only significance measure reported, because hopsy samples
are correlated draws from one polytope, not biological replicates;
per-reaction $p$-values on 2,000-sample pools are overpowered to the
point of meaninglessness. Several findings that are coherent across
every tested representative do not pass the $q$-value gate; those are
reported as *coherence findings* and labelled as such.

**Two contrast families.** Light-regime (F1) contrasts compare the
three regimes at matched hpf (15 contrasts). Developmental (F2)
contrasts compare adjacent stages within a regime (12 contrasts). **F2
uses the matched-floor pooled samples from §2.7.** Sign convention:
positive $\log_{2}\mathrm{FC}$ means the alphabetically earlier
regime (BL before D before LD) in F1, or the earlier hpf in F2, is
higher.

**Direction.** When reported for a subsystem, direction is the sign
of the mean centered $\log_{2}\mathrm{FC}$, cross-checked against the
majority direction of individual tested representatives. Subsystems
whose per-reaction direction count is split closer than 60/40 are
reported with no direction.

**Reaction inventory (presence/absence).** Finding 4 (the peroxisomal
β-oxidation switch at 72 hpf, where LD retains the subsystem one
stage before BL and D) rests on reaction counts per model, not on
sampled flux. For each F1 contrast, `reaction_inventory_F1.csv`
reports which reactions are retained in each cell; a subsystem's
presence difference is a structural extraction artifact independent
of the sampling. Finding 4 is confirmed by this inventory plus the
single-gene check below.

**Single-gene zFPKM check (peroxisomal switch).** For a reaction set
whose presence/absence differs across regimes at a stage, we ask
whether any single gene's zFPKM crosses the inclusion threshold
between those cells. For the LD_72 peroxisomal set versus BL_72 and
D_72, the candidate drivers (acox1, fasn, acaca) are modestly higher
in LD at 72 hpf (zFPKM differences of 0.19 to 0.34), but none crosses
`exp_threshold` on its own. Inclusion is therefore a joint effect of
several genes plus the connectivity rule FASTCORE applies when it
adds reactions to keep the core flux-consistent. The switch is
reported as a topology finding first and a flux finding second.

**BL_96 isoleucine-topology sanity check.** The BL_96 model has
biomass_max roughly 40 % below D_96 and LD_96. Two checks locate the
cause. (i) *Envelope-sum comparison* across the 105 measured
exchanges at 96 hpf. Measured isoleucine (`MAR09039`) uptake is
similar across regimes at 96 hpf, so the difference is not in what
the embryo takes in. (ii) *Binding-constraint dual analysis* on
biomass_max. In BL_96 the dual on `MAR09039` is the largest active
constraint; in D_96 and LD_96 it is small, and the §11.17
perturbation confirms a 10 % isoleucine relaxation moves biomass_max
by +8.2 % at BL_96 and +0.7 % at BL_120. The BL_96 extraction has
retained a network that lacks the internal isoleucine-supplying
reactions the D_96 and LD_96 extractions kept. This one topology
difference drives the chapter's BL_96 thread across §3.1, §3.2 and
§3.3.

**Clock gene panel (§3.4).** The 29-gene panel in
`chapter_assets/tables/table_S8_clock_panel_full.csv` is the
canonical zebrafish core-clock set plus paralogues: *bmal1a*,
*bmal1b*, *bmal2*, *clocka*, *clockb*, *per1a*, *per1b*, *per2*,
*per3*, *cry1a*, *cry1b*, *cry2*, *cry3a*, *cry3b*, *nr1d1*,
*nr1d2a*, *nr1d2b*, *roraa*, *rorab*, *rorca*, *rorcb*, *nfil3*,
*nfil3-4*, *nfil3-5*, *nfil3-6*, *bhlhe40*, *bhlhe41*, *dbpa*,
*dbpb*. Six genes had no zFPKM value in the 120 hpf samples
(*clocka*, *clockb*, *rorca*, *nfil3-4*, *nfil3-6*, *dbpa*) and were
excluded, leaving 23. The 0.5-zFPKM margin rule, LD exceeds both
constant regimes by at least 0.5 units, is an interpretation
threshold chosen before the panel was tabulated; it keeps the panel
small and the gap visible on single-time-point bulk data, which
cannot separate rhythm amplitude from the LD sampling phase.

## 2.11 Visualization (`pipeline/05_chapter/scripts/05, 06, 09`)

**PCA (`05_pca_raw.py`, `06_pca_centered.py`).** Each cell is reduced
to a vector of log mean non-zero flux across the 3,905 shared
representatives, columns standardised, and projected into the first
three principal components. The centered variant subtracts $m_{c}$
from each cell before standardisation. Loadings and explained
variance per PC are exported.

**Subsystem density ridges (`09_density_plots.py`).** Two sets of
per-subsystem density plots are produced via kernel density
estimation on log flux values.

- `figures/subsystem_density/developmental/`: one PDF per subsystem,
  three panels (BL, D, LD), each panel overlaying five KDE curves
  (one per hpf, viridis gradient). Uses raw log flux to preserve the
  magnitude scaling that captures developmental ramp.
- `figures/subsystem_density/light_regime/`: one PDF per subsystem
  that has a cross-regime finding at any hpf, five panels (24, 48,
  72, 96, 120 hpf), each panel overlaying three KDE curves (BL, D,
  LD coloured). Uses per-cell-centered log flux so the regime
  comparison is on composition, not magnitude. A significance box is
  drawn on panels where the subsystem is a cross-regime finding at
  that hpf.

## 2.12 Software and reproducibility

- Python 3.13 with pandas, numpy, scipy, matplotlib, scikit-learn and
  h5py (version pins in `pipeline/05_chapter/requirements.txt`)
- Gurobi (academic licence) is used for all LPs: PolyRound rounding,
  FASTCORE and rFASTCORMICS extraction, CycleFreeFlux loop removal,
  and all capacity LPs
- hopsy (Paul et al. 2024) over HOPS (Jadebeck et al. 2021);
  PolyRound (Theorell et al. 2022); COBRA Toolbox
  (Heirendt et al. 2019); rFASTCORMICS (Pacheco et al. 2019)
- MATLAB R2025b for stages `00_bake` and `01_extract`
- Random seeds are set per hopsy chain; the four chains of each cell
  use the same seed set across cells

**Code and data availability.** All analysis code, the frozen
curation ledger, the metabolomics envelope CSV, the chapter output
tables, and this methods document are at
`github.com/blamp23/CSM_Building_Testing_Analyzing_ZPP`, tag
`chapter-v1`. The large binary artifacts (base and extracted `.mat`
model files, hopsy sample pools of about 1 GB per cell) stay on the
analysis host; the `CSM_DATA_ROOT` environment variable in
`pipeline/05_chapter/scripts/config.py` points at them.
Preregistration and amendment log: `legacy/docs/FAMILY_PREREG.md`.
