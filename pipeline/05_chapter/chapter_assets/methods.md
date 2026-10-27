# 2 Methods

The pipeline is implemented in `CSM_Building_Testing_Analyzing_ZPP/` under
`pipeline/00_bake` through `pipeline/05_chapter`, with shared assets in
`functions/`, `data/` and `rules/`. The six stages are:

| Stage | Does what |
|---|---|
| **00_bake** | Build five OCR-anchored base models from Zebrafish1 |
| **01_extract** | rFASTCORMICS extraction → 15 context-specific models |
| **02_sampling** | PolyRound rounding + hopsy billiard-walk sampling |
| **03_loop_removal** | CycleFreeFlux correction of sampled fluxes |
| **04_matched_floor** | F2 resample at matched biomass floors for 9 cells |
| **05_chapter** | ATP + xenobiotic capacity LPs + per-reaction comparison framework |

Transcriptomics preprocessing, metabolomics processing and the boundary
helper (`functions/apply_context_bounds.py`) feed every stage. Each extracted
model is probed with an independent single-substrate ATP capability sweep
(§2.10). The sections below specify each stage, with math for the central
modelling choices.

## 2.1 Experimental data

*Danio rerio* embryos were raised under three light regimes (constant blue
light, BL; constant dark, D; 12:12 light/dark, LD) and collected at 24, 48,
72, 96 and 120 hpf, with three biological replicates per regime × stage for
RNA-seq and metabolomics. Oxygen-consumption rate (OCR) was measured by
Seahorse at each stage; its stage-wise mean becomes the oxygen uptake bound
in §2.4.

## 2.2 Transcriptomics

Reads were quality-trimmed, aligned to GRCz11 and quantified per gene as
fragments per kilobase of transcript per million mapped reads:

$$\mathrm{FPKM}_{g,s} = \frac{10^{9}\,C_{g,s}}{L_{g}\,N_{s}}$$

for gene $g$ in sample $s$, with read count $C$, gene length $L$ and
mapped-read total $N$. FPKM corrects for gene length and depth but not for
distribution shifts between samples, which matters because extraction needs
a per-gene expressed-or-not call. A fixed cutoff calls different fractions
of the genome expressed depending on library yield; a percentile cutoff
forces every sample to the same expressed count; raw FPKM as a continuous
score (GIMME) conflates abundance with activity. We used zFPKM
(Hart et al. 2013) instead.

Across a sample's genome, $\log_{2}(\mathrm{FPKM})$ is bimodal, with a
noise mode on the left and an approximately Gaussian expressed mode on the
right. For each sample we estimate the density by kernel density
estimation (excluding $\mathrm{FPKM}=0$), locate the expressed-mode peak
$\mu(s)$ and width $\sigma(s)$, and standardise:

$$\mathrm{zFPKM}_{g,s} = \frac{\log_{2}\mathrm{FPKM}_{g,s}-\mu_{s}}{\sigma_{s}}$$

A zFPKM near 0 is a typical expressed gene; below −3 is noise.
rFASTCORMICS uses these scores for a three-level call. Genes at or above
`exp_threshold` (default 0) are expressed and enter the FASTCORE core,
genes at or below `inex_threshold` (default about −3) are removal
candidates, and genes in between are left for flux consistency to decide.
Thresholds are per sample (`zfpkm_thresholds_by_cell.csv`); in our data
`exp_threshold` sits near 0, `inex_threshold` runs from about −1.2 to −2.5,
and expressed-gene counts vary 10–15 % across cells.

Gene calls reach reactions through gene-protein-reaction (GPR) rules,
evaluated as OR (maximum across isoforms) and AND (minimum across
subunits). The symbol-to-Ensembl map
(`data/gene_map_ensdarg_to_model.csv`) is many-to-many through paralogues
and renames, and where several Ensembl IDs match one symbol we take the
maximum zFPKM.

## 2.3 Metabolomics envelope

Metabolomics enters the models as bounds on exchange reactions. For
metabolite $m$ in sample $s$, the measured abundance is referenced to the
fifteen-cell median:

$$\mathrm{fc}_{m,s} = \frac{\tilde{I}_{m,s}}{\operatorname{median}_{s}\tilde{I}_{m,s}}$$

The median is robust to single-sample outliers. The result is relative
abundance rather than absolute concentration, which is why the model's
flux units are relative too. The processed table
`data/phase2_exchange_bounds_long.csv` holds 1,575 rows: 105 exchanges × 15
cells. Those 105 are the extracellular metabolites that were quantified —
sugars, amino acids, nucleosides, short-chain fatty acids, water-soluble
vitamins, choline, carnitine, bile acids, selected acylcarnitines and
ketone bodies. Each maps to exactly one Zebrafish1 boundary exchange. The
remaining ~1,560 boundary exchanges are unmeasured and are handled by the
curation ledger (§2.6).

Each measured exchange $j$ in condition $(c,h)$ is bounded by its
fold-change:

$$\mathrm{lb}_{j} = -\mathrm{BOUND\_EX}\cdot\mathrm{fc}_{j,c,h},\quad \mathrm{ub}_{j} = +\mathrm{BOUND\_EX}\cdot\mathrm{fc}_{j,c,h}$$

with `BOUND_EX = 10` fixed across every cell. The bounds are symmetric
because a concentration measurement does not reveal flux direction. The
magnitude follows the fold-change, so a cell with doubled extracellular
lactate gets twice the lactate-exchange headroom. The absolute scale is a
convention, not a measurement — no single bound means anything on its own.
But because every cell uses the same `BOUND_EX`, the ratio between two
cells' bounds on the same exchange is the ratio of what was measured, and
those ratios are the only quantities here that should be read as physical.

The envelope is applied *after* extraction, not during it. Extraction
decides which reactions a cell keeps; the envelope decides how fast each
one can run. Running extraction with cell-specific bounds would let
metabolomics influence the first decision as well as the second: a cell
with a low measured fold-change on some substrate would carry a tight
bound, and a reaction downstream of it could fail the flux-consistency
check in that cell alone — not because any gene was off, but because the
bound left it too little room. The two models would then differ in
topology, and a comparison between them would report a network difference
whose only cause was a concentration.

## 2.4 Reference GEM and OCR-anchored base (`pipeline/00_bake/`)

Zebrafish1 (Wang et al. 2021) provides $S$, bounds, annotations and GPR
rules — 12,909 reactions and 8,344 metabolites across 9 compartments. Of
those reactions, 7,861 (60.9 %) carry a GPR and 5,048 (39.1 %) carry none;
the GPRs reference 2,697 unique gene symbols.

Three MATLAB/Python scripts in `pipeline/00_bake/` build the per-stage
bases from Zebrafish1:

```
01_bake_ocr.m           Wang GEM + O2/H2O2 leak fixes + OCR cap →  baked_ocr_{hpf}.mat
02_build_ledger.py      base parent + s128c v5 flags             →  rules/v7_ledger.csv
03_apply_v7_boundary.m  bake + ledger + ε                        →  baked_ocr_{hpf}_v7.mat
```

Step 01 corrects twenty curated oxygen/peroxide leak reactions
(`o2_leak_fixes.csv`), sets every exchange reaction uniformly to
`lb = −10, ub = +10` as a dimensionless envelope radius, and caps the
oxygen exchange `MAR09048` at the measured OCR for the stage:

| hpf | OCR (µmol O₂ gDW⁻¹ h⁻¹) |
|---|---|
| 24 | 96.79 |
| 48 | 173.35 |
| 72 | 226.63 |
| 96 | 273.07 |
| 120 | 313.00 |

Step 02 builds the curation ledger (§2.6). Step 03 applies the ledger to
the base — zeroing four reduced-redox substrates (FADH₂, NADH, NADPH,
ubiquinol) in the cofactor pool `MAR00022`, setting the pool
pseudo-metabolite `MAM01602c` coefficient in the biomass reaction
`MAR00021` to $\varepsilon = 10^{-3}$ (keeps the cofactor pool a trace
demand on growth, not a growth ceiling), closing twenty class-2/class-3
exchange reactions, and applying a Part 1 / Part 2 default to the
remaining ~1,500 unmeasured non-ledger exchanges (`ub = +1000` for
secretion to prevent byproduct-sink artefacts; `lb = 0` for uptake to
prevent LP shortcuts through unmeasured substrates).

The five resulting `baked_ocr_{hpf}_v7.mat` files are the input to
extraction.

## 2.5 Context extraction with rFASTCORMICS (`pipeline/01_extract/`)

Extraction starts from a base network and a core set $C$ of reactions
that must stay active, and asks for the smallest flux-consistent network
containing them. FASTCORE (Vlassis, Pacheco & Sauter 2014) finds that
network through a sequence of LPs:

$$\begin{aligned}
\max_{v,z}\;&\sum_{j\in C'} z_{j}\\
\text{s.t.}\;&Sv=0,\quad \mathrm{lb}\leq v\leq \mathrm{ub},\\
&v_{j}\geq \varepsilon\,z_{j},\quad 0\leq z_{j}\leq 1\quad(j\in C'),
\end{aligned}$$

where $C'$ is the working core subset, $\varepsilon = 10^{-3}$ is the
flux threshold separating active from numerically zero, and $z$ is a
continuous relaxation of the active indicator. Each LP maximises how many
core reactions reach $\varepsilon$; iterating over $C'$ and absorbing
newly activated reactions converges in at most $|C|$ LPs. The output is a
retained set $A$ that is flux-consistent by construction. FASTCC
(Vlassis et al. 2014) runs first and removes reactions that cannot carry
flux under any bounds.

FASTCORE needs every gene called either on or off, and genes near the
threshold change answer between replicates. rFASTCORMICS
(Pacheco et al. 2019) gives them a third option. Expressed genes put
their reactions in the core; unexpressed genes mark theirs for removal;
uncertain genes leave the decision to flux consistency. The reaction
stays if the rest of the network can use it. A gene is called by the
majority of its three replicates. If the core comes out disconnected, the
method moves borderline reactions from unexpressed to uncertain one at a
time until it connects, and logs each one. We ran it per sample with
$\varepsilon = 10^{-3}$, giving one model per cell
(`04_extract_v7.m` → 15 × `trans_rfastcormics_{cell}.mat`).

Some boundary reactions are decided by the ledger rather than by
expression. Protected into the core regardless of transcriptomic status
are class-1 essential-nutrient uptakes (vitamins, essential amino acids,
bulk inorganics), yolk-derived lipoproteins not synthesised de novo, the
biomass reaction (`MAR00021`) and ATP maintenance (`MAR03964`). Closed
into the exclusion set are class-2 rows that are synthesisable from a
class-1 precursor and so do not warrant an uptake shortcut. The ledger
(`rules/v7_ledger.csv`, 75 rows) is the audit record for every one of
these decisions.

**Table 1. Extracted model sizes.**

| cell | reactions | metabolites | genes |
|---|---|---|---|
| BL_24 | 5,136 | 3,722 | 1,881 |
| BL_48 | 5,001 | 3,737 | 1,863 |
| BL_72 | 5,245 | 3,903 | 1,897 |
| BL_96 | 5,528 | 4,115 | 1,909 |
| BL_120 | 5,504 | 4,103 | 1,913 |
| D_24 | 5,097 | 3,704 | 1,870 |
| D_48 | 4,841 | 3,690 | 1,878 |
| D_72 | 5,213 | 3,876 | 1,888 |
| D_96 | 5,464 | 4,066 | 1,910 |
| D_120 | 5,559 | 4,116 | 1,917 |
| LD_24 | 4,815 | 3,705 | 1,873 |
| LD_48 | 4,862 | 3,710 | 1,877 |
| LD_72 | 5,327 | 3,956 | 1,886 |
| LD_96 | 5,552 | 4,093 | 1,926 |
| LD_120 | 5,592 | 4,163 | 1,919 |

Two kinds of difference separate the fifteen networks, and they were made
differently. Between stages a regime's network changes by up to 1,310
reactions, and expression decided a third to a half of them. Between
regimes a stage's three networks differ by only 500 to 600 reactions, and
expression decided as few as 4 % of those at 24 and 48 hpf, with
flux consistency settling the rest. Development rewrites the network;
the transcriptome does the writing.

## 2.6 Boundary application and curation ledger (`functions/apply_context_bounds.py` + `rules/`)

Three constraints set the boundary of each cell's polytope. The 105
measured exchanges take the envelope of §2.3, scaled to that cell's
fold-changes. Oxygen uptake is capped at the measured OCR for the stage.
Four nutrients — α- and γ-tocotrienol, β-carotene and inositol — are
limited to what the biomass reaction needs from them rather than being
supplied freely. Those four get special treatment because unlimited
supply of any of them gave the model an energy source the embryo does not
have (§3.6 of the preregistration). A consistency check then confirms
that every closed exchange is closed, every medium component is open, and
every measured exchange was scaled.

The curation ledger (`rules/v7_ledger.csv`, 75 rows) records every
boundary decision not derivable from the data. Six entries recur in the
results:

- **FAD** is allowed in as a nutrient. The two enzymes that make it (RFK,
  FLAD1) drop out of several extractions, so without the uptake those
  cells could not make flavin at all; FAD is measured, so its uptake is
  bounded by the envelope like any other measured exchange.
- **Lipoprotein deliveries** (HDL, LDL, VLDL, chylomicron, apoA1, apoB100)
  are open, since a yolk-fed embryo receives them and does not make them.
- **5-deoxyadenosine**, a by-product of radical-SAM reactions that
  nothing in the network consumes, is given a closed sink so it has a
  mass-balanced place to go.
- The small coefficients in the biomass pool are set to $10^{-3}$, the
  same $\varepsilon$ used in extraction, so biomass is solvable without
  those trace components inflating its yield.
- **β-carotene** uptake, identified as a phantom energy source, is
  bounded to biomass demand (see §11.19 of the preregistration).
- **Carotenoid biosynthesis** gets no boundary action because vitamin-A
  derivatives are internal to the network and have no exchange.

With bounds in place, `apply_context_bounds` closes all uptake, reopens
only the medium, the four limited nutrients, one unit of glucose and
oxygen, and asks for the most ATP each model can make. Every extracted
model returns 31.5 (the textbook complete-oxidation yield), against 351 for
the unextracted base under the same test. Limiting the four nutrients and
extracting on expression together remove the reactions the phantom energy
ran through.

## 2.7 Flux sampling (`pipeline/02_sampling/` + `pipeline/04_matched_floor/`)

Sampling draws from each cell's polytope near its growth ceiling. For
each cell, $\mathrm{biomass}_{\max}$ is found by FBA under the §2.6
bounds, and a floor of $0.9 \cdot \mathrm{biomass}_{\max}$ is imposed on
the biomass reaction. The 0.9 is a convention; rather than defend it by
precedent, its consequence is controlled directly. For developmental (F2)
contrasts between adjacent stages within a regime, both cells are
re-sampled at the lower cell's floor so a difference between them cannot
come from one polytope being larger than the other
(`pipeline/04_matched_floor/`; nine cells × 4 chains).

The polytope is badly shaped for sampling, with some directions millions
of times longer than others, so it is rounded first. PolyRound
(Theorell et al. 2022) removes redundant constraints, projects into the
null space of $S$ (taking the problem from ~5,300 to ~1,500 dimensions),
and applies an affine map that makes the shape close to isotropic.

We sample the rounded polytope with the billiard walk in hopsy
(Paul et al. 2024; Jadebeck et al. 2021). Hit-and-run samplers pick a
random direction and stop at the first wall; in that many dimensions,
almost every random direction meets a wall almost at once, so consecutive
samples are nearly identical and mixing time grows exponentially with
dimension. Both CHRR and ACHR failed to mix on these polytopes, returning
effective sample sizes in the single digits from hundreds of draws. The
billiard walk picks a random direction *and* a random flight distance,
travels in a straight line, and reflects on each face until the distance
is spent. A single step can bounce many times, so a direction that
started out pointing across a thin axis ends up carrying the walker a
long way along a long one. Mixing time grows in proportion to dimension
rather than exponentially.

Each cell gets four chains of 150,000 proposed steps, thinned to 500
retained samples per chain, for 2,000 samples per cell. Chains start from
the Chebyshev centre of the rounded polytope. Rank-$\hat{R}$
(Vehtari et al. 2021) has a median of 1.013–1.023 across cells, with
98.3–100 % of reactions below 1.10, and median bulk ESS is 155–249. That
clears the usual $\hat{R} < 1.05$ and ESS > 100 bars but not the stricter
$\hat{R} < 1.01$ and ESS > 400. Every claim in §3 is therefore made at
the subsystem level, where a distribution built from many reactions
averages out per-reaction noise, and no per-reaction confidence intervals
are reported.

## 2.8 Loop removal and pooling (`pipeline/03_loop_removal/`)

A sampled flux vector satisfies $Sv = 0$ and the bounds, but it can still
carry flux around internal loops that exchange nothing with the outside
and would need no thermodynamic driving force to run. Those loops are not
biology; they are a property of $S$ that the sampler cannot see.
CycleFreeFlux (Desouki et al. 2015) removes them from each sample by
finding the flux vector that keeps the same exchange and biomass fluxes
and the same sign on every internal reaction, with the least total
internal flux:

$$\begin{aligned}
\min_{v}\;&\sum_{j\in \mathrm{internal}} |v_{j}|\\
\text{s.t.}\;&Sv=0,\quad \mathrm{lb}\leq v\leq \mathrm{ub},\\
&v_{j}=v_{j}^{(0)}\quad (j\in \mathrm{exchange}\cup\{\mathrm{biomass}\}),\\
&\operatorname{sign}(v_{j})=\operatorname{sign}(v_{j}^{(0)})\quad (j\in \mathrm{internal})
\end{aligned}$$

The absolute values are handled by splitting each internal flux into
positive and negative parts, which makes the problem a linear program.
Across all sixty chains the correction removes 89.4 ± 0.1 % of internal
flux magnitude, which is how much of what the sampler produced was loop
flux. Every analysis uses the corrected samples.

**Pooling at the representative level.** The ~5,300 reactions in a model
are not independent. Many run in fixed proportion to one another because
the stoichiometry allows no other way, and treating each as a separate
observation would count the same flux several times. We group reactions
whose fluxes are forced into perfect proportion by $S$ — columns with
cosine similarity above 0.9999 on a basis of the null space — and report
one representative per group, chosen as the reaction with the highest
mean absolute flux (ties broken by reaction ID). Grouping depends on $S$
alone, not on bounds, so a given reaction belongs to the same subset in
every cell that retains it. For LD_72 this takes 5,300 reactions down to
2,649 representatives. The four chains are pooled per cell into a
representatives × 2,000 matrix (`pooled_v7/{cell}_pooled.csv`) with the
subset membership map (`subsets_v7/{cell}_subset_membership.csv`). These
are the input to §2.9–§2.10.

## 2.9 Capability probes (`pipeline/05_chapter/scripts/01..03, 08`)

Four linear programs probe each extracted model's capacity under its own
bounds.

**ATP yield sweep (`01_substrate_yields.py`).** For each model and each
of 44 candidate substrates — hexoses, three-carbon intermediates,
acetate, ketone bodies, long-chain fatty acids, nucleosides, the twenty
proteogenic amino acids, choline and inositol — we start from the §2.6
bounds, close every organic uptake, open the 17 inorganic medium
components at `lb = −1000`, hold the four limited nutrients at their
demand bounds, force the substrate in at exactly one unit
(`lb = ub = −1`), leave biomass free, and maximise ATP maintenance
(`MAR03964`), once with oxygen capped at the measured OCR and once with
oxygen closed. A substrate counts as usable if its aerobic yield exceeds
0.5 ATP. Capping oxygen is important: with unlimited oxygen the LP finds
ways to oxidise the inorganic medium components themselves and reports
over 1,000 ATP per mmol; capped at the physiological OCR, oxygen becomes
the limiting reagent and the probe returns textbook yields (31.5 for
glucose, ~106 for palmitate).

**ATP ceiling under physiological OCR (`02_atp_ceiling_physio.py`).** For
each cell, maximum flux through ATP maintenance under the full §2.6
bounds — measured envelope + OCR cap + ledger — under three biomass
regimes: free biomass (primary measurement), per-cell floor (biomass
pinned at $0.9 \cdot \mathrm{biomass}_{\max}$ for that cell), and
matched floor (both cells of a contrast pinned at the lower of their
two floors). Dual values on the first-order KKT conditions identify
which bound sets the ceiling.

**Glucose-only ATP ceiling (`03_atp_ceiling_glucose.py`).** Variant of
the physio-OCR ceiling that closes every carbon exchange except glucose
(`MAR09034` opened to `lb = −1000`), while keeping the OCR cap and
measured secretion bounds in place. Separates the question of whether a
model knows how to make ATP from glucose from the question of what its
ATP throughput is under the measured envelope.

**Xenobiotic clearance (`08_capacity_figure.py`).** A generic xenobiotic
X is added to each extracted model with four reactions: uptake, a
Phase I hydroxylation that consumes NADPH and O₂ (a CYP-like step), a
Phase II conjugation of the hydroxylated product with either
UDP-glucuronate or glutathione, and export of the conjugate. For each
cell and each conjugation route we maximise X uptake with biomass pinned
at $0.9\cdot\mathrm{biomass}_{\max}$ under the §2.6 bounds and OCR cap.
X has no identity and no toxicity. The LP measures the cofactor drain
(NADPH, O₂, UDP-glucuronate or glutathione) each model can divert to
clearance once 90 % of its growth demand is met.

## 2.10 Comparison framework (`pipeline/05_chapter/scripts/04, 07`)

For each contrast of two cells $A$ and $B$ we ask one question of every
reaction both models carry: how much does mean non-zero sampled flux
differ, and does the difference cluster in particular subsystems?

**Per-reaction $\log_{2}$ fold change.** For reaction $r$ carried by
both cells, with mean absolute non-zero flux $\bar{v}_{r,A}$ and
$\bar{v}_{r,B}$:

$$\log_{2}\mathrm{FC}^{\mathrm{raw}}_{r} = \log_{2}\bar{v}_{r,A} - \log_{2}\bar{v}_{r,B}$$

Reactions with $\max|v_{r,c}| > 100$ in either cell are flagged as
loop-carrying and dropped from the comparison. Pairs where $\bar{v}_r < 10^{-10}$ in
either cell are dropped as numerically zero.

**Per-cell median centering.** Each cell's log-flux vector is centered on
its own median across the 3,905 reactions shared by all 15 cells:

$$m_{c} = \mathrm{median}_{r}\,\log_{2}\bar{v}_{r,c}\quad (r \in \mathrm{shared})$$

$$\log_{2}\mathrm{FC}_{r} = \log_{2}\mathrm{FC}^{\mathrm{raw}}_{r} - (m_{A}-m_{B})$$

The shift $m_{A}-m_{B}$ removes the magnitude baseline that scales with
the OCR cap (which rises with hpf by construction) and with each cell's
$\mathrm{biomass}_{\max}$ (which varies by extraction topology). What
remains is pathway composition: a subsystem's reactions shift relative to
the rest of the model, not relative to a common zero.

**Reaction gates.** A reaction is called *notable* if $|\log_{2}\mathrm{FC}_{r}| \geq 1$,
carried by both cells, and not loop-flagged.

**Subsystem rollup.** Each subsystem that has at least five tested
reactions in a contrast gets a row with:
- $n_{\mathrm{tested}}$, $n_{\mathrm{up}}$, $n_{\mathrm{down}}$ — direction counts after centering
- $n_{\mathrm{notable\_up}}$, $n_{\mathrm{notable\_down}}$
- coherence $= |n_{\mathrm{up}} - n_{\mathrm{down}}| / (n_{\mathrm{up}} + n_{\mathrm{down}})$
- mean, median, and max-magnitude $\log_{2}\mathrm{FC}$

**Hypergeometric enrichment.** For each contrast, we test whether notable
reactions are over-represented in each subsystem relative to the rest of
the tested set. For a subsystem with $a$ notable and $b$ non-notable
reactions out of the contrast's $a+c$ notable and $b+d$ non-notable total
tested reactions, Fisher's exact one-sided $p$-value tests
$a \geq a_{\mathrm{expected}}$; Benjamini–Hochberg correction within
contrast gives the reported enrichment $q$-value. This is the *only*
significance measure reported because hopsy samples are correlated draws
from one polytope, not biological replicates; per-reaction $p$-values on
2,000-sample pools are overpowered to the point of meaninglessness.
Several findings that are coherent across every tested reaction do not
pass the $q$-value gate and are reported as *coherence findings* with
that label.

**Two contrast families.** Light-regime (F1) contrasts compare the three
regimes at matched hpf (15 contrasts). Developmental (F2) contrasts
compare adjacent stages within a regime (12 contrasts). F2 uses the
matched-floor pooled samples (§2.7). Sign convention: positive
$\log_{2}\mathrm{FC}$ means the alphabetically earlier regime (BL before D
before LD) in F1, or the earlier hpf in F2, is higher.

Direction, when reported for a subsystem, is the sign of the mean
centered $\log_{2}\mathrm{FC}$, cross-checked against the majority
direction of individual tested reactions. Subsystems whose per-reaction
direction count is split closer than 60/40 are reported with no
direction.

## 2.11 Visualization (`pipeline/05_chapter/scripts/05, 06, 09`)

**PCA (`05_pca_raw.py`, `06_pca_centered.py`).** Each cell is reduced to
a vector of log mean non-zero flux across the 3,905 shared reactions,
columns standardised, and projected into the first three principal
components. The centered variant subtracts $m_{c}$ from each cell before
standardization. Loadings and explained variance per PC are exported.

**Subsystem density ridges (`09_density_plots.py`).** Two sets of
per-subsystem density plots are produced via kernel density estimation
on log flux values.

- `figures/subsystem_density/developmental/` — one PDF per subsystem,
  three panels (BL, D, LD), each panel overlaying five KDE curves (one
  per hpf, viridis gradient). Uses raw log flux to preserve the
  magnitude scaling that captures developmental ramp.
- `figures/subsystem_density/light_regime/` — one PDF per subsystem that
  has a cross-regime finding at any hpf, five panels (24, 48, 72, 96,
  120 hpf), each panel overlaying three KDE curves (BL, D, LD coloured).
  Uses per-cell-centered log flux so the regime comparison is on
  composition, not magnitude. A significance box is drawn on panels
  where the subsystem is a cross-regime finding at that hpf.
