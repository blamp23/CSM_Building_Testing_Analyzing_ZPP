# base_model/

This folder turns the published Wang zebrafish GEM into the five per-hpf parent models that extraction consumes. The output is `baked_ocr_{24,48,72,96,120}_v7.mat`, each one a self-contained COBRA model with a thermodynamically consistent topology, oxygen bounded by measured respiration, every exchange sitting at a boundary rule written down in `rules/v7_ledger.csv`, and the pool coefficient in biomass set so that cofactors stop acting as a growth ceiling.

Three scripts, run in order.

```
01_bake_ocr.m           Wang GEM + O2/H2O2 leak fixes + OCR cap  ->  baked_ocr_{hpf}.mat
02_build_ledger.py      base parent + s128c v5 flags             ->  rules/v7_ledger.csv
03_apply_v7_boundary.m  bake + ledger + epsilon                  ->  baked_ocr_{hpf}_v7.mat
```

The rest of this document explains what each step is doing and why.

## The starting model

We start from `Benji_Wang_Zebrafish-GEM_v2.mat`, Wang's zebrafish adaptation of Human-GEM / Human1. 12,909 reactions, 8,344 metabolites, 2,055 genes. It's a reconstruction, not a context-specific model: every reaction with evidence in any human tissue (and a zebrafish ortholog) is present, every exchange is open at ±1000, there is no transcriptomics layer, no measured growth rate, and the biomass reaction has standard eukaryotic composition at coefficient −1 per substrate.

Published-state Wang won't give you anything useful under flux balance analysis. Four things interact to make it unusable before any of the biology is addressed:

Oxygen and peroxide can be produced from water through thermodynamically-infeasible but mass-balance-feasible loops Human-GEM inherits from its own reconstruction history. The LP will exploit these as free energy to maximise biomass. Specific reactions (around twenty, mostly in the ROS defense subsystems) need to have their directions constrained before any FBA is meaningful.

The cofactor pool reaction `MAR00022` treats sixteen cofactors as consumed one-to-one with biomass. Several are catalytically cycled in vivo and turn over thousands of times per cell cycle. Treating them as 1:1 consumables makes the pool a hard growth ceiling: in a wide-open run, biomass_max is limited not by what the embryo actually needs for growth but by how much cytochrome-C, NADH, or FADH₂ the LP has to "produce" to satisfy a fictional stoichiometric demand. The pool coefficients are a modelling convention, not a measurement.

The exchange envelope at ±1000 means the LP can take up any metabolite at an unphysical rate. An embryo in a yolk-sac environment is a closed system that does not have access to arbitrary carbon sources, and bound tightening against measured metabolomics is a core part of any serious context-specific extraction. We handle this in two stages: a uniform ±10 convention baked at step 01, and a per-cell fold-change multiplier applied downstream in `apply_context_bounds`.

Finally, oxygen uptake has been directly measured for zebrafish embryos at each developmental stage, and respiration is one of the most informative single constraints we can apply. There is no good reason to leave `MAR09048` wide open.

Steps 01–03 address these four problems in order.

## 01_bake_ocr.m — the Wang-to-OCR bake

This step produces one baseline per developmental stage. Each file is byte-identical to the others except for the oxygen bound and the metadata.

The input is Wang zebrafish GEM v2 with its original cofactor pool stoichiometry intact. We do not pre-modify any pool coefficients at this layer — the pool question is answered more cleanly in step 03 by a single parameter at the biomass level (epsilon, below) rather than by hand-picking which cofactor coefficients to deflate.

The script loads Wang v2, applies the O₂ / H₂O₂ leak fixes from `o2_leak_fixes.csv` (around twenty curated bound changes, each a `block_fwd`, `block_rev`, or `block_both` action against a specific reaction ID), sets every exchange reaction uniformly to lb = −10, ub = +10, then caps oxygen at the stage-specific OCR:

```
24 hpf   OCR = 96.79  umol O2 / g DW / hr
48 hpf       = 173.35
72 hpf       = 226.63
96 hpf       = 273.07
120 hpf      = 313.00
```

OCR values are from `Data/Metabolomics/physiology_by_timepoint.csv`, averaged across replicates per hpf. The sign convention follows Wang: `MAR09048` writes O₂ as `MAM02630e → ∅`, so uptake is a negative flux. Capping uptake at the measured value means `m.lb = -OCR, m.ub = 0`. Secretion is blocked; oxygen can only be consumed.

The ±10 convention on the other exchanges is worth a short remark. It is not a measured flux value and should not be read as one. It is a dimensionless envelope radius that the metabolomics layer multiplies, per cell, by the measured fold-change for the 105 exchanges we have composition data for. An exchange at fc = 1.0 sits at ±10. An exchange at fc = 0.5 sits at ±5. Using ±1000 (Wang default) lets the LP find mass-balance-feasible but physiologically nonsense solutions; using the measured absolute flux at this step would conflate the bake with the metabolomics envelope, which belongs post-extraction. ±10 is a working midpoint.

Each output file carries a `baked_meta` struct recording the source file, the number of leak fixes applied, the exchange-bound convention, the per-hpf OCR, and the baseline biomass_max from a sanity FBA. The baseline biomass_max comes out around 9.5 at every hpf, limited by amino acid uptakes at the ±10 ceiling. This value does not matter for extraction; it just confirms the model is feasible.

What `01_bake_ocr.m` does not do: no transcriptomics, no metabolomics fold-changes, no ledger-driven boundary rule, no epsilon. The output is the cleanest version of Wang's zebrafish GEM for this project, with the two pre-conditions (thermodynamic sanity, physiological oxygen) satisfied and nothing else.

## 02_build_ledger.py — the curation ledger

The published Wang GEM has about 1,660 exchange reactions. Some of them are measurable in metabolomics and get a fold-change bound. Many are not. The question is: what should the model do with the ones that aren't measured?

Before v7 we never answered this cleanly. We applied the fc bounds to the 105 and left everything else at ±10 default. That default was the single largest source of modelling artefacts in v6: the LP found growth-sustaining routes through unmeasured exchanges that an embryo cannot physiologically access. A dropped riboflavin kinase in some cells meant the LP routed around it by uptaking FAD directly, which looked like real biology until we audited it. The reviewer's correct response was that FAD was never the issue; the problem was that we had no principle for the boundary at all. The ledger is the principle.

The ledger is a table in `rules/v7_ledger.csv`, one row per (metabolite, exchange MAR id) pair across the union of three sets: substrates of the biomass reaction `MAR00021`, substrates of the cofactor pool reaction `MAR00022`, and the uptake-essential exchanges flagged by `s128c` on the v5 models. Every row carries a class, an action, and a one-line justification. There are only three classes.

Class 1 is a true dietary essential. The embryo cannot synthesise it (either the animal lineage never had the pathway, or Zebrafish-GEM lacks the enzyme), so it must cross the boundary. All measured metabolites in the 105 default to class 1. Essential amino acids, essential fatty acids, true vitamins, bulk inorganic ions, O₂, water, and protons all live here. The action is `keep_exchange`: the exchange stays open at its default bound (either ±10 or the fc-scaled value if measured). This class also contains two subcategories that are logically class 1 but tagged with a different action for provenance: `medium_uptake` for unmeasured vitamin / carotenoid uptakes that get set to lb = −1000 so they do not become default-±10 limiters, and `ysl_lipid_delivery` for the lipoprotein-mediated yolk lipid delivery (apoA1, apoB100, chylomicron, VLDL, LDL, HDL) which stays at the Wang default ±10.

Class 2 is something the embryo can make from a class-1 precursor, where the pathway exists in Zebrafish-GEM with a gene rule. FAD comes from riboflavin via RFK and FLAD1. Calcidiol comes from vitamin D3 via CYP2R1 or CYP27A1. Lipoic acid comes from octanoyl-ACP via LIAS. Cholesterol comes from acetyl-CoA via the mevalonate pathway. The action is `close_exchange`: the model is told to make it internally, not take it up. If the gate check later reveals that the internal route is topologically absent (as happened with inositol — Zebrafish-GEM lacks ISYNA1), the metabolite gets reclassified to class 1 with a documented reason.

Class 3 is a boundary artefact. A pool pseudo-metabolite like `PI pool` or `PG-CL pool` that is produced by an internal assembly reaction and whose exchange exists only as a Wang modelling convenience. A complex assembly like LDL remnant or chylomicron remnant that the LP can treat as a shortcut. A reduced redox state like NADH where the exchange is non-physical. A xenobiotic conjugate like bilirubin-bisglucuronoside. The action is also `close_exchange`.

Metabolites that are not in any of the three class dictionaries get labelled `?` and have to be reviewed by hand before the ledger is used. In v7 as frozen on 2026-10-01 there are zero unclassified rows.

One rule overrides everything else: if an exchange is in the 105 measured set, it is class 1 regardless of what its chemistry would otherwise say. The reasoning is simple. If we measured the metabolite, we have evidence about its abundance across cells, and the right thing to do with that evidence is to let the measurement set the bound. Class-2-by-chemistry becomes class-1-by-measurement, and the fc multiplier handles the per-cell variation. This is why FAD, which looks like a class-2 metabolite on paper, is actually class 1 in the frozen ledger: `MAR01939` is in the 105. The reviewer's "FAD was never special" observation drops out of this rule as a corollary.

The ledger is built by `02_build_ledger.py`. The script reads the base parent, enumerates the biomass and pool substrates and the s128c flags, classifies each species via three dictionaries encoded in the script, applies the measurement override, flags envelope leaks (species where more than one exchange exists and some are measured while others are not), and writes `rules/v7_ledger.csv` plus `rules/v7_pool_table.csv` (the 16-row decision table for MAR00022 pool substrates — see below). The output is deterministic. Rerunning the script against the same base parent and the same s128c v5 CSV produces the same ledger.

The ledger is immutable in v7. If v8 ever happens, the ledger changes and the version tag rolls. This is enforced by the amendment log in `FAMILY_PREREG.md` §10, which lists the ledger among the frozen configuration items.

## 03_apply_v7_boundary.m — the overlay, and what epsilon is

With a baked_ocr file and a ledger in hand, this is the final step: produce `baked_ocr_{hpf}_v7.mat` by applying the ledger's rules to the base parent at every hpf. There are five things it does, in order.

First, it edits the biomass pool reaction `MAR00022`. The pool table in `rules/v7_pool_table.csv` has one row per pool substrate with a `net_consumed` boolean and a `new_coef`. The rule is: reduced redox states are not accumulating demands and get zeroed. In practice this zeroes four coefficients — FADH₂, NADH, NADPH, and ubiquinol — because the right representation of a biomass demand for flavin or NAD(P) is through riboflavin or niacin uptake and internal pool assembly, not through the reduced redox state. The other twelve substrates (lipoyl-lysine, biotin, CoA, cobamide-coenzyme, cytochrome-C, L-carnitine, riboflavin, THF, tetrahydrobiopterin, and vitamins A/D/E) keep their Wang coefficients.

Second, it changes the pool's coefficient in the biomass reaction `MAR00021` itself, and this is where epsilon comes in. `MAM01602c`, the "cofactors and vitamins" pseudo-metabolite that `MAR00022` produces, has a coefficient in `MAR00021` of −1 in the Wang reconstruction. That is a convention, not a measurement: Wang chose to make the pool contribute at the same stoichiometric weight as a single amino acid, because the alternative would have been to invent a measured coefficient they did not have. We follow the same logic in the opposite direction. Since we have no measurement for how much cofactor pool a zebrafish embryo demands per unit biomass, we set the coefficient to something small enough that the pool stops acting as a growth ceiling. We write it as ε (epsilon) in the script and in the amendment, with a default value of ε = 10⁻³.

This single edit is the only parameter that controls how much the cofactor pool contributes to growth stoichiometry. The design intent is: whatever the twelve-substrate pool reaction produces, multiply its contribution to biomass by ε so that the pool is a trace demand on growth rather than an inflated one. An ε sensitivity scan over {10⁻², 10⁻³, 10⁻⁴} on all fifteen extracted models shows biomass_max varying by a factor of 1.06 to 1.23 between the ends of the scan, depending on cell. That is not perfectly flat (meaning: the pool is not entirely irrelevant to biomass_max), and it is not a growth ceiling either (meaning: the pool is not setting biomass_max). It sits in the right physiological range — a trace demand that affects the number by a few percent, which is where a cofactor pool should sit if the model is representing it honestly. If ε were set to 10⁻², the pool would start to limit growth; at 10⁻⁴, the pool would be numerically invisible. 10⁻³ is a defensible middle.

Third, it closes the ledger's twenty `close_exchange` rows by setting lb = ub = 0 on those exchange reactions. These are the class 2 and class 3 rows.

Fourth, it applies the Part 1 / Part 2 boundary rule to every remaining exchange that is not O₂, not in the 105, not in the ledger, and not in a hardcoded medium list of 17 bulk inorganics (water, CO₂, HCO₃⁻, H⁺, Pi, NH₃, sulfate, Zn, K⁺, Na⁺, Ca²⁺, Mg²⁺, Cu²⁺, Fe²⁺, chloride, iodide, carbonate). The medium gets lb = −1000, ub = +1000. The remainder gets lb = 0, ub = +1000. The semantic rule is: an exchange reaction we did not measure, that is not a bulk ion, and that is not on the ledger's keep list, should default to being unable to take up anything but able to secrete anything. Part 1 (ub = 1000) says: if the model generates a byproduct that has nowhere else to go, let it leave. This is what fixed the 5-deoxyadenosine problem the Fritzemeier-style ATP tests surfaced. Part 2 (lb = 0) says: do not let the LP find a growth-sustaining uptake we never measured.

Fifth, it saves the result as `baked_ocr_{hpf}_v7.mat` and records the sanity biomass_max at the extraction-bound convention (around 10 for every cell, limited by amino acid uptakes).

## Where this leaves the model

After step 03, every file in `models/baked_ocr_{hpf}_v7.mat` is a Wang zebrafish GEM at a specific developmental stage. The thermodynamic artefacts are blocked, oxygen is physiologically bounded, the cofactor pool reduces to a trace biomass demand through a single parameter at the biomass level (ε = 10⁻³), four reduced-redox substrates are zeroed inside `MAR00022` on principle (not because of their numerical contribution), the boundary rule for every single exchange reaction is either measured, kept, closed, or set to the Part 1 / Part 2 default, and a complete audit trail from Wang v2 to the baked file lives in the `baked_meta` struct. This is the input to extraction.

Nothing in the base model file reflects the per-cell transcriptomics or the per-cell metabolomics fold-change. Those are applied downstream — at extraction time (transcriptomics, via rFASTCORMICS) and at analysis time (metabolomics envelope, via `functions/apply_context_bounds`). The base model and the ledger are the invariant parts that every one of the 15 extracted models is built from.

## For a methods section

Context-specific extraction was performed on five hpf-specific parent models derived from the Wang zebrafish GEM v2, a zebrafish adaptation of Human1. The published Wang reconstruction was first corrected for thermodynamically infeasible oxygen and hydrogen-peroxide pathways by applying twenty curated bound restrictions (`o2_leak_fixes.csv`), after which every exchange reaction was bounded uniformly at ±10 as a dimensionless envelope radius, with the metabolomics fold-change scaling applied post-extraction. Five baselines were saved, one per developmental timepoint, with the oxygen exchange `MAR09048` capped at the measured oxygen consumption rate at that stage (`physiology_by_timepoint.csv`). A curation ledger (`rules/v7_ledger.csv`) classifies every exchange reaction into one of three classes: dietary essentials and all 105 metabolomics-measured exchanges that must remain open (class 1), metabolites with a known internal synthesis route that are closed to force the model to make them (class 2), and boundary artefacts such as pool pseudo-metabolites and complex lipoprotein assemblies that are closed (class 3). The ledger is applied along with three additional edits to produce the v7 baseline per hpf: zeroing four reduced-redox substrates (FADH₂, NADH, NADPH, ubiquinol) in the cofactor pool `MAR00022`; setting the pool pseudo-metabolite `MAM01602c` coefficient in the biomass reaction `MAR00021` to ε = 10⁻³, which keeps the pool a trace demand on growth stoichiometry rather than a growth ceiling (an ε sensitivity scan over {10⁻², 10⁻³, 10⁻⁴} shows biomass_max varying by a factor of 1.06–1.23 across the scan, placing the pool as a physiologically plausible modest contributor); closing twenty class-2 and class-3 exchange reactions; and applying a Part 1 / Part 2 default to the remaining 1,500 unmeasured non-ledger exchanges (ub = +1000 for secretion to prevent byproduct-sink artefacts, lb = 0 for uptake to prevent LP shortcuts through unmeasured substrates). The resulting five `baked_ocr_{hpf}_v7.mat` files are the input to rFASTCORMICS extraction.
