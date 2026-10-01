# base_model/ — baking the Wang zebrafish GEM into 5 per-hpf parents

This folder documents **everything that happened to the published Wang zebrafish GEM before any v7-specific curation**. The output of this stage is 5 per-hpf baseline `.mat` files (`baked_ocr_{24,48,72,96,120}.mat`) that are the input to `03_apply_v7_boundary.m`.

Scripts in this folder:

| step | script | input | output |
|---|---|---|---|
| 1 | `01_bake_ocr.m` | `models/zebrafishGEM_v2_modcofpool.mat` | `models/baked_ocr_{hpf}.mat` × 5 |
| 2 | `02_build_ledger.py` | base parent + s128c v5 output | `rules/v7_ledger.csv`, `rules/v7_pool_table.csv` |
| 3 | `03_apply_v7_boundary.m` | `baked_ocr_{hpf}.mat` + ledger | `baked_ocr_{hpf}_v7.mat` × 5 |

Steps 2 and 3 are the v7 overlay. **This README is about step 1 and everything upstream of it** — i.e. how the Wang model becomes `zebrafishGEM_v2_modcofpool.mat` and then `baked_ocr_{hpf}.mat`.

---

## 1. Starting point — Wang zebrafish GEM v2

- **File**: `Benji_Wang_Zebrafish-GEM_v2.mat` (upstream of this repo; the model is a zebrafish adaptation of Human1 / Human-GEM by Wang et al.)
- **Dimensions**: 12,909 reactions × 8,344 metabolites × 2,055 genes (Zebrafish-GEM 1.0.0 lineage, Wang adaptation).
- **Compartments**: cytosol (c), mitochondrion (m), nucleus (n), ER (r), peroxisome (x), lysosome (l), Golgi (g), extracellular (e), ∅ (s).
- **Published state**: Wang's model is a *reconstruction*, not a context-specific model. It has every reaction with any gene evidence across human tissues (and zebrafish orthologs), every exchange reaction open at ±1000, no transcriptomics, no proteomics, no measured growth rate. Biomass reaction (`MAR00021`) is a standard eukaryotic composition at coefficient `−1` per substrate.

The published state is **not usable for context-specific modeling as-is**. Four things must be fixed before we can even start thinking about which reactions the zebrafish embryo uses:

1. **Thermodynamically-infeasible shortcuts** that Wang's reconstruction inherits from Human-GEM — specifically O₂ and H₂O₂ cycles that let the LP make ATP from nothing.
2. **The biomass cofactor pool stoichiometry** — Wang treats several catalytically-cycled cofactors (NADH, cytochrome-C, etc.) as if they are consumed 1 : 1 with biomass, which is biologically wrong and numerically pathological.
3. **The exchange envelope** — ±1000 lets the LP use any metabolite as food, which is unphysical for an embryo in a closed yolk-sac environment.
4. **Oxygen uptake** — no OCR is applied, so the LP can respire at any rate; the actual embryo's OCR has been measured and is one of the single most informative constraints we can impose.

Everything in this folder exists to turn Wang's published-state model into something that reflects a developmentally-timed, metabolomics-aware, OCR-anchored embryo — but **before** any of the per-cell transcriptomics or ledger curation applies.

---

## 2. Upstream curation (not in this folder, but part of the bake chain)

These scripts live elsewhere in the repository tree because they're one-shot infrastructure that produced the model artefacts we start from. They are documented here for provenance.

### 2.1 Gene ID mapping — `phase0_audit/`

ENSDARG gene IDs (used in zebrafish transcriptomics) map to the model's internal gene names via `data/gene_map_ensdarg_to_model.csv` (produced by `phase0_audit/s03_map_genes.py`). This mapping is used downstream for context-specific extraction; it does **not** affect the base model.

### 2.2 O₂ / H₂O₂ leak block — `gimme_transcriptomics_pipeline/results/pipeline/o2_leak_fixes.csv`

**Why**: in Wang's reconstruction, several reactions can produce O₂ or H₂O₂ in directions that violate thermodynamics but satisfy mass balance. The LP exploits these as "free energy" loops when maximizing growth with wide-open exchanges.

**How**: ~20 reactions identified by curated review of superoxide / peroxide / oxygen-species producers; each gets its forward or reverse direction blocked (`block_fwd`, `block_rev`, or `block_both`). The actions are tabulated in `o2_leak_fixes.csv` and applied inside `01_bake_ocr.m` (not as a separate step).

**Effect on biomass**: minimal under constrained bounds; essential under wide-open bounds. Blocks the "ATP from H₂O" pathology.

### 2.3 Cofactor pool modification — `phase2_metabolomics/s04_modify_cofactor_pool.m`

**Why**: Wang's `MAR00022` (the cofactor-pool-formation reaction that produces `MAM01602c`, "cofactors and vitamins") has 16 substrates at coefficient `−1` each, meaning growth consumes 1 : 1 of each of them. Several of those (NADH, NADPH, FADH2, CoA, L-carnitine, riboflavin, THF, BH4, ubiquinol, cytochrome-C) are **catalytically cycled** in vivo — they turn over thousands of times per cell cycle and are not consumed stoichiometrically with new biomass. Treating them as 1 : 1 consumables makes the pool a hard growth ceiling that is biologically meaningless.

**How** (minimal change to Wang's stoichiometry): the current pipeline reduces **only the cytochrome-C coefficient** from `−1` to `−0.1` (see `results/qc/phase2_modcofpool.txt`). This is the single most impactful change: in the un-constrained Wang model, cyto-C is the dominant pool substrate (removing it alone lifts biomass_max from 9.04 to 500 — a 55-fold change). The remaining 15 substrates keep their Wang coefficients.

A more aggressive variant (an earlier branch of the pipeline, now superseded) stripped 9 cycled cofactors entirely. **The current pipeline deliberately keeps all 16 pool members at their Wang coefficients, except cyto-C, because**:

- The pool as a whole is a *modelling convention*, not a biological measurement. Changing many coefficients at once introduces unexamined degrees of freedom.
- The subsequent v7 overlay (`03_apply_v7_boundary.m`) addresses the pool-mass issue at a different layer — by setting the pool's coefficient in `MAR00021` to `ε = 1e-3`, which makes the pool itself a trace demand regardless of what's in it (an ε scan confirms biomass_max is flat across ε ∈ {1e-2, 1e-3, 1e-4}).
- Four of the 16 pool substrates are reduced-redox states (FADH2, NADH, NADPH, ubiquinol) that are also set to coefficient 0 by the v7 overlay for the same reason.

So the base-model file `zebrafishGEM_v2_modcofpool.mat` is **Wang's `MAR00022` with cyto-C coef `-1 → -0.1`, and nothing else**. The deeper pool simplification happens later, downstream of the OCR bake.

### 2.4 Phase 2 metabolomics envelope (not applied here)

For completeness: the metabolomics-informed exchange envelope (105 measured exchanges × 15 condition × hpf groups, fold-change × class-share × class-mass-fraction × μ chain → `data/phase2_exchange_bounds_long.csv`) is **not applied to the base model**. The exchange envelope is applied **per cell, post-extraction**, via `apply_context_bounds` in `functions/`. See `data/phase2_exchange_bounds_long.csv` for the fc and flux columns and `FAMILY_PREREG.md` §3 for the derivation.

---

## 3. `01_bake_ocr.m` — the OCR bake (five parents, one per hpf)

This is the only script in this folder that produces one of the five `baked_ocr_{hpf}.mat` outputs. Its job is to take `zebrafishGEM_v2_modcofpool.mat` and freeze it into **one baseline per developmental stage** with:

- O₂ / H₂O₂ leak fixes applied (step 2.2 above).
- Every exchange reaction uniformly bounded at `±10`.
- O₂ exchange `MAR09048` capped at `−OCR(hpf)` on the uptake side, `0` on secretion.
- A reproducibility `baked_meta` struct saved alongside the model, recording the source file, number of leak fixes applied, exchange bound convention, O₂ bound, and the biomass_max observed under the baseline (sanity check that the model is feasible).

### 3.1 Why ±10 for every exchange

**Scale convention**. ±10 is **not** a measured flux value — it is a dimensionless envelope radius that the metabolomics layer (post-extraction, via `apply_context_bounds`) will later multiply by a per-cell fold-change (fc). An exchange at fc = 1.0 gets `lb = -10, ub = +10`; an exchange at fc = 0.5 gets `lb = -5, ub = +5`. This is a **relative-units convention** — only cross-cell ratios at a fixed BOUND_EX are physically interpretable. See `FAMILY_PREREG.md` §10.5 for the full statement.

Using ±1000 (Wang's default) would let the LP find mass-balance-feasible but biologically wrong solutions. Using the measured absolute flux at this step would conflate the bake with the metabolomics envelope, which should be applied per cell. ±10 is the midpoint: tight enough to force the LP to respect the envelope, loose enough that the fc multiplier does the per-cell work.

### 3.2 Why the OCR cap specifically

**O₂ uptake is the single most constrained flux in the embryo**, and it has been measured directly. The published OCR for zebrafish embryos:

| hpf | OCR (µmol O₂ · g DW⁻¹ · hr⁻¹) |
|---:|---:|
| 24  |  96.79 |
| 48  | 173.35 |
| 72  | 226.63 |
| 96  | 273.07 |
| 120 | 313.00 |

(Source: `Data/Metabolomics/physiology_by_timepoint.csv`; measurements averaged per hpf across replicates.)

Capping O₂ at the measured value does three things:

1. **Prevents the LP from "breathing" at an unphysical rate** — Wang's default ±1000 would let the embryo oxidize anything.
2. **Couples oxidative phosphorylation output to a measured anchor** — downstream analyses can trust that ATP yields are tethered to actual gas exchange, not to a modelling default.
3. **Varies the oxidative ceiling with developmental stage** — a 24 hpf embryo doesn't have the same mitochondrial capacity as a 120 hpf embryo. Baking per hpf captures this without any per-condition information.

Note the sign: `MAR09048` in Wang convention is `MAM02630e → ∅` (secretion-positive). "Uptake OCR" of 97 µmol means `lb = -97, ub = 0` — the exchange can only consume oxygen, up to that rate.

### 3.3 Why 5 separate files rather than 1 with a flexible OCR

**Reproducibility**. Each cell line's analysis (BL_24, D_24, LD_24 at 24 hpf; BL_48, …) uses `baked_ocr_24.mat` as its parent. Loading one file per hpf makes the parent fully determined by the file path — no runtime OCR arithmetic. The 5 files are byte-identical except for `MAR09048`'s `lb` and the `baked_meta.ocr_umol_O2_gDW_hr` field.

### 3.4 Reproducibility of this step

Rerunning `01_bake_ocr.m` from `zebrafishGEM_v2_modcofpool.mat` + `o2_leak_fixes.csv` + `physiology_by_timepoint.csv` produces the five output files deterministically. The script sets no random seeds (there's no stochasticity), uses MATLAB sparse ops for the stoichiometry, and relies on COBRA Toolbox + Gurobi for the sanity-check biomass_max solve (any LP solver works; the biomass_max value is reported for provenance, not consumed downstream).

### 3.5 What this does NOT include

- No transcriptomics (that's extraction, step 4).
- No metabolomics envelope on measured exchanges (that's `apply_context_bounds`, applied post-extraction).
- No ledger closes, no Part 1 / Part 2 boundary, no ε pool coefficient (all of that is `03_apply_v7_boundary.m`).
- No gene-expression-driven pruning.

The output of `01_bake_ocr.m` is the **cleanest possible version of Wang's zebrafish GEM**: thermodynamically consistent, physiologically bounded on oxygen, and ready to accept condition- and stage-specific curation on top.

---

## 4. Verification

The script saves a `baked_meta` struct inside each `baked_ocr_{hpf}.mat` with:

```matlab
baked_meta.source               = 'zebrafishGEM_v2_modcofpool.mat'
baked_meta.o2_h2o2_fix_csv      = (path)
baked_meta.n_o2_h2o2_fixes      = N (number of leak fixes applied)
baked_meta.exchange_bound       = 10
baked_meta.n_exchanges          = M (total exchange reactions)
baked_meta.hpf                  = {24, 48, 72, 96, 120}
baked_meta.ocr_umol_O2_gDW_hr   = {96.79, 173.35, 226.63, 273.07, 313.00}
baked_meta.o2_rxn               = 'MAR09048'
baked_meta.biomass_max          = (sanity FBA result)
```

The baseline `biomass_max` from `optimizeCbModel(changeObjective(m, 'MAR00021'), 'max')` is approximately **9.5 across all 5 hpf** (bounded by amino acid uptakes at ±10, not by O₂). This is the number that gets lifted to ~10 by the v7 overlay's ε = 1e-3 pool coefficient change in `03_apply_v7_boundary.m`.

---

## 5. One-paragraph summary for a methods section

> *Context-specific extraction was performed on five hpf-specific parent models derived from the Wang zebrafish GEM v2 (a zebrafish adaptation of Human1 / Human-GEM). The published Wang reconstruction was first corrected for thermodynamically infeasible oxygen and hydrogen-peroxide pathways by applying 18 curated bound restrictions (`o2_leak_fixes.csv`), after which the cofactor pool reaction `MAR00022` had its cytochrome-C coefficient reduced from −1 to −0.1 (justified by cyto-C's role as a 1:1 non-consumed structural pool substrate; `phase2_modcofpool.txt`). Every exchange reaction was then bounded uniformly at ±10 as a dimensionless envelope radius (the metabolomics fold-change applied post-extraction defines the physical scale). Finally, five baselines were saved — one per developmental timepoint — with the O₂ exchange `MAR09048` capped at the measured oxygen consumption rate at that stage (24 hpf: 96.79, 48: 173.35, 72: 226.63, 96: 273.07, 120: 313.00 µmol O₂ g DW⁻¹ hr⁻¹; `physiology_by_timepoint.csv`). These five files, `baked_ocr_{24,48,72,96,120}.mat`, are the input to the v7 boundary rule overlay (ledger closes, Part 1 secretion uncap, Part 2 unmeasured uptake close, medium whitelist, biomass pool coefficient ε = 1 × 10⁻³; see `FAMILY_PREREG.md` §10.3) and ultimately to rFASTCORMICS extraction per cell.*
