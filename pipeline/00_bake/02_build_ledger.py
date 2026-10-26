"""
s130 -- v7 curation ledger builder.

Builds Step 1's biomass-demand ledger for the v7 rebuild. Produces:

  s130_v7_ledger.csv       one row per (MAM species, exchange MAR id).
                           class ∈ {1,2,3}, action ∈ {keep_exchange,
                           close_exchange}. Includes rows with no exchange
                           (substrates internally produced only).
  s130_envelope_leaks.csv  species where >=2 exchanges exist AND at least one
                           is in the metabolomics 105 while others are not.
                           These bypass the measured bound; per Q2 the
                           unscaled duplicates are closed like class 3.

Inputs:
  - ocr_anchored_extraction/models/baked_ocr_72.mat   (base parent)
  - reviewer_packet_v4_slim/s128c_candidates_v5.csv   (v5 flagged exchanges)
  - results/qc/phase2_exchange_bounds_long.csv        (105 measured exchanges)

Classification rules:

  Class 1 -- keep_exchange. True vitamin / yolk-derived essential nutrient
     the embryo cannot synthesise. Also any exchange in the measured 105
     (override, per Q2).
  Class 2 -- close_exchange. Embryo can synthesise from a class-1 precursor;
     pathway exists in Zebrafish1 base model. Includes FAD/FMN chain,
     calcidiol/calcitriol, lipoic acid, NAD(P)(H), CoA, THF, cobamide.
  Class 3 -- close_exchange. Boundary artefact: pool metabolites, complex
     lipoprotein assemblies, acyl-CoAs, pseudo-metabolites, xenobiotics
     that shouldn't cross an embryo boundary.
  \"biomass\" pseudo-metabolite (MAM03971): ignored, not written to ledger.

Cycled cofactors in the pool (NADH/NADPH/CoA/etc.) after the Phase-2 trim
still appear as MAR00021/22 substrates in the base parent (the trim edits
MAR00022 stoichiometry, not membership). net_consumed=False for those; the
biomass reaction just recycles them.
"""
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
from pathlib import Path

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
BASE_PARENT = REPO / "ocr_anchored_extraction/models/baked_ocr_72.mat"
S128C_CANDS = REPO / "reviewer_packet_v4_slim/s128c_candidates_v5.csv"
BOUNDS_CSV  = REPO / "results/qc/phase2_exchange_bounds_long.csv"
OUT_LEDGER  = REPO / "reviewer_packet_v4_slim/s130_v7_ledger.csv"
OUT_LEAKS   = REPO / "reviewer_packet_v4_slim/s130_envelope_leaks.csv"
OUT_POOL    = REPO / "reviewer_packet_v4_slim/s130_pool_table.csv"

# 16-row pool decisions (per user 2026-09-29: net_consumed rule).
# Reduced redox states (NADH, NADPH, FADH2, ubiquinol) get coef=0.
# Everything else (prosthetic groups, cofactor pools, stored vitamins,
# small-molecule dietary essentials, structural cyt-C) stays.
POOL_DECISIONS = {
    "MAM00209c": (True,  "[protein]-N6-(lipoyl)lysine — prosthetic group; per-cell accumulating"),
    "MAM01401c": (True,  "biotin — carboxylase prosthetic; per-cell accumulating"),
    "MAM01597c": (True,  "CoA — coenzyme pool; per-cell needed from pantothenate"),
    "MAM01600m": (True,  "cobamide-coenzyme — B12-derived prosthetic"),
    "MAM01631m": (True,  "cytochrome-C — structural protein; coef already 0.1"),
    "MAM01803c": (False, "FADH2 — reduced redox state; never the demand (see riboflavin)"),
    "MAM02348c": (True,  "L-carnitine — small pool for FA transport, per-cell accumulating"),
    "MAM02553c": (False, "NADH — reduced redox state; net accumulating pool is NAD+ from niacin"),
    "MAM02555c": (False, "NADPH — reduced redox state; net accumulating pool is NADP+ from niacin"),
    "MAM02842c": (True,  "riboflavin — vitamin, feeds flavin pool (FMN/FAD synthesis)"),
    "MAM02978c": (True,  "tetrahydrobiopterin — cofactor pool"),
    "MAM02980c": (True,  "THF — one-carbon carrier cofactor"),
    "MAM03102m": (False, "ubiquinol — reduced redox state; ETC turnover, not accumulating demand"),
    "MAM03139c": (True,  "vitamin A derivatives — stored vitamin"),
    "MAM03140c": (True,  "vitamin D derivatives — stored vitamin"),
    "MAM03143c": (True,  "vitamin E derivatives — stored vitamin"),
}

BIOMASS_RXN = "MAR00021"
POOL_RXN    = "MAR00022"


# ---- HDF5 loader (same shape as s128c) -----------------------------------
def _resolve(f, x):
    if isinstance(x, h5py.Reference):
        return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object:
            return _resolve(f, x[0])
        try:
            return "".join(chr(int(c)) for c in x)
        except Exception:
            return ""
    return str(x)


def load_model(fp, root_key):
    with h5py.File(fp, "r") as f:
        m = f[root_key]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        metNames = [_resolve(f, r) for r in np.asarray(m["metNames"]).ravel()]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        S = sp.csc_matrix((data, ir, jc), shape=(int(ir.max())+1, len(jc)-1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
        try:
            rxnNames = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
        except Exception:
            rxnNames = [""] * len(rxns)
    return dict(rxns=rxns, mets=mets, metNames=metNames, rxnNames=rxnNames,
                S=S, lb=lb, ub=ub)


def strip_comp(mid):
    if len(mid) > 1 and mid[-1] in "cemnrlxsg":
        return mid[:-1]
    return mid


# ---- Classification tables (species-level; extracellular exchanges) ------
# Class 1: cannot be synthesised by animal metabolism, must be dietary
CLASS1_SPECIES = {
    # Vitamins
    "MAM01361": "vitamin B12 (aquacob(III)alamin) — not synthesised by animals",
    "MAM01362": "vitamin B12 variants",
    "MAM01600": "cobamide-coenzyme — B12 derivative, precursor is dietary B12",
    "MAM02842": "riboflavin — vitamin B2, dietary",
    "MAM02834": "retinol — vitamin A, dietary (yolk)",
    "MAM03139": "vitamin A derivatives (pool)",
    "MAM01418": "cholecalciferol (vitamin D3) — no UVB conversion under any of BL/D/LD, yolk-derived",
    "MAM03140": "vitamin D derivatives (pool) — upstream of calcidiol",
    "MAM03143": "vitamin E derivatives",
    "MAM03144": "vitamin K",
    "MAM02980": "THF — folate-derived; dietary folate",
    "MAM02976": "folate",
    "MAM03127": "thiamin (B1)",
    "MAM02585": "niacin (B3)",
    "MAM02685": "pantothenate (B5)",
    "MAM02812": "pyridoxine (B6)",
    "MAM01401": "biotin (B7)",
    "MAM00209": "protein-N6-(lipoyl)lysine — downstream of lipoate but treated as scaffold",
    "MAM02348": "L-carnitine — dietary in fish embryo (yolk)",
    "MAM02978": "tetrahydrobiopterin (BH4)",
    "MAM03102": "ubiquinol / ubiquinone — fat-soluble; keep as dietary until synthesis confirmed",
    # Essential amino acids
    "MAM02125": "histidine — essential AA",
    "MAM02184": "isoleucine — essential AA",
    "MAM02360": "leucine — essential AA",
    "MAM02426": "lysine — essential AA",
    "MAM02471": "methionine — essential AA",
    "MAM02724": "phenylalanine — essential AA",
    "MAM02993": "threonine — essential AA",
    "MAM03089": "tryptophan — essential AA",
    "MAM03135": "valine — essential AA",
    "MAM01365": "arginine — conditionally essential in embryo",
    # O2 / H2O / minerals — boundary specification, not curation
    "MAM02630": "O2",
    "MAM02040": "H2O",
    "MAM02039": "H+",
    "MAM02209": "HCO3-",
    "MAM01596": "CO2",
    "MAM01596c": "CO2",
    "MAM02751": "Pi (inorganic phosphate)",
    "MAM01821": "Fe2+",
    "MAM02467": "Mg2+",
    "MAM01829": "Fe3+",
    "MAM02904": "sodium",
    "MAM02393": "K+",
    "MAM02875": "sulfate",
    "MAM01986": "glycine — non-essential but ambiguous in embryo; kept as class 1 default",
    "MAM01974": "glutamate",
    "MAM01975": "glutamine",
    "MAM02896": "serine",
    "MAM02770": "proline",
    "MAM02986": "tyrosine — conditionally essential from phenylalanine",
    "MAM03101": "tyrosine",
    "MAM01307": "alanine",
    "MAM01369": "asparagine",
    "MAM01370": "aspartate",
    "MAM01628": "cysteine",
}

# Class 2: embryo can synthesise from a class-1 precursor
CLASS2_SPECIES = {
    "MAM01802": "FAD — synthesised from riboflavin via RFK+FLAD1 (retired v6 forcing = special case of this rule)",
    "MAM01803": "FADH2 — FAD redox partner",
    "MAM01828": "FMN — riboflavin → FMN via RFK",
    "MAM02828": "FMN (extracellular variant)",
    "MAM01415": "calcidiol — CYP2R1/CYP27A1 hydroxylation of vitamin D3",
    "MAM01417": "calcitriol — CYP27B1 hydroxylation of calcidiol",
    "MAM02394": "lipoic acid — LIAS synthesis from octanoyl-ACP",
    "MAM02552": "NAD",
    "MAM02554": "NADP",
    "MAM02555": "NADPH",
    "MAM01597": "CoA — synthesised from pantothenate; accumulating pool per cell so pool coef stays, but exchange closes because CoA doesn't cross membranes as CoA (pantothenate uptake covers it)",
    "MAM01722": "DNA-5-methylcytosine — internal",
    "MAM01631": "cytochrome-C — internal",
    "MAM01602": "cofactors and vitamins (pseudo) — internal aggregation",
    "MAM01450": "cholesterol — mevalonate pathway exists; close exchange. Note: yolk is an alternate source in early embryo but modelling as synthesis is cleaner and consistent with the general rule",
}

# NOTE: earlier draft mislabeled MAM02685 as "pantothenate" and added it as
# class 1 keep_exchange. That was wrong — MAM02685 is the PE-LD pool; real
# pantothenate is MAM02680 (exchange MAR09145), which is already in the 105
# measured set. MAM02685 belongs in class 3 (see CLASS3_SPECIES below).

# Manual additions from s137 minimal-medium result (2026-09-29):
# All 13 support members classified per user's decision. Plus LDL moved
# out of class 3 into the lipoprotein group.
MANUAL_ROWS = [
    # class 1, medium_uptake -- vitamin/carotenoid (not synthesised by animals)
    dict(MAM_id="MAM01330", exchange_MAR_id="MAR09152", cls=1, action="medium_uptake",
         justification="alpha-tocotrienol — vitamin E form, yolk-derived, not synthesisable"),
    dict(MAM_id="MAM01938", exchange_MAR_id="MAR09154", cls=1, action="medium_uptake",
         justification="gamma-tocotrienol — vitamin E form, yolk-derived, not synthesisable"),
    dict(MAM_id="MAM01385", exchange_MAR_id="MAR09276", cls=1, action="medium_uptake",
         justification="beta-carotene — vitamin A precursor, yolk-derived, not synthesisable de novo"),
    # class 1, ysl_lipid_delivery -- yolk syncytial layer lipoprotein packaging
    # (Fraher et al. 2016 Cell Rep): apoB / apoA lipoproteins secreted into
    # embryonic circulation. Rate unmeasured, default bound.
    dict(MAM_id="MAM03147", exchange_MAR_id="MAR09049", cls=1, action="ysl_lipid_delivery",
         justification="VLDL — YSL yolk-lipid delivery (Fraher 2016); default bound; unmeasured"),
    dict(MAM_id="MAM01570", exchange_MAR_id="MAR09024", cls=1, action="ysl_lipid_delivery",
         justification="chylomicron — YSL yolk-lipid delivery (Fraher 2016); default bound; unmeasured"),
    dict(MAM_id="MAM02048", exchange_MAR_id="MAR09050", cls=1, action="ysl_lipid_delivery",
         justification="HDL — YSL yolk-lipid delivery (Fraher 2016); default bound; unmeasured"),
    dict(MAM_id="MAM01350", exchange_MAR_id="MAR04110", cls=1, action="ysl_lipid_delivery",
         justification="apoA1 — YSL apolipoprotein, HDL scaffold (Fraher 2016); default bound; unmeasured"),
    dict(MAM_id="MAM01351", exchange_MAR_id="MAR11987", cls=1, action="ysl_lipid_delivery",
         justification="apoB100 — YSL apolipoprotein, VLDL/LDL scaffold (Fraher 2016); default bound; unmeasured"),
    dict(MAM_id="MAM02353", exchange_MAR_id="MAR09051", cls=1, action="ysl_lipid_delivery",
         justification="LDL — reopened to treat all lipoprotein exchanges identically (Fraher 2016 yolk lipid delivery); default bound; unmeasured"),
    # class 2, close_exchange -- synthesisable with GPR present
    dict(MAM_id="MAM02049", exchange_MAR_id="MAR09107", cls=2, action="close_exchange",
         justification="heme — synthesised via ALAS→heme biosynthesis pathway; close uptake"),
    dict(MAM_id="MAM02171", exchange_MAR_id="MAR09361", cls=1, action="medium_uptake",
         justification="myo-inositol reclassified class 1 (2026-09-29) after closure failure; de novo ISYNA1/IMPA route not verified in this base (producibility test pending, non-blocking); dietary requirement documented for several teleosts (NRC 2011)"),
    # class 3, close_exchange -- boundary artefacts (alternates L1 picked)
    dict(MAM_id="MAM02044", exchange_MAR_id="MAR04119", cls=3, action="close_exchange",
         justification="haptoglobin — plasma protein artefact, not physiologically uptaken by embryo"),
    dict(MAM_id="MAM02173", exchange_MAR_id="MAR10263", cls=3, action="close_exchange",
         justification="inositol-1-phosphate — intracellular intermediate, exchange is boundary artefact"),
    dict(MAM_id="MAM03597", exchange_MAR_id="MAR10554", cls=3, action="close_exchange",
         justification="Glutaminyl-Histidyl-Histidine — tripeptide, boundary artefact; His demand covered by MAR09038"),
]
# LDL was originally class 3 in the auto-classification via CLASS3_SPECIES;
# remove it there so MANUAL_ROWS's ysl_lipid_delivery entry wins.
# (Handled by adding MANUAL_MAM_ids set — see below.)
MANUAL_MAM_IDS = {r["MAM_id"] for r in MANUAL_ROWS}

# Class 3: boundary artefact
CLASS3_SPECIES = {
    "MAM02553": "NADH — reduced redox state, exchange is modelling artefact (real accumulating pool is NAD+ from niacin)",
    "MAM02685": "PE-LD pool — phosphatidylethanolamine-lipid-droplet pool, assembly reaction internal; exchange shouldn't exist (2026-09-30 correction; was mislabeled as pantothenate)",
    "MAM02750": "PI pool — assembly reaction internal; exchange shouldn't exist",
    "MAM02715": "PG-CL pool — phosphatidylglycerol/cardiolipin pool",
    "MAM01589": "CL pool",
    "MAM02908": "SM pool — sphingomyelin pool",
    "MAM01451": "cholesterol-ester pool",
    "MAM02733": "phosphatidate-LD-TAG pool",
    "MAM02392": "lipid droplet",
    "MAM03161": "glycogen — synthesised internally",
    # MAM02353 (LDL) reclassified to class 1 ysl_lipid_delivery via MANUAL_ROWS (2026-09-29)
    "MAM02352": "LDL remnant",
    "MAM01396": "bilirubin — heme degradation product, boundary junk in embryo",
    "MAM01397": "bilirubin-bisglucuronoside — phase-II excretion product",
    "MAM03109": "UDP-glucuronate — activated sugar, internal only",
    "MAM02122": "hexanoyl-CoA — internally activated from hexanoate",
    "MAM01721": "DNA — internal biomass component",
    "MAM02847": "RNA — internal biomass component",
}

# Cycled cofactors (Phase-2 trim removed from pool but MAR00021 still lists)
CYCLED_COFACTORS = {"MAM02552", "MAM02553", "MAM02554", "MAM02555",
                    "MAM01597", "MAM02040", "MAM02039", "MAM01596",
                    "MAM02751"}

IGNORE_SPECIES = {"MAM03971"}  # biomass pseudo-metabolite (sink)


def classify(species, measured_in_105_species):
    if species in IGNORE_SPECIES:
        return None, "ignore — biomass sink pseudo-metabolite"
    if species in measured_in_105_species:
        base_class, base_just = 1, "measured in the metabolomics 105 (override)"
        if species in CLASS2_SPECIES:
            base_just += "; would be class 2 (" + CLASS2_SPECIES[species] + ")"
        elif species in CLASS3_SPECIES:
            base_just += "; would be class 3 (" + CLASS3_SPECIES[species] + ")"
        return base_class, base_just
    if species in CLASS1_SPECIES:
        return 1, CLASS1_SPECIES[species]
    if species in CLASS2_SPECIES:
        return 2, CLASS2_SPECIES[species]
    if species in CLASS3_SPECIES:
        return 3, CLASS3_SPECIES[species]
    return None, "unclassified — review manually"


# ---- Build ---------------------------------------------------------------
print("loading base parent …")
base = load_model(BASE_PARENT, "m")
Snnz = np.diff(base["S"].tocsc().indptr)
is_ex = Snnz == 1
# exchange metabolite per ex rxn
ex_of_met = {}   # species -> [(ex_rxn, mm_id)]
mm_of_ex  = {}   # ex_rxn -> mm_id
for j in np.where(is_ex)[0]:
    col = base["S"].getcol(j).toarray().ravel()
    nz = np.where(col != 0)[0]
    if len(nz) != 1:
        continue
    mm = base["mets"][nz[0]]
    sp_ = strip_comp(mm)
    ex_of_met.setdefault(sp_, []).append((base["rxns"][j], mm))
    mm_of_ex[base["rxns"][j]] = mm

def stoich(rxn_id):
    j = base["rxns"].index(rxn_id)
    col = base["S"].getcol(j).toarray().ravel()
    return {base["mets"][i]: float(col[i]) for i in np.where(col != 0)[0]}
bio  = stoich(BIOMASS_RXN)
pool = stoich(POOL_RXN)
bio_subs  = {mm for mm, c in bio.items()  if c < 0}
pool_subs = {mm for mm, c in pool.items() if c < 0}

# metabolomics 105
mb = pd.read_csv(BOUNDS_CSV)
measured_ex_rxns = set(mb.ex_rxn.unique())
measured_species = {strip_comp(mm_of_ex[r]) for r in measured_ex_rxns
                    if r in mm_of_ex}
print(f"metabolomics 105: {len(measured_ex_rxns)} ex_rxn, "
      f"{len(measured_species)} distinct species")

# s128c candidates
cands = pd.read_csv(S128C_CANDS)
cand_species = set(cands.species.unique())

# universe of species to include in ledger
species_union = set()
for mm in bio_subs  | pool_subs:
    species_union.add(strip_comp(mm))
species_union |= cand_species
# Manual additions (rule-level, not demand-driven): precursors that need to
# stay class-1 open because a downstream product's exchange has been closed.
species_union.add("MAM02685")  # pantothenate — CoA precursor, keep open since MAR04256 closed

name_of = dict(zip(base["mets"], base["metNames"]))

rows = []
# Species handled by MANUAL_ROWS are skipped in the auto loop to avoid conflict
for species in sorted(species_union):
    if species in MANUAL_MAM_IDS:
        continue
    # find canonical name
    name = "?"
    for mm in base["mets"]:
        if strip_comp(mm) == species:
            name = name_of.get(mm, "?"); break

    cls, just = classify(species, measured_species)
    if cls is None:
        # write anyway with cls = "?" so the user sees it
        cls = "?"
    action = "keep_exchange" if cls == 1 else "close_exchange"
    if species in IGNORE_SPECIES:
        continue

    # in biomass? in pool? net_consumed?
    in_bio  = species in {strip_comp(mm) for mm in bio_subs}
    in_pool = species in {strip_comp(mm) for mm in pool_subs}
    net_consumed = not (species in CYCLED_COFACTORS)
    from_s128c = species in cand_species

    # exchange rows
    exs = ex_of_met.get(species, [])
    if not exs:
        rows.append(dict(metabolite=name, MAM_id=species,
                         exchange_MAR_id="",
                         in_biomass=in_bio, in_pool=in_pool,
                         from_s128c=from_s128c,
                         net_consumed=net_consumed,
                         **{"class": cls}, action="none_no_exchange",
                         measured_in_105=False,
                         justification=just))
    else:
        for ex_rxn, mm_id in sorted(exs):
            measured = ex_rxn in measured_ex_rxns
            row_action = action
            # per Q2: envelope-leak duplicates get closed
            if any(other != ex_rxn and other in measured_ex_rxns
                   for other, _ in exs) and not measured and cls == 1:
                row_action = "close_exchange"
                row_just = just + "; envelope-leak duplicate: another exchange of this species is measured, this one bypasses the bound"
            else:
                row_just = just
            rows.append(dict(metabolite=name, MAM_id=species,
                             exchange_MAR_id=ex_rxn,
                             in_biomass=in_bio, in_pool=in_pool,
                             from_s128c=from_s128c,
                             net_consumed=net_consumed,
                             **{"class": cls}, action=row_action,
                             measured_in_105=measured,
                             justification=row_just))

# Append MANUAL_ROWS to ledger
for r in MANUAL_ROWS:
    mm_full = None
    for candidate in (r["MAM_id"] + "e", r["MAM_id"] + "c"):
        if candidate in base["mets"]:
            mm_full = candidate; break
    metname = name_of.get(mm_full, "?") if mm_full else "?"
    measured = r["exchange_MAR_id"] in measured_ex_rxns
    rows.append(dict(
        metabolite=metname, MAM_id=r["MAM_id"],
        exchange_MAR_id=r["exchange_MAR_id"],
        in_biomass=False, in_pool=False,
        from_s128c=r["MAM_id"] in cand_species,
        net_consumed=True,
        **{"class": r["cls"]}, action=r["action"],
        measured_in_105=measured,
        justification=r["justification"] + " (manual s137)"))

ledger = pd.DataFrame(rows).sort_values(
    ["class","MAM_id","exchange_MAR_id"], kind="stable")
ledger.to_csv(OUT_LEDGER, index=False)
print(f"wrote {OUT_LEDGER.name} ({len(ledger)} rows)")

# ---- 16-row pool table (MAR00022 decisions) ------------------------------
pool_rows = []
for mm, coef in [(m, c) for m, c in pool.items() if c < 0]:
    keep, rationale = POOL_DECISIONS.get(mm, (True, "unknown — kept by default"))
    pool_rows.append(dict(MAM_id=mm, metabolite=name_of.get(mm, "?"),
                          current_coef=coef,
                          new_coef=coef if keep else 0.0,
                          net_consumed=keep, rationale=rationale))
pool_df = pd.DataFrame(pool_rows).sort_values("MAM_id")
pool_df.to_csv(OUT_POOL, index=False)
print(f"wrote {OUT_POOL.name} ({len(pool_df)} rows; "
      f"{int((~pool_df.net_consumed).sum())} to be zeroed)")

# ---- Envelope leaks ------------------------------------------------------
leaks = []
for species, exs in ex_of_met.items():
    if len(exs) < 2:
        continue
    measured_here = {r for r, _ in exs if r in measured_ex_rxns}
    unmeasured_here = {r for r, _ in exs if r not in measured_ex_rxns}
    if measured_here and unmeasured_here:
        cls, _ = classify(species, measured_species)
        name = name_of.get(species+"c", name_of.get(exs[0][1], "?"))
        for ex_rxn in sorted(unmeasured_here):
            leaks.append(dict(metabolite=name, MAM_id=species,
                              unscaled_exchange=ex_rxn,
                              measured_partner=next(iter(sorted(measured_here))),
                              **{"class": cls}, action="close_exchange"))
leaks_df = pd.DataFrame(leaks)
leaks_df.to_csv(OUT_LEAKS, index=False)
print(f"wrote {OUT_LEAKS.name} ({len(leaks_df)} rows)")

# ---- Summary -------------------------------------------------------------
print("\n=== ledger summary ===")
print(ledger.groupby(["class","action"]).size().to_string())
print("\n=== unclassified species (need manual review) ===")
u = ledger[ledger["class"] == "?"]
if len(u):
    print(u[["MAM_id","metabolite","in_biomass","in_pool","from_s128c",
             "exchange_MAR_id","measured_in_105"]].drop_duplicates().to_string(index=False))
else:
    print("(none)")

print("\n=== envelope leaks ===")
if len(leaks_df):
    print(leaks_df.to_string(index=False))
else:
    print("(none)")

print("\n=== class-2 close_exchange rows (these are the retired-forcing candidates) ===")
c2 = ledger[(ledger["class"] == 2) & (ledger.action == "close_exchange")]
print(c2[["MAM_id","metabolite","exchange_MAR_id","justification"]].to_string(index=False))
