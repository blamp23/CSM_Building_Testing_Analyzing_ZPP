"""
s128c -- generalised uptake-essentiality audit.

Motivation. s128b enumerated only species directly consumed by MAR00021 or
MAR00022 and missed FAD, because the pool substrate is FADH2 (MAM01803) not
FAD (MAM01802). The reviewer's ask is more general: any exchanged metabolite
whose internal synthesis has been pruned is a candidate artefact, whether or
not it appears literally in the biomass stoichiometry.

Test. For every extracted cell, and every exchange reaction in that cell,
knock out that exchange (lb=ub=0) and recompute biomass_max under
sampling-style bounds. Ratio = biomax_ko / baseline.
    ratio < 0.01 -> uptake is essential -> ARTEFACT CANDIDATE
    ratio > 0.95 -> synthesis exists, safe
    between      -> partial

Cross-referenced afterwards with:
    - MAR00021 substrate list (direct biomass precursor)
    - MAR00022 substrate list (pool member)
    - chemical family lookup (e.g. FAD <-> FADH2 <-> FMN <-> riboflavin)
so downstream reporting can distinguish "direct precursor artefact" from
"upstream cofactor artefact".
"""
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.optimize import linprog
from pathlib import Path
import argparse

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
BASE_PARENT = REPO / "ocr_anchored_extraction/models/baked_ocr_72.mat"
BOUNDS_CSV  = REPO / "results/qc/phase2_exchange_bounds_long.csv"
BIOMASS_RXN = "MAR00021"
POOL_RXN    = "MAR00022"
O2_EX       = "MAR09048"
BOUND_EX    = 10.0
OCR = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}
CELLS = [f"{c}_{h}" for c in ("BL","D","LD") for h in (24,48,72,96,120)]

ap = argparse.ArgumentParser()
ap.add_argument("--models_dir", default="ocr_anchored_extraction/models")
ap.add_argument("--suffix", default="_v5")
ap.add_argument("--essential_thresh", type=float, default=0.01)
args = ap.parse_args()


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
        try:
            metNames = [_resolve(f, r) for r in np.asarray(m["metNames"]).ravel()]
        except Exception:
            metNames = [""] * len(mets)
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


def apply_sampling_bounds(m, cond, hpf, mb_df):
    # 2026-09-30 CORRECTION: use apply_context_bounds. Old body reset every
    # exchange to ±BOUND_EX before fc, which discarded the v7 boundary.
    import sys
    sys.path.insert(0, str(REPO / "functions"))
    from apply_context_bounds import apply_context_bounds
    lb = m["lb"].copy(); ub = m["ub"].copy()
    Snnz = np.diff(m["S"].tocsc().indptr); is_ex = Snnz == 1
    lb, ub, idx, is_ex = apply_context_bounds(lb, ub, m["rxns"], is_ex, cond, hpf)
    return lb, ub, is_ex, idx


def biomax(m, lb, ub, idx, A_eq=None, b_eq=None, c=None):
    n_rxn = len(m["rxns"]); n_met = m["S"].shape[0]
    if c is None:
        c = np.zeros(n_rxn); c[idx[BIOMASS_RXN]] = -1.0
    if A_eq is None:
        A_eq = m["S"].toarray()
    if b_eq is None:
        b_eq = np.zeros(n_met)
    res = linprog(c, A_eq=A_eq, b_eq=b_eq,
                  bounds=list(zip(lb, ub)), method="highs")
    return -res.fun if res.success else np.nan


# --- annotation from base parent -----------------------------------------
print(f"loading base parent {BASE_PARENT.name} ...")
base = load_model(BASE_PARENT, "m")
def stoichiometry(mm, rxn_id):
    j = mm["rxns"].index(rxn_id)
    col = mm["S"].getcol(j).toarray().ravel()
    return {mm["mets"][i]: float(col[i]) for i in np.where(col != 0)[0]}
bio_subs  = {mm for mm, cc in stoichiometry(base, BIOMASS_RXN).items() if cc < 0}
pool_subs = {mm for mm, cc in stoichiometry(base, POOL_RXN).items() if cc < 0}
bio_species  = {strip_comp(x) for x in bio_subs}
pool_species = {strip_comp(x) for x in pool_subs}
name_of = dict(zip(base["mets"], base["metNames"]))

# curated chemical families -- redox / phosphorylation partners that map onto
# the same pool substrate. Extend as needed; empty for anything not on the
# reviewer's line.
FAMILIES = {
    "MAM01802": "FAD/FADH2",   "MAM01803": "FAD/FADH2",
    "MAM02828": "FMN",         "MAM01828": "FMN",
    "MAM02842": "riboflavin",
    "MAM02553": "NAD/NADH",    "MAM02552": "NAD/NADH",
    "MAM02555": "NADP/NADPH",  "MAM02554": "NADP/NADPH",
}

# --- exchange KO across all cells ----------------------------------------
mb = pd.read_csv(BOUNDS_CSV)
rows = []
for cell in CELLS:
    cond, hpf = cell.split("_"); hpf = int(hpf)
    fp = REPO / args.models_dir / f"trans_rfastcormics_{cell}.mat"
    if not fp.exists():
        print(f"  {cell}: MISSING"); continue
    m = load_model(fp, "parent")
    lb, ub, is_ex, idx = apply_sampling_bounds(m, cond, hpf, mb)
    # cache dense S once per cell (dominant cost otherwise)
    A_eq = m["S"].toarray()
    b_eq = np.zeros(m["S"].shape[0])
    c_obj = np.zeros(len(m["rxns"])); c_obj[idx[BIOMASS_RXN]] = -1.0
    base_biomax = biomax(m, lb, ub, idx, A_eq=A_eq, b_eq=b_eq, c=c_obj)
    print(f"{cell}: baseline={base_biomax:.4f}  ({int(is_ex.sum())} exchanges)", flush=True)
    if not np.isfinite(base_biomax) or base_biomax < 1e-6:
        continue
    ex_js = np.where(is_ex)[0]
    for j in ex_js:
        # identify the exchanged metabolite
        col = m["S"].getcol(j).toarray().ravel()
        nz = np.where(col != 0)[0]
        if len(nz) != 1:
            continue
        mm_id = m["mets"][nz[0]]
        species = strip_comp(mm_id)
        # 2026-10-01 bookkeeping: separate uptake- vs secretion-essentiality.
        # uptake KO: lb -> 0 only. secretion KO: ub -> 0 only.
        lb_up = lb.copy(); lb_up[j] = 0.0
        bm_up = biomax(m, lb_up, ub, idx, A_eq=A_eq, b_eq=b_eq, c=c_obj)
        if not np.isfinite(bm_up): bm_up = 0.0
        r_up = bm_up / base_biomax
        uptake_essential = r_up < args.essential_thresh

        ub_sec = ub.copy(); ub_sec[j] = 0.0
        bm_sec = biomax(m, lb, ub_sec, idx, A_eq=A_eq, b_eq=b_eq, c=c_obj)
        if not np.isfinite(bm_sec): bm_sec = 0.0
        r_sec = bm_sec / base_biomax
        secretion_essential = r_sec < args.essential_thresh

        # combined (both bounds zero) kept for backwards compat
        lb2 = lb.copy(); ub2 = ub.copy()
        lb2[j] = 0.0; ub2[j] = 0.0
        bm = biomax(m, lb2, ub2, idx, A_eq=A_eq, b_eq=b_eq, c=c_obj)
        if not np.isfinite(bm): bm = 0.0
        ratio = bm / base_biomax
        essential = ratio < args.essential_thresh
        rows.append(dict(
            cell=cell, ex_rxn=m["rxns"][j], ex_name=m["rxnNames"][j],
            metabolite=mm_id, species=species,
            met_name=name_of.get(mm_id, name_of.get(species+"c", "?")),
            baseline_biomax=base_biomax, knockout_biomax=bm, ratio=ratio,
            uptake_only_biomax=bm_up, uptake_ratio=r_up,
            secretion_only_biomax=bm_sec, secretion_ratio=r_sec,
            uptake_essential=uptake_essential,
            secretion_essential=secretion_essential,
            essential=essential,
            direct_biomass_substrate=(mm_id in bio_subs) or (species in bio_species),
            pool_substrate=(mm_id in pool_subs) or (species in pool_species),
            chem_family=FAMILIES.get(species, ""),
        ))

df = pd.DataFrame(rows)
df.to_csv(f"s128c_exchange_essentiality{args.suffix}.csv", index=False)
print(f"\nwrote s128c_exchange_essentiality{args.suffix}.csv ({len(df)} rows)")

# 2026-10-01: split uptake-essential and secretion-essential tables.
cand_up = df[df.uptake_essential & ~df.secretion_essential].copy()
cand_up.to_csv(f"s128c_uptake_essential{args.suffix}.csv", index=False)
print(f"wrote s128c_uptake_essential{args.suffix}.csv ({len(cand_up)} rows) — gate uses this")

cand_sec = df[df.secretion_essential & ~df.uptake_essential].copy()
cand_sec.to_csv(f"s128c_secretion_essential{args.suffix}.csv", index=False)
print(f"wrote s128c_secretion_essential{args.suffix}.csv ({len(cand_sec)} rows) — byproduct sinks, listed separately")

cand_both = df[df.uptake_essential & df.secretion_essential].copy()
cand_both.to_csv(f"s128c_both_essential{args.suffix}.csv", index=False)
print(f"wrote s128c_both_essential{args.suffix}.csv ({len(cand_both)} rows) — direction-ambiguous")

cand = df[df.essential].copy()  # backwards compat (combined KO)
cand.to_csv(f"s128c_candidates{args.suffix}.csv", index=False)
print(f"wrote s128c_candidates{args.suffix}.csv ({len(cand)} rows) [combined KO, legacy column]")

print(f"\n==== ARTEFACT CANDIDATES (ratio < {args.essential_thresh}) ====")
if len(cand):
    show = cand[["cell","ex_rxn","metabolite","met_name",
                 "direct_biomass_substrate","pool_substrate","chem_family",
                 "baseline_biomax","knockout_biomax","ratio"]] \
        .sort_values(["cell","ex_rxn"])
    with pd.option_context("display.max_rows", 200, "display.width", 200):
        print(show.to_string(index=False))

    by_species = cand.groupby(["species","met_name","chem_family",
                               "direct_biomass_substrate","pool_substrate"]) \
        .size().reset_index(name="n_cells") \
        .sort_values("n_cells", ascending=False)
    by_species.to_csv(f"s128c_candidates_by_species{args.suffix}.csv", index=False)
    print(f"\n==== BY SPECIES (# cells uptake-essential) ====")
    print(by_species.to_string(index=False))
else:
    print("(none)")
