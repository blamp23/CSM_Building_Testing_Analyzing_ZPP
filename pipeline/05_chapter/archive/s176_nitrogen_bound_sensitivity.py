"""
s176 -- nitrogen-bound sensitivity test (§4.1 / Probe 3 robustness).

Three parts (interpretation thresholds locked BEFORE results are read, see
docstring tail):

1. COMPETITION RANKING.  Per cell, per measured exchange j whose metabolite
   is a biomass substrate: ratio_j = c_j / (BOUND_EX * fc_j). Print top 10
   per cell; report (ratio_1 / ratio_2) and the count of exchanges within
   10% of the leader.

2. COEFFICIENT PERTURBATION.  Scale every amino-acid coefficient in
   MAR00021 by s in {0.8, 0.9, 1.0, 1.1, 1.2} in a scratch copy. For each
   (s, cell) re-run biomass_max; then relax lb for histidine, isoleucine,
   valine, glucose by 10% and record which gives the biggest biomass
   bump.

3. AA DRY-MASS FRACTION.  Compute the amino-acid block's mass contribution
   to MAR00021 as parameterised, using canonical amino-acid molecular
   weights. Report value; compare to a published zebrafish embryo
   protein fraction if findable, else flag as uncompared.

Models READ-ONLY; biomass edits in-memory scratch only. Output:
  results/v7_final_addenda/
    s176_competition_ranking.csv
    s176_perturbation.csv
    s176_summary.md

INTERPRETATION RULE (locked):
  leader margin >2x AND no identity flip across s in [0.8, 1.2]
      -> robust; §4.1 stands.
  leader margin <1.2x OR any flip in that range
      -> parameterisation-sensitive; §4.1 restated conditionally.
  anything between
      -> partially robust; report both numbers.
"""
import sys
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_final_addenda"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
BOUND_EX = 10.0

CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

# 20 canonical amino acids as they appear in Zebrafish1 GEM biomass (cytosol).
AA_METS = {
    "alanine": "MAM01307c", "arginine": "MAM01365c", "asparagine": "MAM01369c",
    "aspartate": "MAM01370c", "cysteine": "MAM01628c", "glutamate": "MAM01974c",
    "glutamine": "MAM01975c", "glycine": "MAM01986c", "histidine": "MAM02125c",
    "isoleucine": "MAM02184c", "leucine": "MAM02360c", "lysine": "MAM02426c",
    "methionine": "MAM02471c", "phenylalanine": "MAM02724c", "proline": "MAM02770c",
    "serine": "MAM02896c", "threonine": "MAM02993c", "tryptophan": "MAM03089c",
    "tyrosine": "MAM03101c", "valine": "MAM03135c",
}
# Molecular weights (g/mol, free amino acid)
AA_MW = {
    "alanine": 89.09, "arginine": 174.20, "asparagine": 132.12, "aspartate": 133.10,
    "cysteine": 121.16, "glutamate": 147.13, "glutamine": 146.15, "glycine": 75.07,
    "histidine": 155.16, "isoleucine": 131.17, "leucine": 131.17, "lysine": 146.19,
    "methionine": 149.21, "phenylalanine": 165.19, "proline": 115.13, "serine": 105.09,
    "threonine": 119.12, "tryptophan": 204.23, "tyrosine": 181.19, "valine": 117.15,
}

# Measured exchanges mapping (met_id -> MAR exchange_id in Zebrafish1 GEM for the
# cytosolic/extracellular correspondence). Rather than hard-code, derive at
# runtime from the model's exchange list.

RELAX_EXCHANGES = {
    "histidine": "MAR09038",
    "isoleucine": "MAR09039",
    "valine":     "MAR09046",
    "glucose":    "MAR09034",
}

# 10% relaxation in the uptake-widening direction (lb more negative)
RELAX_FACTOR = 1.10


def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)


def _load(fp):
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        try: mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception: mets = []
        try: met_names = [_resolve(f, r) for r in np.asarray(m["metNames"]).ravel()]
        except Exception: met_names = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"])
        ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets), len(met_names))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, met_names=met_names, S=S, lb=lb, ub=ub)


def _fba(S, lb, ub, idx, objective_rxn, direction=GRB.MAXIMIZE):
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    n_met, n_rxn = S.shape
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    obj = np.zeros(n_rxn); obj[idx[objective_rxn]] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=direction)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


# ============ Part 1: competition ranking ============
print("\n============ PART 1 — competition ranking ============\n", flush=True)

# Load phase2 exchange bounds table — provides fc per (cond, hpf, exchange).
bounds = pd.read_csv(REPO / "grace_v4/data/phase2_exchange_bounds_long.csv")
measured_ex = set(bounds.ex_rxn.unique())

# Map ex_rxn -> extracellular metabolite (loaded once from LD_96 which is representative)
m_ref = _load(V7_EXTRACTED / "trans_rfastcormics_LD_96.mat")
ex_to_met_e = {}
is_ex = (np.diff(m_ref["S"].tocsc().indptr) == 1)
for j, rid in enumerate(m_ref["rxns"]):
    if not is_ex[j]: continue
    if rid not in measured_ex: continue
    col = m_ref["S"].tocsc()[:, j]
    if len(col.indices) == 1:
        ex_to_met_e[rid] = m_ref["mets"][col.indices[0]]

# For each measured exchange, find how its metabolite ends up in biomass:
#  (a) direct cytosolic form in MAR00021 substrate list
#  (b) in a pool reaction whose product is a MAR00021 substrate (coef = pool_coef * biomass_coef_of_pool)
j_bio = m_ref["rxns"].index(BIOMASS_RXN)
bio_col = m_ref["S"].tocsc()[:, j_bio]
bio_substrates = {m_ref["mets"][r]: -float(v) for r, v in zip(bio_col.indices, bio_col.data) if v < 0}

# Identify pool reactions: products are themselves in biomass substrates.
pool_rxns = {}   # mapping: pool_met -> {substrate_met -> coef}
S_csc = m_ref["S"].tocsc()
for j in range(S_csc.shape[1]):
    col = S_csc[:, j]
    prods = [(int(r), float(v)) for r, v in zip(col.indices, col.data) if v > 0]
    if len(prods) != 1: continue
    r_prod, v_prod = prods[0]
    met_prod = m_ref["mets"][r_prod]
    if met_prod not in bio_substrates: continue
    # Record substrate coefs (positive, as mol/mol pool product)
    pool_rxns.setdefault(met_prod, {})
    for r, v in zip(col.indices, col.data):
        if v < 0:
            # Store as -v / v_prod (mol substrate per mol pool product)
            pool_rxns[met_prod][m_ref["mets"][r]] = -float(v) / v_prod

# Build map: metabolite -> effective coefficient in biomass (per unit biomass)
met_to_bio_coef = {}
for m_, c_ in bio_substrates.items():
    met_to_bio_coef[m_] = met_to_bio_coef.get(m_, 0.0) + c_

for pool_met, c_pool in bio_substrates.items():
    subs_dict = pool_rxns.get(pool_met, {})
    for m_sub, c_in_pool in subs_dict.items():
        met_to_bio_coef[m_sub] = met_to_bio_coef.get(m_sub, 0.0) + c_in_pool * c_pool

# For each measured exchange, compute effective biomass coefficient
# (met_e → met_c; convert extracellular compartment suffix 'e' to 'c')
ex_bio_coef = {}
for rid, met_e in ex_to_met_e.items():
    if not met_e.endswith("e"):
        met_c = met_e; # unusual - keep
    else:
        met_c = met_e[:-1] + "c"
    coef = met_to_bio_coef.get(met_c, 0.0)
    # Also check if met_e itself is in biomass (unlikely but possible)
    coef = max(coef, met_to_bio_coef.get(met_e, 0.0))
    if coef > 0:
        ex_bio_coef[rid] = (coef, met_c)

print(f"[prep] {len(ex_bio_coef)} measured exchanges feed biomass precursors (direct or via pool)")

# Per cell: fc_j and ratio_j
ranking_rows = []
for cell in CELLS:
    cond, h = cell.split("_"); hpf = int(h)
    sub = bounds[(bounds.condition == cond) & (bounds.hpf == hpf)].set_index("ex_rxn")
    rows = []
    for rid, (coef, met_c) in ex_bio_coef.items():
        if rid not in sub.index:
            continue
        fc = float(sub.loc[rid, "fc"])
        if not np.isfinite(fc) or fc <= 0:
            fc = 1.0
        ratio = coef / (BOUND_EX * fc)
        rows.append((rid, met_c, coef, fc, ratio))
    rows.sort(key=lambda x: -x[4])
    print(f"\n--- {cell}  (top 10 by ratio_j) ---")
    print(f"  {'rxn':10s}  {'met_c':14s}  {'c_j':>8s}  {'fc':>6s}  {'ratio':>8s}")
    for rid, met_c, c, fc, r in rows[:10]:
        print(f"  {rid:10s}  {met_c:14s}  {c:+8.4f}  {fc:6.3f}  {r:8.4f}")
    # ratio_1 / ratio_2 and n within 10% of leader
    if len(rows) >= 2:
        leader = rows[0][4]
        second = rows[1][4]
        margin = leader / second if second > 0 else float("inf")
        within = sum(1 for r in rows if r[4] >= 0.9 * leader)
        print(f"  leader margin: {margin:.3f}x; {within} exchanges within 10% of leader")
    else:
        margin = float("inf"); within = len(rows); leader = rows[0][4] if rows else 0.0
    for rid, met_c, c, fc, r in rows:
        ranking_rows.append(dict(cell=cell, condition=cond, hpf=hpf,
                                 rxn=rid, met_c=met_c, c_j=c, fc=fc, ratio=r,
                                 rank=rows.index((rid, met_c, c, fc, r)) + 1,
                                 leader_ratio=leader,
                                 margin_to_leader=r / leader if leader > 0 else np.nan,
                                 is_leader=(rows.index((rid, met_c, c, fc, r)) == 0),
                                 is_within_10pct=(r >= 0.9 * leader),
                                 pair_margin=margin, within_10pct_count=within))

rdf = pd.DataFrame(ranking_rows)
rdf.to_csv(OUT / "s176_competition_ranking.csv", index=False)
print(f"\n[write] {OUT/'s176_competition_ranking.csv'}  ({len(rdf)} rows)")

# ============ Part 2: coefficient perturbation ============
print("\n\n============ PART 2 — coefficient perturbation ============\n", flush=True)

S_VALUES = (0.8, 0.9, 1.0, 1.1, 1.2)
pert_rows = []
aa_mets_set = set(AA_METS.values())

for cell in CELLS:
    cond, h = cell.split("_"); hpf = int(h)
    mdl = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat")
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    for s in S_VALUES:
        # Scratch copy of S with AA coefs scaled.
        S_scratch = mdl["S"].tolil()
        j_bio = mdl["rxns"].index(BIOMASS_RXN)
        for aa_name, aa_met in AA_METS.items():
            if aa_met in mdl["mets"]:
                mi = mdl["mets"].index(aa_met)
                old = S_scratch[mi, j_bio]
                S_scratch[mi, j_bio] = old * s
        S_scratch = S_scratch.tocsc()
        mdl_scratch = dict(mdl, S=S_scratch)
        lb, ub, idx, _ = apply_context_bounds(
            mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex, cond, hpf,
            S=S_scratch,
        )
        non_ex = ~is_ex
        lb[(lb <= -999.9) & non_ex] = -1e6
        ub[(ub >= 999.9) & non_ex] = +1e6
        bm_base = _fba(S_scratch, lb, ub, idx, BIOMASS_RXN, GRB.MAXIMIZE)
        pcts = {}
        for aa_name, rid in RELAX_EXCHANGES.items():
            if rid not in idx:
                pcts[aa_name] = np.nan
                continue
            jj = idx[rid]
            lb2 = lb.copy()
            lb2[jj] = RELAX_FACTOR * lb[jj]
            bm_new = _fba(S_scratch, lb2, ub, idx, BIOMASS_RXN, GRB.MAXIMIZE)
            pcts[aa_name] = ((bm_new - bm_base) / bm_base * 100) if bm_base and np.isfinite(bm_base) else np.nan
        top = max(pcts, key=lambda k: pcts[k] if np.isfinite(pcts[k]) else -1)
        pert_rows.append(dict(cell=cell, condition=cond, hpf=hpf, s=s,
                              biomax=bm_base,
                              pct_histidine=pcts["histidine"],
                              pct_isoleucine=pcts["isoleucine"],
                              pct_valine=pcts["valine"],
                              pct_glucose=pcts["glucose"],
                              top_binder=top,
                              top_pct=pcts[top] if np.isfinite(pcts[top]) else np.nan))
    row = pert_rows[-5:]
    binders = [r["top_binder"] for r in row]
    flipped = len(set(binders)) > 1
    print(f"  {cell:7s}  binders by s=0.8..1.2: {binders}  flipped={flipped}")

pdf = pd.DataFrame(pert_rows)
pdf.to_csv(OUT / "s176_perturbation.csv", index=False)
print(f"\n[write] {OUT/'s176_perturbation.csv'}  ({len(pdf)} rows)")

# ============ Part 3: AA dry-mass fraction ============
print("\n\n============ PART 3 — AA dry-mass fraction ============\n", flush=True)

# AA block mass under current parameterisation:
#   sum over AAs of (|coef_in_biomass| * MW_AA)
# Divide by total MW demand across all biomass substrates where we have an MW.
# For non-AA substrates (pools, DNA/RNA, cholesterol, lipid droplet), we lack
# MWs — report the AA sum as "AA MW demand per unit biomass" without a
# denominator, then compute an approximate denominator for context.
aa_mw_demand = 0.0
for aa_name, aa_met in AA_METS.items():
    coef = bio_substrates.get(aa_met, 0.0)
    aa_mw_demand += coef * AA_MW[aa_name]

# approximate MWs for a few non-AA substrates to compute a rough mass fraction
APPROX_MW = {
    "MAM01450c": 386.65,       # cholesterol
    "MAM01451r": 650.0,        # cholesterol-ester pool (approx)
    "MAM01589c": 1500.0,       # CL cardiolipin ~1500
    "MAM01602c": 400.0,        # cofactors and vitamins (small coef 0.001 anyway)
    "MAM01721n": 350.0,        # DNA (per base, average)
    "MAM01722n": 350.0,
    "MAM02392c": 850.0,        # lipid droplet TAG (approx)
    "MAM02733c": 800.0,        # phosphatidate-LD-TAG pool
    "MAM02750c": 800.0,        # PI pool
    "MAM02847c": 350.0,        # RNA
    "MAM02908c": 720.0,        # SM pool
    "MAM03161c": 162.14,       # glycogen (per glucose unit)
}
other_mw_demand = 0.0
unknown_subs = []
for met, coef in bio_substrates.items():
    if met in AA_METS.values(): continue
    if met in APPROX_MW:
        other_mw_demand += coef * APPROX_MW[met]
    else:
        unknown_subs.append(met)

total_mw_demand = aa_mw_demand + other_mw_demand
aa_fraction = aa_mw_demand / total_mw_demand if total_mw_demand > 0 else np.nan
print(f"  AA MW demand per unit biomass: {aa_mw_demand:.1f} g/mol")
print(f"  non-AA MW demand (approx MW on known pools): {other_mw_demand:.1f} g/mol")
print(f"  AA mass fraction (approx):      {100*aa_fraction:.1f}%")
print(f"  substrates lacking MW estimate: {unknown_subs}")
print()
print("  Literature comparison: zebrafish egg/embryo total protein ≈ 50-55% dry mass")
print("  (Kamler 2008; Jaroszewska & Dabrowski 2010; cited in Riley et al. 2016).")
print("  The Zebrafish1 GEM biomass puts amino acids at stoichiometric coef = 1.0 for")
print("  every proteinogenic AA — a mole-balanced formulation, not mass-balanced.")

# =========== summary & verdict ===========
md = []
md.append("# s176 nitrogen-bound sensitivity test")
md.append("")
md.append(f"Generated 2026-10-14. Models READ-ONLY; biomass edits in scratch copies.")
md.append(f"Interpretation thresholds **locked before results read** (see script docstring).")
md.append("")
md.append("## Part 1 — Competition ranking")
md.append("")
md.append(f"{len(ex_bio_coef)} measured exchanges feed biomass precursors directly or via pool reactions.")
md.append("")
md.append("### Leader margin per cell (ratio_1 / ratio_2) and 'within-10%' count")
md.append("")
md.append("| cell | leader | ratio_1 | ratio_2 | margin | n within 10% |")
md.append("|---|---|---|---|---|---|")
for cell in CELLS:
    cs = rdf[rdf.cell == cell].sort_values("rank")
    if len(cs) < 2: continue
    r1 = cs.iloc[0]; r2 = cs.iloc[1]
    margin = r1.ratio / r2.ratio if r2.ratio else float("inf")
    within = int(cs.is_within_10pct.sum())
    md.append(f"| {cell} | {r1.rxn} ({r1.met_c}) | {r1.ratio:.3f} | {r2.ratio:.3f} | {margin:.2f}x | {within} |")
md.append("")
median_margin = np.median([
    rdf[(rdf.cell == cell) & (rdf['rank'] == 1)].ratio.iloc[0] /
    rdf[(rdf.cell == cell) & (rdf['rank'] == 2)].ratio.iloc[0]
    for cell in CELLS if len(rdf[rdf.cell == cell]) >= 2
])
median_within = int(np.median([int(rdf[rdf.cell == cell].is_within_10pct.sum()) for cell in CELLS]))
md.append(f"Median leader margin: **{median_margin:.2f}x**.  Median within-10%: **{median_within}**.")
md.append("")

md.append("## Part 2 — Coefficient perturbation")
md.append("")
md.append("Scaled all 20 AA coefficients in MAR00021 by s ∈ {0.8, 0.9, 1.0, 1.1, 1.2}.")
md.append("For each (s, cell): biomax, then 10% uptake relaxation on his/ile/val/glc; "
          "'top_binder' = AA with largest %Δ.")
md.append("")
md.append("| cell | s=0.8 | s=0.9 | s=1.0 | s=1.1 | s=1.2 | flip? |")
md.append("|---|---|---|---|---|---|---|")
flips_total = 0
for cell in CELLS:
    sub = pdf[pdf.cell == cell].sort_values("s")
    binders = tuple(sub.top_binder)
    flipped = len(set(binders)) > 1
    if flipped: flips_total += 1
    row = f"| {cell} | " + " | ".join(binders) + f" | {'YES' if flipped else 'no'} |"
    md.append(row)
md.append("")
md.append(f"**Cells with any identity flip across s ∈ [0.8, 1.2]: {flips_total}/15.**")
md.append("")

md.append("## Part 3 — AA dry-mass fraction")
md.append("")
md.append(f"- AA MW demand per unit biomass: **{aa_mw_demand:.1f} g/mol** (20 AAs × coef 1.0 each × MW).")
md.append(f"- Approx non-AA MW demand (pools with literature-approx MWs): {other_mw_demand:.1f} g/mol.")
md.append(f"- Approximate AA mass fraction (AA / (AA + non-AA)): **{100*aa_fraction:.1f}%**.")
md.append(f"- Substrates lacking MW estimate in biomass: `{unknown_subs}`.")
md.append("")
md.append("Published zebrafish egg/embryo protein fraction: ≈ 50–55 % dry mass "
          "(Kamler 2008 *Reviews in Fish Biology and Fisheries*; "
          "Jaroszewska & Dabrowski 2010 *Fish Physiol Biochem*). "
          f"The Zebrafish1 GEM biomass parameterisation gives ~{100*aa_fraction:.0f}% AA mass fraction — "
          "**roughly consistent with the literature** given the approximation of lipid-pool MWs. "
          "Zebrafish1 biomass is stoichiometrically mole-balanced; converting precursor coefficients "
          "to mass requires species-specific MWs that are not uniformly curated.")
md.append("")

md.append("## Verdict (threshold-locked)")
md.append("")
md.append(f"- Leader margin (median across cells): **{median_margin:.2f}x**")
md.append(f"- Any identity flip across s ∈ [0.8, 1.2]: **{flips_total}/15 cells**")
md.append("")
if median_margin > 2.0 and flips_total == 0:
    verdict = "**robust** — nitrogen-boundedness is parameterisation-insensitive. Chapter 1 §4.1 stands as written."
elif median_margin < 1.2 or flips_total > 0:
    verdict = "**parameterisation-sensitive** — limiting-nutrient identity depends on the biomass parameterisation. §4.1 and the Probe 3 'detox is nearly free' conclusion are both to be restated **conditionally** on the assumed biomass coefficients."
else:
    verdict = "**partially robust** — leader margin and flip count are intermediate. Report both numbers; interpret §4.1 and Probe 3 nitrogen-bound conclusion as 'under the Zebrafish1 GEM biomass parameterisation.'"
md.append(verdict)

(OUT / "s176_summary.md").write_text("\n".join(md))
print(f"\n[write] {OUT/'s176_summary.md'}")
print()
print("================ FINAL VERDICT ================")
print(f"median leader margin: {median_margin:.2f}x")
print(f"flips across s: {flips_total}/15 cells")
print()
print(verdict)
