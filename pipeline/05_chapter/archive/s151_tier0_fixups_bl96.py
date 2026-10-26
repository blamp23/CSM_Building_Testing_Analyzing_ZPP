"""
s151 -- Tier 0 fixups + BL_96 drill-down.

Rewrites:
  tier0_2x2_v7.csv            4 corners + top/env/interaction/gap in absolute
                              biomass units; dominance per pair; shares only
                              where |gap| > 5% of smaller corner.
  tier0_parent_child_v7.csv   Replaces single top binding with the top-3
                              exchange constraints by dual × bound per cell
                              (lb=0 irreversibility rows excluded).

New files:
  bl96_lacking_reactions.csv       54 transporters + 43 exchanges in LD_96 but
                                   not BL_96, with metabolite names + compartments.
  bl96_greedy_addback.csv          Ordered greedy add-back from 8.535 until ≥ 95%
                                   of the 12.255 ceiling. Each row = one pick.
  bl96_greedy_gene_zfpkm.csv       For each picked reaction's genes, zFPKM in
                                   BL_96 vs BL_96's exp/inex thresholds.
"""
import sys, re
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
V7_BASE_DIR  = REPO / "ocr_anchored_extraction/models"
FC_TABLE     = REPO / "grace_v4/data/phase2_exchange_bounds_long.csv"
ZFPKM_FP     = REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv"
ZFPKM_TH_FP  = REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv"
OUT          = REPO / "grace_v4/results/tier0_v7"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
CONDS = ("BL","D","LD")
HPFS  = (24,48,72,96,120)
F1_PAIRS = [(f"{a}_{h}", f"{b}_{h}") for h in HPFS
            for a,b in [("BL","D"),("BL","LD"),("D","LD")]]

MB = pd.read_csv(FC_TABLE)
MEASURED_SET = set(MB.ex_rxn.unique())


# ---------- MATLAB v7.3 helpers --------------------------------------------
def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)


def _str_list(f, m, key):
    try:
        arr = np.asarray(m[key]).ravel()
        return [_resolve(f, r) for r in arr]
    except Exception:
        return []


def _load_full(fp, root_key):
    with h5py.File(fp, "r") as f:
        m = f[root_key]
        rxns = _str_list(f, m, "rxns")
        try: subs = _str_list(f, m, "subSystems")
        except Exception: subs = [""] * len(rxns)
        try: metNames = _str_list(f, m, "metNames")
        except Exception: metNames = []
        try: mets = _str_list(f, m, "mets")
        except Exception: mets = []
        try: metComps = np.asarray(m["metComps"]).ravel().astype(int)
        except Exception: metComps = np.zeros(len(mets), dtype=int)
        try: compNames = _str_list(f, m, "compNames")
        except Exception: compNames = []
        try: comps = _str_list(f, m, "comps")
        except Exception: comps = []
        try: grRules = _str_list(f, m, "grRules")
        except Exception: grRules = [""] * len(rxns)
        try: genes = _str_list(f, m, "genes")
        except Exception: genes = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        try: n_met_hdf = int(np.asarray(m["mets"]).shape[0])
        except Exception: n_met_hdf = 0
        n_met = max(int(ir.max()) + 1, n_met_hdf, len(mets), len(metNames))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc)-1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, S=S, lb=lb, ub=ub,
                subs=subs, mets=mets, metNames=metNames, metComps=metComps,
                comps=comps, compNames=compNames, grRules=grRules, genes=genes)


# ---------- FBA via Gurobi + dual decomposition ----------------------------
def fba(model, lb, ub, idx, with_duals=False):
    S = model["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("fba", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c_obj = np.zeros(n_rxn); c_obj[idx[BIOMASS_RXN]] = 1.0
    gm.setMObjective(None, c_obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    if gm.Status != GRB.OPTIMAL:
        gm.dispose(); env.dispose()
        return (np.nan, np.zeros(n_rxn)) if with_duals else np.nan
    bm = gm.ObjVal
    if with_duals:
        rcs = np.array([v[j].RC for j in range(n_rxn)])
        gm.dispose(); env.dispose()
        return bm, rcs
    gm.dispose(); env.dispose()
    return bm


def top3_exchange_duals(model, lb, ub, idx, is_ex, n_top=3):
    """Top n exchange constraints by |rc × bound|, excluding lb = 0 (irreversible)."""
    bm, rcs = fba(model, lb, ub, idx, with_duals=True)
    hits = []
    for j in range(len(model["rxns"])):
        if not is_ex[j]: continue
        if model["rxns"][j] == BIOMASS_RXN: continue
        rc = rcs[j]
        if abs(rc) < 1e-9: continue
        if rc > 0:  # lower-bound binding
            if lb[j] == 0.0: continue   # exclude irreversibility
            if abs(lb[j]) < 1e-12: continue
            hits.append((abs(rc * lb[j]), model["rxns"][j], "lb", lb[j], rc))
        else:  # upper-bound binding
            if abs(ub[j]) < 1e-12: continue
            hits.append((abs(rc * ub[j]), model["rxns"][j], "ub", ub[j], rc))
    hits.sort(reverse=True)
    return bm, hits[:n_top]


# ===========================================================================
# PART 1 — rewrite tier0_parent_child_v7.csv with top-3 exchange duals
# ===========================================================================
print("=== PART 1 — top-3 exchange duals per cell ===", flush=True)
pc_rows = []
base_cache = None; last_hpf = -1
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]
for cell in CELLS:
    cond, hpf = cell.split("_"); hpf = int(hpf)
    if hpf != last_hpf:
        base_cache = _load_full(V7_BASE_DIR / f"baked_ocr_{hpf}_v7.mat", "m")
        last_hpf = hpf
    is_ex_p = (np.diff(base_cache["S"].tocsc().indptr) == 1)
    lbp, ubp, idxp, _ = apply_context_bounds(
        base_cache["lb"].copy(), base_cache["ub"].copy(),
        base_cache["rxns"], is_ex_p, cond, hpf)
    p_bm, p_hits = top3_exchange_duals(base_cache, lbp, ubp, idxp, is_ex_p)

    ch = _load_full(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    is_ex_c = (np.diff(ch["S"].tocsc().indptr) == 1)
    lbc, ubc, idxc, _ = apply_context_bounds(
        ch["lb"].copy(), ch["ub"].copy(),
        ch["rxns"], is_ex_c, cond, hpf)
    c_bm, c_hits = top3_exchange_duals(ch, lbc, ubc, idxc, is_ex_c)

    rec = dict(cell=cell, condition=cond, hpf=hpf,
               parent_biomax=p_bm, child_biomax=c_bm,
               parent_over_child=p_bm/c_bm if c_bm else np.nan)
    for k in range(3):
        if k < len(p_hits):
            c0,r0,s0,b0,rc0 = p_hits[k]
            rec[f"parent_top{k+1}_rxn"]=r0; rec[f"parent_top{k+1}_side"]=s0
            rec[f"parent_top{k+1}_bound"]=b0; rec[f"parent_top{k+1}_dual"]=rc0
            rec[f"parent_top{k+1}_contrib"]=c0
        else:
            for suf in ["rxn","side","bound","dual","contrib"]: rec[f"parent_top{k+1}_{suf}"]=""
        if k < len(c_hits):
            c0,r0,s0,b0,rc0 = c_hits[k]
            rec[f"child_top{k+1}_rxn"]=r0; rec[f"child_top{k+1}_side"]=s0
            rec[f"child_top{k+1}_bound"]=b0; rec[f"child_top{k+1}_dual"]=rc0
            rec[f"child_top{k+1}_contrib"]=c0
        else:
            for suf in ["rxn","side","bound","dual","contrib"]: rec[f"child_top{k+1}_{suf}"]=""
    pc_rows.append(rec)
    print(f"  {cell}: p={p_bm:.3f} top3={[h[1] for h in p_hits]}  "
          f"c={c_bm:.3f} top3={[h[1] for h in c_hits]}", flush=True)

pc_df = pd.DataFrame(pc_rows)
pc_df.to_csv(OUT / "tier0_parent_child_v7.csv", index=False)
print(f"[write] tier0_parent_child_v7.csv rows={len(pc_df)}", flush=True)


# ===========================================================================
# PART 2 — 2x2 in absolute biomass units, per-pair dominance
# ===========================================================================
print("\n=== PART 2 — 2×2 in absolute units, per-pair dominance ===", flush=True)
child_cache = {}
def get_child(cell):
    if cell not in child_cache:
        child_cache[cell] = _load_full(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat","parent")
    return child_cache[cell]

tt_rows = []
for cA, cB in F1_PAIRS:
    condA, hpf = cA.split("_"); hpf = int(hpf)
    condB = cB.split("_")[0]
    mA = get_child(cA); mB = get_child(cB)
    is_ex_A = (np.diff(mA["S"].tocsc().indptr) == 1)
    is_ex_B = (np.diff(mB["S"].tocsc().indptr) == 1)
    def corner(mX, is_ex_X, cond_env):
        lb_, ub_, idx_, _ = apply_context_bounds(
            mX["lb"].copy(), mX["ub"].copy(), mX["rxns"], is_ex_X, cond_env, hpf)
        return fba(mX, lb_, ub_, idx_)
    aa = corner(mA, is_ex_A, condA); ab = corner(mA, is_ex_A, condB)
    ba = corner(mB, is_ex_B, condA); bb = corner(mB, is_ex_B, condB)
    total_gap = aa - bb
    top_eff = ((aa-ba) + (ab-bb)) / 2
    env_eff = ((aa-ab) + (ba-bb)) / 2
    interaction = (aa - ab) - (ba - bb)
    small = min(aa, bb)
    material = abs(total_gap) > 0.05 * small
    if material:
        top_share  = top_eff / total_gap * 100 if abs(total_gap)>1e-9 else np.nan
        env_share  = env_eff / total_gap * 100 if abs(total_gap)>1e-9 else np.nan
        int_share  = interaction / total_gap * 100 if abs(total_gap)>1e-9 else np.nan
        dominance  = "topology" if abs(top_eff) >= abs(env_eff) else "envelope"
        verdict    = f"{dominance}-dominated"
    else:
        top_share=env_share=int_share=np.nan
        verdict = "no material gap"
    rxns_A = set(mA["rxns"]); rxns_B = set(mB["rxns"])
    unswap = (MEASURED_SET & (rxns_A ^ rxns_B))
    tt_rows.append(dict(
        hpf=hpf, pair=f"{cA}_vs_{cB}", cell_A=cA, cell_B=cB,
        bmax_A_top_A_env=aa, bmax_A_top_B_env=ab,
        bmax_B_top_A_env=ba, bmax_B_top_B_env=bb,
        total_gap_diag=total_gap,
        topology_effect=top_eff, envelope_effect=env_eff, interaction=interaction,
        topology_share_pct=top_share, envelope_share_pct=env_share,
        interaction_share_pct=int_share,
        gap_pct_of_smaller_corner=abs(total_gap)/small*100,
        dominance=verdict,
        n_unswappable_exchanges=len(unswap),
        unswappable_exchanges=";".join(sorted(unswap)),
    ))
    print(f"  {cA} vs {cB}  AA={aa:.3f} AB={ab:.3f} BA={ba:.3f} BB={bb:.3f}  "
          f"gap={total_gap:+.3f}  top={top_eff:+.3f}  env={env_eff:+.3f}  int={interaction:+.3f}  "
          f"→ {verdict}", flush=True)

tt_df = pd.DataFrame(tt_rows)
tt_df.to_csv(OUT / "tier0_2x2_v7.csv", index=False)
print(f"[write] tier0_2x2_v7.csv rows={len(tt_df)}", flush=True)


# ===========================================================================
# PART 3 — BL_96 drill-down
# ===========================================================================
print("\n=== PART 3 — BL_96 drill-down ===", flush=True)
bl96 = get_child("BL_96")
ld96 = get_child("LD_96")
bl_rxns = set(bl96["rxns"])
ld_rxns_arr = np.asarray(ld96["rxns"])

# ---------- 3a: list the 54 transporters + 43 exchanges ----------
# metComps -> comp code mapping
ld_metComps = np.asarray(ld96["metComps"]) if len(ld96["metComps"]) else None
ld_compNames = ld96["compNames"] if len(ld96["compNames"]) else ["" for _ in range(99)]
def comp_of_met(mi):
    if ld_metComps is None or mi >= len(ld_metComps): return ""
    c = ld_metComps[mi]
    try:
        return ld_compNames[c-1] if c-1 < len(ld_compNames) else str(c)
    except Exception:
        return str(c)

ld_mets = ld96["mets"]; ld_metNames = ld96["metNames"]
ld_S_csc = ld96["S"].tocsc()
def rxn_metinfo(j):
    col = ld_S_csc[:, j]
    rows = col.indices; vals = col.data
    subs = []; prods = []
    for r_, v_ in zip(rows, vals):
        name = ld_metNames[r_] if r_ < len(ld_metNames) else ld_mets[r_]
        cmp  = comp_of_met(r_)
        tag = f"{name} [{cmp}]"
        (subs if v_ < 0 else prods).append(tag)
    return "; ".join(subs), "; ".join(prods)

missing_rxns = []
for j, r in enumerate(ld_rxns_arr):
    if r in bl_rxns: continue
    sub, prod = rxn_metinfo(j)
    missing_rxns.append(dict(
        rxn=r,
        subsystem=ld96["subs"][j] if j < len(ld96["subs"]) else "",
        substrates=sub, products=prod,
        gpr=ld96["grRules"][j] if j < len(ld96["grRules"]) else "",
        is_exchange=(np.diff(ld_S_csc.indptr)[j] == 1),
    ))
miss_df = pd.DataFrame(missing_rxns)
# flag the 54 transporters and 43 exchanges from the Part 2 subsystem buckets
is_T = miss_df.subsystem == "Transport reactions"
is_E = miss_df.subsystem == "Exchange/demand reactions"
print(f"  Transport reactions missing: {int(is_T.sum())}", flush=True)
print(f"  Exchange/demand reactions missing: {int(is_E.sum())}", flush=True)
miss_df[is_T | is_E].to_csv(OUT / "bl96_lacking_reactions.csv", index=False)
print(f"[write] bl96_lacking_reactions.csv rows={int((is_T|is_E).sum())}", flush=True)


# ---------- 3b: greedy add-back from BL_96 baseline toward 12.255 ceiling ----------
print("\n--- greedy add-back ---", flush=True)
base96 = _load_full(V7_BASE_DIR / "baked_ocr_96_v7.mat", "m")
base_rxns_arr = np.asarray(base96["rxns"])
base_idx_map = {r: i for i, r in enumerate(base_rxns_arr)}
is_ex_base = (np.diff(base96["S"].tocsc().indptr) == 1)
# Pre-apply BL,96 context bounds on the base
lb_b, ub_b, idx_b, _ = apply_context_bounds(
    base96["lb"].copy(), base96["ub"].copy(),
    base96["rxns"], is_ex_base, "BL", 96)

# Initial mask: BL_96's reactions
def set_mask(mask_set, lb, ub):
    for i, r in enumerate(base_rxns_arr):
        if r not in mask_set:
            lb[i] = 0.0; ub[i] = 0.0
    return lb, ub

allowed = set(bl_rxns)
lb_cur = lb_b.copy(); ub_cur = ub_b.copy()
lb_cur, ub_cur = set_mask(allowed, lb_cur, ub_cur)
baseline = fba(base96, lb_cur, ub_cur, idx_b)
print(f"  BL_96 baseline (base-sized LP, BL env): {baseline:.4f}", flush=True)

# Ceiling: BL ∪ all 173 missing
all_missing_rxns = set(miss_df.rxn)
lb_ceil = lb_b.copy(); ub_ceil = ub_b.copy()
lb_ceil, ub_ceil = set_mask(bl_rxns.union(all_missing_rxns), lb_ceil, ub_ceil)
ceiling = fba(base96, lb_ceil, ub_ceil, idx_b)
print(f"  ceiling (BL ∪ all 173 extras): {ceiling:.4f}", flush=True)
target = baseline + 0.95 * (ceiling - baseline)
print(f"  target (95% of ceiling): {target:.4f}", flush=True)

# Build persistent Gurobi model: all 173 extras start at lb=ub=0; toggle to native values
# For efficiency, keep a single model and use setAttr to flip bounds per iteration.
# base96 model has n_rxn ~ 12909; LP solves in <1s each after warm start.

candidate_pool = sorted(all_missing_rxns)
cand_idx = [base_idx_map[r] for r in candidate_pool]
# native bounds for candidates (after apply_context_bounds)
native_lb = {r: lb_b[base_idx_map[r]] for r in candidate_pool}
native_ub = {r: ub_b[base_idx_map[r]] for r in candidate_pool}

env_g = gp.Env(empty=True); env_g.setParam("OutputFlag",0); env_g.start()
gm = gp.Model("addback", env=env_g)
S_csr = base96["S"].tocsr(); n_met, n_rxn = S_csr.shape
v = gm.addMVar(n_rxn, lb=lb_cur, ub=ub_cur)
gm.addMConstr(S_csr, v, "=", np.zeros(n_met))
c_obj = np.zeros(n_rxn); c_obj[idx_b[BIOMASS_RXN]] = 1.0
gm.setMObjective(None, c_obj, 0.0, sense=GRB.MAXIMIZE)
gm.setParam("Method", 1)   # dual simplex (good for warm restart under bound changes)
gm.optimize()
cur_bm = gm.ObjVal
print(f"  (gurobi persistent) baseline={cur_bm:.4f}", flush=True)

picks = []
remaining = set(candidate_pool)
step = 0
while cur_bm < target and remaining:
    step += 1
    best_rxn = None; best_bm = cur_bm; best_lb_new = 0.0; best_ub_new = 0.0
    # evaluate every remaining candidate
    for r in remaining:
        j = base_idx_map[r]
        old_lb, old_ub = v[j].lb, v[j].ub
        v[j].lb = native_lb[r]; v[j].ub = native_ub[r]
        gm.optimize()
        if gm.Status == GRB.OPTIMAL and gm.ObjVal > best_bm + 1e-6:
            best_bm = gm.ObjVal; best_rxn = r
            best_lb_new = native_lb[r]; best_ub_new = native_ub[r]
        v[j].lb = old_lb; v[j].ub = old_ub
    if best_rxn is None:
        print(f"  step {step}: no further gain found; stopping.", flush=True)
        break
    j = base_idx_map[best_rxn]
    v[j].lb = native_lb[best_rxn]; v[j].ub = native_ub[best_rxn]
    gm.optimize(); cur_bm = gm.ObjVal
    sub = miss_df.loc[miss_df.rxn == best_rxn, "subsystem"].iloc[0]
    picks.append(dict(step=step, rxn=best_rxn, subsystem=sub,
                      biomax_after=cur_bm,
                      gain_this_step=cur_bm - (picks[-1]["biomax_after"] if picks else baseline),
                      cum_gain=cur_bm - baseline,
                      frac_of_full_recovery=(cur_bm - baseline) / max(ceiling - baseline, 1e-9)))
    remaining.discard(best_rxn)
    print(f"  step {step}: +{best_rxn} ({sub[:40]})  bm→{cur_bm:.4f}  frac={picks[-1]['frac_of_full_recovery']:.3f}",
          flush=True)
    if step > 60:
        print("  safety break at 60 steps", flush=True); break

greedy_df = pd.DataFrame(picks)
greedy_df.to_csv(OUT / "bl96_greedy_addback.csv", index=False)
print(f"[write] bl96_greedy_addback.csv rows={len(greedy_df)}", flush=True)
gm.dispose(); env_g.dispose()


# ---------- 3c: zFPKM check on the picked reactions' genes ----------
print("\n--- zFPKM check on picked genes ---", flush=True)
zf = pd.read_csv(ZFPKM_FP)
zth = pd.read_csv(ZFPKM_TH_FP).set_index("cell")
bl96_zf = zf[["ensembl_gene_id","BL_96"]].rename(columns={"BL_96":"zfpkm_BL_96"})
bl96_zf = bl96_zf.set_index("ensembl_gene_id")
bl96_exp_thr   = float(zth.loc["BL_96","exp_threshold_zFPKM"])
bl96_inex_thr  = float(zth.loc["BL_96","inex_threshold_zFPKM"])
print(f"  BL_96 thresholds: exp={bl96_exp_thr:.3f}, inex={bl96_inex_thr:.3f}", flush=True)

def parse_genes(gpr):
    return re.findall(r"ENSDARG\d+", gpr or "")

zf_rows = []
for _, r in greedy_df.iterrows():
    gpr = miss_df.loc[miss_df.rxn == r.rxn, "gpr"].iloc[0]
    genes = parse_genes(gpr)
    verdicts = []
    for g in genes:
        if g in bl96_zf.index:
            v_ = bl96_zf.loc[g, "zfpkm_BL_96"]
            if pd.isna(v_): cls = "NA"
            elif v_ >= bl96_exp_thr:  cls = "expressed"
            elif v_ < bl96_inex_thr:  cls = "not_expressed"
            else:                      cls = "unknown"
            verdicts.append(f"{g}:{v_ if pd.isna(v_) else f'{v_:.2f}'}:{cls}")
        else:
            verdicts.append(f"{g}:NA:absent_from_zfpkm")
    zf_rows.append(dict(step=r.step, rxn=r.rxn, subsystem=r.subsystem,
                        gpr=gpr, n_genes=len(genes),
                        gene_zfpkm_class="; ".join(verdicts),
                        any_expressed=any("expressed" in v and "not_expressed" not in v and "NA_expressed" not in v for v in verdicts),
                        all_not_expressed=(len(verdicts)>0 and all(":not_expressed" in v for v in verdicts))))
zf_df = pd.DataFrame(zf_rows)
zf_df.to_csv(OUT / "bl96_greedy_gene_zfpkm.csv", index=False)
print(f"[write] bl96_greedy_gene_zfpkm.csv rows={len(zf_df)}", flush=True)

n_expr   = int(zf_df.any_expressed.sum())
n_notexp = int(zf_df.all_not_expressed.sum())
n_mixed  = len(zf_df) - n_expr - n_notexp
print(f"  genes expressed above threshold (any): {n_expr}", flush=True)
print(f"  all genes below not-expressed threshold: {n_notexp}", flush=True)
print(f"  mixed / unknown / no GPR:               {n_mixed}", flush=True)
print(f"  verdict: biology={n_notexp} / FASTCORE-pruning={n_expr}", flush=True)

print("\n=== DONE ===", flush=True)
