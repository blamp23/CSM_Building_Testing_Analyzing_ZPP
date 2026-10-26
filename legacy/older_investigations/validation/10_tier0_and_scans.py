"""
s146 -- Tier 0 constraint decomposition + ε scan.

Runs on the v7 base (baked_ocr_{hpf}_v7.mat) and, if present, the v7
extracted models (v7_models/). Uses apply_context_bounds (flux column).

(A) Tier 0 constraint decomposition
    For each cell (v7 base + 15 extracted):
      FBA under apply_context_bounds.
      For each binding exchange constraint: contribution = shadow × |bound|.
      Report top-5 exchanges by contribution and their share of biomass_max
      (they sum to biomass_max in the LP dual sense).

(B) ε scan {1e-2, 1e-3, 1e-4}
    In-memory edit of MAR00021's MAM01602c coefficient. FBA. Report a
    15-cell × 3-ε table of biomass_max under flux bounds.
"""
import sys, h5py, numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.optimize import linprog
from pathlib import Path

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds

V7_MODELS = REPO / "reviewer_packet_v4_slim/v7_models"
BASE_DIR  = REPO / "ocr_anchored_extraction/models"
BIOMASS_RXN = "MAR00021"
POOL_PSEUDO_MET = "MAM01602c"
CELLS = [f"{c}_{h}" for c in ("BL","D","LD") for h in (24,48,72,96,120)]
EPS_SCAN = [1e-2, 1e-3, 1e-4]
BOUND_EX_SCAN = [3.0, 10.0, 30.0]


def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)


def load_model(fp, root_key):
    with h5py.File(fp, "r") as f:
        m = f[root_key]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        try:
            metNames = [_resolve(f, r) for r in np.asarray(m["metNames"]).ravel()]
        except Exception: metNames = [""] * len(mets)
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        S = sp.csc_matrix((data, ir, jc), shape=(int(ir.max())+1, len(jc)-1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, metNames=metNames, S=S, lb=lb, ub=ub)


def biomax_with_res(m, lb, ub, idx):
    n_rxn = len(m["rxns"]); n_met = m["S"].shape[0]
    c = np.zeros(n_rxn); c[idx[BIOMASS_RXN]] = -1.0
    A_eq = m["S"].toarray(); b_eq = np.zeros(n_met)
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=list(zip(lb, ub)), method="highs")
    return ((-res.fun) if res.success else np.nan), res


def constraint_decomposition(res, m, is_ex, lb, ub, n_top=5):
    """Return top-n exchange constraints by contribution to biomass_max
    (contribution = shadow × |bound|). Sum of contributions ≈ biomass_max
    by LP duality (only exchanges with nonzero shadow contribute)."""
    if res is None: return [], np.nan
    lb_m = res.lower.marginals; ub_m = res.upper.marginals
    hits = []
    total = 0.0
    for j in range(len(m["rxns"])):
        if not is_ex[j]: continue
        if m["rxns"][j] == BIOMASS_RXN: continue
        s_lb = abs(lb_m[j]); s_ub = abs(ub_m[j])
        if s_lb > 1e-9 and abs(lb[j]) > 1e-12:
            contrib = s_lb * abs(lb[j])
            hits.append((contrib, m["rxns"][j], "lb", lb[j], s_lb, res.x[j]))
            total += contrib
        if s_ub > 1e-9 and abs(ub[j]) > 1e-12:
            contrib = s_ub * abs(ub[j])
            hits.append((contrib, m["rxns"][j], "ub", ub[j], s_ub, res.x[j]))
            total += contrib
    hits.sort(reverse=True)
    return hits[:n_top], total


def _set_pool_coef(m, eps):
    """Return (lb_bak, ub_bak, S_val_bak) to restore afterwards."""
    if POOL_PSEUDO_MET not in m["mets"]: return None
    if BIOMASS_RXN not in m["rxns"]: return None
    i = m["mets"].index(POOL_PSEUDO_MET)
    j = m["rxns"].index(BIOMASS_RXN)
    old = m["S"][i, j]
    S_new = m["S"].tolil()
    S_new[i, j] = -eps
    m["S"] = S_new.tocsc()
    return old


def _restore_pool_coef(m, old):
    if old is None: return
    i = m["mets"].index(POOL_PSEUDO_MET); j = m["rxns"].index(BIOMASS_RXN)
    S_new = m["S"].tolil(); S_new[i, j] = old; m["S"] = S_new.tocsc()


# ==== (A) Tier 0 constraint decomposition ================================
print("=" * 76)
print("(A) TIER 0 constraint decomposition (flux bounds, ε=1e-3)")
print("=" * 76)
rows = []; last_hpf = -1; base_cache = None
for cell in CELLS:
    cond, hpf = cell.split("_"); hpf = int(hpf)
    if hpf != last_hpf:
        base_cache = load_model(BASE_DIR / f"baked_ocr_{hpf}_v7.mat", "m")
        last_hpf = hpf

    # v7 base
    lb0 = base_cache["lb"].copy(); ub0 = base_cache["ub"].copy()
    Snnz = np.diff(base_cache["S"].tocsc().indptr); is_ex_b = Snnz == 1
    lb, ub, idx, _ = apply_context_bounds(lb0, ub0, base_cache["rxns"], is_ex_b, cond, hpf)
    p_bm, p_res = biomax_with_res(base_cache, lb, ub, idx)
    p_hits, p_sum = constraint_decomposition(p_res, base_cache, is_ex_b, lb, ub, 5)

    # v7 extracted (if exists)
    child_fp = V7_MODELS / f"trans_rfastcormics_{cell}.mat"
    if child_fp.exists():
        child = load_model(child_fp, "parent")
        lb0c = child["lb"].copy(); ub0c = child["ub"].copy()
        Snnz_c = np.diff(child["S"].tocsc().indptr); is_ex_c = Snnz_c == 1
        lb2, ub2, idx2, _ = apply_context_bounds(lb0c, ub0c, child["rxns"], is_ex_c, cond, hpf)
        c_bm, c_res = biomax_with_res(child, lb2, ub2, idx2)
        c_hits, c_sum = constraint_decomposition(c_res, child, is_ex_c, lb2, ub2, 5)
    else:
        c_bm = np.nan; c_hits = []; c_sum = np.nan

    print(f"\n{cell:6s}  parent={p_bm:.4e}  child={c_bm:.4e}  (Δ={c_bm-p_bm:+.2e})")
    print(f"  parent top-5 (sum={p_sum:.4e}):")
    for contrib, rxn, side, bnd, sh, v in p_hits:
        met_name = ""
        col = base_cache["S"].getcol(base_cache["rxns"].index(rxn)).toarray().ravel()
        nz = np.where(col != 0)[0]
        if len(nz):
            mm = base_cache["mets"][nz[0]]
            met_name = str(base_cache["metNames"][nz[0]])[:24]
        share = 100 * contrib / p_bm if p_bm > 0 else np.nan
        print(f"    {rxn:10s} {side}={bnd:+.3e} shadow={sh:.3e} contrib={contrib:.3e} ({share:5.1f}%)  {met_name}")

    if child_fp.exists():
        print(f"  child top-5 (sum={c_sum:.4e}):")
        for contrib, rxn, side, bnd, sh, v in c_hits:
            met_name = ""
            col = child["S"].getcol(child["rxns"].index(rxn)).toarray().ravel()
            nz = np.where(col != 0)[0]
            if len(nz):
                mm = child["mets"][nz[0]]
                met_name = str(child["metNames"][nz[0]])[:24]
            share = 100 * contrib / c_bm if c_bm > 0 else np.nan
            print(f"    {rxn:10s} {side}={bnd:+.3e} shadow={sh:.3e} contrib={contrib:.3e} ({share:5.1f}%)  {met_name}")

    rows.append(dict(cell=cell, parent_biomax=p_bm, child_biomax=c_bm,
                     parent_top1=(p_hits[0][1] if p_hits else ""),
                     child_top1=(c_hits[0][1] if c_hits else "")))

df = pd.DataFrame(rows)
df.to_csv(REPO/"reviewer_packet_v4_slim/s146_tier0_decomposition.csv", index=False)

# ==== (B) ε scan =========================================================
print("\n" + "=" * 76)
print("(B) ε SCAN on biomass_max under flux bounds")
print("=" * 76)
scan_rows = []
for cell in CELLS:
    cond, hpf = cell.split("_"); hpf = int(hpf)
    child_fp = V7_MODELS / f"trans_rfastcormics_{cell}.mat"
    if not child_fp.exists():
        scan_rows.append(dict(cell=cell, **{f"eps_{e}": np.nan for e in EPS_SCAN}))
        continue
    child = load_model(child_fp, "parent")
    lb0 = child["lb"].copy(); ub0 = child["ub"].copy()
    Snnz_c = np.diff(child["S"].tocsc().indptr); is_ex_c = Snnz_c == 1
    lb, ub, idx, _ = apply_context_bounds(lb0, ub0, child["rxns"], is_ex_c, cond, hpf)

    row = dict(cell=cell)
    for eps in EPS_SCAN:
        old = _set_pool_coef(child, eps)
        bm, _ = biomax_with_res(child, lb, ub, idx)
        row[f"eps_{eps}"] = bm
        _restore_pool_coef(child, old)
    scan_rows.append(row)
    row_str = "  ".join(f"ε={e}: {row[f'eps_{e}']:.4e}" for e in EPS_SCAN)
    # flat if max/min < 1.01
    vals = [row[f"eps_{e}"] for e in EPS_SCAN]
    ratio = max(vals) / min(vals) if min(vals) > 0 else np.nan
    print(f"  {cell:6s}  {row_str}  ratio={ratio:.4f}")

df2 = pd.DataFrame(scan_rows)
df2.to_csv(REPO/"reviewer_packet_v4_slim/s146_epsilon_scan.csv", index=False)
print(f"\nwrote s146_epsilon_scan.csv")

# ==== (C) BOUND_EX scan ==================================================
print("\n" + "=" * 76)
print("(C) BOUND_EX SCAN on biomass_max (fc*BOUND_EX envelope)")
print("=" * 76)
bex_rows = []
for cell in CELLS:
    cond, hpf = cell.split("_"); hpf = int(hpf)
    child_fp = V7_MODELS / f"trans_rfastcormics_{cell}.mat"
    if not child_fp.exists():
        bex_rows.append(dict(cell=cell, **{f"bex_{b}": np.nan for b in BOUND_EX_SCAN}))
        continue
    child = load_model(child_fp, "parent")
    Snnz_c = np.diff(child["S"].tocsc().indptr); is_ex_c = Snnz_c == 1

    row = dict(cell=cell); tops = {}
    for bex in BOUND_EX_SCAN:
        lb0 = child["lb"].copy(); ub0 = child["ub"].copy()
        lb, ub, idx, _ = apply_context_bounds(lb0, ub0, child["rxns"], is_ex_c,
                                                cond, hpf, bound_ex=bex)
        bm, res = biomax_with_res(child, lb, ub, idx)
        row[f"bex_{bex}"] = bm
        # capture top binding for reporting
        hits, _ = constraint_decomposition(res, child, is_ex_c, lb, ub, 1)
        tops[bex] = hits[0][1] if hits else ""
    # cross-cell ratio: biomass_max / biomass_max at bex=10 (reference)
    # (ratios between cells are the interpretable quantity)
    ref = row["bex_10.0"]
    ratio_str = "  ".join(f"bex={b}: {row[f'bex_{b}']:.4e}" for b in BOUND_EX_SCAN)
    scale_ratio = row["bex_30.0"] / row["bex_3.0"] if (row["bex_3.0"] and row["bex_3.0"] > 0) else np.nan
    top10 = tops.get(10.0, "")
    print(f"  {cell:6s}  {ratio_str}  (30/3 ratio={scale_ratio:.2f})  top@bex=10: {top10}")
    bex_rows.append(row)

df3 = pd.DataFrame(bex_rows)
df3.to_csv(REPO/"reviewer_packet_v4_slim/s146_bound_ex_scan.csv", index=False)
print(f"\nwrote s146_bound_ex_scan.csv")

# ==== (D) cross-cell ratio summary at bex=10 =============================
print("\n" + "=" * 76)
print("(D) CROSS-CELL RATIOS at BOUND_EX=10 (the interpretable quantity)")
print("=" * 76)
# reference = LD_72 (median-ish); print each cell's ratio to LD
if "bex_10.0" in df3.columns:
    ref = df3.loc[df3.cell == "LD_72", "bex_10.0"].iloc[0]
    if ref and ref > 0:
        print(f"reference: LD_72 biomax = {ref:.4e}")
        for _, r in df3.iterrows():
            v = r["bex_10.0"]
            print(f"  {r.cell:6s}  biomax={v:.4e}  ratio_vs_LD_72={v/ref:.3f}")
