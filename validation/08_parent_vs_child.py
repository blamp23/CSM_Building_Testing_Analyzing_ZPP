"""
s143 -- two pre-launch sanity checks.

(1) Parent >= child. For each of 15 cells, load the hpf-appropriate v7
    baseline (baked_ocr_{hpf}_v7.mat) and apply that cell's fc envelope
    to measured exchanges. Compute biomass_max + dominant binding
    constraint. Compare to the extracted model's biomass_max under the
    same fc envelope (apply_sampling_bounds convention). Assert
    extracted <= base.

(2) MAR00463 — formula / subsystem / GPR / flux consistency on v7 base.
"""
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.optimize import linprog
from pathlib import Path

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
V7_MODELS = REPO / "reviewer_packet_v4_slim/v7_models"
BASE_DIR  = REPO / "ocr_anchored_extraction/models"
BOUNDS    = REPO / "results/qc/phase2_exchange_bounds_long.csv"
BIOMASS_RXN = "MAR00021"
O2_EX = "MAR09048"
CELLS = [f"{c}_{h}" for c in ("BL","D","LD") for h in (24,48,72,96,120)]


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
        except Exception:
            metNames = [""] * len(mets)
        try:
            rxnNames = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
        except Exception:
            rxnNames = [""] * len(rxns)
        try:
            grRules = [_resolve(f, r) for r in np.asarray(m["grRules"]).ravel()]
        except Exception:
            grRules = [""] * len(rxns)
        try:
            subSystems = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception:
            subSystems = [""] * len(rxns)
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        S = sp.csc_matrix((data, ir, jc), shape=(int(ir.max())+1, len(jc)-1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, metNames=metNames, rxnNames=rxnNames,
                grRules=grRules, subSystems=subSystems, S=S, lb=lb, ub=ub)


def apply_fc_envelope(m, cond, hpf, mb, base_bounds=True):
    """Apply cell-specific fc envelope to measured exchanges.
    base_bounds=True: preserve everything else (base curation applies).
    """
    lb = m["lb"].copy(); ub = m["ub"].copy()
    idx = {r: i for i, r in enumerate(m["rxns"])}
    Snnz = np.diff(m["S"].tocsc().indptr); is_ex = Snnz == 1
    sel = (mb.condition == cond) & (mb.hpf == hpf)
    BOUND_EX = 10.0
    for rid, fc in zip(mb.loc[sel, "ex_rxn"], mb.loc[sel, "fc"]):
        if rid in idx:
            j = idx[rid]
            fc = float(fc) if np.isfinite(fc) and fc > 0 else 1.0
            # apply fc as a multiplier on ±BOUND_EX for measured exchanges
            lb[j] = -BOUND_EX * fc
            ub[j] = +BOUND_EX * fc
    return lb, ub, idx, is_ex


def biomax_with_dual(m, lb, ub, idx):
    n_rxn = len(m["rxns"]); n_met = m["S"].shape[0]
    c = np.zeros(n_rxn); c[idx[BIOMASS_RXN]] = -1.0
    A_eq = m["S"].toarray(); b_eq = np.zeros(n_met)
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=list(zip(lb, ub)), method="highs")
    if not res.success: return np.nan, None
    return -res.fun, res


def top_binding(res, m, is_ex):
    """Return the top-shadow-price non-biomass constraint."""
    lb_marg = res.lower.marginals
    ub_marg = res.upper.marginals
    idx = {r: i for i, r in enumerate(m["rxns"])}
    hits = []
    for j in range(len(m["rxns"])):
        if m["rxns"][j] == BIOMASS_RXN: continue
        s = max(abs(lb_marg[j]), abs(ub_marg[j]))
        if s > 1e-6:
            hits.append((s, m["rxns"][j], m["rxnNames"][j][:40] if isinstance(m["rxnNames"][j], str) else str(m["rxnNames"][j])[:40],
                        m["lb"][j], m["ub"][j], res.x[j], bool(is_ex[j])))
    hits.sort(reverse=True)
    return hits[:5]  # top 5


mb = pd.read_csv(BOUNDS)

# ---- (1) parent >= child ------------------------------------------------
print("=" * 70)
print("(1) PARENT (v7 base + fc envelope) vs CHILD (v7 extracted + fc envelope)")
print("=" * 70)
rows = []
last_hpf = -1; base_cache = None
for cell in CELLS:
    cond, hpf = cell.split("_"); hpf = int(hpf)
    if hpf != last_hpf:
        base_fp = BASE_DIR / f"baked_ocr_{hpf}_v7.mat"
        print(f"\n[loading v7 base {base_fp.name} …]")
        base_cache = load_model(base_fp, "m")
        last_hpf = hpf

    # parent
    lb, ub, idx_p, is_ex_p = apply_fc_envelope(base_cache, cond, hpf, mb)
    parent_bm, parent_res = biomax_with_dual(base_cache, lb, ub, idx_p)
    parent_binding = top_binding(parent_res, base_cache, is_ex_p) if parent_res else []

    # child
    child_fp = V7_MODELS / f"trans_rfastcormics_{cell}.mat"
    child = load_model(child_fp, "parent")
    lb2, ub2, idx_c, is_ex_c = apply_fc_envelope(child, cond, hpf, mb)
    child_bm, child_res = biomax_with_dual(child, lb2, ub2, idx_c)

    ok = child_bm <= parent_bm + 1e-4
    top1 = f"{parent_binding[0][1]} shadow={parent_binding[0][0]:.3f}" if parent_binding else "(none)"
    rows.append(dict(cell=cell, parent=parent_bm, child=child_bm,
                     delta=child_bm - parent_bm, parent_ok=ok,
                     parent_top_binding=top1))
    mark = "OK" if ok else "VIOLATE"
    print(f"  {cell:6s}  parent={parent_bm:.4f}  child={child_bm:.4f}  "
          f"Δ={child_bm - parent_bm:+.4f}  [{mark}]  top-binding: {top1}")
    if parent_binding:
        for s, rxn, name, lb_v, ub_v, v, ex in parent_binding[:3]:
            print(f"    dual  {rxn:10s}  shadow={s:.3f}  lb={lb_v:+.2f} ub={ub_v:+.2f} v={v:+.3f}  {'EX' if ex else 'IN'}  {name}")

df = pd.DataFrame(rows)
df.to_csv(REPO / "reviewer_packet_v4_slim/s143_parent_vs_child.csv", index=False)
print(f"\nwrote s143_parent_vs_child.csv")
n_ok = df.parent_ok.sum()
print(f"parent_ok gate: {n_ok}/15")

# ---- (2) MAR00463 -------------------------------------------------------
print("\n" + "=" * 70)
print("(2) MAR00463 diagnostic")
print("=" * 70)
m = base_cache  # last hpf loaded; MAR00463 is topology-only, hpf-invariant
if "MAR00463" not in m["rxns"]:
    print("MAR00463 NOT in base"); raise SystemExit(1)
j = m["rxns"].index("MAR00463")
col = m["S"].getcol(j).toarray().ravel()
nz = np.where(col != 0)[0]
subs = " + ".join(f"{-col[i]:g} {m['mets'][i]} ({m['metNames'][i]})" for i in nz if col[i] < 0)
prods = " + ".join(f"{col[i]:g} {m['mets'][i]} ({m['metNames'][i]})" for i in nz if col[i] > 0)
print(f"formula:  {subs}  ->  {prods}")
print(f"bounds:   lb={m['lb'][j]}, ub={m['ub'][j]}")
subsys = m["subSystems"][j] if isinstance(m["subSystems"][j], str) else str(m["subSystems"][j])
print(f"subsystem: {subsys}")
grrule = m["grRules"][j] if isinstance(m["grRules"][j], str) else str(m["grRules"][j])
print(f"GPR:       {grrule if grrule else '(no gene rule)'}")

# flux consistency: max and min v[MAR00463] subject to S @ v = 0, lb<=v<=ub
n_rxn = len(m["rxns"]); n_met = m["S"].shape[0]
c_pos = np.zeros(n_rxn); c_pos[j] = -1.0  # max v_j
A_eq = m["S"].toarray(); b_eq = np.zeros(n_met)
res_max = linprog(c_pos, A_eq=A_eq, b_eq=b_eq, bounds=list(zip(m["lb"], m["ub"])), method="highs")
c_neg = np.zeros(n_rxn); c_neg[j] = 1.0   # min v_j
res_min = linprog(c_neg, A_eq=A_eq, b_eq=b_eq, bounds=list(zip(m["lb"], m["ub"])), method="highs")
v_max = -res_max.fun if res_max.success else None
v_min = res_min.fun if res_min.success else None
print(f"flux range: v ∈ [{v_min}, {v_max}]")
inconsistent = (v_max is not None and v_min is not None and abs(v_max) < 1e-6 and abs(v_min) < 1e-6)
print(f"flux-consistent? {'NO — BLOCKED' if inconsistent else 'yes'}")
