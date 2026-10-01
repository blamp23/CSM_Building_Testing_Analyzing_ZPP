"""
s132 -- Step 3 closure check on v7 base.

Compares biomass_max under extraction bounds:
   v6 base = baked_ocr_72.mat  (unchanged Wang + OCR bake)
   v7 base = baked_ocr_72_v7.mat  (pool coef edits + 17 exchange closures)

Pass = v7 biomass_max >= v6 biomass_max within 1e-4.
Fail = per-close_exchange reopen sweep: reopen one closed exchange at a
time and record biomass_max; anything that recovers biomass identifies a
gap in the internal synthesis path for that metabolite. Report the gap;
do NOT force reactions.
"""
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.optimize import linprog
from pathlib import Path

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
V6_MAT = REPO / "ocr_anchored_extraction/models/baked_ocr_72.mat"
V7_MAT = REPO / "ocr_anchored_extraction/models/baked_ocr_72_v7.mat"
LEDGER = REPO / "reviewer_packet_v4_slim/s130_v7_ledger.csv"
BIOMASS_RXN = "MAR00021"
TOL = 1e-4


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
    return dict(rxns=rxns, mets=mets, metNames=metNames, S=S, lb=lb, ub=ub)


def biomax(m, lb=None, ub=None):
    if lb is None: lb = m["lb"]
    if ub is None: ub = m["ub"]
    idx = {r: i for i, r in enumerate(m["rxns"])}
    n_rxn = len(m["rxns"]); n_met = m["S"].shape[0]
    c = np.zeros(n_rxn); c[idx[BIOMASS_RXN]] = -1.0
    A_eq = m["S"].toarray(); b_eq = np.zeros(n_met)
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=list(zip(lb, ub)), method="highs")
    return -res.fun if res.success else np.nan


print("=== v6 base FBA ===")
v6 = load_model(V6_MAT, "m")
v6_max = biomax(v6)
print(f"v6 base biomass_max: {v6_max:.6f}")

if not V7_MAT.exists():
    print(f"\nv7 base not built yet ({V7_MAT.name}). Run s135 first.")
    raise SystemExit(1)

print("\n=== v7 base FBA (Part 1 + Part 2 as built by s135) ===")
v7 = load_model(V7_MAT, "m")
v7_max = biomax(v7)
print(f"v7 base biomass_max (Part 1 + Part 2): {v7_max:.6f}")

# --- Part 1 alone: rebuild bounds in memory ------------------------------
# Reload base + apply only pool edits + close + Part 1 (ub=1000 for
# unmeasured non-close non-keep non-medium). Skip Part 2 (no lb=0 uptake close).
print("\n=== v7 Part-1-only variant (in-memory, for comparison) ===")
mb = pd.read_csv(REPO/"results/qc/phase2_exchange_bounds_long.csv")
measured = set(mb.ex_rxn.unique())
ledger_df = pd.read_csv(LEDGER)
close_set = set(ledger_df.loc[ledger_df.action == "close_exchange", "exchange_MAR_id"].dropna())
keep_set  = set(ledger_df.loc[ledger_df.action == "keep_exchange", "exchange_MAR_id"].dropna())
MEDIUM = {"MAR09047","MAR09058","MAR09072","MAR09073","MAR09074","MAR09076",
          "MAR09077","MAR09078","MAR09079","MAR09080","MAR09081","MAR09082",
          "MAR09148","MAR09150","MAR13066","MAR13072","MAR13073"}
# start from v6 base + pool edits + close (i.e. v7 Step 1+2 minus Part 1/2).
# Easier: start from v7 (which already has pool edits + close + Part 1 + Part 2),
# revert Part 2 by relaxing lb=0 back to -10 for the rows Part 2 closed.
import scipy.sparse as sp
v7_p1 = dict(v7)
lb1 = v7["lb"].copy(); ub1 = v7["ub"].copy()
Snnz = np.diff(v7["S"].tocsc().indptr); is_ex = Snnz == 1
for j in np.where(is_ex)[0]:
    rid = v7["rxns"][j]
    if rid == "MAR09048": continue
    if rid in close_set: continue
    if rid in measured: continue
    if rid in MEDIUM: continue
    if rid in keep_set: continue
    # this exchange got Part 1 + Part 2 (lb=0, ub=1000); Part 1 alone means
    # ub=1000 but lb reset to -10 (the previous default).
    lb1[j] = -10.0
p1_max = biomax(v7_p1, lb=lb1, ub=ub1)
print(f"v7 base biomass_max (Part 1 only): {p1_max:.6f}")

delta = v7_max - v6_max
print(f"\nΔ (Part1+Part2 vs v6) = {delta:+.6f}  ({100*delta/v6_max:+.2f}%)")

# --- audit: every exchange at its bound in the final base optimum -------
def biomax_with_duals(m, lb, ub):
    n_rxn = len(m["rxns"]); n_met = m["S"].shape[0]
    idx = {r: i for i, r in enumerate(m["rxns"])}
    c = np.zeros(n_rxn); c[idx[BIOMASS_RXN]] = -1.0
    A_eq = m["S"].toarray(); b_eq = np.zeros(n_met)
    return linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=list(zip(lb, ub)), method="highs")

print("\n=== every exchange at its bound in v7 (Part 1 + Part 2) optimum ===")
res = biomax_with_duals(v7, v7["lb"], v7["ub"])
Snnz = np.diff(v7["S"].tocsc().indptr); is_ex_v7 = Snnz == 1
rows_at_bound = []
for j in np.where(is_ex_v7)[0]:
    v = res.x[j]; lb = v7["lb"][j]; ub = v7["ub"][j]
    if abs(v) < 1e-6: continue
    at_ub = abs(v - ub) < 1e-4 and ub < 999
    at_lb = abs(v - lb) < 1e-4 and lb > -999
    if at_ub or at_lb:
        col = v7["S"].getcol(j).toarray().ravel()
        mm = v7["rxns"][j] # fallback
        for i in np.where(col != 0)[0]:
            mm = v7["mets"][i]; break
        rows_at_bound.append(dict(rxn=v7["rxns"][j], metabolite=mm,
                                   lb=lb, ub=ub, v=v,
                                   at="ub" if at_ub else "lb"))
if rows_at_bound:
    df = pd.DataFrame(rows_at_bound)
    df.to_csv(REPO/"reviewer_packet_v4_slim/s132_bounds_at_optimum.csv", index=False)
    print(df.to_string(index=False))
else:
    print("(no exchanges at their bound with nonzero flux)")

if v7_max >= v6_max - TOL:
    print("\nPASS — v7 base is feasible at ≥ v6 biomass_max.")
    raise SystemExit(0)

print("\nFAIL — v7 base biomass_max dropped. Running per-exchange reopen sweep …")
ledger = pd.read_csv(LEDGER)
close_rows = ledger[ledger.action == "close_exchange"].copy()
idx = {r: i for i, r in enumerate(v7["rxns"])}
rows = []
for _, row in close_rows.iterrows():
    rid = row.exchange_MAR_id
    if not isinstance(rid, str) or rid == "":
        continue
    if rid not in idx:
        continue
    j = idx[rid]
    lb2 = v7["lb"].copy(); ub2 = v7["ub"].copy()
    # reopen this exchange only
    lb2[j] = -10.0; ub2[j] = 10.0
    bm = biomax(v7, lb2, ub2)
    rec = dict(exchange=rid, metabolite=row.metabolite,
               **{"class": row["class"]}, biomax=bm, delta=bm - v7_max)
    rows.append(rec)
    print(f"  reopen {rid} ({row.metabolite:32s}) → biomax = {bm:.4f} "
          f"(Δ = {bm - v7_max:+.4f})")
sweep = pd.DataFrame(rows).sort_values("biomax", ascending=False)
sweep.to_csv(REPO/"reviewer_packet_v4_slim/s132_reopen_sweep.csv", index=False)
print("\nwrote s132_reopen_sweep.csv")
print("\nAny row with a large positive delta identifies a pathway gap in v7 —")
print("stop and report which metabolite; do NOT force reactions.")
