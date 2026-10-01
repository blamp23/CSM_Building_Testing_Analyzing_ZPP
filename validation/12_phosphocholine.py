"""
s148 -- phosphocholine report (read-only).
Checks: in 105 measured? which cells flagged by s128c? route present?
"""
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
from pathlib import Path

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
mb = pd.read_csv(REPO / "results/qc/phase2_exchange_bounds_long.csv")
cands = pd.read_csv(REPO / "reviewer_packet_v4_slim/s128c_candidates_v7.csv")

print("=== phosphocholine (MAR09845 / MAM02738) ===")
in105 = "MAR09845" in set(mb.ex_rxn.unique())
print(f"in 105 measured? {in105}")

# cells flagged
pc = cands[cands.species == "MAM02738"] if "species" in cands.columns else pd.DataFrame()
print(f"\ns128c flagged in {len(pc)} cells:")
if len(pc):
    print(pc[["cell","ex_rxn","met_name","baseline_biomax","knockout_biomax","ratio"]].to_string(index=False))

# route check on v7 base at hpf 48 (where some cells were flagged)
def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)

with h5py.File(REPO / "ocr_anchored_extraction/models/baked_ocr_48_v7.mat", "r") as f:
    m = f["m"]
    rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
    mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
    metNames = [_resolve(f, r) for r in np.asarray(m["metNames"]).ravel()]
    rxnNames = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
    S_grp = m["S"]
    S = sp.csc_matrix((np.asarray(S_grp["data"]), np.asarray(S_grp["ir"]), np.asarray(S_grp["jc"])),
                      shape=(int(np.asarray(S_grp["ir"]).max())+1, len(np.asarray(S_grp["jc"]))-1))
    lb = np.asarray(m["lb"]).ravel()
    ub = np.asarray(m["ub"]).ravel()

print(f"\n=== reactions touching MAM02738 (phosphocholine) in v7 base ===")
targets = [mm for mm in mets if "MAM02738" in mm]
print(f"compartments: {targets}")
for mm in targets:
    i = mets.index(mm)
    row = S.getrow(i).toarray().ravel()
    for j in np.where(row != 0)[0]:
        col = S.getcol(j).toarray().ravel()
        nz = np.where(col != 0)[0]
        subs = " + ".join(f"{-col[k]:g} {mets[k]}" for k in nz if col[k] < 0)
        prods = " + ".join(f"{col[k]:g} {mets[k]}" for k in nz if col[k] > 0)
        role = "produces" if row[j] > 0 else "consumes"
        print(f"  {rxns[j]:10s} ({role})  [{lb[j]:+g},{ub[j]:+g}]  {rxnNames[j][:40]}")
        print(f"    {subs}  ->  {prods}")
