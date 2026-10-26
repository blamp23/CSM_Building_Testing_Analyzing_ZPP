"""
s82 -- §11.15 positive control: knock out Steroid metabolism in LD_72.

Reads the v7 LD_72 extracted model, zeros out lb and ub for every reaction
whose subsystem is "Steroid metabolism", writes the KO model as
v7_models_ko_steroid/trans_rfastcormics_LD_72_KO_steroid.mat with the same
structure as the other v7 models so s41d_hopsy_sample_direct.py can
preprocess + sample it without further changes.

Reports: n_reactions in Steroid metabolism, n_reactions zeroed, biomass_max
before and after (via Gurobi), and the KO cell id to use downstream.
"""
import numpy as np, scipy.sparse as sp, h5py, shutil
from pathlib import Path
import gurobipy as gp
from gurobipy import GRB

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
SRC_DIR = REPO/"reviewer_packet_v4_slim/v7_models"
KO_DIR  = REPO/"reviewer_packet_v4_slim/v7_models_ko_steroid"
KO_DIR.mkdir(parents=True, exist_ok=True)

CELL = "LD_72"
TARGET_SUB = "Steroid metabolism"
BIOMASS_RXN = "MAR00021"

src = SRC_DIR/f"trans_rfastcormics_{CELL}.mat"
dst = KO_DIR/f"trans_rfastcormics_{CELL}_KO_steroid.mat"
shutil.copyfile(src, dst)

def resolve_str(f, obj):
    arr = np.asarray(obj).ravel()
    try: return "".join(chr(int(c)) for c in arr)
    except: return ""

# Identify reactions in the target subsystem
with h5py.File(dst, "r") as f:
    m = f["parent"]
    S_grp = m["S"]
    data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
    n_rxn = len(jc) - 1; n_met = int(ir.max()) + 1
    S = sp.csc_matrix((data, ir, jc), shape=(n_met, n_rxn))
    rxns = [resolve_str(f, f[r]) for r in np.asarray(m["rxns"]).ravel()]
    sub_refs = np.asarray(m["subSystems"]).ravel()
    subs = []
    for ref in sub_refs:
        obj = f[ref]; arr = np.asarray(obj).ravel()
        if arr.dtype == object and len(arr) > 0:
            subs.append(resolve_str(f, f[arr[0]]))
        else:
            subs.append(resolve_str(f, obj))
    lb = np.asarray(m["lb"]).ravel().astype(float)
    ub = np.asarray(m["ub"]).ravel().astype(float)
    biomass_idx = rxns.index(BIOMASS_RXN)

ko_idx = [j for j, s in enumerate(subs) if s == TARGET_SUB]
print(f"Target subsystem: {TARGET_SUB}")
print(f"n reactions in {TARGET_SUB}: {len(ko_idx)}")
print(f"n total reactions in {CELL}: {n_rxn}")
if len(ko_idx) == 0:
    print(f"[FATAL] no reactions found in {TARGET_SUB} for {CELL}; aborting")
    raise SystemExit(1)

# biomass_max BEFORE KO (sanity)
env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
mo = gp.Model("bm_before", env=env)
mo.setParam("Threads", 2); mo.setParam("OutputFlag", 0)
v = mo.addMVar(n_rxn, lb=lb, ub=ub)
mo.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
mo.setObjective(v[biomass_idx], GRB.MAXIMIZE)
mo.optimize()
bm_before = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
print(f"biomass_max BEFORE KO: {bm_before:.4f}")

# Apply KO: lb = ub = 0 for all Steroid metabolism reactions
new_lb = lb.copy(); new_ub = ub.copy()
new_lb[ko_idx] = 0.0; new_ub[ko_idx] = 0.0

# biomass_max AFTER KO
mo = gp.Model("bm_after", env=env)
mo.setParam("Threads", 2); mo.setParam("OutputFlag", 0)
v = mo.addMVar(n_rxn, lb=new_lb, ub=new_ub)
mo.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
mo.setObjective(v[biomass_idx], GRB.MAXIMIZE)
mo.optimize()
bm_after = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
print(f"biomass_max AFTER  KO: {bm_after:.4f}")
if not np.isfinite(bm_after) or bm_after < 0.1:
    print(f"[WARN] biomass_max collapsed under KO ({bm_after}). Model may be unbuildable.")

# Write the modified lb/ub back
with h5py.File(dst, "r+") as f:
    del f["parent/lb"]; f["parent/lb"] = new_lb.reshape(-1, 1)
    del f["parent/ub"]; f["parent/ub"] = new_ub.reshape(-1, 1)

print(f"wrote {dst}")
print(f"next step: preprocess + sample this KO model with s41d, outdir_name=hopsy_v7_ko_steroid")
print(f"  cell token = '{CELL}_KO_steroid' (requires passing the KO cell name to s41d)")
