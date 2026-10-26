"""
s82b (§11.16 redesign) -- percentile-based KD for Steroid metabolism in LD_72.

Instead of capping at 10% of FVA max (which never binds for interior
polytope samples), caps at 10% of the *actual sampled* max |v| per reaction
from the original LD_72 pooled samples.

Starts from Fix A tightened LD_72, applies new caps to Steroid reactions,
re-FVA tightens everything, writes v7_models_kd2_steroid/.
"""
import sys, h5py, shutil
import numpy as np, pandas as pd, scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB
from pathlib import Path

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
SRC = REPO/"reviewer_packet_v4_slim/v7_models_tight/trans_rfastcormics_LD_72.mat"
POOLED = REPO/"grace_v4/results/pooled_v7/LD_72_pooled.csv"
KD2_DIR = REPO/"reviewer_packet_v4_slim/v7_models_kd2_steroid"
KD2_DIR.mkdir(exist_ok=True, parents=True)

TARGET_SUB = "Steroid metabolism"
BIOMASS_RXN = "MAR00021"
KD_FRACTION = 0.10
MARGIN = 1e-4
ZERO_TOL = 1e-10


def resolve_str(f, obj):
    arr = np.asarray(obj).ravel()
    try: return "".join(chr(int(c)) for c in arr)
    except: return ""


def load_model(mat_fp):
    with h5py.File(mat_fp, "r") as f:
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
    return S, rxns, subs, lb, ub, n_rxn, n_met


def fva_tighten(S, lb, ub, n_rxn, n_met, biomass_idx):
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    mo = gp.Model("fva", env=env)
    mo.setParam("Threads", 4); mo.setParam("OutputFlag", 0)
    v = mo.addMVar(n_rxn, lb=lb, ub=ub)
    mo.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    mo.setObjective(v[biomass_idx], GRB.MAXIMIZE); mo.optimize()
    if mo.Status != GRB.OPTIMAL:
        raise RuntimeError(f"biomass_max LP failed: status={mo.Status}")
    bmax = mo.ObjVal
    v[biomass_idx].LB = max(0.9 * bmax, lb[biomass_idx])
    print(f"  biomass_max = {bmax:.4f}, floor = {max(0.9*bmax, lb[biomass_idx]):.4f}", flush=True)
    mn = np.zeros(n_rxn); mx = np.zeros(n_rxn)
    for j in range(n_rxn):
        if j % 500 == 0:
            print(f"    FVA {j}/{n_rxn}", flush=True)
        mo.setObjective(v[j], GRB.MINIMIZE); mo.optimize()
        mn[j] = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
        mo.setObjective(v[j], GRB.MAXIMIZE); mo.optimize()
        mx[j] = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
    new_lb = mn - MARGIN; new_ub = mx + MARGIN
    new_lb[biomass_idx] = max(0.9 * bmax, lb[biomass_idx])
    new_ub[biomass_idx] = ub[biomass_idx]
    return new_lb, new_ub, bmax


print(f"[source] Fix A tightened LD_72: {SRC}", flush=True)
S, rxns, subs, lb, ub, n_rxn, n_met = load_model(SRC)
biomass_idx = rxns.index(BIOMASS_RXN)
steroid_idx = [j for j, s in enumerate(subs) if s == TARGET_SUB]
print(f"[info] {len(steroid_idx)} reactions in '{TARGET_SUB}'", flush=True)

# Read original pooled samples to get per-reaction empirical max |v|
print(f"[info] loading original LD_72 pooled samples for percentile caps...", flush=True)
XA = pd.read_csv(POOLED, index_col=0)
rxn_to_pos = {r: j for j, r in enumerate(rxns)}

print(f"\n=== Knockdown (B2): 10% of EMPIRICAL max |v| per rxn ===", flush=True)
kd_lb = lb.copy(); kd_ub = ub.copy()
diag_rows = []
for j in steroid_idx:
    rxn_id = rxns[j]
    if rxn_id not in XA.index:
        continue
    samples = XA.loc[rxn_id].values
    pos_samples = samples[samples > ZERO_TOL]
    neg_samples = samples[samples < -ZERO_TOL]
    max_pos = pos_samples.max() if len(pos_samples) else 0.0
    min_neg = neg_samples.min() if len(neg_samples) else 0.0
    # Preserve sign: cap on each end at 10% of the empirical extreme (if any)
    new_ub_val = KD_FRACTION * max_pos if max_pos > 0 else min(0.0, ub[j])
    new_lb_val = KD_FRACTION * min_neg if min_neg < 0 else max(0.0, lb[j])
    kd_ub[j] = min(new_ub_val, ub[j])
    kd_lb[j] = max(new_lb_val, lb[j])
    diag_rows.append(dict(rxn=rxn_id,
                           orig_ub=ub[j], orig_lb=lb[j],
                           emp_max=max_pos, emp_min=min_neg,
                           kd_ub=kd_ub[j], kd_lb=kd_lb[j]))
diag = pd.DataFrame(diag_rows)
print(diag.to_string(index=False))
diag.to_csv(KD2_DIR/"steroid_kd_bounds.csv", index=False)
print(f"[KD2] applied empirical-cap on {len(diag)} steroid reactions", flush=True)

print(f"[KD2] FVA-tightening full KD2 model...", flush=True)
new_lb_kd, new_ub_kd, bmax_kd = fva_tighten(S, kd_lb, kd_ub, n_rxn, n_met, biomass_idx)
# Preserve explicit KD2 caps on steroid
for j in steroid_idx:
    new_lb_kd[j] = max(new_lb_kd[j], kd_lb[j])
    new_ub_kd[j] = min(new_ub_kd[j], kd_ub[j])

dst = KD2_DIR/"trans_rfastcormics_LD_72.mat"
shutil.copyfile(SRC, dst)
with h5py.File(dst, "r+") as f:
    del f["parent/lb"]; f["parent/lb"] = new_lb_kd.reshape(-1, 1)
    del f["parent/ub"]; f["parent/ub"] = new_ub_kd.reshape(-1, 1)
print(f"[KD2] biomass_max = {bmax_kd:.4f}", flush=True)
print(f"[KD2] wrote {dst}", flush=True)
print(f"\nNext: preprocess + sample")
print(f"  python scripts/s41d_hopsy_sample_direct.py --cell LD_72 \\")
print(f"      --models_dir ../reviewer_packet_v4_slim/v7_models_kd2_steroid \\")
print(f"      --outdir_name hopsy_v7_kd2_steroid --preprocess_only --skip_apply_bounds")
