"""
s82 (corrected 2026-10-05 per §11.15 amendment) --
positive control model construction for LD_72 Steroid metabolism.

Produces two modified LD_72 models from the Fix A tightened base:

  (A) KO  — all Steroid metabolism reactions forced to lb = ub = 0.
            Tier 1 (zero-fraction Fisher) validator.
  (B) KD  — each Steroid metabolism reaction capped at 10 % of its
            Frame A FVA range, preserving sign.
            Tier 2 (W1) validator.

Both are built from the Fix A tightened LD_72 so MVE rounding stays
tractable. After the perturbation, FVA tightening is re-run on the
full model (the KO/KD may make prior tight bounds of neighbouring
reactions infeasible, so bounds must be recomputed).

Writes:
  reviewer_packet_v4_slim/v7_models_ko_steroid/trans_rfastcormics_LD_72.mat
  reviewer_packet_v4_slim/v7_models_kd_steroid/trans_rfastcormics_LD_72.mat

So s41d can address both as --cell LD_72 with different --models_dir.
"""
import sys, h5py, shutil
import numpy as np, scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB
from pathlib import Path

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
SRC = REPO/"reviewer_packet_v4_slim/v7_models_tight/trans_rfastcormics_LD_72.mat"
KO_DIR = REPO/"reviewer_packet_v4_slim/v7_models_ko_steroid"
KD_DIR = REPO/"reviewer_packet_v4_slim/v7_models_kd_steroid"
KO_DIR.mkdir(exist_ok=True, parents=True)
KD_DIR.mkdir(exist_ok=True, parents=True)

TARGET_SUB = "Steroid metabolism"
BIOMASS_RXN = "MAR00021"
KD_FRACTION = 0.10
MARGIN = 1e-4


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
    """Return (new_lb, new_ub) via FVA on every reaction under biomass
    floor = 0.9 × biomass_max of the input bounds."""
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
    print(f"  biomass_max = {bmax:.4f}, floor = {0.9*bmax:.4f}", flush=True)

    mn = np.zeros(n_rxn); mx = np.zeros(n_rxn)
    for j in range(n_rxn):
        if j % 500 == 0:
            print(f"    FVA {j}/{n_rxn}", flush=True)
        mo.setObjective(v[j], GRB.MINIMIZE); mo.optimize()
        mn[j] = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
        mo.setObjective(v[j], GRB.MAXIMIZE); mo.optimize()
        mx[j] = mo.ObjVal if mo.Status == GRB.OPTIMAL else np.nan
    new_lb = mn - MARGIN
    new_ub = mx + MARGIN
    # keep biomass bounds at FVA floor / loaded ceiling (not margin-shrunk)
    new_lb[biomass_idx] = max(0.9 * bmax, lb[biomass_idx])
    new_ub[biomass_idx] = ub[biomass_idx]
    return new_lb, new_ub, bmax


def write_model(src_fp, dst_fp, new_lb, new_ub):
    shutil.copyfile(src_fp, dst_fp)
    with h5py.File(dst_fp, "r+") as f:
        del f["parent/lb"]; f["parent/lb"] = new_lb.reshape(-1, 1)
        del f["parent/ub"]; f["parent/ub"] = new_ub.reshape(-1, 1)
    print(f"  wrote {dst_fp}", flush=True)


# ---------------------------------------------------------------- load base

print(f"[source] Fix A tightened LD_72: {SRC}", flush=True)
S, rxns, subs, lb, ub, n_rxn, n_met = load_model(SRC)
biomass_idx = rxns.index(BIOMASS_RXN)
steroid_idx = [j for j, s in enumerate(subs) if s == TARGET_SUB]
print(f"[info] {len(steroid_idx)} reactions in '{TARGET_SUB}' (of {n_rxn} total)", flush=True)
if len(steroid_idx) == 0:
    sys.exit("[FATAL] no reactions found in target subsystem")

# ---------------------------------------------------------------- KD build

print(f"\n=== Knockdown ({int(KD_FRACTION*100)}% cap, preserve sign) ===", flush=True)
# Current Fix A tight bounds for steroid reactions are already ≈ FVA ranges
# under biomass_lb = 11.41, so cap at 10 % of those bounds directly.
kd_lb = lb.copy(); kd_ub = ub.copy()
for j in steroid_idx:
    kd_ub[j] = KD_FRACTION * max(ub[j], 0.0)
    kd_lb[j] = KD_FRACTION * min(lb[j], 0.0)
print(f"[KD] applied cap on {len(steroid_idx)} steroid reactions", flush=True)
print(f"[KD] FVA-tightening full KD model (recomputing all bounds)...", flush=True)
new_lb_kd, new_ub_kd, bmax_kd = fva_tighten(S, kd_lb, kd_ub, n_rxn, n_met, biomass_idx)
# Preserve the explicit KD caps even if FVA says a steroid rxn could flux
# more (it can't under the LD_72 flux space, but belt-and-braces)
for j in steroid_idx:
    new_lb_kd[j] = max(new_lb_kd[j], kd_lb[j])
    new_ub_kd[j] = min(new_ub_kd[j], kd_ub[j])
write_model(SRC, KD_DIR/"trans_rfastcormics_LD_72.mat", new_lb_kd, new_ub_kd)
print(f"[KD] biomass_max = {bmax_kd:.4f}", flush=True)

# ---------------------------------------------------------------- KO build

print(f"\n=== Knockout (lb = ub = 0 on all steroid reactions) ===", flush=True)
ko_lb = lb.copy(); ko_ub = ub.copy()
ko_lb[steroid_idx] = 0.0; ko_ub[steroid_idx] = 0.0
print(f"[KO] applied full KO on {len(steroid_idx)} steroid reactions", flush=True)
print(f"[KO] FVA-tightening full KO model (recomputing all bounds)...", flush=True)
new_lb_ko, new_ub_ko, bmax_ko = fva_tighten(S, ko_lb, ko_ub, n_rxn, n_met, biomass_idx)
# Preserve KO zeros for steroid
for j in steroid_idx:
    new_lb_ko[j] = 0.0
    new_ub_ko[j] = 0.0
write_model(SRC, KO_DIR/"trans_rfastcormics_LD_72.mat", new_lb_ko, new_ub_ko)
print(f"[KO] biomass_max = {bmax_ko:.4f}", flush=True)

print(f"\nDONE. Next: preprocess both models with s41d.")
print(f"  KD cmd: python scripts/s41d_hopsy_sample_direct.py --cell LD_72 \\")
print(f"           --models_dir ../reviewer_packet_v4_slim/v7_models_kd_steroid \\")
print(f"           --outdir_name hopsy_v7_kd_steroid --preprocess_only --skip_apply_bounds")
print(f"  KO cmd: same, --models_dir v7_models_ko_steroid --outdir_name hopsy_v7_ko_steroid")
