"""
s169 -- Energy-gate validation on masked base-model corners.

The attribution (s168) uses the shared v7 base model + a per-cell reaction
mask to build 2×2 corners (mask_A × env_B, etc.). apply_context_bounds was
called with skip_energy_gate=True to pass the base-model LP (which fails
the gate because it has free-energy routes rFASTCORMICS prunes).

This script runs the §11.19 energy gate explicitly on EACH corner
(mask × env combination), with the mask applied. The gate passes iff the
masked LP is structurally equivalent to an extracted model (which does pass).

If every corner passes the gate, the s168 attribution values are clean.
If any corner fails, we have a phantom route and the attribution is
contaminated.

Output: v7_susceptibility/energy_gate_mask_check.csv  (one row per corner)
"""
import sys
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds, MEDIUM, OCR, _load_demand_bounds

V7_BASES = REPO / "ocr_anchored_extraction/models"
V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"

BIOMASS_RXN = "MAR00021"
ATP_RXN = "MAR03964"
O2_EX = "MAR09048"

TRIGGERS = [
    ("BL_24",  "LD_24",  "GSH"),
    ("BL_24",  "D_24",   "GSH"),
    ("D_120",  "LD_120", "GSH"),
    ("BL_120", "D_120",  "GSH"),
    ("D_96",   "LD_96",  "NADPH"),
    ("BL_96",  "D_96",   "NADPH"),
    ("BL_96",  "D_96",   "UDP_GlcA"),
    ("BL_72",  "D_72",   "UDP_GlcA"),
    ("D_72",   "LD_72",  "UDP_GlcA"),
    ("D_120",  "LD_120", "UDP_GlcA"),
    ("BL_72",  "D_72",   "PAPS"),
    ("D_72",   "LD_72",  "PAPS"),
    ("BL_96",  "LD_96",  "PAPS"),
    ("BL_120", "LD_120", "PAPS"),
]


def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)


def _load(fp, root_key):
    with h5py.File(fp, "r") as f:
        m = f[root_key]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = int(ir.max()) + 1 if len(ir) else 0
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, S=S, lb=lb, ub=ub)


def _atp_max_null_uptake(mdl, lb, ub, idx, cond, hpf):
    """Replicate §11.19 energy gate: close every exchange for uptake except
    MEDIUM, glucose at -1, O2 at -1000, free-supplied at demand bounds.
    Max MAR03964."""
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb2 = lb.copy(); ub2 = ub.copy()
    for j in np.where(is_ex)[0]:
        lb2[j] = 0.0
    for rid in MEDIUM:
        if rid in idx: lb2[idx[rid]] = -1000.0
    demand = _load_demand_bounds()
    if demand is not None:
        for (c, h, rid), lbv in demand.items():
            if c == cond and h == hpf and rid in idx:
                lb2[idx[rid]] = lbv
    if "MAR09034" in idx:
        lb2[idx["MAR09034"]] = -1.0; ub2[idx["MAR09034"]] = -1.0
    if O2_EX in idx:
        lb2[idx[O2_EX]] = -1000.0; ub2[idx[O2_EX]] = 0.0
    if BIOMASS_RXN in idx:
        j = idx[BIOMASS_RXN]; lb2[j] = 0.0; ub2[j] = max(ub2[j], 1000.0)
    if ATP_RXN in idx:
        ub2[idx[ATP_RXN]] = 1e9
    S = mdl["S"].tocsr()
    n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(n_rxn, lb=lb2, ub=ub2)
    gm.addMConstr(S, v, "=", np.zeros(n_met))
    obj = np.zeros(n_rxn); obj[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


# Pre-load rxns per cell + bases per hpf
cells_needed = set()
for a, b, _ in TRIGGERS: cells_needed.update((a, b))

cell_rxns = {}
for cell in cells_needed:
    fp = V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat"
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        cell_rxns[cell] = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]

base_by_hpf = {}
for hpf in (24, 48, 72, 96, 120):
    base_by_hpf[hpf] = _load(V7_BASES / f"baked_ocr_{hpf}_v7.mat", "m")

# For each corner, apply context bounds + widen + mask, then run gate
rows = []
for a, b, cof in TRIGGERS:
    cond_a, hpf_a = a.split("_"); hpf_a = int(hpf_a)
    cond_b, hpf_b = b.split("_"); hpf_b = int(hpf_b)
    assert hpf_a == hpf_b
    base = base_by_hpf[hpf_a]
    is_ex = (np.diff(base["S"].tocsc().indptr) == 1)
    idx = {r: i for i, r in enumerate(base["rxns"])}

    for corner, mask_cell, env_cell in [
        ("AA", a, a), ("AB", a, b), ("BA", b, a), ("BB", b, b)
    ]:
        env_cond, env_hpf = env_cell.split("_"); env_hpf = int(env_hpf)
        lb, ub, _, _ = apply_context_bounds(
            base["lb"].copy(), base["ub"].copy(), base["rxns"], is_ex,
            env_cond, env_hpf, S=base["S"], skip_energy_gate=True)
        non_ex = ~is_ex
        lb[(lb <= -999.9) & non_ex] = -1e6
        ub[(ub >= 999.9) & non_ex] = +1e6
        mask_rxns = set(cell_rxns[mask_cell])
        for j, rid in enumerate(base["rxns"]):
            if rid not in mask_rxns:
                lb[j] = 0.0; ub[j] = 0.0
        atp = _atp_max_null_uptake(base, lb, ub, idx, env_cond, env_hpf)
        verdict = "PASS" if (np.isfinite(atp) and atp - 31.5 < 0.5) else "FAIL"
        rows.append(dict(cell_A=a, cell_B=b, cofactor=cof, corner=corner,
                         mask_cell=mask_cell, env_cell=env_cell,
                         atp_null_uptake=atp, verdict=verdict))
        print(f"  {a} vs {b}  {cof:9s}  {corner}  mask={mask_cell}  env={env_cell}  "
              f"ATP={atp:.2f}  {verdict}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "energy_gate_mask_check.csv", index=False)
print(f"\n[write] {OUT/'energy_gate_mask_check.csv'}")
fails = df[df.verdict == "FAIL"]
print(f"\n--- VERDICT ---")
print(f"corners checked: {len(df)}")
print(f"PASS: {(df.verdict=='PASS').sum()}")
print(f"FAIL: {(df.verdict=='FAIL').sum()}")
if len(fails):
    print("\nFailed corners (phantom-route risk):")
    print(fails[['cell_A','cell_B','cofactor','corner','mask_cell','env_cell','atp_null_uptake']].to_string(index=False))
