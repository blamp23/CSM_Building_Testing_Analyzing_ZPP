"""
s168 -- 2×2 topology × envelope attribution for Probe 1 sink gaps.

Scaffold modelled on s150 (Tier 0 2×2) but applied to a cofactor sink.

For each specified pair (cell_A, cell_B, cofactor):
  Corner (mask_A, env_A) = sink_A at own matched floor      [= Probe 1 result]
  Corner (mask_A, env_B) = cell_A reactions kept × cell_B envelope bounds
  Corner (mask_B, env_A) = cell_B reactions kept × cell_A envelope bounds
  Corner (mask_B, env_B) = sink_B at own matched floor      [= Probe 1 result]

The reaction mask is built from each cell's extracted-model rxns list; applied
on the v7 base model (parent of same hpf). apply_context_bounds supplies the
envelope under (cond, hpf) of whichever envelope corner is selected. For the
cross corners, biomass floor = matched_floor of the pair (same across both
cross corners).

Attribution: for a pair A vs B with |sink_A - sink_B| > 10%,
  topology_effect  = (mask_A, env_B) - (mask_B, env_B)          # same env; mask swap
  envelope_effect  = (mask_A, env_A) - (mask_A, env_B)          # same mask; env swap
Report which effect carries the gap.

Model source: v7 BASE models (ocr_anchored_extraction/models/baked_ocr_{hpf}_v7.mat).
All biomass floors and bounds are reapplied; cell-specific reaction masks
set lb=ub=0 for reactions NOT retained by that cell's extraction.

Internal ±1000 Wang defaults widened to ±1e6 (matches patched s165).

Output: probe1_attribution_2x2.csv
"""
import sys, time, json
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds, compute_biomax

V7_BASES = REPO / "ocr_anchored_extraction/models"
V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
COFACTORS = {
    "NADPH":    "MAM02555c",
    "GSH":      "MAM02026c",
    "UDP_GlcA": "MAM03109c",
    "PAPS":     "MAM02681c",
}

# Triggers: top-6 clear envelope cases + UDP-GlcA/PAPS pairs
TRIGGERS = [
    # (cell_A, cell_B, cofactor)
    ("BL_24",  "LD_24",  "GSH"),
    ("BL_24",  "D_24",   "GSH"),
    ("D_120",  "LD_120", "GSH"),
    ("BL_120", "D_120",  "GSH"),
    ("D_96",   "LD_96",  "NADPH"),
    ("BL_96",  "D_96",   "NADPH"),
    # UDP-GlcA
    ("BL_96",  "D_96",   "UDP_GlcA"),
    ("BL_72",  "D_72",   "UDP_GlcA"),
    ("D_72",   "LD_72",  "UDP_GlcA"),
    ("D_120",  "LD_120", "UDP_GlcA"),
    # PAPS
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
        try: mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception: mets = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, S=S, lb=lb, ub=ub)


def _max_sink(base, lb, ub, idx, met_id, biomass_lb):
    """Augmented LP: sink of met_id with biomass floored at biomass_lb.
    Returns sink_max."""
    S = base["S"].tocsc()
    n_met, n_rxn = S.shape
    if met_id not in base["mets"]:
        return np.nan
    mi = base["mets"].index(met_id)
    extra = sp.csc_matrix(([-1.0], ([mi], [0])), shape=(n_met, 1))
    S_aug = sp.hstack([S, extra]).tocsr()
    sink_j = n_rxn
    lb_aug = np.concatenate([lb, [0.0]])
    ub_aug = np.concatenate([ub, [1000.0]])
    lb_aug[idx[BIOMASS_RXN]] = biomass_lb
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(S_aug.shape[1], lb=lb_aug, ub=ub_aug)
    gm.addMConstr(S_aug, v, "=", np.zeros(S_aug.shape[0]))
    obj = np.zeros(S_aug.shape[1]); obj[sink_j] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


def _corner(base_mdl, mask_cell, env_cell, cof_met, biomass_lb):
    """Build v7 base at env_cell's (cond, hpf), apply mask (block rxns not in
    mask_cell's extraction), widen internal ±1000, then max sink."""
    env_cond, env_hpf = env_cell.split("_"); env_hpf = int(env_hpf)
    base = dict(base_mdl)
    lb0 = base["lb"].copy(); ub0 = base["ub"].copy()
    is_ex = (np.diff(base["S"].tocsc().indptr) == 1)
    # Energy gate must be skipped on the base model (demand bounds were
    # calibrated against the EXTRACTED cells; the base violates the gate
    # until the extraction mask is applied). After masking, the resulting
    # LP is structurally the same as the extracted cell under env_cell's
    # bounds, which does pass.
    lb, ub, idx, _ = apply_context_bounds(lb0, ub0, base["rxns"], is_ex,
                                           env_cond, env_hpf, S=base["S"],
                                           skip_energy_gate=True)
    non_ex = ~is_ex
    lb[(lb <= -999.9) & non_ex] = -1e6
    ub[(ub >= 999.9) & non_ex] = +1e6
    # Mask: block reactions not in mask_cell's extraction
    mask_rxns = set(cell_rxns[mask_cell])
    for j, rid in enumerate(base["rxns"]):
        if rid not in mask_rxns:
            lb[j] = 0.0; ub[j] = 0.0
    return _max_sink(base, lb, ub, idx, cof_met, biomass_lb)


# ---------- pre-load each cell's rxns + base per hpf + matched floors ----------
cell_rxns = {}
for cell in [c for c, _, _ in TRIGGERS] + [c for _, c, _ in TRIGGERS]:
    fp = V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat"
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
    cell_rxns[cell] = rxns

base_by_hpf = {}
for hpf in (24, 48, 72, 96, 120):
    base_by_hpf[hpf] = _load(V7_BASES / f"baked_ocr_{hpf}_v7.mat", "m")
    print(f"  loaded base hpf={hpf}  rxns={len(base_by_hpf[hpf]['rxns'])}", flush=True)

# biomass_max per cell under final v7 bounds (needed for matched floor)
biomax = {}
for cell in set(c for c, _, _ in TRIGGERS) | set(c for _, c, _ in TRIGGERS):
    cond, h = cell.split("_"); hpf = int(h)
    mdl = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(mdl["lb"].copy(), mdl["ub"].copy(),
                                           mdl["rxns"], is_ex, cond, hpf, S=mdl["S"])
    bm, _ = compute_biomax(mdl, lb, ub, idx)
    biomax[cell] = bm
print(f"  biomax: {biomax}", flush=True)

# ---------- 2×2 swap per trigger ----------
rows = []
for a, b, cof_name in TRIGGERS:
    cond_a, hpf_a = a.split("_"); hpf_a = int(hpf_a)
    cond_b, hpf_b = b.split("_"); hpf_b = int(hpf_b)
    assert hpf_a == hpf_b, f"triggers must be same hpf: {a} vs {b}"
    cof_met = COFACTORS[cof_name]
    base_mdl = base_by_hpf[hpf_a]
    floor = 0.9 * min(biomax[a], biomax[b])

    sAA = _corner(base_mdl, mask_cell=a, env_cell=a, cof_met=cof_met, biomass_lb=floor)
    sAB = _corner(base_mdl, mask_cell=a, env_cell=b, cof_met=cof_met, biomass_lb=floor)
    sBA = _corner(base_mdl, mask_cell=b, env_cell=a, cof_met=cof_met, biomass_lb=floor)
    sBB = _corner(base_mdl, mask_cell=b, env_cell=b, cof_met=cof_met, biomass_lb=floor)

    # Attribution:
    #   envelope effect for A's mask: sAA - sAB  (A env vs B env, both A's topology)
    #   topology effect under B env: sAB - sBB   (A mask vs B mask, both B env)
    #   total gap: sAA - sBB
    env_eff   = sAA - sAB
    topo_eff  = sAB - sBB
    total_gap = sAA - sBB
    denom = total_gap if abs(total_gap) > 1e-6 else np.nan
    env_frac  = env_eff / denom if np.isfinite(denom) else np.nan
    topo_frac = topo_eff / denom if np.isfinite(denom) else np.nan
    rows.append(dict(
        cell_A=a, cell_B=b, cofactor=cof_name, hpf=hpf_a, floor=floor,
        sAA_massA_envA=sAA, sAB_massA_envB=sAB,
        sBA_massB_envA=sBA, sBB_massB_envB=sBB,
        total_gap=total_gap, env_effect=env_eff, topo_effect=topo_eff,
        env_frac=env_frac, topo_frac=topo_frac,
    ))
    print(f"  {a:7s} vs {b:7s}  {cof_name:9s}  "
          f"AA={sAA:7.2f}  AB={sAB:7.2f}  BA={sBA:7.2f}  BB={sBB:7.2f}  "
          f"gap={total_gap:+.2f}  env_frac={env_frac:+.2f}  topo_frac={topo_frac:+.2f}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "probe1_attribution_2x2.csv", index=False)
print(f"\n[write] {OUT/'probe1_attribution_2x2.csv'}")
print("\n--- summary ---")
print("For |env_frac| > |topo_frac|: envelope dominates.")
print("For |topo_frac| > |env_frac|: topology dominates.")
print()
df["dominant"] = np.where(df.env_frac.abs() > df.topo_frac.abs(), "envelope", "topology")
print(df.groupby(["cofactor","dominant"]).size().unstack(fill_value=0))
