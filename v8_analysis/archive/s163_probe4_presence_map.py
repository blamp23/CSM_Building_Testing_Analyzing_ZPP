"""
s163 -- Probe 4: FVA presence map for the stress-response subsystems.

Target subsystems:
  - "Drug metabolism"
  - "Xenobiotics metabolism"
  - "Glutathione metabolism"
  - "ROS detoxification"
  - "Metabolism of xenobiotics by cytochrome P450"   (synonym variants kept)

For each of the 15 extracted cells, under `apply_context_bounds` (Frame A
unconstrained -- no biomass floor), run FVA (min, max) on every target
reaction and classify:
  - BLOCKED     if max - min < 1e-6 and |max| < 1e-6
  - FORCED      if min > 1e-6 or max < -1e-6 (same sign at both ends)
  - OPTIONAL    otherwise (can carry, not required)

Also records whether the reaction has a non-empty GPR.

Read-only on model files. Writes to v7_susceptibility/ (never v7_final/).

Outputs:
  probe4_presence.csv            rxn, subsystem, cell, class, min, max, has_gpr
  probe4_presence_matrix.csv     rxn × 15 cells {blocked|forced|optional|absent}
  probe4_summary.md              per-subsystem counts per cell
"""
import sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds, compute_biomax

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

TARGET_SUBS = {
    "Drug metabolism",
    "Xenobiotics metabolism",
    "Glutathione metabolism",
    "ROS detoxification",
    "Metabolism of xenobiotics by cytochrome P450",
}

FVA_TOL = 1e-6


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


def _load(fp, root_key):
    with h5py.File(fp, "r") as f:
        m = f[root_key]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        try:
            rxn_names = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
        except Exception:
            rxn_names = [""] * len(rxns)
        try:
            subs = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception:
            subs = [""] * len(rxns)
        try:
            grrules = [_resolve(f, r) for r in np.asarray(m["grRules"]).ravel()]
        except Exception:
            grrules = [""] * len(rxns)
        S_grp = m["S"]
        data = np.asarray(S_grp["data"])
        ir = np.asarray(S_grp["ir"])
        jc = np.asarray(S_grp["jc"])
        try:
            n_met_hdf = int(np.asarray(m["mets"]).shape[0])
        except Exception:
            n_met_hdf = 0
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, n_met_hdf)
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, rxn_names=rxn_names, subs=subs, grrules=grrules,
                S=S, lb=lb, ub=ub)


def _fva(gm, v, idx_rxns, biomass_j, biomax, floor_frac=0.9):
    """FVA with biomass floored at floor_frac × biomax. Returns dict rxn_j -> (min, max)."""
    # Pin biomass lb to floor
    v[biomass_j].lb = floor_frac * biomax
    results = {}
    for j in idx_rxns:
        # min
        gm.setObjective(v[j], GRB.MINIMIZE)
        gm.optimize()
        v_min = v[j].X if gm.Status == GRB.OPTIMAL else np.nan
        # max
        gm.setObjective(v[j], GRB.MAXIMIZE)
        gm.optimize()
        v_max = v[j].X if gm.Status == GRB.OPTIMAL else np.nan
        results[j] = (v_min, v_max)
    return results


def _classify(vmin, vmax):
    if not (np.isfinite(vmin) and np.isfinite(vmax)):
        return "infeasible"
    if abs(vmax - vmin) < FVA_TOL and abs(vmax) < FVA_TOL:
        return "blocked"
    if vmin > FVA_TOL or vmax < -FVA_TOL:
        return "forced"
    return "optional"


rows = []
per_cell_rxns = {}  # rxn -> {cell: class}

for cell in CELLS:
    cond, h = cell.split("_")
    hpf = int(h)
    print(f"\n[cell] {cell}", flush=True)
    mdl = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)

    # §11.19-gated bounds
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex, cond, hpf,
        S=mdl["S"],
    )
    biomax, _ = compute_biomax(mdl, lb, ub, idx)
    print(f"  biomax = {biomax:.3f}", flush=True)

    # Build Gurobi model once per cell, reused for all FVA LPs
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    n_met, n_rxn = mdl["S"].shape
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(mdl["S"].tocsr(), v, "=", np.zeros(n_met))
    # pin biomass floor via MVar lb update
    bio_j = idx[BIOMASS_RXN]
    new_lb = lb.copy()
    new_lb[bio_j] = 0.9 * biomax
    v.setAttr("LB", new_lb.tolist())
    gm.update()

    # Target reactions in this model
    target_js = [j for j, s in enumerate(mdl["subs"]) if s in TARGET_SUBS]
    print(f"  target rxns in model: {len(target_js)}", flush=True)

    for j in target_js:
        obj = np.zeros(n_rxn); obj[j] = 1.0
        gm.setMObjective(None, obj, 0.0, sense=GRB.MINIMIZE)
        gm.optimize()
        vmin = v.X[j] if gm.Status == GRB.OPTIMAL else np.nan
        gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
        gm.optimize()
        vmax = v.X[j] if gm.Status == GRB.OPTIMAL else np.nan
        cls = _classify(vmin, vmax)
        rid = mdl["rxns"][j]
        rows.append(dict(
            cell=cell, condition=cond, hpf=hpf,
            rxn=rid,
            rxn_name=mdl["rxn_names"][j] if j < len(mdl["rxn_names"]) else "",
            subsystem=mdl["subs"][j] if j < len(mdl["subs"]) else "",
            grrule=mdl["grrules"][j] if j < len(mdl["grrules"]) else "",
            has_gpr=bool((mdl["grrules"][j] if j < len(mdl["grrules"]) else "").strip()),
            fva_min=vmin, fva_max=vmax, cls=cls,
        ))
        per_cell_rxns.setdefault(rid, {})[cell] = cls

    gm.dispose(); env.dispose()

df = pd.DataFrame(rows)
df.to_csv(OUT / "probe4_presence.csv", index=False)
print(f"\n[write] {OUT/'probe4_presence.csv'}  ({len(df)} rows)", flush=True)

# Reaction x cell matrix
all_rxns = sorted(per_cell_rxns.keys())
mat = pd.DataFrame(index=all_rxns, columns=CELLS)
for rid in all_rxns:
    for c in CELLS:
        mat.loc[rid, c] = per_cell_rxns[rid].get(c, "absent")
# Attach reaction metadata (take first non-null across cells)
meta = df.drop_duplicates("rxn").set_index("rxn")[["rxn_name", "subsystem", "has_gpr"]]
mat = meta.join(mat, how="right")
mat.to_csv(OUT / "probe4_presence_matrix.csv")
print(f"[write] {OUT/'probe4_presence_matrix.csv'}", flush=True)

# Summary
md = ["# Probe 4 presence map — stress-response subsystems",
      "",
      "FVA at Frame A (biomass ≥ 0.9 × biomass_max per cell) under final v7",
      "bounds (incl. §11.19 demand bounds). Read-only on model files.",
      "",
      "Classes: `blocked` (fva_min=fva_max=0), `forced` (same-sign range),",
      "`optional` (crosses 0), `absent` (reaction not in extracted model).",
      "",
      "## Per-subsystem counts per cell",
      ""]
for sub in sorted(TARGET_SUBS):
    sub_rxns = df[df.subsystem == sub].rxn.unique()
    if len(sub_rxns) == 0:
        continue
    md.append(f"### {sub} ({len(sub_rxns)} reactions in union)")
    md.append("")
    md.append("| cell | forced | optional | blocked | absent |")
    md.append("|---|---|---|---|---|")
    for cell in CELLS:
        counts = {"forced": 0, "optional": 0, "blocked": 0, "absent": 0}
        for rid in sub_rxns:
            counts[per_cell_rxns.get(rid, {}).get(cell, "absent")] += 1
        md.append(f"| {cell} | {counts['forced']} | {counts['optional']} | {counts['blocked']} | {counts['absent']} |")
    md.append("")
(OUT / "probe4_summary.md").write_text("\n".join(md))
print(f"[write] {OUT/'probe4_summary.md'}", flush=True)

print("\n--- verdict ---")
tot = len(df)
print(f"total (rxn × cell) rows: {tot}")
for cls in ("forced", "optional", "blocked", "infeasible"):
    n = (df.cls == cls).sum()
    print(f"  {cls:12s}: {n} ({100*n/tot:.1f}%)")
with_gpr = df[df.has_gpr].shape[0]
print(f"  with GPR: {with_gpr}/{tot} ({100*with_gpr/tot:.1f}%)")
