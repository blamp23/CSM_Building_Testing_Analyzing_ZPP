"""
C4: Loops-allowed FVA on LD_72, D_72, BL_96, LD_48 for every reaction in
7 subsystems, under OLD vs NEW bounds.

Subsystems (as named in Human-GEM):
  - Glycolysis / Gluconeogenesis
  - Pentose phosphate pathway
  - Pentose and glucuronate interconversions
  - Inositol phosphate metabolism
  - Vitamin A metabolism
  - Beta oxidation of (saturated/unsaturated/odd-chain/branched) fatty acids
  - Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism

Per reaction: old [min, max], new [min, max], fractional change on each.
Per subsystem per cell: count and fraction of reactions whose range changes
> 20% (either endpoint).
"""
import sys, re
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
import importlib
import apply_context_bounds as acb; importlib.reload(acb)
from s250_atp_sweep import load_model, V7_EXTRACTED

OUT = REPO / "grace_v4/results/atp_sweep_v7"
BIOMASS = "MAR00021"
FREE = ["MAR09152","MAR09154","MAR09276","MAR09361"]

SUBSYSTEMS = {
    "Glycolysis": r"^Glycolysis",
    "Pentose phosphate pathway": r"^Pentose phosphate",
    "Pentose/glucuronate interconversions": r"^Pentose and glucuronate",
    "Inositol phosphate metabolism": r"^Inositol phosphate",
    "Vitamin A / retinol metabolism": r"[Vv]itamin A|[Rr]etinol",
    "Fatty acid oxidation": r"[Bb]eta oxidation|^Fatty acid",
    "TCA": r"[Tt]ricarboxylic|TCA",
}
import os
CELLS_TO_RUN = os.environ.get("CELLS", "LD_72,D_72,BL_96,LD_48").split(",")


def apply_old(model, cond, hpf):
    is_ex = (np.diff(model["S"].indptr) == 1)
    lb, ub, idx, _ = acb.apply_context_bounds(model["lb"].copy(), model["ub"].copy(),
                                              model["rxns"], is_ex, cond, hpf,
                                              skip_energy_gate=True)
    for rid in FREE:
        if rid in idx: lb[idx[rid]] = -1000.0
    return lb, ub, idx

def apply_new(model, cond, hpf):
    is_ex = (np.diff(model["S"].indptr) == 1)
    lb, ub, idx, _ = acb.apply_context_bounds(model["lb"].copy(), model["ub"].copy(),
                                              model["rxns"], is_ex, cond, hpf,
                                              skip_energy_gate=True)
    return lb, ub, idx


def fva_range(S_csr, lb, ub, n_rxn, n_met, j, env, gm):
    """FVA on reaction j: returns (vmin, vmax). Reuses env, builds new model
    with warm-start context."""
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S_csr, v, "=", np.zeros(n_met))
    c_mn = np.zeros(n_rxn); c_mn[j] = 1.0
    gm.setMObjective(None, c_mn, 0.0, sense=GRB.MINIMIZE); gm.optimize()
    vmin = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.setMObjective(None, c_mn, 0.0, sense=GRB.MAXIMIZE); gm.optimize()
    vmax = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    return vmin, vmax


def run_fva(model, lb, ub, idx, rxns_to_fva):
    """Returns dict rxn -> (vmin, vmax)."""
    S = model["S"]; n_met, n_rxn = S.shape
    S_csr = S.tocsr()
    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.setParam("Method", 1)      # dual simplex (best warm restart)
    env.setParam("Threads", 1)     # single-threaded per LP — better for a
                                   # sequence of thousands of small LPs
                                   # (overhead of multi-threading dominates)
    env.setParam("Presolve", 1)
    env.start()
    out = {}
    # One persistent model, change objective
    gm = gp.Model("fva", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S_csr, v, "=", np.zeros(n_met))
    for r in rxns_to_fva:
        if r not in idx: continue
        j = idx[r]
        c = np.zeros(n_rxn); c[j] = 1.0
        gm.setMObjective(None, c, 0.0, sense=GRB.MINIMIZE); gm.optimize()
        vmin = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
        gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm.optimize()
        vmax = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
        out[r] = (vmin, vmax)
    gm.dispose(); env.dispose()
    return out


all_rows = []
summary = []
for cell in CELLS_TO_RUN:
    cond, h = cell.split("_"); h = int(h)
    print(f"\n[cell] {cell}", flush=True)
    m = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    idx = {r:i for i,r in enumerate(m["rxns"])}
    subs = m["subs"]
    # Build subsystem → rxn list from the model
    sub2rxn = {}
    for name, pat in SUBSYSTEMS.items():
        rxns = [m["rxns"][j] for j in range(len(m["rxns"])) if re.search(pat, subs[j], re.IGNORECASE)]
        sub2rxn[name] = rxns
        print(f"  {name:40s}  n_rxns={len(rxns)}", flush=True)
    all_fva_rxns = sorted(set().union(*sub2rxn.values()))
    print(f"  total unique rxns across 7 subsystems: {len(all_fva_rxns)}", flush=True)

    lb_o, ub_o, _ = apply_old(m, cond, h)
    lb_n, ub_n, _ = apply_new(m, cond, h)
    print(f"  running OLD FVA...", flush=True)
    out_o = run_fva(m, lb_o, ub_o, idx, all_fva_rxns)
    print(f"  running NEW FVA...", flush=True)
    out_n = run_fva(m, lb_n, ub_n, idx, all_fva_rxns)

    # Per-rxn report
    for sub_name, rxns in sub2rxn.items():
        n_total = 0; n_changed = 0
        for r in rxns:
            if r not in idx: continue
            n_total += 1
            vmo, vxo = out_o.get(r, (np.nan, np.nan))
            vmn, vxn = out_n.get(r, (np.nan, np.nan))
            # fractional change: |new - old| / max(|old|, 1e-6)
            def fc(a, b):
                if not (np.isfinite(a) and np.isfinite(b)): return np.nan
                scale = max(abs(a), 1e-6)
                return (b - a) / scale
            fc_min = fc(vmo, vmn)
            fc_max = fc(vxo, vxn)
            changed = (abs(fc_min) > 0.2 and abs(vmo-vmn) > 1e-4) or \
                      (abs(fc_max) > 0.2 and abs(vxo-vxn) > 1e-4)
            if changed: n_changed += 1
            all_rows.append(dict(
                cell=cell, subsystem=sub_name, rxn=r,
                old_min=vmo, old_max=vxo, new_min=vmn, new_max=vxn,
                fc_min=fc_min, fc_max=fc_max, changed=changed))
        frac = n_changed / n_total if n_total else 0
        summary.append(dict(cell=cell, subsystem=sub_name,
                            n_total=n_total, n_changed=n_changed, fraction_changed=frac))
        print(f"    {sub_name:40s}  n_rxns={n_total:3d}  >20%_changed={n_changed:3d} ({frac*100:4.1f}%)", flush=True)

per_rxn = pd.DataFrame(all_rows)
suffix = "_" + "_".join(CELLS_TO_RUN) if len(CELLS_TO_RUN) < 4 else ""
per_rxn.to_csv(OUT / f"fva_C4_per_rxn{suffix}.csv", index=False)
per_sub = pd.DataFrame(summary)
per_sub.to_csv(OUT / f"fva_C4_per_subsystem_summary{suffix}.csv", index=False)
print(f"\n[write] fva_C4_per_rxn.csv rows={len(per_rxn)}", flush=True)
print(f"[write] fva_C4_per_subsystem_summary.csv rows={len(per_sub)}", flush=True)
print("\n=== Per-subsystem per-cell: n_changed (>20% in min or max) ===")
pv = per_sub.pivot(index="subsystem", columns="cell", values="n_changed")
print(pv.to_string())
print("\n=== Fraction changed ===")
pv2 = per_sub.pivot(index="subsystem", columns="cell", values="fraction_changed")
print(pv2.to_string())
