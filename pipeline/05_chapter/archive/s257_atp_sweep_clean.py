"""
Clean ATP-sweep rerun.

Changes vs s250:
  (1) Class-1 vitamin exchanges (α-tocotrienol, γ-tocotrienol, β-carotene,
      inositol) are closed for uptake (lb=0) during each probe UNLESS they
      are the probed substrate. The s255/s256 diagnosis showed that the
      2 107-ATP base_120 glucose result was almost entirely driven by
      β-carotene (C40H56) and inositol (C6H12O6) feeding central carbon
      metabolism at lb=-1000; closing them gives textbook glucose = 31.50.
  (2) O2 is capped at the per-hpf physiological OCR (24→97, 48→173,
      72→227, 96→273, 120→313) rather than at the stoichiometric
      -(C + H/4 - O/2) per substrate. With class-1 closed the O2 cap is
      no longer "hiding" an inflation; using the physiological OCR means
      the heatmap reflects the network's capacity under in-vivo oxygen,
      not an artefact of the cap choice.

Outputs: atp_sweep_v7_clean.csv + the summary tables, written alongside
the original files so both can be compared.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from s250_atp_sweep import (
    load_model, class1_vitamin_exchanges, prep_bounds,
    V7_BASE_DIR, V7_EXTRACTED, ATP_RXN, O2_EX, BIOMASS_RXN,
    HPFS, CELLS, CONDS,
)

OUT = REPO / "grace_v4/results/atp_sweep_v7"
OCR = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}


def prep_bounds_clean(model, class1_set, probed_mar):
    """Variant of prep_bounds that closes class-1 vitamins EXCEPT the probed one."""
    lb, ub = prep_bounds(model, class1_set)
    for rid in class1_set:
        if rid == probed_mar: continue
        if rid in model["rxns"]:
            j = model["rxns"].index(rid)
            lb[j] = 0.0; ub[j] = 1000.0
    return lb, ub


def max_atp(model, lb, ub, idx):
    S = model["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val


panel = pd.read_csv(OUT / "atp_sweep_substrates.csv")
class1 = class1_vitamin_exchanges()

models_to_run = []
for h in HPFS:
    models_to_run.append((f"base_{h}", h, V7_BASE_DIR / f"baked_ocr_{h}_v7.mat", "m"))
for cell in CELLS:
    _, h = cell.split("_"); h = int(h)
    models_to_run.append((cell, h, V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent"))

rows = []
for name, h, fp, root in models_to_run:
    print(f"[model] {name}", flush=True)
    model = load_model(fp, root)
    idx = {r: i for i, r in enumerate(model["rxns"])}
    if ATP_RXN not in idx:
        print(f"  [warn] ATP rxn missing, skip"); continue
    o2_cap = -float(OCR[h])
    for _, sub in panel.iterrows():
        if sub.mar not in idx:
            rows.append(dict(model=name, hpf=h, substrate=sub.label,
                             mar=sub.mar, carbon_count=int(sub.carbon),
                             atp_aerobic=np.nan, atp_anaerobic=np.nan,
                             note="substrate_mar_not_in_model"))
            continue
        lb0, ub0 = prep_bounds_clean(model, class1, sub.mar)
        sj = idx[sub.mar]
        # Aerobic
        lb, ub = lb0.copy(), ub0.copy()
        lb[sj] = -1.0; ub[sj] = -1.0
        if O2_EX in idx:
            lb[idx[O2_EX]] = o2_cap; ub[idx[O2_EX]] = 0.0
        atp_aer = max_atp(model, lb, ub, idx)
        # Anaerobic
        lb, ub = lb0.copy(), ub0.copy()
        lb[sj] = -1.0; ub[sj] = -1.0
        if O2_EX in idx:
            lb[idx[O2_EX]] = 0.0; ub[idx[O2_EX]] = 0.0
        atp_ana = max_atp(model, lb, ub, idx)
        rows.append(dict(model=name, hpf=h, substrate=sub.label,
                         mar=sub.mar, carbon_count=int(sub.carbon),
                         atp_aerobic=atp_aer, atp_anaerobic=atp_ana, note=""))

df = pd.DataFrame(rows)
df["atp_per_C_aerobic"] = df["atp_aerobic"] / df["carbon_count"]
df["usable"] = df["atp_aerobic"] > 0.5
df.to_csv(OUT / "atp_sweep_v7_clean.csv", index=False)
print(f"[write] atp_sweep_v7_clean.csv rows={len(df)}", flush=True)

# ---- Positive controls ------------------------------------------------------
print("\n=== Glucose control ===", flush=True)
for _, r in df[df.substrate == "glucose"].iterrows():
    print(f"  {r.model:12s} aer={r.atp_aerobic:.3f}  ana={r.atp_anaerobic:.3f}")

print("\n=== Palmitate aerobic ===", flush=True)
for _, r in df[df.substrate == "palmitate"].iterrows():
    print(f"  {r.model:12s} aer={r.atp_aerobic:.3f}")

# ---- Summaries --------------------------------------------------------------
cnt = df.groupby("model")["usable"].sum().rename("n_usable").reset_index()
cnt.to_csv(OUT / "atp_sweep_v7_clean_usable_count.csv", index=False)
print(f"\n[write] atp_sweep_v7_clean_usable_count.csv")

base_usable = {h: set(df[(df.model == f"base_{h}") & df.usable].substrate) for h in HPFS}
bmc_rows = []
for cell in CELLS:
    cond, h = cell.split("_"); h = int(h)
    cu = set(df[(df.model == cell) & df.usable].substrate)
    for s in sorted(base_usable[h] - cu):
        bmc_rows.append(dict(cell=cell, condition=cond, hpf=h, substrate=s, direction="base_only"))
    for s in sorted(cu - base_usable[h]):
        bmc_rows.append(dict(cell=cell, condition=cond, hpf=h, substrate=s, direction="cell_only"))
pd.DataFrame(bmc_rows).to_csv(OUT / "atp_sweep_v7_clean_base_minus_cell.csv", index=False)

cond_rows = []
for h in HPFS:
    sb = {c: set(df[(df.model == f"{c}_{h}") & df.usable].substrate) for c in CONDS}
    for s in sorted(set().union(*sb.values())):
        flags = {c: s in sb[c] for c in CONDS}
        if len(set(flags.values())) > 1:
            cond_rows.append(dict(hpf=h, substrate=s, **{c: flags[c] for c in CONDS}))
pd.DataFrame(cond_rows).to_csv(OUT / "atp_sweep_v7_clean_condition_differences.csv", index=False)
print(f"[write] atp_sweep_v7_clean_condition_differences.csv rows={len(cond_rows)}")
