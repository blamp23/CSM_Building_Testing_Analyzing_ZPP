"""
Diagnose the 0.75 anaerobic glucose yield at BL_24/D_24/LD_24 and D_48/LD_48.

Expected: 2.0 ATP (glucose → 2 lactate via EM glycolysis + LDH).
Observed in those cells: 0.75.

Protocol: in BL_24, run the clean-sweep anaerobic probe (glucose -1, O2=0,
class-1 closed), then greedy add-back from the base's pFBA-missing-in-cell
to find the single reaction (or set) that restores ATP to 2.0.
"""
import sys, re
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from s250_atp_sweep import (
    load_model, class1_vitamin_exchanges, prep_bounds,
    V7_BASE_DIR, V7_EXTRACTED, ATP_RXN, O2_EX,
)

OUT = REPO / "grace_v4/results/atp_sweep_v7"
OCR = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}

def prep_bounds_clean(model, class1_set, probed_mar):
    lb, ub = prep_bounds(model, class1_set)
    for rid in class1_set:
        if rid == probed_mar: continue
        if rid in model["rxns"]:
            j = model["rxns"].index(rid)
            lb[j] = 0.0; ub[j] = 1000.0
    return lb, ub

def parse_gpr_genes(gpr):
    toks = re.findall(r"[A-Za-z0-9_.\-]+", gpr)
    return [t for t in toks if t.lower() not in ("and","or","not")]

zfp = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv").set_index("ensembl_gene_id")
zthr = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv")
thr_map = dict(zip(zthr.iloc[:,0], zthr["exp_threshold_zFPKM"]))
gmap = pd.read_csv(REPO / "results/qc/gene_map_ensdarg_to_model.csv")
SYM_TO_ENS = (gmap.groupby("model_gene")["ensembl_gene_id"]
              .apply(lambda s: sorted(set(s.dropna()))).to_dict())
class1 = class1_vitamin_exchanges()
GLU_EX = "MAR09034"

PROBES = [("BL_24", 24), ("D_24", 24), ("LD_24", 24), ("D_48", 48), ("LD_48", 48)]

rows_summary, rows_detail = [], []
for cell, h in PROBES:
    print(f"\n[probe] {cell} glucose anaerobic", flush=True)
    base = load_model(V7_BASE_DIR / f"baked_ocr_{h}_v7.mat", "m")
    cm = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    cell_rxns = set(cm["rxns"])
    lb0, ub0 = prep_bounds_clean(base, class1, GLU_EX)
    idx = {r:i for i,r in enumerate(base["rxns"])}
    lb0[idx[GLU_EX]] = -1.0; ub0[idx[GLU_EX]] = -1.0
    lb0[idx[O2_EX]] = 0.0;    ub0[idx[O2_EX]] = 0.0

    # Simulated cell anaerobic: SKO base-only rxns
    missing_all = [r for r in base["rxns"] if r not in cell_rxns]
    lb_cur, ub_cur = lb0.copy(), ub0.copy()
    for r in missing_all:
        j = idx[r]; lb_cur[j] = 0.0; ub_cur[j] = 0.0

    S = base["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag",0); env.start()
    gm = gp.Model("fba", env=env); vv = gm.addMVar(n_rxn, lb=lb_cur, ub=ub_cur)
    gm.addMConstr(S.tocsr(), vv, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm.optimize()
    atp0 = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
    print(f"  cell-sim anaerobic ATP = {atp0:.4f} (target 2.0)", flush=True)

    # Base pFBA to size the catabolic path
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    for r in missing_all:
        j = idx[r]; vv[j].lb = lb0[idx[r]]; vv[j].ub = ub0[idx[r]]
    gm.optimize()
    atp_base = gm.ObjVal
    gm.dispose()
    # pFBA
    gm2 = gp.Model("pfba", env=env)
    vp = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(ub0>0, ub0, 0.0))
    vn = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(lb0<0, -lb0, 0.0))
    SS = sp.hstack([S, -S], format="csr")
    gm2.addMConstr(SS, gp.MVar.fromlist(list(vp.tolist())+list(vn.tolist())),
                   "=", np.zeros(n_met))
    gm2.addConstr(vp[idx[ATP_RXN]] - vn[idx[ATP_RXN]] == atp_base)
    gm2.setMObjective(None, np.ones(2*n_rxn), 0.0, sense=GRB.MINIMIZE); gm2.optimize()
    flux = np.array([vp[j].X - vn[j].X for j in range(n_rxn)])
    sparse_set = {base["rxns"][j] for j in range(n_rxn) if abs(flux[j]) > 1e-6}
    gm2.dispose()
    missing = [r for r in missing_all if r in sparse_set]
    print(f"  base anaerobic ATP*={atp_base:.3f}  sparse route={len(sparse_set)}  missing-in-cell∩sparse={len(missing)}", flush=True)

    # Greedy add-back
    native_lb = {r: lb0[idx[r]] for r in missing}
    native_ub = {r: ub0[idx[r]] for r in missing}
    lb_cur, ub_cur = lb0.copy(), ub0.copy()
    for r in missing_all:
        j = idx[r]; lb_cur[j] = 0.0; ub_cur[j] = 0.0
    gm = gp.Model("addback", env=env); vg = gm.addMVar(n_rxn, lb=lb_cur, ub=ub_cur)
    gm.addMConstr(S.tocsr(), vg, "=", np.zeros(n_met))
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm.setParam("Method",1); gm.optimize()
    cur = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
    picks = []
    remaining = set(missing)
    target = 1.9       # ~textbook 2.0
    step = 0
    while cur < target and remaining and step < 20:
        step += 1
        best_rxn, best_atp = None, cur
        for r in remaining:
            j = idx[r]
            ol, ou = vg[j].lb, vg[j].ub
            vg[j].lb = native_lb[r]; vg[j].ub = native_ub[r]
            gm.optimize()
            if gm.Status == GRB.OPTIMAL and gm.ObjVal > best_atp + 1e-9:
                best_atp = gm.ObjVal; best_rxn = r
            vg[j].lb = ol; vg[j].ub = ou
        if best_rxn is None:
            print(f"  step {step}: no single gain available"); break
        j = idx[best_rxn]
        vg[j].lb = native_lb[best_rxn]; vg[j].ub = native_ub[best_rxn]
        gm.optimize(); cur = gm.ObjVal
        gpr_raw = base["grRules"][j]
        gpr = str(gpr_raw).replace("\x00","").strip() if gpr_raw is not None else ""
        subs = base["subs"][j] if j < len(base["subs"]) else ""
        syms = parse_gpr_genes(gpr) if gpr else []
        thr = thr_map.get(cell, 0)
        zinfo = []
        for s in syms:
            ens = SYM_TO_ENS.get(s, [])
            zs = [float(zfp.loc[e,cell]) for e in ens if e in zfp.index and cell in zfp.columns]
            if zs: zinfo.append(f"{s}={max(zs):.2f}({'exp' if max(zs)>=thr else 'noexp'})")
            else: zinfo.append(f"{s}=NA")
        picks.append(dict(cell=cell, step=step, rxn=best_rxn,
                          atp_after=cur, subsystem=subs, gpr=gpr,
                          zfpkm_check=";".join(zinfo)))
        print(f"  step {step}: +{best_rxn}  ATP→{cur:.3f}  subs={subs[:40]}  gpr={gpr!r}  {';'.join(zinfo)}",
              flush=True)
    rows_summary.append(dict(cell=cell, baseline_atp=atp0, final_atp=cur,
                             n_steps=len(picks), base_atp=atp_base,
                             restoring_set=";".join(p["rxn"] for p in picks)))
    rows_detail.extend(picks)
    gm.dispose(); env.dispose()

pd.DataFrame(rows_detail).to_csv(OUT / "anaerobic_glucose_ldh_detail.csv", index=False)
pd.DataFrame(rows_summary).to_csv(OUT / "anaerobic_glucose_ldh_summary.csv", index=False)
print(f"\n[write] anaerobic_glucose_ldh_summary.csv rows={len(rows_summary)}", flush=True)
print(pd.DataFrame(rows_summary).to_string())
