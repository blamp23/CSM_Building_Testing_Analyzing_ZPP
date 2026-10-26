"""
Greedy add-back for the 9 substrates lost in every extracted model.
Run on LD_72 as a representative cell; report minimal restoring sets with
GPR + per-gene zFPKM. If the restoring set is dominated by gene-less
transport reactions, the methods chapter must say "GPR-less transporter
pruning is a known cost of expression-gated extraction".
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
MIN_ATP = 0.5

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

# Data
sweep = pd.read_csv(OUT / "atp_sweep_v7_clean.csv")
panel = pd.read_csv(OUT / "atp_sweep_substrates.csv")
bmc = pd.read_csv(OUT / "atp_sweep_v7_clean_base_minus_cell.csv")
# universally lost = substrate 'base_only' in all 15 cells
cnt = bmc[bmc.direction=="base_only"].groupby("substrate").cell.count()
universal = sorted(cnt[cnt == 15].index.tolist())
print(f"[info] universally-lost substrates: {universal}", flush=True)

zfp = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv")
zthr = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv")
thr_map = dict(zip(zthr.iloc[:,0], zthr["exp_threshold_zFPKM"]))
zfp = zfp.set_index(zfp.columns[0])
gmap = pd.read_csv(REPO / "results/qc/gene_map_ensdarg_to_model.csv")
SYM_TO_ENS = (gmap.groupby("model_gene")["ensembl_gene_id"]
              .apply(lambda s: sorted(set(s.dropna()))).to_dict())
class1 = class1_vitamin_exchanges()

cell = "LD_72"; h = 72
base = load_model(V7_BASE_DIR / f"baked_ocr_{h}_v7.mat", "m")
cm = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
cell_rxns = set(cm["rxns"])

def classify(gpr, cell, thr):
    if gpr is None or (isinstance(gpr, float) and np.isnan(gpr)):
        return "no_gpr", "", ""
    g = str(gpr).replace("\x00","").strip()
    if not g: return "no_gpr", "", ""
    syms = parse_gpr_genes(g)
    vals, states = [], []
    for s in syms:
        ens = SYM_TO_ENS.get(s, [])
        if not ens: vals.append(f"{s}=NOMAP"); states.append("nomap"); continue
        zs = [float(zfp.loc[e,cell]) for e in ens if e in zfp.index and cell in zfp.columns]
        if not zs: vals.append(f"{s}=NA"); states.append("na"); continue
        zm = max(zs); vals.append(f"{s}={zm:.2f}")
        states.append("exp" if zm >= thr else "noexp")
    exp = [s for s in states if s in ("exp","noexp")]
    if not exp: tag = "no_mapping"
    elif all(s=="exp" for s in exp): tag = "all_expressed"
    elif all(s=="noexp" for s in exp): tag = "none_expressed"
    else: tag = "mixed"
    return tag, ";".join(vals), ";".join(states)

all_detail, summary = [], []
for sub_label in universal:
    sub = panel[panel.label == sub_label].iloc[0]
    sub_mar = sub.mar
    print(f"\n[probe] {sub_label} ({sub_mar}) failing in {cell}", flush=True)
    lb0, ub0 = prep_bounds_clean(base, class1, sub_mar)
    idx = {r:i for i,r in enumerate(base["rxns"])}
    if sub_mar not in idx:
        print(f"  substrate exchange {sub_mar} not in base — skipping"); continue
    lb0[idx[sub_mar]] = -1.0; ub0[idx[sub_mar]] = -1.0
    lb0[idx[O2_EX]]  = -float(OCR[h]); ub0[idx[O2_EX]] = 0.0

    missing_all = [r for r in base["rxns"] if r not in cell_rxns]
    # pFBA sparse route in base under probe
    S = base["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag",0); env.start()
    gm1 = gp.Model("fba", env=env); vv = gm1.addMVar(n_rxn, lb=lb0, ub=ub0)
    gm1.addMConstr(S.tocsr(), vv, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm1.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm1.optimize()
    atp_star = gm1.ObjVal if gm1.Status == GRB.OPTIMAL else 0.0
    if atp_star < MIN_ATP:
        print(f"  base LP ATP={atp_star:.2f} < {MIN_ATP}, base can't use this substrate under clean policy; skipping");
        gm1.dispose(); env.dispose()
        summary.append(dict(substrate=sub_label, cell=cell, base_atp=atp_star,
                            n_pFBA=0, n_candidates=0, final_atp=atp_star,
                            restoring_set="", verdict="base_also_infeasible"))
        continue
    gm1.dispose()
    gm2 = gp.Model("pfba", env=env)
    vp = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(ub0>0, ub0, 0.0))
    vn = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(lb0<0, -lb0, 0.0))
    SS = sp.hstack([S, -S], format="csr")
    gm2.addMConstr(SS, gp.MVar.fromlist(list(vp.tolist())+list(vn.tolist())),
                   "=", np.zeros(n_met))
    gm2.addConstr(vp[idx[ATP_RXN]] - vn[idx[ATP_RXN]] == atp_star)
    gm2.setMObjective(None, np.ones(2*n_rxn), 0.0, sense=GRB.MINIMIZE); gm2.optimize()
    flux = np.array([vp[j].X - vn[j].X for j in range(n_rxn)])
    sparse_set = {base["rxns"][j] for j in range(n_rxn) if abs(flux[j]) > 1e-6}
    gm2.dispose()
    missing = [r for r in missing_all if r in sparse_set]
    print(f"  base ATP*={atp_star:.2f}  pFBA_sparse={len(sparse_set)}  missing-in-cell∩sparse={len(missing)}", flush=True)
    native_lb = {r: lb0[idx[r]] for r in missing}
    native_ub = {r: ub0[idx[r]] for r in missing}
    lb_cur = lb0.copy(); ub_cur = ub0.copy()
    for r in missing_all:
        j = idx[r]; lb_cur[j] = 0.0; ub_cur[j] = 0.0
    gm = gp.Model("addback", env=env); vg = gm.addMVar(n_rxn, lb=lb_cur, ub=ub_cur)
    gm.addMConstr(S.tocsr(), vg, "=", np.zeros(n_met))
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE); gm.setParam("Method",1); gm.optimize()
    atp0 = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
    # Unlock all pFBA-missing at once
    for r in missing:
        j = idx[r]; vg[j].lb = native_lb[r]; vg[j].ub = native_ub[r]
    gm.optimize()
    atp_all = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
    print(f"  SKO baseline={atp0:.3f}  all-pFBA-unlocked={atp_all:.3f}", flush=True)
    # If all-unlock still < 0.5, probe is harder; report best we can
    picks = []
    if atp_all >= MIN_ATP:
        # Trim: pFBA on the all-unlocked LP, keep only active rxns in candidate pool
        lb_now = np.array([vg[j].lb for j in range(n_rxn)])
        ub_now = np.array([vg[j].ub for j in range(n_rxn)])
        gm3 = gp.Model("pfba2", env=env)
        p2 = gm3.addMVar(n_rxn, lb=0.0, ub=np.where(ub_now>0, ub_now, 0.0))
        n2 = gm3.addMVar(n_rxn, lb=0.0, ub=np.where(lb_now<0, -lb_now, 0.0))
        gm3.addMConstr(SS, gp.MVar.fromlist(list(p2.tolist())+list(n2.tolist())),
                       "=", np.zeros(n_met))
        gm3.addConstr(p2[idx[ATP_RXN]] - n2[idx[ATP_RXN]] == atp_all)
        gm3.setMObjective(None, np.ones(2*n_rxn), 0.0, sense=GRB.MINIMIZE); gm3.optimize()
        fv = np.array([p2[j].X - n2[j].X for j in range(n_rxn)])
        active = [r for r in missing if abs(fv[idx[r]]) > 1e-6]
        gm3.dispose()
        print(f"  pFBA on all-unlocked: {len(active)} active restoring rxns", flush=True)
        for r in active:
            j = idx[r]
            gpr = base["grRules"][j]
            subs_str = base["subs"][j] if j < len(base["subs"]) else ""
            tag, zstr, sstr = classify(gpr, cell, thr_map.get(cell, np.nan))
            flux_j = fv[j]
            picks.append(dict(substrate=sub_label, rxn=r, atp_after=atp_all,
                              gpr=str(gpr).replace("\x00","").strip(),
                              subsystem=subs_str, flux=float(flux_j),
                              gene_zfpkms=zstr, gene_states=sstr, rxn_class=tag))
            print(f"    +{r:12s} v={flux_j:+7.3f}  cls={tag:14s}  gpr={gpr!r} subs={subs_str[:35]}", flush=True)
    gm.dispose(); env.dispose()
    tags = [p["rxn_class"] for p in picks]
    if not picks:
        verdict = "no_restoring"; reason = f"all-pFBA-unlocked ATP={atp_all:.2f}"
    elif any(t == "none_expressed" for t in tags):
        verdict = "transcriptional"
        reason = f"{sum(1 for t in tags if t=='none_expressed')} restoring rxn(s) with all genes below threshold"
    elif all(t == "no_gpr" for t in tags):
        verdict = "artefact_gpr_less"
        reason = "ALL restoring reactions have no GPR"
    elif all(t in ("no_gpr","all_expressed") for t in tags):
        verdict = "artefact"
        reason = "mix of gene-less + above-threshold; model dropped a reaction whose genes were available"
    else:
        verdict = "mixed"
        reason = f"tags: {dict((t, tags.count(t)) for t in set(tags))}"
    summary.append(dict(substrate=sub_label, cell=cell, base_atp=atp_star,
                        n_pFBA=len(sparse_set), n_candidates=len(missing),
                        final_atp=atp_all, restoring_set=";".join(p["rxn"] for p in picks),
                        n_rxns_in_set=len(picks),
                        n_no_gpr=sum(1 for t in tags if t=="no_gpr"),
                        n_none_expressed=sum(1 for t in tags if t=="none_expressed"),
                        n_all_expressed=sum(1 for t in tags if t=="all_expressed"),
                        verdict=verdict, reason=reason))
    all_detail.extend(picks)

det = pd.DataFrame(all_detail)
det.to_csv(OUT / "universal_lost_addback_detail.csv", index=False)
pd.DataFrame(summary).to_csv(OUT / "universal_lost_addback_summary.csv", index=False)
print(f"\n[write] universal_lost_addback_detail.csv rows={len(det)}", flush=True)
print(f"[write] universal_lost_addback_summary.csv rows={len(summary)}", flush=True)
print("\n=== Verdicts ===")
print(pd.DataFrame(summary)[["substrate","final_atp","n_rxns_in_set","n_no_gpr",
                              "n_none_expressed","n_all_expressed","verdict","reason"]].to_string())
