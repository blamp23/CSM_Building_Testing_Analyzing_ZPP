"""
pFBA + single-reaction knockout essentiality for the condition-specific
ATP-sweep failures.

For each (hpf, substrate, failing_cell) where base is usable + 2/3
conditions usable + 1 condition infeasible:
  1. pFBA in the BASE at that hpf with substrate forced to lb=ub=-1 and
     O2 capped at the stoichiometric demand. Phase 1: max ATP = v*.
     Phase 2: min sum(|v_j|) subject to v[ATP] = v*. Sparse route = {j:
     |v_j| > 1e-6}.
  2. Candidate missing = sparse_route ∩ (base.rxns - cell.rxns).
  3. For each candidate, re-solve ATP in the base with that reaction
     forced to v = 0 (SKO). If max ATP < 0.5, it's essential for this
     substrate in the base.
  4. Classify essentials by GPR + zFPKM in the failing cell.

Writes: atp_sweep_pfba_essentials.csv
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from s250_atp_sweep import (
    load_model, class1_vitamin_exchanges, prep_bounds,
    V7_BASE_DIR, V7_EXTRACTED, ATP_RXN, O2_EX,
)

OUT = REPO / "grace_v4/results/atp_sweep_v7"
CONDS = ("BL", "D", "LD")
TOL = 1e-6
MIN_ATP = 0.5

# -------- helpers ------------------------------------------------------------
def parse_gpr_genes(gpr):
    toks = re.findall(r"[A-Za-z0-9_.\-]+", gpr)
    return [t for t in toks if t.lower() not in ("and", "or", "not")]


def pfba_route(model, lb, ub, idx):
    """Phase 1: max ATP. Phase 2: min sum|v| s.t. ATP = ATP*. Returns (ATP*, flux, env, S)."""
    S = model["S"]
    n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()

    # Phase 1
    gm1 = gp.Model("max_atp", env=env)
    v = gm1.addMVar(n_rxn, lb=lb, ub=ub)
    gm1.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm1.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm1.optimize()
    if gm1.Status != GRB.OPTIMAL:
        gm1.dispose(); env.dispose()
        return np.nan, None
    atp_star = gm1.ObjVal
    gm1.dispose()

    # Phase 2: split v = vp - vn, both >= 0, min sum(vp + vn)
    # with lb <= v <= ub and v[ATP] = atp_star
    gm2 = gp.Model("pfba", env=env)
    vp = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(ub > 0, ub, 0.0))
    vn = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(lb < 0, -lb, 0.0))
    # v = vp - vn
    # lb <= vp - vn <= ub
    I = sp.eye(n_rxn, format="csr")
    # Equality: S (vp - vn) = 0  <==>  [S, -S][vp;vn] = 0
    SS = sp.hstack([S, -S], format="csr")
    gm2.addMConstr(SS, gp.MVar.fromlist(list(vp.tolist()) + list(vn.tolist())), "=", np.zeros(n_met))
    # Lower bound on v
    gm2.addMConstr(sp.hstack([I, -I], format="csr"),
                   gp.MVar.fromlist(list(vp.tolist()) + list(vn.tolist())),
                   ">=", lb)
    gm2.addMConstr(sp.hstack([I, -I], format="csr"),
                   gp.MVar.fromlist(list(vp.tolist()) + list(vn.tolist())),
                   "<=", ub)
    # ATP equality
    atp_row = np.zeros(2 * n_rxn)
    atp_row[idx[ATP_RXN]] = 1.0
    atp_row[n_rxn + idx[ATP_RXN]] = -1.0
    gm2.addConstr(gp.quicksum(float(atp_row[k]) * (vp[k] if k < n_rxn else vn[k - n_rxn])
                              for k in range(2 * n_rxn) if atp_row[k] != 0) == atp_star)
    # Min sum
    one = np.ones(2 * n_rxn)
    gm2.setMObjective(None, one, 0.0, sense=GRB.MINIMIZE)
    gm2.optimize()
    if gm2.Status != GRB.OPTIMAL:
        gm2.dispose(); env.dispose()
        return atp_star, None
    flux = np.array([vp[j].X - vn[j].X for j in range(n_rxn)])
    gm2.dispose(); env.dispose()
    return atp_star, flux


def max_atp_with_sko(model, lb, ub, idx, sko_j):
    """Max ATP with sko_j forced to v=0."""
    S = model["S"]
    n_met, n_rxn = S.shape
    lb2 = lb.copy(); ub2 = ub.copy()
    lb2[sko_j] = 0.0; ub2[sko_j] = 0.0
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("sko", env=env)
    v = gm.addMVar(n_rxn, lb=lb2, ub=ub2)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
    gm.dispose(); env.dispose()
    return val


# -------- load data ----------------------------------------------------------
sweep = pd.read_csv(OUT / "atp_sweep_v7.csv")
panel = pd.read_csv(OUT / "atp_sweep_substrates.csv")
diffs = pd.read_csv(OUT / "atp_sweep_condition_differences.csv")

zfp = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv")
zthr = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv")
thr_map = dict(zip(zthr.iloc[:, 0], zthr["exp_threshold_zFPKM"]))
zfp = zfp.set_index(zfp.columns[0])

gmap = pd.read_csv(REPO / "results/qc/gene_map_ensdarg_to_model.csv")
SYM_TO_ENS = (gmap.groupby("model_gene")["ensembl_gene_id"]
              .apply(lambda s: sorted(set(s.dropna()))).to_dict())

class1 = class1_vitamin_exchanges()

base_cache = {}
def get_base(h):
    if h not in base_cache:
        base_cache[h] = load_model(V7_BASE_DIR / f"baked_ocr_{h}_v7.mat", "m")
    return base_cache[h]

cell_cache = {}
def get_cell(cell):
    if cell not in cell_cache:
        cell_cache[cell] = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    return cell_cache[cell]

# -------- iterate ------------------------------------------------------------
rows = []
for _, row in diffs.iterrows():
    h = int(row.hpf); sub_label = row.substrate
    flags = {c: bool(row[c]) for c in CONDS}
    if sum(flags.values()) != 2: continue
    base_row = sweep[(sweep.model == f"base_{h}") & (sweep.substrate == sub_label)]
    if base_row.empty or not bool(base_row.iloc[0].usable): continue
    failing_cond = next(c for c, v in flags.items() if not v)
    cell = f"{failing_cond}_{h}"
    sub_meta = panel[panel.label == sub_label].iloc[0]
    sub_mar = sub_meta.mar; o2_cap = -float(sub_meta.o2_demand)
    print(f"[probe] hpf={h} sub={sub_label} failing={cell}", flush=True)

    base = get_base(h)
    bidx = {r: i for i, r in enumerate(base["rxns"])}
    lb, ub = prep_bounds(base, class1)
    if sub_mar not in bidx or O2_EX not in bidx: continue
    lb[bidx[sub_mar]] = -1.0; ub[bidx[sub_mar]] = -1.0
    lb[bidx[O2_EX]]   = o2_cap; ub[bidx[O2_EX]] = 0.0

    atp_star, flux = pfba_route(base, lb, ub, bidx)
    if flux is None:
        print(f"  [warn] pFBA failed, skip"); continue
    sparse_route = [base["rxns"][j] for j in range(len(flux)) if abs(flux[j]) > TOL]
    print(f"  pFBA ATP*={atp_star:.2f}  sparse_route={len(sparse_route)}", flush=True)

    cm = get_cell(cell)
    cell_rxns = set(cm["rxns"])
    candidates = [r for r in sparse_route if r not in cell_rxns]
    print(f"  candidates (sparse ∩ missing-in-cell): {len(candidates)}", flush=True)

    # pFBA sparse-route membership IS the "on the catabolic route" filter;
    # single-rxn SKO is too strict given base redundancy (never drops ATP
    # below 0.5 for any single rxn, because alternative routes exist).
    # The cell lacks the ENSEMBLE of missing rxns, not any single one.
    # So we classify every candidate and report its pFBA flux magnitude
    # (bigger |v| in base => more load it was carrying in the parsimonious
    # path). Also still probe single-SKO atp_after_sko for reference.
    cell_thr = thr_map.get(cell, np.nan)
    for r in candidates:
        j = bidx[r]
        atp_sko = max_atp_with_sko(base, lb, ub, bidx, j)
        gpr_raw = base["grRules"][j]
        gpr = "" if (gpr_raw is None or (isinstance(gpr_raw, float) and np.isnan(gpr_raw))) else str(gpr_raw).strip()
        subs = base["subs"][j] if j < len(base["subs"]) else ""
        v_pfba = float(flux[j])
        if not gpr:
            cls = "extraction_artefact_no_gpr"
            gene_zfpkms = ""; gene_states = ""
        else:
            symbols = parse_gpr_genes(gpr)
            vals = []; states = []
            for sym in symbols:
                ens_list = SYM_TO_ENS.get(sym, [])
                if not ens_list:
                    vals.append(f"{sym}=NOMAP"); states.append("nomap"); continue
                zs = [float(zfp.loc[e, cell]) for e in ens_list
                      if e in zfp.index and cell in zfp.columns]
                if not zs:
                    vals.append(f"{sym}=NA"); states.append("na"); continue
                zmax = max(zs)
                vals.append(f"{sym}={zmax:.2f}")
                states.append("exp" if zmax >= cell_thr else "noexp")
            gene_zfpkms = ";".join(vals); gene_states = ";".join(states)
            if any(s == "exp" for s in states):
                cls = "extraction_artefact_despite_expression"
            elif any(s == "noexp" for s in states):
                cls = "transcriptional"
            else:
                cls = "no_gene_data"
        rows.append(dict(hpf=h, substrate=sub_label, failing_cell=cell,
                         base_atp=atp_star, n_sparse=len(sparse_route),
                         n_candidates=len(candidates),
                         rxn=r, pfba_flux=v_pfba, atp_after_sko=atp_sko,
                         gpr=gpr, subsys=subs,
                         gene_zfpkms=gene_zfpkms, gene_states=gene_states,
                         classification=cls))
    print(f"  classified: {len(candidates)}", flush=True)

out = pd.DataFrame(rows)
out.to_csv(OUT / "atp_sweep_pfba_essentials.csv", index=False)
print(f"\n[write] atp_sweep_pfba_essentials.csv rows={len(out)}", flush=True)
if len(out):
    print(out.classification.value_counts())
    print("\n=== per-probe essentials ===")
    agg = (out.groupby(["hpf", "substrate", "failing_cell", "classification"])
             .size().unstack(fill_value=0).reset_index())
    agg.to_csv(OUT / "atp_sweep_pfba_essentials_summary.csv", index=False)
    print(agg.to_string())
