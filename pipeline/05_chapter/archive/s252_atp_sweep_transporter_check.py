"""
Gene-less-transporter check for condition-specific ATP-sweep failures.

For every (hpf, substrate) where the substrate is usable in the base and in
>=2 of 3 conditions but not the third, run FBA in the base with that
substrate forced (same policy as s250) and identify reactions carrying
flux. From that catabolic route, list reactions present in the base but
missing from the failing extracted model. For each missing reaction,
record its GPR and the per-gene zFPKM in the failing cell. Classify:
  - no GPR                           -> extraction_artefact (gap)
  - GPR, all genes expressed         -> extraction_artefact (model drop)
  - GPR, >=1 gene not expressed      -> transcriptional
Writes atp_sweep_transporter_check.csv
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from s250_atp_sweep import (
    load_model, class1_vitamin_exchanges, prep_bounds,
    V7_BASE_DIR, V7_EXTRACTED, ATP_RXN, O2_EX, BIOMASS_RXN, MEDIUM,
    o2_demand_for_complete_oxidation,
)

OUT = REPO / "grace_v4/results/atp_sweep_v7"
CONDS = ("BL", "D", "LD")

# -------- helpers ------------------------------------------------------------
GENE_RE = re.compile(r"ENS[A-Z]*\w+|[A-Za-z0-9_-]+")

def parse_gpr_genes(gpr):
    if not gpr: return []
    # strip operators + parens
    toks = re.findall(r"[A-Za-z0-9_.\-]+", gpr)
    return [t for t in toks if t.lower() not in ("and","or","not")]

def fba_solution(model, lb, ub, idx):
    import scipy.sparse as sp
    S = model["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("atp", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    if gm.Status != GRB.OPTIMAL:
        gm.dispose(); env.dispose()
        return np.nan, None
    val = gm.ObjVal
    flux = np.array([v[j].X for j in range(n_rxn)])
    gm.dispose(); env.dispose()
    return val, flux

def active_route(model, lb, ub, idx, tol=1e-6):
    val, flux = fba_solution(model, lb, ub, idx)
    if flux is None: return val, set()
    active = {model["rxns"][j] for j in range(len(flux)) if abs(flux[j]) > tol}
    return val, active

# -------- load data ----------------------------------------------------------
sweep = pd.read_csv(OUT / "atp_sweep_v7.csv")
panel = pd.read_csv(OUT / "atp_sweep_substrates.csv")
diffs = pd.read_csv(OUT / "atp_sweep_condition_differences.csv")

zfp = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv")
zthr = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv")
# Thresholds: index by cell -> exp_threshold_zFPKM
thr_map = dict(zip(zthr.iloc[:,0], zthr["exp_threshold_zFPKM"]))
# zFPKM: expects first col gene (Ensembl ID), rest cells
gene_col = zfp.columns[0]
zfp = zfp.set_index(gene_col)

# GPR tokens use gene symbols; zFPKM uses Ensembl IDs. Build mapping.
gmap = pd.read_csv(REPO / "results/qc/gene_map_ensdarg_to_model.csv")
# model_gene (symbol) -> list of ensembl IDs
SYM_TO_ENS = (gmap.groupby("model_gene")["ensembl_gene_id"]
              .apply(lambda s: sorted(set(s.dropna()))).to_dict())

class1 = class1_vitamin_exchanges()

# -------- cache base models --------------------------------------------------
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

# -------- iterate failing (hpf, substrate, failing_cond) --------------------
rows = []
for _, row in diffs.iterrows():
    h = int(row.hpf); sub_label = row.substrate
    # Pick the failing condition: the lone False among {BL,D,LD}
    flags = {c: bool(row[c]) for c in CONDS}
    n_true = sum(flags.values())
    if n_true != 2:                 # want 2/3 usable, 1 missing
        continue
    # Also require base is usable at this hpf for this substrate (filter)
    base_row = sweep[(sweep.model == f"base_{h}") & (sweep.substrate == sub_label)]
    if base_row.empty or not bool(base_row.iloc[0].usable):
        continue
    failing_cond = next(c for c, v in flags.items() if not v)
    cell = f"{failing_cond}_{h}"
    print(f"[probe] hpf={h} sub={sub_label} failing={cell}", flush=True)

    sub_meta = panel[panel.label == sub_label].iloc[0]
    sub_mar = sub_meta.mar
    o2_cap = -float(sub_meta.o2_demand)

    # Base FBA
    base = get_base(h)
    bidx = {r:i for i,r in enumerate(base["rxns"])}
    lb, ub = prep_bounds(base, class1)
    if sub_mar not in bidx or O2_EX not in bidx:
        continue
    lb[bidx[sub_mar]] = -1.0; ub[bidx[sub_mar]] = -1.0
    lb[bidx[O2_EX]] = o2_cap;  ub[bidx[O2_EX]] = 0.0
    base_atp, base_route = active_route(base, lb, ub, bidx)
    if not base_route:
        print(f"  [warn] base route empty, skip")
        continue

    # Failing cell: reactions present
    cm = get_cell(cell)
    cell_rxns = set(cm["rxns"])
    missing = sorted(base_route - cell_rxns)
    if not missing:
        # catabolic reactions all present; maybe something subtler (bound diff)
        # record a single row noting this
        rows.append(dict(hpf=h, substrate=sub_label, failing_cell=cell,
                         base_atp=base_atp, n_route=len(base_route),
                         n_missing=0, rxn="(none)", gpr="", subsys="",
                         classification="no_missing_rxns",
                         note="all route rxns present in cell"))
        continue

    # For each missing reaction, record GPR + zFPKM in failing cell
    cell_thr = thr_map.get(cell, np.nan)
    for mr in missing:
        j = bidx[mr]
        gpr_raw = base["grRules"][j]
        gpr = "" if (gpr_raw is None or (isinstance(gpr_raw, float) and np.isnan(gpr_raw))) else str(gpr_raw).strip()
        subs = base["subs"][j] if j < len(base["subs"]) else ""
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
                # take the max zFPKM across isoform ensembls for this symbol
                zs = [float(zfp.loc[e, cell]) for e in ens_list if e in zfp.index and cell in zfp.columns]
                if not zs:
                    vals.append(f"{sym}=NA"); states.append("na"); continue
                zmax = max(zs)
                vals.append(f"{sym}={zmax:.2f}")
                states.append("exp" if zmax >= cell_thr else "noexp")
            gene_zfpkms = ";".join(vals)
            gene_states = ";".join(states)
            # Classification: GPR is often "A or B" (any gene suffices)
            # or "A and B" (both needed). Human-GEM mostly uses OR for isoforms
            # of transporters. Simplification: expressed if ANY gene is "exp".
            if any(s == "exp" for s in states):
                cls = "extraction_artefact_despite_expression"
            elif any(s == "noexp" for s in states):
                cls = "transcriptional"
            else:
                cls = "no_gene_data"
        rows.append(dict(hpf=h, substrate=sub_label, failing_cell=cell,
                         base_atp=base_atp, n_route=len(base_route),
                         n_missing=len(missing), rxn=mr,
                         gpr=gpr, subsys=subs,
                         gene_zfpkms=gene_zfpkms, gene_states=gene_states,
                         classification=cls))

out = pd.DataFrame(rows)
out.to_csv(OUT / "atp_sweep_transporter_check.csv", index=False)
print(f"\n[write] atp_sweep_transporter_check.csv rows={len(out)}", flush=True)
print(out.classification.value_counts())

# Condensed per-probe summary + narrow focus on transport reactions
transport = out[out.subsys.str.contains("Transport", na=False, case=False)].copy()
transport.to_csv(OUT / "atp_sweep_transporter_check_transport_only.csv", index=False)
print(f"[write] atp_sweep_transporter_check_transport_only.csv rows={len(transport)}", flush=True)

summary = (out.groupby(["hpf","substrate","failing_cell","classification"])
             .size().unstack(fill_value=0).reset_index())
summary.to_csv(OUT / "atp_sweep_transporter_check_summary.csv", index=False)
print(f"[write] atp_sweep_transporter_check_summary.csv", flush=True)
print(summary.to_string())
