"""
ATP-yield sweep over single carbon sources — all 15 v7 extracted models + the
5 v7 base per-hpf models. Read-only capability probe (no changes to models
or ledger).

Per-model protocol:
  Start from the model's stored bounds, then re-set exchange bounds:
    - Every ORGANIC uptake (metabolite formula contains C) : lb = 0, ub = 1000
    - Every INORGANIC exchange (incl. water, protons, salts): lb = -1000, ub = 1000
    - Class-1 vitamins (ledger action=medium_uptake)        : lb = -1000, ub = 1000
    - O2 (MAR09048): aerobic  lb=-1000, anaerobic lb = 0
    - Biomass (MAR00021): lb = 0, ub = 1000 (no floor)
  For each probe substrate: that exchange is OVERRIDDEN to lb = ub = -1.
  Objective: maximise MAR03964 (ATP maintenance / ATP generation).

Outputs under grace_v4/results/atp_sweep_v7/:
  atp_sweep_substrates.csv
  atp_sweep_v7.csv
  atp_sweep_summary_usable_count.csv
  atp_sweep_summary_base_minus_cell.csv
  atp_sweep_condition_differences.csv
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import V7_MODELS as V7_EXTRACTED, BASE_MODELS as V7_BASE_DIR
from config import LEDGER, ATP_SWEEP as OUT, CONDS, HPFS, CELLS
OUT.mkdir(parents=True, exist_ok=True)
CELLS = list(CELLS)

BIOMASS_RXN = "MAR00021"
ATP_RXN     = "MAR03964"
O2_EX       = "MAR09048"

# 17 bulk-inorganic exchanges (water, protons, bulk salts, O2, CO2). Same set
# used by functions/apply_context_bounds.MEDIUM. Opened to lb=-1000.
MEDIUM = frozenset({
    "MAR09047","MAR09058","MAR09072","MAR09073","MAR09074","MAR09076",
    "MAR09077","MAR09078","MAR09079","MAR09080","MAR09081","MAR09082",
    "MAR09148","MAR09150","MAR13066","MAR13072","MAR13073",
})

# Plausible-yolk-nutrient substrate panel. Keys are label; regex is anchored on metName.
PATTERNS = {
    "glucose":     r"^glucose$|^D-glucose$",
    "fructose":    r"^fructose$|^D-fructose$",
    "galactose":   r"^galactose$|^D-galactose$",
    "mannose":     r"^mannose$|^D-mannose$",
    "glycerol":    r"^glycerol$",
    "L-lactate":   r"^L-lactate$|^lactate$|^\(S\)-lactate$",
    "pyruvate":    r"^pyruvate$",
    "acetate":     r"^acetate$",
    "beta-hydroxybutyrate": r"beta.?hydroxybutyrate|\(R\)-3-hydroxybutanoate|R-3-hydroxybutyrate",
    "acetoacetate": r"^acetoacetate$",
    "palmitate":   r"^palmitate$|^hexadecanoate$",
    "stearate":    r"^stearate$|^octadecanoate$",
    "oleate":      r"^oleate$|^cis-9-octadecenoate$",
    "linoleate":   r"^linoleate$",
    "linolenate":  r"^linolenate$|^alpha-linolenate$",
    "DHA":         r"^DHA$|docosahexaenoate",
    "EPA":         r"^EPA$|eicosapentaenoate",
    "choline":     r"^choline$",
    "inositol":    r"^inositol$|^myo-inositol$",
    "adenosine":   r"^adenosine$",
    "guanosine":   r"^guanosine$",
    "cytidine":    r"^cytidine$",
    "uridine":     r"^uridine$",
    "thymidine":   r"^thymidine$",
}
for aa in ["alanine","arginine","asparagine","aspartate","cysteine","glutamate",
           "glutamine","glycine","histidine","isoleucine","leucine","lysine",
           "methionine","phenylalanine","proline","serine","threonine",
           "tryptophan","tyrosine","valine"]:
    PATTERNS[f"L-{aa}"] = rf"^L-{aa}$|^{aa}$"


# ---------- MATLAB v7.3 loader (minimal fields we need) ---------------------
def _resolve(f, x):
    if isinstance(x, h5py.Reference):
        return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)

def _str_list(f, m, key):
    try:
        arr = np.asarray(m[key]).ravel()
        return [_resolve(f, r) for r in arr]
    except Exception:
        return []

def load_model(fp, root_key):
    with h5py.File(fp, "r") as f:
        m = f[root_key]
        rxns     = _str_list(f, m, "rxns")
        rxnNames = _str_list(f, m, "rxnNames")
        mets     = _str_list(f, m, "mets")
        metNames = _str_list(f, m, "metNames")
        metFormulas = _str_list(f, m, "metFormulas")
        try: grRules = _str_list(f, m, "grRules")
        except Exception: grRules = [""] * len(rxns)
        try: subs = _str_list(f, m, "subSystems")
        except Exception: subs = [""] * len(rxns)
        try: comps = _str_list(f, m, "comps")
        except Exception: comps = []
        try: metComps = np.asarray(m["metComps"]).ravel().astype(int)
        except Exception: metComps = np.zeros(len(mets), dtype=int)
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        n_met_hdf = int(np.asarray(m["mets"]).shape[0])
        n_met = max(int(ir.max()) + 1, n_met_hdf, len(mets), len(metNames))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc)-1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, rxnNames=rxnNames, mets=mets, metNames=metNames,
                metFormulas=metFormulas, metComps=metComps, comps=comps,
                S=S, lb=lb, ub=ub, grRules=grRules, subs=subs)


# ---------- substrate panel discovery (base model only) ---------------------
def count_carbon(formula):
    if not formula: return 0
    m = re.search(r"C(\d*)(?![a-z])", formula)
    if not m: return 0
    return int(m.group(1)) if m.group(1) else 1

def parse_formula(formula):
    """Return (C, H, O, N, S, charge-ignored) for a Hill-form formula."""
    if not formula: return {}
    atoms = {}
    for sym, num in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if not sym: continue
        atoms[sym] = atoms.get(sym, 0) + (int(num) if num else 1)
    return atoms

def o2_demand_for_complete_oxidation(formula):
    """
    O2 needed to fully oxidize CcHhNnOoSs -> cCO2 + (h/2)H2O + (s)SO4 + (n/2)N2.
    Standard formula: O2 = C + H/4 - O/2 (ignoring N, S for simplicity;
    adequate for the yolk-nutrient panel which is C/H/O-dominated).
    """
    a = parse_formula(formula)
    C = a.get("C", 0); H = a.get("H", 0); O = a.get("O", 0)
    return max(0.0, C + H/4.0 - O/2.0)

def exchange_table(model):
    """Return DataFrame of boundary exchanges: (mar, mam, met_name, formula, carbon, ex_idx)."""
    S = model["S"]
    is_ex = (np.diff(S.indptr) == 1)
    rows = []
    for j in np.where(is_ex)[0]:
        col = S.getcol(j).tocoo()
        met_i = int(col.row[0])
        rows.append(dict(mar=model["rxns"][j], mam=model["mets"][met_i],
                         met_name=model["metNames"][met_i],
                         formula=model["metFormulas"][met_i],
                         carbon=count_carbon(model["metFormulas"][met_i]),
                         ex_idx=j))
    return pd.DataFrame(rows)

def build_substrate_panel(base_model):
    exdf = exchange_table(base_model)
    rows = []
    for label, pat in PATTERNS.items():
        hits = exdf[exdf.met_name.str.match(pat, case=False, na=False)]
        if not len(hits):
            print(f"  [warn] {label}: no match in base model")
            continue
        if len(hits) > 1:
            # prefer the one that's a direct boundary of the e-compartment metabolite
            print(f"  [warn] {label}: {len(hits)} matches, taking first: {hits.iloc[0].mar}")
        r = hits.iloc[0]
        rows.append(dict(label=label, mar=r.mar, mam=r.mam,
                         met_name=r.met_name, formula=r.formula,
                         carbon=int(r.carbon),
                         o2_demand=o2_demand_for_complete_oxidation(r.formula)))
    return pd.DataFrame(rows)


# ---------- ledger: class-1 vitamins (medium_uptake) -----------------------
def class1_vitamin_exchanges():
    led = pd.read_csv(LEDGER)
    return set(led.loc[led.action == "medium_uptake", "exchange_MAR_id"].dropna())


# ---------- LP: ATP max with substrate constraints -------------------------
def prep_bounds(model, class1_set):
    """
    Return (lb, ub) with the universal probe policy applied.
      - organic uptakes closed (lb=0)
      - inorganic exchanges opened for uptake (lb=-1000)
      - class-1 vitamin exchanges opened for uptake (lb=-1000)
      - biomass unconstrained (lb=0, ub=1000)
      - secretion uncapped (ub=1000) on every exchange
    O2 is handled by the caller (aerobic/anaerobic).
    Substrate override (lb=ub=-1) applied by caller.
    """
    S = model["S"]
    is_ex = (np.diff(S.indptr) == 1)
    lb = model["lb"].copy()
    ub = model["ub"].copy()
    rxns = model["rxns"]; mets = model["mets"]; metFormulas = model["metFormulas"]

    for j in np.where(is_ex)[0]:
        col = S.getcol(j).tocoo()
        met_i = int(col.row[0])
        rid = rxns[j]
        if rid in MEDIUM:
            # bulk inorganics / water / protons / O2 / CO2
            lb[j] = -1000.0; ub[j] = 1000.0
        elif rid in class1_set:
            # alpha/gamma-tocotrienol, beta-carotene, inositol
            lb[j] = -1000.0; ub[j] = 1000.0
        else:
            # everything else closed for uptake, open for secretion
            lb[j] = 0.0; ub[j] = 1000.0

    if BIOMASS_RXN in rxns:
        j = rxns.index(BIOMASS_RXN)
        lb[j] = 0.0; ub[j] = 1000.0
    return lb, ub


def max_atp(model, lb, ub, idx):
    S = model["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("atp", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c_obj = np.zeros(n_rxn); c_obj[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c_obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    status = gm.Status
    val = gm.ObjVal if status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val, status


# ---------- main ------------------------------------------------------------
def main():
    class1 = class1_vitamin_exchanges()
    print(f"[info] class-1 vitamin exchanges: {sorted(class1)}", flush=True)

    print("[step] loading base_120 for substrate panel", flush=True)
    base_120 = load_model(V7_BASE_DIR / "baked_ocr_120_v7.mat", "m")
    panel = build_substrate_panel(base_120)
    panel.to_csv(OUT / "atp_sweep_substrates.csv", index=False)
    print(f"[write] atp_sweep_substrates.csv n={len(panel)}", flush=True)

    # Also probe each hpf base model + 15 extracted
    models_to_run = []
    for h in HPFS:
        models_to_run.append((f"base_{h}", V7_BASE_DIR / f"baked_ocr_{h}_v7.mat", "m"))
    for cell in CELLS:
        models_to_run.append((cell, V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent"))

    rows = []
    for name, fp, root in models_to_run:
        print(f"[model] {name}", flush=True)
        model = load_model(fp, root)
        idx = {r: i for i, r in enumerate(model["rxns"])}

        # ATP rxn present?
        if ATP_RXN not in idx:
            print(f"  [warn] {name}: {ATP_RXN} missing, skipping", flush=True)
            continue

        lb0, ub0 = prep_bounds(model, class1)

        for _, sub in panel.iterrows():
            if sub.mar not in idx:
                rows.append(dict(model=name, substrate=sub.label,
                                 mar=sub.mar, carbon_count=int(sub.carbon),
                                 atp_aerobic=np.nan, atp_anaerobic=np.nan,
                                 note="substrate_mar_not_in_model"))
                continue
            sj = idx[sub.mar]
            # O2 is bounded at the exact complete-oxidation stoichiometric
            # demand (C + H/4 - O/2 per unit substrate) so O2 is the
            # limiter — this gives textbook ATP for glucose (6 O2 -> 31.5)
            # and palmitate (23 O2 -> ~106) without side-oxidation inflation.
            o2_cap = -float(sub.o2_demand)
            # Aerobic
            lb = lb0.copy(); ub = ub0.copy()
            lb[sj] = -1.0; ub[sj] = -1.0
            if O2_EX in idx:
                lb[idx[O2_EX]] = float(o2_cap); ub[idx[O2_EX]] = 0.0
            atp_aer, st_aer = max_atp(model, lb, ub, idx)
            # Anaerobic
            lb = lb0.copy(); ub = ub0.copy()
            lb[sj] = -1.0; ub[sj] = -1.0
            if O2_EX in idx:
                lb[idx[O2_EX]] = 0.0; ub[idx[O2_EX]] = 0.0
            atp_ana, st_ana = max_atp(model, lb, ub, idx)
            note = ""
            if st_aer != GRB.OPTIMAL: note = f"aer_status={st_aer};"
            if st_ana != GRB.OPTIMAL: note += f"ana_status={st_ana};"
            rows.append(dict(model=name, substrate=sub.label,
                             mar=sub.mar, carbon_count=int(sub.carbon),
                             atp_aerobic=atp_aer, atp_anaerobic=atp_ana,
                             note=note))

    df = pd.DataFrame(rows)
    df["atp_per_C_aerobic"] = df["atp_aerobic"] / df["carbon_count"]
    df["usable"] = df["atp_aerobic"] > 0.5
    df.to_csv(OUT / "atp_sweep_v7.csv", index=False)
    print(f"[write] atp_sweep_v7.csv rows={len(df)}", flush=True)

    # ------- Positive control: glucose + palmitate ---------------
    print("\n=== Positive control: glucose should be ≈31.5 aer / ≈2.0 ana ===", flush=True)
    glc = df[df.substrate == "glucose"]
    for _, r in glc.iterrows():
        aer_ok = (not np.isnan(r.atp_aerobic)) and abs(r.atp_aerobic - 31.5) < 1.0
        ana_ok = (not np.isnan(r.atp_anaerobic)) and abs(r.atp_anaerobic - 2.0) < 0.5
        tag = "OK" if (aer_ok and ana_ok) else "FAIL"
        print(f"  [{tag}] {r.model:12s} aer={r.atp_aerobic:.3f} ana={r.atp_anaerobic:.3f}", flush=True)

    print("\n=== Palmitate aerobic ===", flush=True)
    palm = df[df.substrate == "palmitate"]
    for _, r in palm.iterrows():
        print(f"  {r.model:12s} aer={r.atp_aerobic:.3f}  ana={r.atp_anaerobic:.3f}", flush=True)

    # ------- Summary 1: usable count per model -----------------------
    cnt = df.groupby("model").apply(lambda d: int(d.usable.sum())).rename("n_usable")
    cnt_df = cnt.reset_index()
    cnt_df.to_csv(OUT / "atp_sweep_summary_usable_count.csv", index=False)

    # ------- Summary 2: base-usable but cell-infeasible -----------------
    # For each cell_hpf, compare against base_hpf.
    base_usable = {h: set(df[(df.model == f"base_{h}") & df.usable].substrate) for h in HPFS}
    bmc_rows = []
    for cell in CELLS:
        cond, h = cell.split("_"); h = int(h)
        cell_usable = set(df[(df.model == cell) & df.usable].substrate)
        lost = base_usable[h] - cell_usable
        gained = cell_usable - base_usable[h]
        for s in sorted(lost):
            bmc_rows.append(dict(cell=cell, condition=cond, hpf=h,
                                 substrate=s, direction="base_only"))
        for s in sorted(gained):
            bmc_rows.append(dict(cell=cell, condition=cond, hpf=h,
                                 substrate=s, direction="cell_only"))
    pd.DataFrame(bmc_rows).to_csv(OUT / "atp_sweep_summary_base_minus_cell.csv", index=False)

    # ------- Summary 3: LD vs BL vs D at same hpf -----------------------
    cond_rows = []
    for h in HPFS:
        s_by_cond = {c: set(df[(df.model == f"{c}_{h}") & df.usable].substrate) for c in CONDS}
        for s in sorted(set().union(*s_by_cond.values())):
            flags = {c: (s in s_by_cond[c]) for c in CONDS}
            # Record only if NOT identical across all 3
            if len(set(flags.values())) > 1:
                cond_rows.append(dict(hpf=h, substrate=s,
                                      LD=flags["LD"], BL=flags["BL"], D=flags["D"]))
    pd.DataFrame(cond_rows).to_csv(OUT / "atp_sweep_condition_differences.csv", index=False)

    print(f"\n[write] atp_sweep_summary_usable_count.csv", flush=True)
    print(f"[write] atp_sweep_summary_base_minus_cell.csv rows={len(bmc_rows)}", flush=True)
    print(f"[write] atp_sweep_condition_differences.csv rows={len(cond_rows)}", flush=True)


if __name__ == "__main__":
    main()
