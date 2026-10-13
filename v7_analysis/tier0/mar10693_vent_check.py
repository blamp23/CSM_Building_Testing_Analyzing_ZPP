"""
s160 -- Tier 0 addendum: close the MAR10693 (Trp-Gly-Leu) secretion stub.

For each of 15 child models (extracted v7), FBA biomass_max with
MAR10693 ub = 0; compare to the baseline biomass_max from
tier0_parent_child_v7.csv.

For LD_24 only: list producer reactions of Trp-Gly-Leu (any rxn with it
as a product) + the transporter formula (both intracellular & boundary
MAR10693). Flag whether H+ or Na+ appears in the transporter's
stoichiometry.

Output: results/v7_final_addenda/mar10693.md (iff needed) + a stdout
summary table.
"""
import sys
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
T0_PC = REPO / "grace_v4/results/tier0_v7/tier0_parent_child_v7.csv"
OUT_DIR = REPO / "grace_v4/results/v7_final_addenda"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
STUB = "MAR10693"
CONDS = ("BL","D","LD")
HPFS = (24,48,72,96,120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]


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
        try: rxnNames = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
        except Exception: rxnNames = [""]*len(rxns)
        try: subs = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception: subs = [""]*len(rxns)
        try: metNames = [_resolve(f, r) for r in np.asarray(m["metNames"]).ravel()]
        except Exception: metNames = []
        try: mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception: mets = []
        try: metComps = np.asarray(m["metComps"]).ravel().astype(int)
        except Exception: metComps = np.zeros(len(mets), dtype=int)
        try: compNames = [_resolve(f, r) for r in np.asarray(m["compNames"]).ravel()]
        except Exception: compNames = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        try: n_met_hdf = int(np.asarray(m["mets"]).shape[0])
        except Exception: n_met_hdf = 0
        n_met = max(int(ir.max()) + 1, n_met_hdf, len(mets), len(metNames))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc)-1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, rxnNames=rxnNames, subs=subs, mets=mets, metNames=metNames,
                metComps=metComps, compNames=compNames, S=S, lb=lb, ub=ub)


def fba(mdl, lb, ub, idx):
    S = mdl["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c_obj = np.zeros(n_rxn); c_obj[idx[BIOMASS_RXN]] = 1.0
    gm.setMObjective(None, c_obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    bm = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return bm


# -------- baselines from Tier 0 parent-child file -------
pc = pd.read_csv(T0_PC).set_index("cell")
baselines = pc["child_biomax"].to_dict()
print(f"[baseline] loaded {len(baselines)} child biomass_max values from Tier 0 parent-child", flush=True)

# -------- 15 LPs: close MAR10693.ub -------
rows = []
for cell in CELLS:
    cond, hpf = cell.split("_"); hpf = int(hpf)
    ch = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    is_ex = (np.diff(ch["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(ch["lb"].copy(), ch["ub"].copy(),
                                           ch["rxns"], is_ex, cond, hpf)
    base = baselines[cell]
    # close MAR10693 secretion
    if STUB in idx:
        j = idx[STUB]
        ub[j] = 0.0
        bm_new = fba(ch, lb, ub, idx)
    else:
        bm_new = np.nan
    delta = bm_new - base
    pct = (delta / base) * 100 if (base and np.isfinite(base) and base != 0) else np.nan
    rows.append(dict(cell=cell, hpf=hpf, condition=cond,
                     biomax_baseline=base, biomax_closed=bm_new,
                     delta=delta, pct_change=pct))
    print(f"  {cell}: baseline={base:.3f}  closed={bm_new:.3f}  Δ={delta:+.3f}  ({pct:+.2f}%)",
          flush=True)

df = pd.DataFrame(rows)
max_abs_pct = df["pct_change"].abs().max()
print(f"\nmax |% change| across cells: {max_abs_pct:.2f}%")
threshold_met = (df["pct_change"].abs() >= 5.0).any()

# -------- LD_24 producer + transporter inspection -------
print("\n--- LD_24 producers of Trp-Gly-Leu + transporter MAR10693 ---", flush=True)
ld24 = _load(V7_EXTRACTED / "trans_rfastcormics_LD_24.mat", "parent")
rxn_names = ld24["rxnNames"]; subs = ld24["subs"]
mets = ld24["mets"]; metNames = ld24["metNames"]
metComps = ld24["metComps"]; compNames = ld24["compNames"]
def comp_tag(mi):
    if mi >= len(metComps): return ""
    c = int(metComps[mi]) - 1
    return compNames[c] if 0 <= c < len(compNames) else ""

# locate TrpGlyLeu metabolite(s)
target_subs = {i for i, nm in enumerate(metNames)
               if "trpglyleu" in nm.lower().replace("-","").replace(" ","")
               or "tryptophanyl-glycyl-leucine" in nm.lower()}
print(f"metabolite(s) matching Trp-Gly-Leu: {sorted(target_subs)[:10]}")
for mi in sorted(target_subs):
    print(f"  met {mi}: id={mets[mi] if mi<len(mets) else '?'}  name={metNames[mi]}  comp={comp_tag(mi)}")

# find reactions producing or transporting Trp-Gly-Leu
S_csc = ld24["S"].tocsc()
producer_info = []
for j in range(S_csc.shape[1]):
    col = S_csc[:, j]
    hits = [(int(r), float(v)) for r, v in zip(col.indices, col.data) if r in target_subs]
    if not hits: continue
    rxn_id = ld24["rxns"][j]; rxn_name = rxn_names[j] if j<len(rxn_names) else ""
    sub = subs[j] if j<len(subs) else ""
    # express as formula
    def formula_of(j):
        col = S_csc[:, j]
        subs_ = []; prods_ = []
        for r_, v_ in zip(col.indices, col.data):
            nm = metNames[r_] if r_<len(metNames) else mets[r_] if r_<len(mets) else f"met{r_}"
            cmp = comp_tag(r_)
            tag = f"{nm}[{cmp}]"
            (subs_ if v_<0 else prods_).append((abs(float(v_)), tag))
        s_str = " + ".join([f"{(c if c!=1 else '')}{t}" for c,t in subs_]) or "∅"
        p_str = " + ".join([f"{(c if c!=1 else '')}{t}" for c,t in prods_]) or "∅"
        rev = "<=>" if ld24['lb'][j]<0 else "->"
        return f"{s_str} {rev} {p_str}"
    producer_info.append(dict(rxn=rxn_id, name=rxn_name, subsystem=sub,
                              formula=formula_of(j)))

for p in producer_info:
    print(f"  {p['rxn']}  {p['name']!r}  [sub={p['subsystem']!r}]")
    print(f"    {p['formula']}")

# transporter formula + H+/Na+ scan
transporter_text = ""
for p in producer_info:
    f = p["formula"].lower()
    if "h+" in f or "h2o" not in f and "sodium" in f:
        pass
    if p["rxn"] == STUB or "transport" in p.get("subsystem","").lower():
        transporter_text = p["formula"]

h_present = any("h+[" in p["formula"].lower() for p in producer_info)
na_present = any("na+[" in p["formula"].lower() or "sodium" in p["formula"].lower() for p in producer_info)
print(f"\nH+ in any producer/transporter formula? {h_present}")
print(f"Na+ or 'sodium' in any producer/transporter formula? {na_present}")

# -------- output decision -------
df.to_csv(OUT_DIR/"mar10693_biomax_table.csv", index=False)
print(f"\n[write] {OUT_DIR/'mar10693_biomax_table.csv'}")

if threshold_met:
    md = []
    md.append("# MAR10693 (Trp-Gly-Leu secretion stub) — Tier 0 addendum")
    md.append("")
    md.append(f"**Generated 2026-10-13.** No resample; FBA only on frozen v7 extracted models.")
    md.append("")
    md.append("## Rule applied")
    md.append("")
    md.append("`MAR10693.ub := 0` (close the gene-less tripeptide secretion stub). All other bounds via `apply_context_bounds`. FBA for biomass_max.")
    md.append("")
    md.append("## Per-cell result")
    md.append("")
    md.append("| cell | hpf | baseline biomass_max | closed biomass_max | Δ | % change |")
    md.append("|---|---|---|---|---|---|")
    for _, r in df.iterrows():
        md.append(f"| {r.cell} | {r.hpf} | {r.biomax_baseline:.3f} | {r.biomax_closed:.3f} | {r.delta:+.3f} | {r.pct_change:+.2f}% |")
    md.append("")
    md.append(f"**max |% change|:** {max_abs_pct:.2f}%")
    md.append("")
    md.append("## LD_24: producers of Trp-Gly-Leu + transporter formula")
    md.append("")
    for p in producer_info:
        md.append(f"- **{p['rxn']}** {p['name']!r} [sub={p['subsystem']!r}]")
        md.append(f"  - `{p['formula']}`")
    md.append("")
    md.append(f"- H⁺ present in any formula: **{h_present}**")
    md.append(f"- Na⁺ / sodium present in any formula: **{na_present}**")
    (OUT_DIR/"mar10693.md").write_text("\n".join(md))
    print(f"[write] {OUT_DIR/'mar10693.md'}")
else:
    print(f"\n|%change| < 5% in every cell → MAR10693 is a passive vent; one footnote in §Tier 0 only, no mar10693.md written.")
