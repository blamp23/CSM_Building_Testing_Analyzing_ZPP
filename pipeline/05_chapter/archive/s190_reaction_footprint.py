"""
s190 -- reaction-set footprint analysis.

Part A — stage comparison (per hpf, 3-way across BL/D/LD):
  Shared by all 3 / shared by exactly 2 / unique to 1.
  For every non-shared reaction, record:
    - is_flux_carrying: mean(|v|) > 1e-6 in the pooled samples of
      every cell that has it (bool)
    - gpr_class: GENELESS / UNMAPPABLE / EXPRESSED / UNKNOWN_ZONE / OFF
      per cell that has the reaction, aggregated to one class per rxn
      (worst-case: majority class across holders if they disagree)
  Writes 5 CSVs (one per hpf) + summary.

Part B — regime developmental sequence (per regime, 4 adjacent-stage
2-set comparisons: 24↔48, 48↔72, 72↔96, 96↔120):
  Shared / A-only / B-only. Same flux + GPR classification.

Output:
  grace_v4/results/v7_final_addenda/footprint/
    stage_{hpf}_three_way.csv
    stage_{hpf}_summary.csv
    regime_{BL|D|LD}_adjacent.csv
    footprint_summary.csv
"""
import sys, re, json
from pathlib import Path
import numpy as np, pandas as pd
import h5py

REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
OUT = R / "v7_final_addenda/footprint"
OUT.mkdir(parents=True, exist_ok=True)
V7 = REPO / "reviewer_packet_v4_slim/v7_models"

HPFS = [24, 48, 72, 96, 120]
CONDS = ("BL", "D", "LD")
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

ZFPKM = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv").set_index("ensembl_gene_id")
THRESH = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv").set_index("cell")
_GMAP = pd.read_csv(REPO / "results/qc/gene_map_ensdarg_to_model.csv")
SYM2ENS = {}
for _, r in _GMAP.iterrows():
    SYM2ENS.setdefault(str(r.model_gene).lower(), set()).add(str(r.ensembl_gene_id))


def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)


def _load_cell(cell):
    fp = V7 / f"trans_rfastcormics_{cell}.mat"
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        try: rxn_names = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
        except Exception: rxn_names = [""] * len(rxns)
        try: subs = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception: subs = [""] * len(rxns)
        try: grrules = [_resolve(f, r) for r in np.asarray(m["grRules"]).ravel()]
        except Exception: grrules = [""] * len(rxns)
    pool = pd.read_csv(R / f"pooled_v7/{cell}_pooled.csv").set_index("rxn")
    return dict(rxns=rxns, rxn_names=rxn_names, subs=subs,
                grrules=grrules, pool=pool)


def _parse_gpr(gpr):
    gpr = (gpr or "").strip()
    if not gpr: return []
    toks = re.findall(r"[A-Za-z][A-Za-z0-9._\-]*", gpr)
    toks = [t for t in toks if t.lower() not in ("and", "or")]
    return [(t, SYM2ENS.get(t.lower(), set())) for t in toks]


def _classify_cell(gpr, cell):
    tokens = _parse_gpr(gpr)
    if not tokens: return "GENELESS"
    exp_t  = THRESH.loc[cell, "exp_threshold_zFPKM"] if cell in THRESH.index else 0.0
    inex_t = THRESH.loc[cell, "inex_threshold_zFPKM"] if cell in THRESH.index else -1.0
    n_map = n_exp = n_off = n_unk = 0
    for sym, ens_set in tokens:
        if not ens_set: continue
        n_map += 1
        vals = [ZFPKM.loc[e, cell] for e in ens_set
                if e in ZFPKM.index and cell in ZFPKM.columns and pd.notna(ZFPKM.loc[e, cell])]
        if not vals: continue
        z = max(vals)
        if z >= exp_t: n_exp += 1
        elif z <= inex_t: n_off += 1
        else: n_unk += 1
    if n_map == 0: return "UNMAPPABLE"
    if n_exp > 0: return "EXPRESSED"
    if n_unk > 0 and n_off == 0: return "UNKNOWN_ZONE"
    if n_off > 0 and n_unk == 0: return "OFF"
    return "UNKNOWN_ZONE_mixed"


# Pre-load all 15 cells
print("[load] 15 extracted cells + pooled samples", flush=True)
cells_data = {c: _load_cell(c) for c in CELLS}


# =============== PART A: stage comparison (3-way per hpf) ===============
def _flux_carrying(cell, rxn, tol=1e-6):
    pool = cells_data[cell]["pool"]
    if rxn not in pool.index: return False
    return bool(np.abs(pool.loc[rxn].values).mean() > tol)


stage_summary = []
for hpf in HPFS:
    bl, d, ld = f"BL_{hpf}", f"D_{hpf}", f"LD_{hpf}"
    sets = {c: set(cells_data[c]["rxns"]) for c in (bl, d, ld)}
    all_three = sets[bl] & sets[d] & sets[ld]
    only_two = (
        (sets[bl] & sets[d]) - sets[ld]
    ) | (
        (sets[bl] & sets[ld]) - sets[d]
    ) | (
        (sets[d] & sets[ld]) - sets[bl]
    )
    only_one = (
        sets[bl] - sets[d] - sets[ld]
    ) | (
        sets[d] - sets[bl] - sets[ld]
    ) | (
        sets[ld] - sets[bl] - sets[d]
    )

    rows = []
    for rid in sorted(all_three | only_two | only_one):
        holders = [c for c in (bl, d, ld) if rid in sets[c]]
        tier = ("shared_by_all" if len(holders) == 3
                else "shared_by_two" if len(holders) == 2
                else "unique_to_one")
        if tier == "shared_by_all":
            # No per-rxn flux/GPR detail; just count below
            rows.append(dict(rxn=rid, tier=tier, holders=";".join(holders),
                             flux_carrying=np.nan, gpr_class=""))
            continue
        # For non-shared reactions, collect flux + GPR class per holder
        flux_frac = sum(_flux_carrying(c, rid) for c in holders) / len(holders)
        gpr_classes = []
        rxn_name = ""
        subsystem = ""
        gpr_text = ""
        for c in holders:
            md = cells_data[c]
            try:
                j = md["rxns"].index(rid)
                if not rxn_name: rxn_name = md["rxn_names"][j]
                if not subsystem: subsystem = md["subs"][j]
                if not gpr_text: gpr_text = md["grrules"][j]
                gpr_classes.append(_classify_cell(md["grrules"][j], c))
            except ValueError:
                pass
        # Aggregate GPR class across holders: majority; tie → first
        if gpr_classes:
            # Prefer stricter class in order: GENELESS, UNMAPPABLE, OFF, UNKNOWN_ZONE_mixed, UNKNOWN_ZONE, EXPRESSED
            # (expressed hardest to achieve; use frequency)
            from collections import Counter
            agg = Counter(gpr_classes).most_common(1)[0][0]
        else:
            agg = "n/a"
        rows.append(dict(rxn=rid, tier=tier, holders=";".join(holders),
                         rxn_name=rxn_name, subsystem=subsystem,
                         gpr=gpr_text, flux_carrying=flux_frac, gpr_class=agg))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"stage_{hpf}_three_way.csv", index=False)
    # Summary row
    stage_summary.append(dict(
        hpf=hpf,
        shared_by_all=len(all_three),
        shared_by_two=len(only_two),
        unique_to_one=len(only_one),
        total=len(all_three) + len(only_two) + len(only_one),
    ))
    print(f"  hpf {hpf}:  shared_all={len(all_three)}  shared_two={len(only_two)}  "
          f"unique_one={len(only_one)}", flush=True)

ss = pd.DataFrame(stage_summary)
ss.to_csv(OUT / "stage_summary.csv", index=False)
print(f"\n[write] {OUT}/stage_{{hpf}}_three_way.csv  +  stage_summary.csv")


# =============== PART B: regime adjacent-stage pairs ===============
regime_summary = []
for cond in CONDS:
    rows_all = []
    for h1, h2 in [(24,48), (48,72), (72,96), (96,120)]:
        a, b = f"{cond}_{h1}", f"{cond}_{h2}"
        set_a = set(cells_data[a]["rxns"])
        set_b = set(cells_data[b]["rxns"])
        shared = set_a & set_b
        a_only = set_a - set_b
        b_only = set_b - set_a
        rows_all.append(dict(
            regime=cond, pair=f"{h1}→{h2}",
            shared=len(shared), a_only=len(a_only), b_only=len(b_only),
            total_union=len(set_a | set_b),
        ))
        # Non-shared reactions, per-cell flux + GPR
        for rid in sorted(a_only | b_only):
            holder = a if rid in a_only else b
            md = cells_data[holder]
            j = md["rxns"].index(rid)
            gpr_class = _classify_cell(md["grrules"][j], holder)
            flux = _flux_carrying(holder, rid)
            rows_all.append(dict(
                regime=cond, pair=f"{h1}→{h2}", rxn=rid,
                holder=holder,
                rxn_name=md["rxn_names"][j],
                subsystem=md["subs"][j],
                gpr=md["grrules"][j],
                gpr_class=gpr_class,
                flux_carrying=flux,
                direction="a_only" if rid in a_only else "b_only",
            ))
    df_reg = pd.DataFrame(rows_all)
    df_reg.to_csv(OUT / f"regime_{cond}_adjacent.csv", index=False)
    # summary rows
    for pair_h in ((24,48), (48,72), (72,96), (96,120)):
        row = df_reg[(df_reg.regime==cond) & (df_reg.pair==f"{pair_h[0]}→{pair_h[1]}") & df_reg.get("rxn").isna()]
        if len(row):
            regime_summary.append({**row.iloc[0].to_dict()})

rsum = pd.DataFrame(regime_summary)
rsum.to_csv(OUT / "regime_summary.csv", index=False)
print(f"[write] {OUT}/regime_{{BL|D|LD}}_adjacent.csv + regime_summary.csv")
print()
print("--- Part A: stage three-way counts ---")
print(ss.to_string(index=False))
print()
print("--- Part B: regime adjacent-stage counts ---")
print(rsum.to_string(index=False))
