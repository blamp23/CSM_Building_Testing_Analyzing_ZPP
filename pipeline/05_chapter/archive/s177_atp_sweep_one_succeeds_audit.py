"""
s177 -- Audit the six 'one-succeeds' cases from the ATP sweep.

For each (substrate, hpf, winner_regime): find reactions that touch the
substrate metabolite (any compartment) in the WINNER's extracted model
but are absent from BOTH loser cells. For each such reaction run the
same GPR + zFPKM audit used for the ten single-failure cases (and for
gss / UDP-GlcA / PAPS):
  - GENELESS                       (empty GPR; extraction artefact)
  - TRANSCRIPTIONAL_PARTIAL        (any gene below inex in losers, none above exp)
  - EXTRACTION_UNKNOWN_ZONE        (all mappable genes in [inex, exp])
  - ALGORITHM_UNKNOWN_RETAINED     (some gene above exp in losers but
                                    reaction still pruned)
  - UNMAPPABLE                     (symbol absent from the model-gene map)

Six cases (from atp_sweep_v7_clean.csv one-succeeds):
  cytidine   48  LD
  cytidine   96  LD
  cytidine  120  LD
  thymidine  72  LD
  fructose   96  BL
  L-lysine  120  D

Output: results/v7_final_addenda/
  s177_atp_one_succeeds_audit.csv
  s177_atp_one_succeeds_audit.md
"""
import sys, json, re
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_final_addenda"
OUT.mkdir(parents=True, exist_ok=True)

ZFPKM = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv").set_index("ensembl_gene_id")
THRESH = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv").set_index("cell")
_GMAP = pd.read_csv(REPO / "results/qc/gene_map_ensdarg_to_model.csv")
SYM2ENS = {}
for _, r in _GMAP.iterrows():
    SYM2ENS.setdefault(str(r.model_gene).lower(), set()).add(str(r.ensembl_gene_id))

CASES = [
    # (substrate_name, met_prefix, catabolic_subsystems, hpf, winner, [losers])
    # Scope expanded from substrate-touching only to the relevant catabolic
    # subsystem(s) because 3 of 6 cases had identical first-shell reaction sets.
    ("cytidine",  "MAM01630",
     {"Pyrimidine metabolism", "Pyrimidine catabolism"},  48,  "LD", ["BL", "D"]),
    ("cytidine",  "MAM01630",
     {"Pyrimidine metabolism", "Pyrimidine catabolism"},  96,  "LD", ["BL", "D"]),
    ("cytidine",  "MAM01630",
     {"Pyrimidine metabolism", "Pyrimidine catabolism"}, 120,  "LD", ["BL", "D"]),
    ("thymidine", "MAM02996",
     {"Pyrimidine metabolism", "Pyrimidine catabolism"},  72,  "LD", ["BL", "D"]),
    ("fructose",  "MAM01840",
     {"Fructose and mannose metabolism",
      "Starch and sucrose metabolism",
      "Glycolysis / Gluconeogenesis"},                     96,  "BL", ["D", "LD"]),
    ("L-lysine",  "MAM02426",
     {"Lysine metabolism",
      "Methionine and cysteine metabolism",   # some lysine routes land here
      "Cysteine and methionine metabolism",
      "Tryptophan metabolism"},                           120,  "D",  ["BL", "LD"]),
]


def _resolve(f, x):
    if isinstance(x, h5py.Reference): return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object: return _resolve(f, x[0])
        try: return "".join(chr(int(c)) for c in x)
        except Exception: return ""
    return str(x)


def _load(fp):
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        try: rxn_names = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
        except Exception: rxn_names = [""] * len(rxns)
        try: subs = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception: subs = [""] * len(rxns)
        try: grrules = [_resolve(f, r) for r in np.asarray(m["grRules"]).ravel()]
        except Exception: grrules = [""] * len(rxns)
        try: mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception: mets = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
    return dict(rxns=rxns, rxn_names=rxn_names, subs=subs, grrules=grrules,
                mets=mets, S=S)


def _reactions_touching(mdl, met_prefix):
    """Reaction IDs that touch any metabolite with id starting with met_prefix."""
    S = mdl["S"].tocsc()
    rows_of_interest = {i for i, mid in enumerate(mdl["mets"]) if mid.startswith(met_prefix)}
    hits = set()
    for j in range(S.shape[1]):
        col = S[:, j]
        if any(int(r) in rows_of_interest for r in col.indices):
            hits.add(mdl["rxns"][j])
    return hits


def _parse_gpr(gpr):
    gpr = (gpr or "").strip()
    if not gpr: return []
    toks = re.findall(r"[A-Za-z][A-Za-z0-9._\-]*", gpr)
    toks = [t for t in toks if t.lower() not in ("and", "or")]
    out = []
    for t in toks:
        out.append((t, SYM2ENS.get(t.lower(), set())))
    return out


def _classify(gpr_tokens, loser_cells):
    """Return (verdict, detail list) given list of (sym, {ensIDs}) and the
    two loser-cell names (for zFPKM state lookup)."""
    if len(gpr_tokens) == 0:
        return "GENELESS", []
    gene_detail = []
    n_mappable = n_above_in_any_loser = n_below_inex_in_all_losers = n_unknown = 0
    for sym, ens_set in gpr_tokens:
        if not ens_set:
            gene_detail.append(dict(symbol=sym, states={c: None for c in loser_cells}))
            continue
        n_mappable += 1
        states = {}
        for c in loser_cells:
            zc_vals = [ZFPKM.loc[e, c] for e in ens_set
                       if e in ZFPKM.index and c in ZFPKM.columns and pd.notna(ZFPKM.loc[e, c])]
            if not zc_vals:
                states[c] = ("missing", None); continue
            zc = max(zc_vals)
            exp_t = THRESH.loc[c, "exp_threshold_zFPKM"] if c in THRESH.index else np.nan
            inex_t = THRESH.loc[c, "inex_threshold_zFPKM"] if c in THRESH.index else np.nan
            if pd.notna(exp_t) and zc >= exp_t:   states[c] = ("expressed", float(zc))
            elif pd.notna(inex_t) and zc <= inex_t: states[c] = ("off", float(zc))
            else:                                   states[c] = ("unknown", float(zc))
        gene_detail.append(dict(symbol=sym, states={c: states[c] for c in loser_cells}))
        state_vals = [states[c][0] for c in loser_cells]
        if any(s == "expressed" for s in state_vals):
            n_above_in_any_loser += 1
        elif all(s == "off" for s in state_vals):
            n_below_inex_in_all_losers += 1
        else:
            n_unknown += 1

    if n_mappable == 0:
        return "UNMAPPABLE", gene_detail
    if n_below_inex_in_all_losers > 0 and n_above_in_any_loser == 0:
        return "TRANSCRIPTIONAL_PARTIAL", gene_detail
    if n_above_in_any_loser > 0:
        return f"ALGORITHM_UNKNOWN_RETAINED ({n_above_in_any_loser}/{n_mappable} loser-expressed)", gene_detail
    return f"EXTRACTION_UNKNOWN_ZONE ({n_unknown}/{n_mappable} unknown in losers)", gene_detail


all_rows = []
case_summaries = []

def _reactions_in_subs(mdl, target_subs):
    return {mdl["rxns"][j] for j, s in enumerate(mdl["subs"]) if s in target_subs}


for substrate, met_prefix, catabolic_subs, hpf, winner, losers in CASES:
    winner_cell = f"{winner}_{hpf}"
    loser_cells = [f"{c}_{hpf}" for c in losers]
    mdl_w = _load(V7_EXTRACTED / f"trans_rfastcormics_{winner_cell}.mat")
    mdl_l = [_load(V7_EXTRACTED / f"trans_rfastcormics_{c}.mat") for c in loser_cells]
    # Combine substrate-touching AND catabolic-subsystem scope
    rxns_w = _reactions_touching(mdl_w, met_prefix) | _reactions_in_subs(mdl_w, catabolic_subs)
    rxns_l_union = set().union(
        *(_reactions_touching(m, met_prefix) | _reactions_in_subs(m, catabolic_subs)
          for m in mdl_l))
    diff = rxns_w - rxns_l_union
    print(f"\n=== {substrate}  hpf={hpf}  winner={winner_cell} vs {loser_cells} ===")
    print(f"  reactions touching {met_prefix}* in winner: {len(rxns_w)}")
    print(f"  reactions touching {met_prefix}* in loser-union: {len(rxns_l_union)}")
    print(f"  winner \\ (loser union): {len(diff)} reactions")

    verdicts_count = {}
    rows = []
    for rid in sorted(diff):
        j = mdl_w["rxns"].index(rid)
        gpr = mdl_w["grrules"][j] if j < len(mdl_w["grrules"]) else ""
        gpr_tokens = _parse_gpr(gpr)
        verdict, gene_detail = _classify(gpr_tokens, loser_cells)
        key = verdict.split(" ")[0]
        verdicts_count[key] = verdicts_count.get(key, 0) + 1
        row = dict(substrate=substrate, hpf=hpf, winner=winner_cell,
                   losers=";".join(loser_cells),
                   rxn=rid, rxn_name=mdl_w["rxn_names"][j] if j < len(mdl_w["rxn_names"]) else "",
                   subsystem=mdl_w["subs"][j] if j < len(mdl_w["subs"]) else "",
                   gpr=gpr, verdict=verdict,
                   gene_detail=json.dumps(gene_detail))
        rows.append(row)
        all_rows.append(row)
    print("  verdicts:")
    for k, v in sorted(verdicts_count.items()):
        print(f"    {v:3d}  {k}")
    for r in rows:
        print(f"    {r['rxn']:10s}  [{r['subsystem'][:35]:35s}]  gpr=`{r['gpr'][:45]:45s}`  {r['verdict']}")
    case_summaries.append(dict(substrate=substrate, hpf=hpf, winner=winner_cell,
                               n_diff=len(diff), **verdicts_count))

df = pd.DataFrame(all_rows)
df.to_csv(OUT / "s177_atp_one_succeeds_audit.csv", index=False)
print(f"\n[write] {OUT/'s177_atp_one_succeeds_audit.csv'}  ({len(df)} rows)")

# Summary markdown
md = ["# s177 — ATP sweep one-succeeds audit", ""]
md.append("Six cases from `atp_sweep_v7_clean.csv` where exactly one regime")
md.append("could catabolise the substrate for ATP, and the other two could not.")
md.append("Mirror of the ten single-failure audit; question flipped to: what")
md.append("does the winner retain that both losers lack, and does the")
md.append("retention correspond to a transcriptional difference?")
md.append("")
md.append("### Per-case verdict counts")
md.append("")
md.append("| substrate | hpf | winner | n_diff | geneless | transc. | unknown-zone | algorithm | unmappable |")
md.append("|---|---|---|---|---|---|---|---|---|")
for s in case_summaries:
    row = [s.get("substrate"), str(s.get("hpf")), s.get("winner"), str(s.get("n_diff")),
           str(s.get("GENELESS", 0)), str(s.get("TRANSCRIPTIONAL_PARTIAL", 0)),
           str(s.get("EXTRACTION_UNKNOWN_ZONE", 0)),
           str(s.get("ALGORITHM_UNKNOWN_RETAINED", 0)),
           str(s.get("UNMAPPABLE", 0))]
    md.append("| " + " | ".join(row) + " |")
md.append("")
# Aggregated verdict across all 6 cases
agg = {}
for r in all_rows:
    key = r["verdict"].split(" ")[0]
    agg[key] = agg.get(key, 0) + 1
md.append("### Aggregate across all six cases")
md.append("")
md.append(f"- Total winner-only reactions: **{len(all_rows)}**")
for k in sorted(agg.keys()):
    md.append(f"- {k}: **{agg[k]}** ({100*agg[k]/len(all_rows):.0f}%)")
md.append("")
md.append("### LD asymmetry")
md.append("")
ld_cases = [s for s in case_summaries if s["winner"].startswith("LD_")]
ld_total = sum(s.get("n_diff", 0) for s in ld_cases)
md.append(f"- 4 of 6 one-succeeds cases are LD. Winner-only reaction count in those four: **{ld_total}**")
md.append("  (if dominated by unknown-zone / geneless, this is consistent with LD's extraction retaining")
md.append("  more unknown-zone reactions by extraction-algorithm chance, not with LD-specific biology).")
(OUT / "s177_atp_one_succeeds_audit.md").write_text("\n".join(md))
print(f"[write] {OUT/'s177_atp_one_succeeds_audit.md'}")

print("\n=== AGGREGATE VERDICT ===")
for k in sorted(agg.keys()):
    print(f"  {agg[k]:3d}  {k}  ({100*agg[k]/len(all_rows):.0f}%)")
