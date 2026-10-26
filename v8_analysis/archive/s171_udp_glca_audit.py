"""
s171 -- UDP-GlcA pathway GPR audit for the two topology-dominant pairs.

Same scheme as s170 (gss). For each reaction in A\\B among UDP-GlcA–relevant
subsystems, classify as:
  - GENELESS                 (empty GPR)
  - TRANSCRIPTIONAL_PARTIAL  (any gene below inex threshold, none above exp)
  - EXTRACTION_UNKNOWN_ZONE  (all genes in [inex, exp] — rFASTCORMICS edge case)
  - ALGORITHM_UNKNOWN_RETAINED  (some gene above exp but reaction pruned)

Pairs:
  BL_96  vs D_96   UDP_GlcA  (topo_frac = 0.74)
  D_120  vs LD_120 UDP_GlcA  (topo_frac = 1.26)

Target subsystems (UDP-GlcA synthesis + nucleotide-sugar + UDP-glucose feeders):
  "Pentose and glucuronate interconversions"
  "Nucleotide sugar metabolism" / "Nucleotide metabolism"
  "Starch and sucrose metabolism"
  "Glycolysis / Gluconeogenesis"  (UDP-glucose feed)
  "Pentose phosphate pathway"
"""
import sys, json, re
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"

ZFPKM = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv").set_index("ensembl_gene_id")
THRESH = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv").set_index("cell")
_GMAP = pd.read_csv(REPO / "results/qc/gene_map_ensdarg_to_model.csv")
SYM2ENS = {}
for _, r in _GMAP.iterrows():
    SYM2ENS.setdefault(str(r.model_gene).lower(), set()).add(str(r.ensembl_gene_id))

PAIRS = [
    ("BL_96",  "D_96",   "UDP_GlcA"),
    ("D_120",  "LD_120", "UDP_GlcA"),
]

TARGET_SUBS = {
    "Pentose and glucuronate interconversions",
    "Nucleotide sugar metabolism",
    "Nucleotide metabolism",
    "Starch and sucrose metabolism",
    "Glycolysis / Gluconeogenesis",
    "Pentose phosphate pathway",
    "Fructose and mannose metabolism",
    "Galactose metabolism",
    "Amino sugar metabolism",
    "Keratan sulfate biosynthesis",
    "Keratan sulfate degradation",
    "Chondroitin / heparan sulfate biosynthesis",
    "Chondroitin sulfate degradation",
    "Heparan sulfate degradation",
    "Glycerolipid metabolism",   # UDP-glucose can feed via DAG-glucose
}


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
    return dict(rxns=rxns, rxn_names=rxn_names, subs=subs, grrules=grrules)


def _parse_gpr(gpr):
    gpr = (gpr or "").strip()
    if not gpr: return []
    toks = re.findall(r"[A-Za-z][A-Za-z0-9._\-]*", gpr)
    toks = [t for t in toks if t.lower() not in ("and", "or")]
    out = []
    for t in toks:
        out.append((t, SYM2ENS.get(t.lower(), set())))
    return out


def _audit_pair(a, b):
    mA = _load(V7_EXTRACTED / f"trans_rfastcormics_{a}.mat")
    mB = _load(V7_EXTRACTED / f"trans_rfastcormics_{b}.mat")
    set_B = set(mB["rxns"])
    rows = []
    for j, rid in enumerate(mA["rxns"]):
        sub = mA["subs"][j] if j < len(mA["subs"]) else ""
        if sub not in TARGET_SUBS: continue
        if rid in set_B: continue
        gpr = mA["grrules"][j] if j < len(mA["grrules"]) else ""
        gpr_tokens = _parse_gpr(gpr)
        gene_detail = []
        n_above_exp = n_below_inex = n_unknown = n_mappable = 0
        thresh_B_exp = THRESH.loc[b, "exp_threshold_zFPKM"] if b in THRESH.index else np.nan
        thresh_B_inex = THRESH.loc[b, "inex_threshold_zFPKM"] if b in THRESH.index else np.nan
        for sym, ens_set in gpr_tokens:
            if not ens_set:
                gene_detail.append(dict(symbol=sym, ensembl=None, zFPKM_A=None, zFPKM_B=None, state_B=None))
                continue
            n_mappable += 1
            zA_vals, zB_vals = [], []
            for e in ens_set:
                if e in ZFPKM.index:
                    if a in ZFPKM.columns: zA_vals.append(ZFPKM.loc[e, a])
                    if b in ZFPKM.columns: zB_vals.append(ZFPKM.loc[e, b])
            zA = max([v for v in zA_vals if pd.notna(v)], default=np.nan)
            zB = max([v for v in zB_vals if pd.notna(v)], default=np.nan)
            if pd.notna(zB) and pd.notna(thresh_B_exp) and zB >= thresh_B_exp:
                st = "expressed"; n_above_exp += 1
            elif pd.notna(zB) and pd.notna(thresh_B_inex) and zB <= thresh_B_inex:
                st = "off"; n_below_inex += 1
            else:
                st = "unknown"; n_unknown += 1
            gene_detail.append(dict(symbol=sym, zFPKM_A=float(zA) if pd.notna(zA) else None,
                                    zFPKM_B=float(zB) if pd.notna(zB) else None, state_B=st))

        if len(gpr_tokens) == 0:          verdict = "GENELESS"
        elif n_mappable == 0:              verdict = "UNMAPPABLE"
        elif n_below_inex > 0 and n_above_exp == 0:  verdict = "TRANSCRIPTIONAL_PARTIAL"
        elif n_above_exp > 0:              verdict = f"ALGORITHM_UNKNOWN_RETAINED ({n_above_exp}/{n_mappable} expressed, {n_unknown}/{n_mappable} unknown)"
        else:                              verdict = f"EXTRACTION_UNKNOWN_ZONE ({n_unknown}/{n_mappable} unknown, {n_below_inex}/{n_mappable} off)"

        rows.append(dict(cell_A=a, cell_B=b, rxn=rid,
                         rxn_name=mA["rxn_names"][j] if j < len(mA["rxn_names"]) else "",
                         subsystem=sub, gpr=gpr, n_tokens=len(gpr_tokens),
                         n_mappable=n_mappable, n_above_exp=n_above_exp,
                         n_below_inex=n_below_inex, n_unknown=n_unknown,
                         verdict=verdict, gene_detail=json.dumps(gene_detail),
                         thresh_B_exp=float(thresh_B_exp) if pd.notna(thresh_B_exp) else None,
                         thresh_B_inex=float(thresh_B_inex) if pd.notna(thresh_B_inex) else None))
    return rows


all_rows = []
for a, b, cof in PAIRS:
    rows = _audit_pair(a, b)
    all_rows.extend(rows)
    print(f"\n=== {a} vs {b}  ({cof}) :: {len(rows)} reactions in A\\B among UDP-GlcA subsystems ===")
    cnt = {}
    for r in rows:
        key = r['verdict'].split(' ')[0]
        cnt[key] = cnt.get(key, 0) + 1
    for k, v in sorted(cnt.items()):
        print(f"  {v:3d}  {k}")
    for r in rows:
        print(f"  {r['rxn']:10s}  [{r['subsystem'][:40]:40s}]  gpr=`{r['gpr'][:50]:50s}`  verdict={r['verdict']}")

df = pd.DataFrame(all_rows)
df.to_csv(OUT / "udp_glca_pathway_audit.csv", index=False)
print(f"\n[write] {OUT/'udp_glca_pathway_audit.csv'}  ({len(df)} rows)")
