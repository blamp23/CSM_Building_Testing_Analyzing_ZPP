"""
s170 -- GPR audit: BL_24 vs LD_24 (and D_120 vs LD_120, BL_120 vs D_120) GSH
topology gaps.

Lists trans-sulfuration / Methionine-cysteine / Glutathione metabolism
reactions present in cell_A but absent in cell_B, with GPR and per-gene
per-cell zFPKM score + per-cell expression threshold.

For each reaction in A \ B:
  - GPR empty            -> extraction artefact (gene-less; same class as MAR00466)
  - GPR non-empty, every gene below LD threshold in LD -> transcriptional
                             (real finding; rFASTCORMICS pruned because off)
  - GPR non-empty, at least one gene above threshold in LD -> extraction
                             did not include it despite expression (algorithm
                             artefact or weighted-core issue)

Output:
  gsh_pathway_audit.csv              rxn, present_A/B, subsystem, formula, GPR,
                                     n_genes_above_LD_thresh, gene_detail (JSON)
  gsh_pathway_audit.md               human-readable per-reaction detail
"""
import sys, json, re
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"

ZFPKM = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv")
ZFPKM = ZFPKM.set_index("ensembl_gene_id")
THRESH = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv").set_index("cell")
# Model-gene symbol → ENSDARG mapping (one gene may map to multiple ensemblIDs)
_GMAP = pd.read_csv(REPO / "results/qc/gene_map_ensdarg_to_model.csv")
SYM2ENS = {}
for _, r in _GMAP.iterrows():
    SYM2ENS.setdefault(str(r.model_gene).lower(), set()).add(str(r.ensembl_gene_id))

PAIRS = [
    ("BL_24",  "LD_24",  "GSH"),
    ("BL_24",  "D_24",   "GSH"),
    ("D_120",  "LD_120", "GSH"),
    ("BL_120", "D_120",  "GSH"),
]

# Subsystems whose reactions feed GSH (precursors + synthesis)
TARGET_SUBS = {
    "Methionine and cysteine metabolism",
    "Cysteine and methionine metabolism",
    "Glutathione metabolism",
    "Transsulfuration",
    "Sulfur metabolism",
    "Glycine, serine and threonine metabolism",
}


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
        try: rxn_names = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
        except Exception: rxn_names = [""] * len(rxns)
        try: subs = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception: subs = [""] * len(rxns)
        try: grrules = [_resolve(f, r) for r in np.asarray(m["grRules"]).ravel()]
        except Exception: grrules = [""] * len(rxns)
        try: mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception: mets = []
        try: met_names = [_resolve(f, r) for r in np.asarray(m["metNames"]).ravel()]
        except Exception: met_names = []
        try: met_comps = np.asarray(m["metComps"]).ravel().astype(int)
        except Exception: met_comps = np.zeros(len(mets), dtype=int)
        try: comp_names = [_resolve(f, r) for r in np.asarray(m["compNames"]).ravel()]
        except Exception: comp_names = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets), len(met_names))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
    return dict(rxns=rxns, rxn_names=rxn_names, subs=subs, grrules=grrules,
                mets=mets, met_names=met_names, met_comps=met_comps,
                comp_names=comp_names, S=S)


def _formula(m, j):
    S = m["S"].tocsc(); col = S[:, j]
    subs_, prods_ = [], []
    for r_, v_ in zip(col.indices, col.data):
        nm = m["met_names"][r_] if r_ < len(m["met_names"]) else m["mets"][r_]
        c = int(m["met_comps"][r_]) - 1 if r_ < len(m["met_comps"]) else -1
        cmp = m["comp_names"][c] if 0 <= c < len(m["comp_names"]) else ""
        tag = f"{nm}[{cmp}]"
        (subs_ if v_ < 0 else prods_).append((abs(float(v_)), tag))
    s_str = " + ".join(f"{(c if c!=1 else '')}{t}" for c, t in subs_) or "∅"
    p_str = " + ".join(f"{(c if c!=1 else '')}{t}" for c, t in prods_) or "∅"
    return f"{s_str} -> {p_str}"


# Load all 15 cells
CELLS_ALL = [f"{c}_{h}" for c in ("BL","D","LD") for h in (24,48,72,96,120)]
models = {}
for cell in CELLS_ALL:
    models[cell] = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")


def _parse_gpr(gpr):
    """Extract model gene symbols (lowercase alphanumeric + '.') from a GPR string,
    resolve each to a set of ENSDARG IDs via SYM2ENS. Returns list of (symbol, ensembl_set)."""
    gpr = (gpr or "").strip()
    if not gpr:
        return []
    # Tokens are symbols separated by " and ", " or ", "(", ")"
    toks = re.findall(r"[A-Za-z][A-Za-z0-9._\-]*", gpr)
    # Filter out boolean operators
    toks = [t for t in toks if t.lower() not in ("and", "or")]
    out = []
    for t in toks:
        ens = SYM2ENS.get(t.lower(), set())
        out.append((t, ens))
    return out


def _audit_pair(cell_A, cell_B):
    """Reactions in A \ B among TARGET_SUBS."""
    mA, mB = models[cell_A], models[cell_B]
    set_B = set(mB["rxns"])
    rows = []
    for j, rid in enumerate(mA["rxns"]):
        sub = mA["subs"][j] if j < len(mA["subs"]) else ""
        if sub not in TARGET_SUBS: continue
        if rid in set_B: continue  # present in both
        gpr = mA["grrules"][j] if j < len(mA["grrules"]) else ""
        formula = _formula(mA, j)
        gpr_tokens = _parse_gpr(gpr)   # [(symbol, {ensIDs})]
        gene_detail = []
        n_above_exp = 0        # zFPKM > exp_threshold (actively expressed)
        n_below_inex = 0       # zFPKM < inex_threshold (actively OFF)
        n_unknown = 0          # between -- rFASTCORMICS unknown zone
        n_mappable = 0
        thresh_B_exp = THRESH.loc[cell_B, "exp_threshold_zFPKM"] if cell_B in THRESH.index else np.nan
        thresh_B_inex = THRESH.loc[cell_B, "inex_threshold_zFPKM"] if cell_B in THRESH.index else np.nan
        thresh_A_exp = THRESH.loc[cell_A, "exp_threshold_zFPKM"] if cell_A in THRESH.index else np.nan
        thresh_A_inex = THRESH.loc[cell_A, "inex_threshold_zFPKM"] if cell_A in THRESH.index else np.nan
        for sym, ens_set in gpr_tokens:
            if not ens_set:
                gene_detail.append(dict(symbol=sym, ensembl=None, zFPKM_A=None, zFPKM_B=None, state_B=None))
                continue
            n_mappable += 1
            zA_vals = []; zB_vals = []
            for e in ens_set:
                if e in ZFPKM.index:
                    if cell_A in ZFPKM.columns:
                        zA_vals.append(ZFPKM.loc[e, cell_A])
                    if cell_B in ZFPKM.columns:
                        zB_vals.append(ZFPKM.loc[e, cell_B])
            zA = max([v for v in zA_vals if pd.notna(v)], default=np.nan)
            zB = max([v for v in zB_vals if pd.notna(v)], default=np.nan)
            if pd.notna(zB) and pd.notna(thresh_B_exp) and zB >= thresh_B_exp:
                state_B = "expressed"
                n_above_exp += 1
            elif pd.notna(zB) and pd.notna(thresh_B_inex) and zB <= thresh_B_inex:
                state_B = "off"
                n_below_inex += 1
            else:
                state_B = "unknown"
                n_unknown += 1
            gene_detail.append(dict(symbol=sym,
                                    ensembl=list(ens_set),
                                    zFPKM_A=float(zA) if pd.notna(zA) else None,
                                    zFPKM_B=float(zB) if pd.notna(zB) else None,
                                    state_B=state_B))

        # Verdict policy:
        #   Any gene missing (unmappable or no zFPKM) treated as "unknown".
        #   For an AND (all needed): if ANY required gene is OFF, FASTCORE can prune.
        #   For an OR (any works): if ALL alternatives are OFF, FASTCORE can prune.
        #   Without parsing full AND/OR, we conservatively use these per-reaction labels:
        if len(gpr_tokens) == 0:
            verdict = "GENELESS_EXTRACTION_ARTEFACT"
        elif n_mappable == 0:
            verdict = "UNMAPPABLE_GPR"
        elif n_below_inex > 0 and n_above_exp == 0:
            verdict = "TRANSCRIPTIONAL_PARTIAL"   # at least one gene explicitly OFF, none clearly ON
        elif n_above_exp > 0:
            verdict = f"ALGORITHM_UNKNOWN_ZONE_or_RETAINED ({n_above_exp}/{n_mappable} expressed, {n_unknown}/{n_mappable} unknown)"
        else:
            # All in "unknown" zone (between inex and exp thresholds)
            verdict = f"EXTRACTION_UNKNOWN_ZONE ({n_unknown}/{n_mappable} unknown, {n_below_inex}/{n_mappable} off)"
        n_above = n_above_exp  # keep for backward-compat

        rows.append(dict(
            cell_A=cell_A, cell_B=cell_B, rxn=rid,
            rxn_name=mA["rxn_names"][j] if j < len(mA["rxn_names"]) else "",
            subsystem=sub, formula=formula, gpr=gpr,
            n_tokens=len(gpr_tokens), n_mappable=n_mappable,
            n_above_thresh_in_B=n_above,
            verdict=verdict, gene_detail=json.dumps(gene_detail),
            thresh_B_expressed=float(thresh_B_exp) if pd.notna(thresh_B_exp) else None,
            thresh_B_inex=float(thresh_B_inex) if pd.notna(thresh_B_inex) else None,
            thresh_A_expressed=float(thresh_A_exp) if pd.notna(thresh_A_exp) else None,
            thresh_A_inex=float(thresh_A_inex) if pd.notna(thresh_A_inex) else None,
            n_above_exp=n_above_exp, n_below_inex=n_below_inex, n_unknown=n_unknown,
        ))
    return rows


all_rows = []
for a, b, cof in PAIRS:
    print(f"\n=== {a}  vs  {b}  ({cof})  :: reactions in A \\ B among target subsystems ===")
    rows = _audit_pair(a, b)
    all_rows.extend(rows)
    v_counts = {}
    for r in rows:
        v_counts[r["verdict"]] = v_counts.get(r["verdict"], 0) + 1
    print(f"  total: {len(rows)}")
    for k, v in sorted(v_counts.items()):
        print(f"    {v:3d}  {k}")
    for r in rows:
        print(f"    {r['rxn']:10s}  sub={r['subsystem']:30s}  gpr=`{r['gpr'][:60]}`  verdict={r['verdict']}")

df = pd.DataFrame(all_rows)
df.to_csv(OUT / "gsh_pathway_audit.csv", index=False)
print(f"\n[write] {OUT/'gsh_pathway_audit.csv'}  ({len(df)} rows)")

# Markdown
md = ["# GSH-pathway GPR audit (trans-sulfuration + glutathione / Me-Cys / sulfur / Gly-Ser-Thr)", ""]
md.append(f"Scope: 4 pairs where Probe 1 attribution gave topology-dominant GSH gaps.")
md.append(f"Subsystems scanned: {sorted(TARGET_SUBS)}")
md.append("")
md.append("Each row below is a reaction present in cell A but absent in cell B.")
md.append("")
md.append("### Verdict legend")
md.append("")
md.append("- **GENELESS_EXTRACTION_ARTEFACT** — reaction has no GPR; rFASTCORMICS prunes gene-less reactions; limitation to disclose, not biology.")
md.append("- **TRANSCRIPTIONAL** — every gene in GPR is below the LD_24 expression threshold; genuine condition-specific retention.")
md.append("- **ALGORITHM_RETAINED_N_of_M_above_thresh** — N of M genes are above threshold in cell B yet the reaction was pruned by rFASTCORMICS; worth noting; weighted-core / algorithmic edge case.")
md.append("")
for a, b, cof in PAIRS:
    sub_df = df[(df.cell_A==a)&(df.cell_B==b)]
    md.append(f"## {a} vs {b} ({cof}) — {len(sub_df)} reactions in A \\ B")
    md.append("")
    for _, r in sub_df.iterrows():
        md.append(f"- **{r.rxn}** ({r.subsystem}) `{r.rxn_name}` — **{r.verdict}**")
        md.append(f"    - formula: `{r.formula}`")
        md.append(f"    - GPR: `{r.gpr or '(empty)'}`")
        if r.n_tokens:
            md.append(f"    - {r.n_above_thresh_in_B}/{r.n_mappable} genes above {b} expressed threshold ({r.thresh_B_expressed:.3f})")
            md.append(f"    - gene detail: `{r.gene_detail[:200]}...`" if len(r.gene_detail) > 200 else f"    - gene detail: `{r.gene_detail}`")
    md.append("")

(OUT / "gsh_pathway_audit.md").write_text("\n".join(md))
print(f"[write] {OUT/'gsh_pathway_audit.md'}")
