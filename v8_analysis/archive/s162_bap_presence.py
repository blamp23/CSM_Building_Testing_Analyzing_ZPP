"""
s162 -- §1.3 BaP presence check across v7 extracted + base models.

Searches metNames for polyaromatic xenobiotic substrings; lists every
reaction consuming/producing those metabolites, their GPRs, subsystems,
and per-cell presence across the 15 extracted models and 5 base models.

Read-only. No bound changes. No LP.

Outputs (grace_v4/results/v7_susceptibility/):
  bap_metabolites.csv         met_id, metName, comp, matched_substring, in_each_model
  bap_reactions.csv           rxn, rxnName, subsystem, formula, GPR, in_each_cell, in_each_base
  bap_presence_matrix.csv     rxn × 20 models binary
  bap_presence.md             human-readable summary
"""
import sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
BASES = REPO / "ocr_anchored_extraction/models"
OUT = REPO / "grace_v4/results/v7_susceptibility"
OUT.mkdir(parents=True, exist_ok=True)

CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

SUBSTR = ("pyrene", "naphth", "phenanthr", "benzo", "anthrac")


def _resolve(f, x):
    if isinstance(x, h5py.Reference):
        return _resolve(f, np.asarray(f[x]).ravel())
    if isinstance(x, np.ndarray):
        if x.dtype == object:
            return _resolve(f, x[0])
        try:
            return "".join(chr(int(c)) for c in x)
        except Exception:
            return ""
    return str(x)


def _load(fp, root_key):
    with h5py.File(fp, "r") as f:
        m = f[root_key]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        try:
            rxn_names = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
        except Exception:
            rxn_names = [""] * len(rxns)
        try:
            subs = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception:
            subs = [""] * len(rxns)
        try:
            grrules = [_resolve(f, r) for r in np.asarray(m["grRules"]).ravel()]
        except Exception:
            grrules = [""] * len(rxns)
        try:
            mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception:
            mets = []
        try:
            met_names = [_resolve(f, r) for r in np.asarray(m["metNames"]).ravel()]
        except Exception:
            met_names = []
        try:
            met_comps = np.asarray(m["metComps"]).ravel().astype(int)
        except Exception:
            met_comps = np.zeros(len(mets), dtype=int)
        try:
            comp_names = [_resolve(f, r) for r in np.asarray(m["compNames"]).ravel()]
        except Exception:
            comp_names = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"])
        ir = np.asarray(S_grp["ir"])
        jc = np.asarray(S_grp["jc"])
        n_met_hdf = len(mets)
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, n_met_hdf, len(met_names))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
    return dict(rxns=rxns, rxn_names=rxn_names, subs=subs, grrules=grrules,
                mets=mets, met_names=met_names, met_comps=met_comps,
                comp_names=comp_names, S=S)


def _comp_of(m, mi):
    c = int(m["met_comps"][mi]) - 1 if mi < len(m["met_comps"]) else -1
    return m["comp_names"][c] if 0 <= c < len(m["comp_names"]) else ""


def _match(name):
    low = name.lower()
    hits = [s for s in SUBSTR if s in low]
    return hits


def _formula(m, j):
    S = m["S"]
    col = S.tocsc()[:, j]
    subs_, prods_ = [], []
    for r_, v_ in zip(col.indices, col.data):
        nm = (m["met_names"][r_] if r_ < len(m["met_names"]) else
              m["mets"][r_] if r_ < len(m["mets"]) else f"met{r_}")
        cmp = _comp_of(m, r_)
        tag = f"{nm}[{cmp}]" if cmp else nm
        (subs_ if v_ < 0 else prods_).append((abs(float(v_)), tag))
    s_str = " + ".join(f"{('' if c == 1 else f'{c:g} ')}{t}" for c, t in subs_) or "∅"
    p_str = " + ".join(f"{('' if c == 1 else f'{c:g} ')}{t}" for c, t in prods_) or "∅"
    return f"{s_str} -> {p_str}"


# Load all models
print("[load] 15 extracted + 5 base models", flush=True)
models = {}
for cell in CELLS:
    fp = V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat"
    models[("ex", cell)] = _load(fp, "parent")
for hpf in HPFS:
    fp = BASES / f"baked_ocr_{hpf}_v7.mat"
    if fp.exists():
        models[("base", str(hpf))] = _load(fp, "m")
    else:
        print(f"  [skip] base baked_ocr_{hpf}_v7.mat not found", flush=True)
print(f"[load] got {sum(1 for k in models if k[0]=='ex')} extracted, "
      f"{sum(1 for k in models if k[0]=='base')} base", flush=True)

# ----------- metabolite scan on union of base models (richest met list) -----
# Build canonical met table from LD_24 extracted (always present) and any base
canon_src = models.get(("base", "24")) or models[("ex", "LD_24")]
met_rows = []
for mi, nm in enumerate(canon_src["met_names"]):
    hits = _match(nm)
    if not hits:
        continue
    met_rows.append(dict(
        met_idx=mi,
        met_id=canon_src["mets"][mi] if mi < len(canon_src["mets"]) else f"met{mi}",
        met_name=nm,
        comp=_comp_of(canon_src, mi),
        matched=",".join(hits),
    ))

# Also scan each model's met table for things the canonical missed
extra_mets = {}
for key, m in models.items():
    for mi, nm in enumerate(m["met_names"]):
        hits = _match(nm)
        if not hits:
            continue
        mid = m["mets"][mi] if mi < len(m["mets"]) else f"met{mi}"
        if mid in {r["met_id"] for r in met_rows}:
            continue
        if mid in extra_mets:
            continue
        extra_mets[mid] = dict(
            met_idx=mi,
            met_id=mid,
            met_name=nm,
            comp=_comp_of(m, mi),
            matched=",".join(hits),
            source=f"{key[0]}:{key[1]}",
        )
met_rows.extend(extra_mets.values())

# For each metabolite, record which models contain it
all_met_ids = sorted({r["met_id"] for r in met_rows})
met_presence = {mid: [] for mid in all_met_ids}
for key, m in models.items():
    tag = f"{key[0]}:{key[1]}"
    for mid in all_met_ids:
        if mid in m["mets"]:
            met_presence[mid].append(tag)
for r in met_rows:
    r["n_models"] = len(met_presence[r["met_id"]])
    r["present_in"] = ";".join(met_presence[r["met_id"]])

met_df = pd.DataFrame(met_rows).sort_values(["matched", "met_name"]).reset_index(drop=True)
met_df.to_csv(OUT / "bap_metabolites.csv", index=False)
print(f"[write] {OUT/'bap_metabolites.csv'}  ({len(met_df)} metabolites)", flush=True)

# ----------- reaction scan: any reaction touching those met_ids -------------
# Walk each model; for each matched met, enumerate reactions touching it.
rxn_hits = {}  # rxn_id -> row dict
hit_mets_lower = {r["met_id"] for r in met_rows}

for key, m in models.items():
    tag = f"{key[0]}:{key[1]}"
    mid_to_idx = {mid: i for i, mid in enumerate(m["mets"])}
    target_idx = {mid_to_idx[mid] for mid in hit_mets_lower if mid in mid_to_idx}
    if not target_idx:
        continue
    S_csc = m["S"].tocsc()
    for j in range(S_csc.shape[1]):
        col = S_csc[:, j]
        if not any(int(r) in target_idx for r in col.indices):
            continue
        rid = m["rxns"][j]
        if rid not in rxn_hits:
            rxn_hits[rid] = dict(
                rxn=rid,
                rxn_name=m["rxn_names"][j] if j < len(m["rxn_names"]) else "",
                subsystem=m["subs"][j] if j < len(m["subs"]) else "",
                grrule=m["grrules"][j] if j < len(m["grrules"]) else "",
                formula=_formula(m, j),
                present_in=[],
            )
        rxn_hits[rid]["present_in"].append(tag)

# Flatten
rxn_rows = []
model_keys_ordered = [("base", str(h)) for h in HPFS] + [("ex", c) for c in CELLS]
model_tags = [f"{k[0]}:{k[1]}" for k in model_keys_ordered if k in models]
for rid, info in rxn_hits.items():
    present_list = info["present_in"]  # a list
    row = {k: v for k, v in info.items() if k != "present_in"}
    row["present_in"] = ";".join(present_list)
    for tag in model_tags:
        row[tag] = int(tag in present_list)
    row["n_models"] = sum(row[tag] for tag in model_tags)
    rxn_rows.append(row)

rxn_df = pd.DataFrame(rxn_rows).sort_values(["subsystem", "rxn"]).reset_index(drop=True)
rxn_df.to_csv(OUT / "bap_reactions.csv", index=False)
print(f"[write] {OUT/'bap_reactions.csv'}  ({len(rxn_df)} reactions)", flush=True)

# Presence matrix only
pm = rxn_df[["rxn", "rxn_name", "subsystem"] + model_tags].copy()
pm.to_csv(OUT / "bap_presence_matrix.csv", index=False)
print(f"[write] {OUT/'bap_presence_matrix.csv'}", flush=True)

# ----------- Markdown summary ------------------------------------------------
md = ["# BaP / polyaromatic xenobiotic presence check (v7)",
      "",
      f"Generated by `s162_bap_presence.py`.",
      f"Scanned {sum(1 for k in models if k[0]=='base')} base models + "
      f"{sum(1 for k in models if k[0]=='ex')} extracted cells.",
      f"Substrings searched (case-insensitive): {', '.join(SUBSTR)}.",
      "",
      "## Metabolites matched",
      "",
      f"**{len(met_df)}** metabolite entries matched.",
      "",
      "| met_id | metName | comp | matched | #models | present in |",
      "|---|---|---|---|---|---|"]
for _, r in met_df.iterrows():
    md.append(f"| {r.met_id} | {r.met_name} | {r.comp} | {r.matched} | {r.n_models} | {r.present_in} |")
md += ["", "## Reactions touching those metabolites", "",
       f"**{len(rxn_df)}** reactions total.",
       "",
       "### Count per model (how many of the matched reactions are present)",
       ""]
counts = {tag: int(rxn_df[tag].sum()) if tag in rxn_df.columns else 0 for tag in model_tags}
md.append("| model | rxns present |")
md.append("|---|---|")
for tag in model_tags:
    md.append(f"| {tag} | {counts[tag]} |")
md += ["", "### All matched reactions", ""]
for _, r in rxn_df.iterrows():
    md.append(f"- **{r.rxn}** ({r.subsystem}) `{r.rxn_name}` — present in **{r.n_models}** models")
    md.append(f"    - formula: `{r.formula}`")
    md.append(f"    - GPR: `{r.grrule or '(no GPR)'}`")
(OUT / "bap_presence.md").write_text("\n".join(md))
print(f"[write] {OUT/'bap_presence.md'}", flush=True)

# Verdict
print("\n--- verdict ---")
if rxn_df.empty:
    print("No BaP-related reactions found in any extracted or base model.")
elif rxn_df[[f"ex:{c}" for c in CELLS]].sum().sum() == 0:
    print("BaP-related reactions exist in base models but ZERO carried into any extracted cell.")
    print("→ Probe 5 'as-extracted' is uniformly infeasible; forced-in version is the only interpretable one.")
else:
    carry = (rxn_df[[f"ex:{c}" for c in CELLS]].sum(axis=1) > 0).sum()
    print(f"{carry}/{len(rxn_df)} BaP reactions are present in at least one extracted cell.")
    per_cell = rxn_df[[f"ex:{c}" for c in CELLS]].sum()
    print(per_cell.to_string())
