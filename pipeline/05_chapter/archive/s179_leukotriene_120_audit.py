"""
s179 -- Leukotriene metabolism GPR/zFPKM audit at 120 hpf.

For each enzyme-subset REPRESENTATIVE in the Leukotriene metabolism
subsystem, in each of the three 120-hpf cells (BL_120, D_120, LD_120),
classify the representative's GPR state in that cell:

  GENELESS        — representative's GPR is empty (template-inherited;
                    retained by FASTCORE without gene input)
  UNMAPPABLE      — GPR has symbols but none map to Ensembl in our gene map
  EXPRESSED       — ≥ 1 mappable gene with zFPKM ≥ exp threshold
  UNKNOWN_ZONE    — no gene above exp; ≥ 1 gene in the inex..exp band
  OFF             — all mappable genes with zFPKM ≤ inex threshold
                    (reaction retained despite its genes being OFF —
                    algorithmic / FASTCORE-connectivity rescue)

Each representative is weighted by whether it carries non-zero flux in the
pooled sample (> 2% of its own column variance; strict ∃ k with v_rk != 0).
"Share of flux-carrying reps in each class" = flux-carrying count per
class / total flux-carrying count.

Covers the two Tier A F1 rows:
    BL_120_vs_LD_120  Leukotriene metabolism  W1 = 0.613 direction = +1
    D_120_vs_LD_120   Leukotriene metabolism  W1 = 0.843 direction = +1

Output: grace_v4/results/v7_final_addenda/s179_leukotriene_120_audit.{csv,md}
"""
import sys, json, re
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import h5py

REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
OUT = R / "v7_final_addenda"
V7 = REPO / "reviewer_packet_v4_slim/v7_models"

CELLS = ["BL_120", "D_120", "LD_120"]
SUBSYSTEM = "Leukotriene metabolism"

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


def _load_model(fp):
    with h5py.File(fp, "r") as f:
        m = f["parent"]
        rxns = [_resolve(f, r) for r in np.asarray(m["rxns"]).ravel()]
        try: rxn_names = [_resolve(f, r) for r in np.asarray(m["rxnNames"]).ravel()]
        except Exception: rxn_names = [""] * len(rxns)
        try: grrules = [_resolve(f, r) for r in np.asarray(m["grRules"]).ravel()]
        except Exception: grrules = [""] * len(rxns)
    return dict(rxns=rxns, rxn_names=rxn_names, grrules=grrules)


def _parse_gpr(gpr):
    gpr = (gpr or "").strip()
    if not gpr: return []
    toks = re.findall(r"[A-Za-z][A-Za-z0-9._\-]*", gpr)
    toks = [t for t in toks if t.lower() not in ("and", "or")]
    out = []
    for t in toks:
        out.append((t, SYM2ENS.get(t.lower(), set())))
    return out


def _classify_cell(gpr_tokens, cell):
    """Classify a reaction's GPR state in one cell."""
    if len(gpr_tokens) == 0:
        return "GENELESS", []
    exp_t  = THRESH.loc[cell, "exp_threshold_zFPKM"]
    inex_t = THRESH.loc[cell, "inex_threshold_zFPKM"]
    gene_states = []
    n_mappable = n_exp = n_off = n_unk = 0
    for sym, ens_set in gpr_tokens:
        if not ens_set:
            gene_states.append(dict(symbol=sym, ensembl=None, z=None, state="unmappable"))
            continue
        n_mappable += 1
        zc_vals = [ZFPKM.loc[e, cell] for e in ens_set
                   if e in ZFPKM.index and cell in ZFPKM.columns and pd.notna(ZFPKM.loc[e, cell])]
        if not zc_vals:
            gene_states.append(dict(symbol=sym, ensembl=list(ens_set), z=None, state="missing"))
            continue
        zc = max(zc_vals)
        if zc >= exp_t:
            st = "expressed"; n_exp += 1
        elif zc <= inex_t:
            st = "off"; n_off += 1
        else:
            st = "unknown"; n_unk += 1
        gene_states.append(dict(symbol=sym, ensembl=list(ens_set), z=float(zc), state=st))
    if n_mappable == 0:
        return "UNMAPPABLE", gene_states
    if n_exp > 0:
        return "EXPRESSED", gene_states
    if n_unk > 0 and n_off == 0:
        return "UNKNOWN_ZONE", gene_states
    if n_off > 0 and n_unk == 0:
        return "OFF", gene_states
    return "UNKNOWN_ZONE_mixed", gene_states


# ---------- Collect data per cell ----------
rows = []
summary = {c: {"GENELESS":0, "UNMAPPABLE":0, "EXPRESSED":0, "UNKNOWN_ZONE":0,
               "OFF":0, "UNKNOWN_ZONE_mixed":0} for c in CELLS}
summary_flux = {c: {"GENELESS":0, "UNMAPPABLE":0, "EXPRESSED":0, "UNKNOWN_ZONE":0,
                    "OFF":0, "UNKNOWN_ZONE_mixed":0} for c in CELLS}

for cell in CELLS:
    sub_fp = R / f"subsets_v7/{cell}_subset_membership.csv"
    pool_fp = R / f"pooled_v7/{cell}_pooled.csv"
    sub = pd.read_csv(sub_fp)
    reps_df = sub[(sub.subsystem == SUBSYSTEM) & (sub.rxn == sub.representative)]
    reps = reps_df.representative.unique().tolist()
    pool = pd.read_csv(pool_fp).set_index("rxn")
    mdl = _load_model(V7 / f"trans_rfastcormics_{cell}.mat")
    rxn_to_j = {r: i for i, r in enumerate(mdl["rxns"])}

    print(f"\n=== {cell}  {SUBSYSTEM}  reps={len(reps)} ===")
    for rid in reps:
        if rid not in rxn_to_j:
            # extracted model may call it differently in rare edge cases
            continue
        j = rxn_to_j[rid]
        gpr = mdl["grrules"][j]
        gpr_tokens = _parse_gpr(gpr)
        verdict, gene_states = _classify_cell(gpr_tokens, cell)
        # Flux-carrying: nonzero in any sample column
        if rid in pool.index:
            flux_vals = pool.loc[rid].values
            carries_flux = bool((flux_vals != 0).any())
            frac_nonzero = float((flux_vals != 0).mean())
            mean_abs_flux = float(np.abs(flux_vals[flux_vals != 0]).mean()) if carries_flux else 0.0
        else:
            carries_flux = False; frac_nonzero = 0.0; mean_abs_flux = 0.0

        rows.append(dict(
            cell=cell, rxn=rid,
            rxn_name=mdl["rxn_names"][j] if j < len(mdl["rxn_names"]) else "",
            gpr=gpr, verdict=verdict,
            n_tokens=len(gpr_tokens),
            carries_flux=carries_flux,
            frac_nonzero=frac_nonzero,
            mean_abs_flux=mean_abs_flux,
            gene_states=json.dumps(gene_states),
        ))
        summary[cell][verdict] += 1
        if carries_flux:
            summary_flux[cell][verdict] += 1

df = pd.DataFrame(rows)
df.to_csv(OUT / "s179_leukotriene_120_audit.csv", index=False)
print(f"\n[write] {OUT/'s179_leukotriene_120_audit.csv'}")

# ---------- Summary tables ----------
def _pct_row(cnts):
    tot = sum(cnts.values())
    if tot == 0: return "n = 0"
    return "  ".join(f"{k}={v} ({100*v/tot:.0f}%)" for k, v in cnts.items() if v > 0)

print("\n--- all representatives in Leukotriene metabolism, per cell ---")
for c in CELLS:
    print(f"  {c:7s}  tot={sum(summary[c].values()):3d}  :: {_pct_row(summary[c])}")

print("\n--- flux-carrying representatives only (non-zero in ≥1 pooled sample) ---")
for c in CELLS:
    print(f"  {c:7s}  tot={sum(summary_flux[c].values()):3d}  :: {_pct_row(summary_flux[c])}")

# ---------- Markdown ----------
md = ["# s179 — Leukotriene metabolism GPR/zFPKM audit at 120 hpf",
      "",
      f"Scope: representatives of the enzyme-subset collapse for subsystem",
      f"**{SUBSYSTEM}** in each of **BL_120, D_120, LD_120**. Classified per",
      "cell by GPR state. Mirrors the audits used in §5.1 (gss / UDP-GlcA /",
      "PAPS) and §8.6 (ten single-failure + six one-succeeds).",
      "",
      "Two Tier A F1 rows motivated this audit:",
      "- `BL_120_vs_LD_120  Leukotriene metabolism`  W1 = 0.613  direction = +1  (BL > LD)",
      "- `D_120_vs_LD_120   Leukotriene metabolism`  W1 = 0.843  direction = +1  (D > LD)",
      "",
      "## Why the audit was run",
      "",
      "Direct zFPKM lookup at 120 hpf (see Claude Code message of 2026-10-14):",
      "`lta4h` is flat and EXPRESSED in all three regimes; `alox5a/5b.1-3/5ap`,",
      "`ltc4s`, `ltb4r` are all below detection. The flux is not travelling",
      "through the canonical synthesis step. This audit answers: of the",
      "representatives that *do* carry flux in the subsystem, how many are",
      "backed by expressed genes vs sitting in the four confound classes?",
      "",
      "## Per-cell GPR-state counts over all representatives",
      ""]
md.append("| cell | total | GENELESS | UNMAPPABLE | EXPRESSED | UNKNOWN_ZONE | OFF | UNKNOWN_ZONE_mixed |")
md.append("|---|---|---|---|---|---|---|---|")
for c in CELLS:
    row = summary[c]; tot = sum(row.values())
    def _p(k): return f"{row[k]} ({100*row[k]/tot:.0f} %)" if tot else "0"
    md.append(f"| {c} | {tot} | {_p('GENELESS')} | {_p('UNMAPPABLE')} | "
              f"{_p('EXPRESSED')} | {_p('UNKNOWN_ZONE')} | {_p('OFF')} | {_p('UNKNOWN_ZONE_mixed')} |")

md.append("")
md.append("## Flux-carrying representatives only")
md.append("")
md.append("A representative is 'flux-carrying' if its row in `pooled_v7/{cell}_pooled.csv`")
md.append("contains ≥ 1 non-zero sample. This is the set that actually contributes to")
md.append("the Tier A W1 statistic.")
md.append("")
md.append("| cell | total | GENELESS | UNMAPPABLE | EXPRESSED | UNKNOWN_ZONE | OFF | UNKNOWN_ZONE_mixed |")
md.append("|---|---|---|---|---|---|---|---|")
for c in CELLS:
    row = summary_flux[c]; tot = sum(row.values())
    def _p(k): return f"{row[k]} ({100*row[k]/tot:.0f} %)" if tot else "0"
    md.append(f"| {c} | {tot} | {_p('GENELESS')} | {_p('UNMAPPABLE')} | "
              f"{_p('EXPRESSED')} | {_p('UNKNOWN_ZONE')} | {_p('OFF')} | {_p('UNKNOWN_ZONE_mixed')} |")

# Verdict rule for the 120 hpf Leukotriene prediction
md.append("")
md.append("## Verdict rule (fixed before this audit was run)")
md.append("")
md.append("A pure flux result at 120 hpf is 'confound-free' only if, in all three")
md.append("cells, the share of **flux-carrying reps classified EXPRESSED** is ≥ 50 %")
md.append("and the share classified GENELESS + UNKNOWN_ZONE + UNMAPPABLE + OFF is < 50 %.")
md.append("")
md.append("If the share of flux-carrying reps that are NOT EXPRESSED exceeds 50 %")
md.append("in any of the three cells, the 120 hpf Leukotriene prediction comes out")
md.append("of the headline. The 48–96 hpf F2 rows (developmental program) are NOT")
md.append("touched unless this audit comes back confound; then they are audited next.")
md.append("")
# compute automatic verdict
confound_detected = False
notes = []
for c in CELLS:
    sf = summary_flux[c]; tot = sum(sf.values())
    if tot == 0:
        notes.append(f"  {c}: no flux-carrying reps")
        continue
    exp_pct = 100 * sf["EXPRESSED"] / tot
    confound_pct = 100 - exp_pct
    notes.append(f"  {c}: EXPRESSED = {exp_pct:.0f} %, confound (GENELESS+UNK+UNMAP+OFF) = {confound_pct:.0f} %")
    if confound_pct > 50:
        confound_detected = True
md += ["## Automatic verdict", ""] + notes + [""]
if confound_detected:
    md.append("**CONFOUND detected in at least one cell. The 120 hpf Leukotriene prediction")
    md.append("comes out of the headline. §4.3 should strike the leukotriene bullet from the")
    md.append("convergence narrative. Next action: run the same audit on the four F2 F2-F2")
    md.append("Leukotriene rows at 48→72 and 72→96 hpf.**")
else:
    md.append("**No confound: EXPRESSED is the majority class in every cell. The 120 hpf")
    md.append("Leukotriene prediction stays; it is driven by downstream metabolism, not")
    md.append("by the canonical synthesis step (lta4h is flat and alox5 is undetected).**")

(OUT / "s179_leukotriene_120_audit.md").write_text("\n".join(md))
print(f"[write] {OUT/'s179_leukotriene_120_audit.md'}")
print()
print("--- automatic verdict ---")
for n in notes: print(n)
print("confound_detected:", confound_detected)
