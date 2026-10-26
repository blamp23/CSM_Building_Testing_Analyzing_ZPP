"""
Medium energy-source gate test.

For every ledger row with action in {medium_uptake, keep_exchange,
ysl_lipid_delivery} that is marked class 1 (free-supplied essentials):
on the v7 base, open ONLY that exchange + O2 + water + protons + Pi +
inorganic ions, maximise MAR03964. A free-supplied essential should not
generate ATP by itself. If max ATP > 0.5 the row fails.

Writes:
  medium_energy_gate.csv — one row per class-1 ledger entry: {exchange,
    metabolite, action, carbon, max_atp, passes_gate}
  methods_section_11_19.md — the §11.19 log entry.

Does NOT modify apply_context_bounds.py or the ledger. Reports what must
change.
"""
import sys, re
from pathlib import Path
import numpy as np, pandas as pd
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from s250_atp_sweep import load_model, V7_BASE_DIR, ATP_RXN, O2_EX, MEDIUM

OUT = REPO / "grace_v4/results/atp_sweep_v7"
LEDGER = REPO / "reviewer_packet_v4_slim/s130_v7_ledger.csv"

led = pd.read_csv(LEDGER)
class1 = led[led["class"] == 1].copy()
print(f"[info] class-1 ledger rows: {len(class1)}", flush=True)
print(class1[["metabolite","MAM_id","exchange_MAR_id","action"]].to_string())

# Load base_120 (any base works, biochemistry is identical)
base = load_model(V7_BASE_DIR / "baked_ocr_120_v7.mat", "m")
idx = {r: i for i, r in enumerate(base["rxns"])}

def gate_atp(target_mar):
    """ATP max when only {target_mar + MEDIUM} are open for uptake."""
    S = base["S"]; n_met, n_rxn = S.shape
    is_ex = (np.diff(S.indptr) == 1)
    lb = np.zeros(n_rxn)
    ub = base["ub"].copy()
    # Close every exchange, then selectively open:
    for j in np.where(is_ex)[0]:
        lb[j] = 0.0; ub[j] = 1000.0  # secretion open, uptake closed
    # Open MEDIUM for uptake
    for rid in MEDIUM:
        if rid in idx:
            lb[idx[rid]] = -1000.0
    # Open target
    if target_mar in idx:
        lb[idx[target_mar]] = -1000.0
    # Biomass unconstrained
    if "MAR00021" in idx:
        j = idx["MAR00021"]; lb[j] = 0.0; ub[j] = 1000.0
    # ATP ub lifted to detect any inflation
    ub[idx[ATP_RXN]] = 1e9
    # Also open internal reaction bounds as stored
    for j in range(n_rxn):
        if not is_ex[j]:
            lb[j] = base["lb"][j]; ub[j] = base["ub"][j]
    # Solve
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    import scipy.sparse as sp
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val

def count_C(formula):
    if not isinstance(formula, str) or not formula: return 0
    m = re.search(r"C(\d*)(?![a-z])", formula)
    if not m: return 0
    return int(m.group(1)) if m.group(1) else 1

# Control: nothing opened beyond MEDIUM
print(f"\n[control] gate ATP with no extra uptake: {gate_atp(''):.4f}", flush=True)

# Get metabolite formula map for each MAM
mam_to_formula = dict(zip(base["mets"], base["metFormulas"]))
mam_to_name    = dict(zip(base["mets"], base["metNames"]))

def gate_atp_with_glucose_primer(target_mar):
    """Variant: like gate_atp but also open glucose at lb=-1 (forced uptake)
    to prime central metabolism. Catches class-1 cofactor+primer loops."""
    atp = gate_atp(target_mar)  # fallback if target not in idx
    S = base["S"]; n_met, n_rxn = S.shape
    is_ex = (np.diff(S.indptr) == 1)
    lb = np.zeros(n_rxn); ub = base["ub"].copy()
    for j in np.where(is_ex)[0]:
        lb[j] = 0.0; ub[j] = 1000.0
    for rid in MEDIUM:
        if rid in idx: lb[idx[rid]] = -1000.0
    if target_mar in idx:
        lb[idx[target_mar]] = -1000.0
    # Glucose primer at exactly 1 unit
    if "MAR09034" in idx:
        lb[idx["MAR09034"]] = -1.0; ub[idx["MAR09034"]] = -1.0
    if "MAR00021" in idx:
        j = idx["MAR00021"]; lb[j] = 0.0; ub[j] = 1000.0
    ub[idx[ATP_RXN]] = 1e9
    for j in range(n_rxn):
        if not is_ex[j]:
            lb[j] = base["lb"][j]; ub[j] = base["ub"][j]
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("", env=env)
    v = gm.addMVar(n_rxn, lb=lb, ub=ub)
    import scipy.sparse as sp
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    val = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return val

# Primer baseline: glucose -1 + MEDIUM, no extra uptake. Textbook ~31.5.
primer_baseline = gate_atp_with_glucose_primer("")
print(f"[primer baseline] gate ATP with glucose primer only: {primer_baseline:.4f}  "
      f"(threshold for 'this row inflates past glucose textbook' = primer + 0.5)", flush=True)

rows = []
for _, r in class1.iterrows():
    mar = r.exchange_MAR_id
    mam = r.MAM_id
    if not isinstance(mar, str) or not mar or mar not in idx:
        rows.append(dict(exchange=mar, metabolite=r.metabolite,
                         action=r.action, formula="", carbon=0,
                         max_atp=np.nan, max_atp_glucose_primer=np.nan,
                         passes_gate=False, passes_gate_primer=False,
                         note="exchange not in base model"))
        continue
    mam_e = mam + "e" if not mam.endswith("e") else mam
    formula = mam_to_formula.get(mam_e, "")
    carbon = count_C(formula)
    atp_alone  = gate_atp(mar)
    atp_primed = gate_atp_with_glucose_primer(mar)
    passes       = (atp_alone  <= 0.5)
    # Fail primed gate if adding this row produces more ATP than glucose-alone
    # plus a 0.5 slop.
    passes_primed = (atp_primed <= primer_baseline + 0.5)
    print(f"  {r.metabolite:25s} {mar} {r.action:20s} C={carbon:3d} "
          f"ATP_alone={atp_alone:9.3f} {'PASS' if passes else 'FAIL'}   "
          f"ATP_primed={atp_primed:9.3f} {'PASS' if passes_primed else 'FAIL'}",
          flush=True)
    rows.append(dict(exchange=mar, metabolite=r.metabolite,
                     action=r.action, formula=formula, carbon=carbon,
                     max_atp=atp_alone, passes_gate=passes,
                     max_atp_glucose_primer=atp_primed,
                     passes_gate_primer=passes_primed,
                     excess_over_primer=atp_primed - primer_baseline))

df = pd.DataFrame(rows)
df.to_csv(OUT / "medium_energy_gate.csv", index=False)
print(f"\n[write] medium_energy_gate.csv  rows={len(df)}", flush=True)
print(f"  passes: {int(df.passes_gate.sum())}  fails: {int((~df.passes_gate).sum())}")

# Write §11.19 methods log
fails = df[~df.passes_gate].sort_values("max_atp", ascending=False)
md = ["# §11.19  Medium energy-source gate",
      "",
      "**Date:** 2026-10-12",
      "",
      "## Test",
      "",
      "On the v7 base at 120 hpf, open only one target exchange plus the 17 "
      "MEDIUM inorganics (O₂, H₂O, H⁺, Pi, Na⁺, K⁺, Cl⁻, HCO₃⁻, SO₄²⁻, "
      "Ca²⁺, Mg²⁺, Zn²⁺, Cu²⁺, Fe²⁺, iodide, CO₂, carbonate). All other "
      "exchanges closed for uptake. Biomass unconstrained. Maximise "
      "MAR03964 with ATP upper bound lifted to 10⁹. Any ledger row (any "
      "class-1 free-supplied essential) that produces max ATP > 0.5 fails.",
      "",
      "## Rationale",
      "",
      "A free-supplied essential should not be a carbon/energy donor. "
      "`apply_context_bounds.medium_uptake_set` opens these exchanges at "
      "`lb = -1000` so that unmeasured essentials don't become bound "
      "limiters, but if an exchange's metabolite is organic (C-containing) "
      "the LP will use it as a cheap carbon source. The ATP-sweep "
      "diagnosis (s255/s256, 2026-10-12) traced the base_120 glucose "
      "2 107-ATP result to β-carotene (C₄₀H₅₆) and inositol (C₆H₁₂O₆) "
      "being used in exactly this way under `lb = -1000`.",
      "",
      "## Result",
      f"",
      f"| exchange | metabolite | action | formula | C | max_ATP | passes |",
      f"|---|---|---|---|---:|---:|:---:|",
      ]
for _, r in df.sort_values("max_atp", ascending=False).iterrows():
    md.append(f"| {r.exchange} | {r.metabolite} | {r.action} | {r.formula} | {r.carbon} | {r.max_atp:.2f} | {'✓' if r.passes_gate else '✗'} |")

md += [
    "",
    "## Fix (not yet applied — awaiting review of impact in C3/C4)",
    "",
    "For each failing row:",
    "",
    "1. In `functions/apply_context_bounds.py` and "
    "`grace_upload_f2_floor_resample/staging/functions/apply_context_bounds.py` "
    "(both twins), change the `medium_uptake` free-supply "
    "default from `lb = -1000` to `lb = -10`. Rationale: `-10` matches "
    "`BOUND_EX`, the default upper on measured exchanges — trace supply "
    "available but not an unlimited pool.",
    "2. In `reviewer_packet_v4_slim/s130_v7_ledger.csv`, add a column "
    "`default_lb` = -10 for the failing medium_uptake rows (-1000 for the "
    "passing ones remains valid).",
    "3. Add a permanent post-application assertion in "
    "`apply_context_bounds`: for each row in `medium_uptake_set` and "
    "`keep_exchange_set`, max ATP under the gate test above must be ≤ 0.5.",
    "",
    "## Known impact (to confirm in C3/C4/C5)",
    "",
    "In the current `LD_72` hopsy samples, inositol uptake (MAR09361) is "
    "~−11.4 mmol/gDW/hr at every draw (comparable to measured glucose). "
    "β-carotene (MAR09276) median uptake is ~−0.13. Tocotrienols are "
    "trace (<0.01). This likely means the sampled polytope has been using "
    "class-1 vitamin carbon in central metabolism; subsystem comparisons "
    "across LD/BL/D could be partly reflecting differences in how much "
    "each cell dips into this pool. The §11.19 fix would prevent that.",
]
md_fp = OUT / "methods_section_11_19.md"
md_fp.write_text("\n".join(md))
print(f"[write] methods_section_11_19.md")
