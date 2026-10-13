"""
apply_context_bounds -- canonical per-cell bound application.

Starts from the model's stored bounds (v7 boundary: Part 1 secretion
uncapped, Part 2 unmeasured uptake closed, ledger closes, medium open),
NEVER resets exchanges to a uniform default. Applies:

  (1) lb = -BOUND_EX × fc, ub = +BOUND_EX × fc on the 105 measured
      exchanges (fc from phase2_exchange_bounds_long.csv, fc=1.0 fallback
      when missing/non-positive)
  (2) OCR cap on MAR09048
  (3) Frame A: biomass_lb = 0.9 × biomass_max (or 'matched' for matched-floor)

BOUND_EX = 10 by convention; envelope is a relative scale centred on
fc = 1.0; only cross-cell ratios are interpretable.

Asserts (post-application):
  - every ledger close_exchange row has lb == ub == 0
  - every unmeasured non-medium exchange has lb == 0 AND ub == 1000
  - every medium row has lb < 0
  - exactly 105 exchanges carry fc scaling (in the ones present in the
    extracted model — pruned ones don't count)

Raises AssertionError if any check fails.
"""
import numpy as np
import pandas as pd
from pathlib import Path

BOUND_EX  = 10.0
O2_EX     = "MAR09048"
OCR = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}

# 17 medium exchanges (bulk inorganics + water/gases). Hardcoded, matches s135.
MEDIUM = frozenset({
    "MAR09047","MAR09058","MAR09072","MAR09073","MAR09074","MAR09076",
    "MAR09077","MAR09078","MAR09079","MAR09080","MAR09081","MAR09082",
    "MAR09148","MAR09150","MAR13066","MAR13072","MAR13073",
})

_REPO_DEFAULT = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")


def _find_file(candidates):
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(f"none of these paths exist: {candidates}")


def _load_tables(repo=_REPO_DEFAULT):
    import os
    # Allow override via env var (handy on Grace: V7_REPO=$SCRATCH/grace_v4)
    repo = Path(os.environ.get("V7_REPO", repo))
    # Phase 2 bounds CSV: try repo-local then Grace-local layouts
    bounds_csv = _find_file([
        repo / "results/qc/phase2_exchange_bounds_long.csv",
        repo / "data/phase2_exchange_bounds_long.csv",
        repo / "grace_v4/data/phase2_exchange_bounds_long.csv",
    ])
    # Ledger: try repo-local then grace_v4-local layouts
    ledger_csv = _find_file([
        repo / "reviewer_packet_v4_slim/s130_v7_ledger.csv",
        repo / "s130_v7_ledger.csv",
    ])
    mb = pd.read_csv(bounds_csv)
    led = pd.read_csv(ledger_csv)
    close_set = set(led.loc[led.action == "close_exchange", "exchange_MAR_id"].dropna())
    # non-close ledger rows (keep_exchange, medium_uptake, ysl_lipid_delivery)
    # stay at stored default (±10 usually); NOT Part 1+2 rows.
    keep_set = set(led.loc[led.action.isin(
        ["keep_exchange","medium_uptake","ysl_lipid_delivery"]),
        "exchange_MAR_id"].dropna())
    # medium_uptake rows: lb comes from per-cell demand bounds (see
    # _load_demand_bounds), not a fixed -1000. See §11.19.
    medium_uptake_set = set(led.loc[led.action == "medium_uptake",
                                    "exchange_MAR_id"].dropna())
    return mb, close_set, keep_set, medium_uptake_set


def _load_demand_bounds(repo=_REPO_DEFAULT):
    """
    Per-cell demand bounds for the 4 free-supplied class-1 rows (inositol,
    β-carotene, α/γ-tocotrienols). Written by s266_compute_demand_bounds.py.
    Returns {(cond, hpf, mar): lb_corrected}.  See §11.19.
    """
    import os
    repo = Path(os.environ.get("V7_REPO", repo))
    for candidate in [
        repo / "grace_v4/results/qc/medium_uptake_demand_bounds.csv",
        repo / "results/qc/medium_uptake_demand_bounds.csv",
        repo / "medium_uptake_demand_bounds.csv",
    ]:
        if candidate.exists():
            df = pd.read_csv(candidate)
            out = {}
            for _, r in df.iterrows():
                if np.isfinite(r.lb_corrected):
                    out[(str(r.cond), int(r.hpf), str(r.mar))] = float(r.lb_corrected)
            return out
    # If the file is missing, callers get a loud failure in apply_context_bounds
    return None


def _run_energy_gate(S, lb, ub, idx, rxns, cond, hpf):
    """
    §11.19 energy-source gate.
    All 4 free-supplied rows (medium_uptake) are kept at their demand bounds
    (already in lb). MEDIUM at -1000. Glucose at -1. O2 at -1000.
    All other exchanges closed for uptake. Max ATP (MAR03964). Assert the
    result stays within 31.5 + 0.5 of textbook.
    """
    from scipy.optimize import linprog
    import scipy.sparse as sp
    n_met, n_rxn = S.shape
    is_ex = (np.diff(S.tocsc().indptr) == 1) if sp.issparse(S) else None
    lb2 = lb.copy(); ub2 = ub.copy()
    # Close every exchange for uptake
    for j in range(n_rxn):
        if is_ex is not None and not is_ex[j]: continue
        lb2[j] = 0.0
    # Open MEDIUM
    for rid in MEDIUM:
        if rid in idx: lb2[idx[rid]] = -1000.0
    # Open free-supplied at their demand bounds (preserve from lb)
    demand = _load_demand_bounds()
    if demand is not None:
        for (c, h, rid), lbv in demand.items():
            if c == cond and h == hpf and rid in idx:
                lb2[idx[rid]] = lbv
    # Glucose primer
    if "MAR09034" in idx:
        lb2[idx["MAR09034"]] = -1.0; ub2[idx["MAR09034"]] = -1.0
    # O2 at -1000
    if O2_EX in idx:
        lb2[idx[O2_EX]] = -1000.0; ub2[idx[O2_EX]] = 0.0
    # Biomass free, ATP objective ub large
    if "MAR00021" in idx:
        j = idx["MAR00021"]; lb2[j] = 0.0; ub2[j] = max(ub2[j], 1000.0)
    if "MAR03964" in idx:
        ub2[idx["MAR03964"]] = 1e9
    c_obj = np.zeros(n_rxn); c_obj[idx["MAR03964"]] = -1.0  # minimise -ATP
    A_eq = S.toarray() if sp.issparse(S) else np.asarray(S)
    res = linprog(c_obj, A_eq=A_eq, b_eq=np.zeros(n_met),
                  bounds=list(zip(lb2, ub2)), method="highs")
    atp = -res.fun if res.success else float("inf")
    assert atp - 31.5 < 0.5, (
        f"ENERGY_GATE FAIL for cond={cond} hpf={hpf}: ATP = {atp:.3f} > 32.0 "
        f"(31.5 + 0.5 slack). Free-supplied rows are producing >textbook "
        f"ATP. See §11.19."
    )
    return atp


def apply_context_bounds(lb, ub, rxns, is_ex, cond, hpf,
                         frame="A", matched_biomax=None, repo=_REPO_DEFAULT,
                         bound_ex=None, S=None, skip_energy_gate=False):
    """
    lb, ub, rxns, is_ex: numpy arrays / list from the loaded model.
    cond in {'BL','D','LD'}; hpf in {24,48,72,96,120}.
    frame='A' -> biomass_lb = 0.9 * biomass_max (caller applies via
       separate LP; this function returns lb/ub *without* the biomass
       floor set); pass frame='matched' with matched_biomax to set the
       matched-floor biomass_lb.
    S: stoichiometric matrix (scipy.sparse or dense ndarray). If provided,
       the §11.19 energy-source gate is run as a post-condition assertion.
       If None, gate is skipped silently (legacy callers).
    skip_energy_gate: emergency bypass (set True only for diagnostic runs
       that intentionally violate the gate; e.g. the medium-energy audit).
    Returns (lb_new, ub_new, idx, is_ex).
    """
    mb, close_set, keep_set, medium_uptake_set = _load_tables(repo)
    demand_lb = _load_demand_bounds(repo)
    lb = np.asarray(lb, dtype=float).copy()
    ub = np.asarray(ub, dtype=float).copy()
    idx = {r: i for i, r in enumerate(rxns)}

    # -------- (1) lb=-BOUND_EX×fc, ub=+BOUND_EX×fc on the 105 measured -----
    bex = bound_ex if bound_ex is not None else BOUND_EX
    n_scaled = 0
    sel = (mb.condition == cond) & (mb.hpf == hpf)
    measured_ex_all = set(mb.ex_rxn.unique())
    for rid, fc in zip(mb.loc[sel, "ex_rxn"], mb.loc[sel, "fc"]):
        if rid in idx:
            j = idx[rid]
            fc = float(fc) if np.isfinite(fc) and fc > 0 else 1.0
            lb[j] = -bex * fc
            ub[j] = +bex * fc
            n_scaled += 1

    # -------- (2) OCR cap on MAR09048 --------------------------------------
    if O2_EX in idx:
        j = idx[O2_EX]
        lb[j] = -OCR[hpf]; ub[j] = 0.0

    # -------- (2b) Per-cell demand bounds on medium_uptake rows -----------
    # §11.19 (2026-10-12): medium_uptake rows (α/γ-tocotrienol, β-carotene,
    # inositol) previously opened at lb=-1000. The energy-source gate
    # (s261/s264) showed β-carotene MAR09276 generated ~971 ATP excess over
    # textbook glucose in every extracted cell, and inositol a further
    # ~1104 in the base. The fix: lb_i = -1.1 × |FVA_min(v_i | biomass_max)|,
    # floor at -1e-3. See grace_v4/results/qc/medium_uptake_demand_bounds.csv.
    measured_set = set(mb.ex_rxn.unique())
    for rid in medium_uptake_set:
        if rid in idx and rid not in measured_set:
            key = (cond, hpf, rid)
            if demand_lb is not None and key in demand_lb:
                lb[idx[rid]] = demand_lb[key]
            elif demand_lb is None:
                # Table missing — raise clearly rather than silently default
                raise FileNotFoundError(
                    "medium_uptake_demand_bounds.csv not found; run "
                    "grace_v4/scripts/s266_compute_demand_bounds.py first "
                    "or set V7_REPO to the correct path."
                )
            else:
                # Table exists but this (cond, hpf, mar) not computed —
                # fall back to the pre-§11.19 default so legacy cells still
                # work. Warn so callers can see.
                import warnings
                warnings.warn(f"no demand bound for ({cond}, {hpf}, {rid}); "
                              "falling back to lb=-1000", RuntimeWarning)
                lb[idx[rid]] = -1000.0

    # -------- assertions ---------------------------------------------------
    # a) ledger closes are still closed (or pruned from model)
    for rid in close_set:
        if rid in idx:
            j = idx[rid]
            assert lb[j] == 0 and ub[j] == 0, \
                f"ASSERTION: ledger close {rid} not closed: lb={lb[j]} ub={ub[j]}"
    # b) unmeasured non-medium non-ledger exchanges (Part 1+2 rows) at lb=0, ub=1000
    n_p12 = 0
    for j in np.where(is_ex)[0]:
        rid = rxns[j]
        if rid == O2_EX: continue
        if rid in close_set: continue
        if rid in keep_set: continue   # ledger keep_exchange/medium_uptake/ysl_lipid_delivery
        if rid in MEDIUM: continue
        if rid in measured_ex_all: continue
        # Genuine Part 1+2 row: unmeasured, not on any ledger row.
        assert lb[j] == 0.0 and ub[j] == 1000.0, \
            f"ASSERTION: Part 1+2 row {rid} bounds violated: lb={lb[j]} ub={ub[j]}"
        n_p12 += 1
    # c) medium rows have lb < 0
    for rid in MEDIUM:
        if rid in idx:
            j = idx[rid]
            assert lb[j] < 0, f"ASSERTION: medium {rid} lb not < 0: {lb[j]}"
    # d) exactly (# measured present in model) exchanges got flux scaling
    n_measured_present = sum(1 for rid in measured_ex_all if rid in idx)
    assert n_scaled == n_measured_present, \
        f"ASSERTION: flux scaling applied to {n_scaled}, expected {n_measured_present} present in model (of 105 total)"

    # e) §11.19 energy-source gate (post-condition). Only run if S provided.
    if S is not None and not skip_energy_gate:
        _run_energy_gate(S, lb, ub, idx, rxns, cond, hpf)

    return lb, ub, idx, is_ex


def compute_biomax(m, lb, ub, idx, biomass_rxn="MAR00021"):
    """FBA for biomass, returns (biomax, res)."""
    from scipy.optimize import linprog
    import scipy.sparse as sp
    n_rxn = len(m["rxns"]); n_met = m["S"].shape[0]
    c = np.zeros(n_rxn); c[idx[biomass_rxn]] = -1.0
    A_eq = m["S"].toarray() if sp.issparse(m["S"]) else np.asarray(m["S"])
    b_eq = np.zeros(n_met)
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=list(zip(lb, ub)), method="highs")
    return (-res.fun if res.success else np.nan), res


def set_biomass_floor(lb, ub, idx, biomax, frac=0.9, biomass_rxn="MAR00021"):
    """Set biomass_lb = frac * biomax (Frame A)."""
    j = idx[biomass_rxn]
    lb = lb.copy()
    lb[j] = frac * biomax
    return lb, ub
