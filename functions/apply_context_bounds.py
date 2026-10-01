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
    # medium_uptake (class-1 vitamins/carotenoid, not in 105): lb=-1000.
    medium_uptake_set = set(led.loc[led.action == "medium_uptake",
                                    "exchange_MAR_id"].dropna())
    return mb, close_set, keep_set, medium_uptake_set


def apply_context_bounds(lb, ub, rxns, is_ex, cond, hpf,
                         frame="A", matched_biomax=None, repo=_REPO_DEFAULT,
                         bound_ex=None):
    """
    lb, ub, rxns, is_ex: numpy arrays / list from the loaded model.
    cond in {'BL','D','LD'}; hpf in {24,48,72,96,120}.
    frame='A' -> biomass_lb = 0.9 * biomass_max (caller applies via
       separate LP; this function returns lb/ub *without* the biomass
       floor set); pass frame='matched' with matched_biomax to set the
       matched-floor biomass_lb.
    Returns (lb_new, ub_new, idx, is_ex).
    """
    mb, close_set, keep_set, medium_uptake_set = _load_tables(repo)
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

    # -------- (2b) Free supply for unmeasured class-1 medium_uptake --------
    # 2026-10-01: medium_uptake rows (α/γ-tocotrienol, β-carotene, inositol)
    # not in the 105 measured set get lb=-1000 so unmeasured essentials
    # don't become default-±10 limiters.
    measured_set = set(mb.ex_rxn.unique())
    for rid in medium_uptake_set:
        if rid in idx and rid not in measured_set:
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
