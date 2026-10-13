"""
s150 -- Tier 0: topology × envelope decomposition under final v7 bounds.

Scope (plain FBA, Gurobi, no sampling):
  Part 1  parent (v7 base) vs child (extracted) biomass_max per cell under
          apply_context_bounds (frame A off; envelope + OCR + medium + ledger).
          15 + 15 = 30 LPs.
  Part 2  2x2 topology × envelope swap for each F1 pair at each hpf. 15 pairs
          × 4 corners = 60 LPs. Topology = extracted reaction set; envelope =
          105 fc-scaled bounds applied via apply_context_bounds(cond=env, hpf).
  Part 3  BL_96 subsystem add-back if the gap at hpf=96 is topology-dominated.

Outputs:
  grace_v4/results/tier0_v7/tier0_parent_child_v7.csv
  grace_v4/results/tier0_v7/tier0_2x2_v7.csv
  grace_v4/results/tier0_v7/dominance_per_hpf.txt
  grace_v4/results/tier0_v7/tier0_bl96_addback.csv       (iff topology-dominated)
"""
import sys, re
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
V7_BASE_DIR  = REPO / "ocr_anchored_extraction/models"
FC_TABLE     = REPO / "grace_v4/data/phase2_exchange_bounds_long.csv"
OUT          = REPO / "grace_v4/results/tier0_v7"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
CONDS = ("BL", "D", "LD")
HPFS  = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]
F1_PAIRS = [
    (f"{a}_{h}", f"{b}_{h}") for h in HPFS
    for a, b in [("BL","D"), ("BL","LD"), ("D","LD")]
]

MB = pd.read_csv(FC_TABLE)
MEASURED_SET = set(MB.ex_rxn.unique())
print(f"[init] {len(MEASURED_SET)} measured exchanges in fc table", flush=True)


# ---------- model loaders --------------------------------------------------
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
        try:
            subs = [_resolve(f, r) for r in np.asarray(m["subSystems"]).ravel()]
        except Exception:
            subs = [""] * len(rxns)
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        try: n_met_hdf = int(np.asarray(m["mets"]).shape[0])
        except Exception: n_met_hdf = 0
        n_met = max(int(ir.max()) + 1, n_met_hdf)
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc)-1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, S=S, lb=lb, ub=ub, subs=subs)


def load_base(hpf):
    return _load(V7_BASE_DIR / f"baked_ocr_{hpf}_v7.mat", "m")


def load_child(cell):
    return _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")


# ---------- FBA via Gurobi -------------------------------------------------
def fba_biomass(m, lb, ub, idx, get_duals=True):
    """Return (biomass_max, top_binding_rxn, top_dual_contribution)."""
    S = m["S"]
    n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    model = gp.Model("fba", env=env)
    v = model.addMVar(n_rxn, lb=lb, ub=ub, name="v")
    S_csr = S.tocsr()
    model.addMConstr(S_csr, v, "=", np.zeros(n_met))
    c_obj = np.zeros(n_rxn); c_obj[idx[BIOMASS_RXN]] = 1.0
    model.setMObjective(None, c_obj, 0.0, sense=GRB.MAXIMIZE)
    model.optimize()
    if model.Status != GRB.OPTIMAL:
        return np.nan, None, np.nan
    biomax = model.ObjVal
    if not get_duals:
        return biomax, None, np.nan

    # Dual contributions on variable bounds: for each var, RC × bound side.
    top_rxn = None; top_contrib = 0.0
    for j in range(n_rxn):
        if m["rxns"][j] == BIOMASS_RXN: continue
        rc = v[j].RC
        if abs(rc) < 1e-9: continue
        bound = lb[j] if rc > 0 else ub[j]
        contrib = abs(rc * bound)
        if contrib > abs(top_contrib):
            top_contrib = contrib
            top_rxn = m["rxns"][j]
    model.dispose(); env.dispose()
    return biomax, top_rxn, top_contrib


# ===========================================================================
# PART 1 — parent vs child, 15 + 15 LPs
# ===========================================================================
print("\n=== PART 1 — parent (v7 base) vs child (extracted) per cell ===", flush=True)
pc_rows = []
base_cache = None; last_hpf = -1
for cell in CELLS:
    cond, hpf = cell.split("_"); hpf = int(hpf)
    if hpf != last_hpf:
        base_cache = load_base(hpf); last_hpf = hpf

    # ---- parent ----
    lb0 = base_cache["lb"].copy(); ub0 = base_cache["ub"].copy()
    Snnz = np.diff(base_cache["S"].tocsc().indptr); is_ex_p = Snnz == 1
    lb, ub, idx_p, _ = apply_context_bounds(lb0, ub0, base_cache["rxns"], is_ex_p, cond, hpf)
    p_bm, p_rxn, p_ctrb = fba_biomass(base_cache, lb, ub, idx_p)

    # ---- child ----
    ch = load_child(cell)
    lb0c = ch["lb"].copy(); ub0c = ch["ub"].copy()
    Snnz_c = np.diff(ch["S"].tocsc().indptr); is_ex_c = Snnz_c == 1
    lbc, ubc, idx_c, _ = apply_context_bounds(lb0c, ub0c, ch["rxns"], is_ex_c, cond, hpf)
    c_bm, c_rxn, c_ctrb = fba_biomass(ch, lbc, ubc, idx_c)

    pc_rows.append(dict(
        cell=cell, hpf=hpf, condition=cond,
        parent_biomax=p_bm, parent_top_rxn=p_rxn, parent_top_contrib=p_ctrb,
        child_biomax=c_bm, child_top_rxn=c_rxn, child_top_contrib=c_ctrb,
        parent_over_child=p_bm / c_bm if c_bm else np.nan,
    ))
    print(f"  {cell}: parent={p_bm:.3f} ({p_rxn}, c={p_ctrb:.2f})  "
          f"child={c_bm:.3f} ({c_rxn}, c={c_ctrb:.2f})  ratio={p_bm/c_bm:.3f}",
          flush=True)

pc_df = pd.DataFrame(pc_rows)
pc_df.to_csv(OUT / "tier0_parent_child_v7.csv", index=False)
print(f"[write] tier0_parent_child_v7.csv  rows={len(pc_df)}", flush=True)


# ===========================================================================
# PART 2 — 2×2 topology × envelope, 15 pairs × 4 corners = 60 LPs
# ===========================================================================
print("\n=== PART 2 — 2×2 topology × envelope swap, 15 pairs × 4 corners ===", flush=True)
child_cache = {}
def get_child(cell):
    if cell not in child_cache:
        child_cache[cell] = load_child(cell)
    return child_cache[cell]

tt_rows = []
for cA, cB in F1_PAIRS:
    condA, hpf = cA.split("_"); hpf = int(hpf)
    condB = cB.split("_")[0]
    mA = get_child(cA); mB = get_child(cB)

    Snnz_A = np.diff(mA["S"].tocsc().indptr); is_ex_A = Snnz_A == 1
    Snnz_B = np.diff(mB["S"].tocsc().indptr); is_ex_B = Snnz_B == 1

    rxns_A = set(mA["rxns"]); rxns_B = set(mB["rxns"])
    unswap_A_only = (MEASURED_SET & rxns_A) - rxns_B
    unswap_B_only = (MEASURED_SET & rxns_B) - rxns_A
    n_unswap = len(unswap_A_only) + len(unswap_B_only)

    def corner(mX, is_ex_X, cond_env):
        lbX = mX["lb"].copy(); ubX = mX["ub"].copy()
        lb_, ub_, idx_, _ = apply_context_bounds(lbX, ubX, mX["rxns"], is_ex_X, cond_env, hpf)
        return fba_biomass(mX, lb_, ub_, idx_, get_duals=False)[0]

    aa = corner(mA, is_ex_A, condA)   # A top, A env
    ab = corner(mA, is_ex_A, condB)   # A top, B env
    ba = corner(mB, is_ex_B, condA)   # B top, A env
    bb = corner(mB, is_ex_B, condB)   # B top, B env

    Δ = aa - bb
    # Avg marginal topology effect = ((aa-ba) + (ab-bb)) / 2
    # Avg marginal envelope effect = ((aa-ab) + (ba-bb)) / 2
    # Non-additivity (interaction) = (aa - ab) - (ba - bb) = aa - ab - ba + bb
    top_eff = ((aa - ba) + (ab - bb)) / 2
    env_eff = ((aa - ab) + (ba - bb)) / 2
    nonadd  = (aa - ab) - (ba - bb)
    # Shares of |Δ| — top + env sum to Δ exactly in 2×2 (algebraic identity).
    denom = abs(Δ) if abs(Δ) > 1e-9 else 1.0
    top_share = top_eff / denom * 100
    env_share = env_eff / denom * 100
    nonadd_share = nonadd / denom * 100

    tt_rows.append(dict(
        hpf=hpf, pair=f"{cA}_vs_{cB}", cell_A=cA, cell_B=cB,
        bmax_A_top_A_env=aa, bmax_A_top_B_env=ab,
        bmax_B_top_A_env=ba, bmax_B_top_B_env=bb,
        gap_diag=Δ,
        topology_effect=top_eff, envelope_effect=env_eff, non_additivity=nonadd,
        topology_share_pct=top_share, envelope_share_pct=env_share,
        non_additivity_share_pct=nonadd_share,
        n_unswappable_exchanges=n_unswap,
        unswappable_A_only=";".join(sorted(unswap_A_only)),
        unswappable_B_only=";".join(sorted(unswap_B_only)),
    ))
    dominance = "topology" if abs(top_eff) >= abs(env_eff) else "envelope"
    print(f"  hpf={hpf:3d}  {cA} vs {cB}: AA={aa:.2f} AB={ab:.2f} BA={ba:.2f} BB={bb:.2f}  "
          f"Δ={Δ:+.2f}  top={top_eff:+.2f} ({top_share:+.0f}%)  env={env_eff:+.2f} ({env_share:+.0f}%)  "
          f"nonadd={nonadd:+.2f}  unswap={n_unswap}  → {dominance}", flush=True)

tt_df = pd.DataFrame(tt_rows)
tt_df.to_csv(OUT / "tier0_2x2_v7.csv", index=False)
print(f"[write] tier0_2x2_v7.csv  rows={len(tt_df)}", flush=True)


# ===========================================================================
# Dominance per hpf — one sentence per hpf
# ===========================================================================
dom_lines = []
for hpf in HPFS:
    g = tt_df[tt_df.hpf == hpf]
    top_sh = g.topology_share_pct.abs().mean()
    env_sh = g.envelope_share_pct.abs().mean()
    nadd_sh = g.non_additivity_share_pct.abs().mean()
    dom = "topology" if top_sh > env_sh else "envelope"
    dom_lines.append(
        f"hpf {hpf}: topology {top_sh:.0f}%, envelope {env_sh:.0f}%, non-additivity {nadd_sh:.0f}% — {dom}-dominated."
    )
with open(OUT / "dominance_per_hpf.txt", "w") as fh:
    fh.write("\n".join(dom_lines) + "\n")
for ln in dom_lines: print("  " + ln, flush=True)


# ===========================================================================
# PART 3 — BL_96 subsystem add-back (only if BL_96 pairs are topology-dominated)
# ===========================================================================
bl96 = tt_df[(tt_df.hpf == 96) & tt_df.pair.str.startswith("BL_96")]
bl96_topology_driven = bool((bl96.topology_effect.abs() > bl96.envelope_effect.abs()).all())
print(f"\n=== PART 3 — BL_96 add-back? topology_dominated={bl96_topology_driven} ===", flush=True)

if bl96_topology_driven:
    bl96_model = get_child("BL_96")
    ld96_model = get_child("LD_96")
    bl_rxns = set(bl96_model["rxns"])
    missing = [(r, s) for r, s in zip(ld96_model["rxns"], ld96_model["subs"]) if r not in bl_rxns]
    print(f"  reactions in LD_96 missing from BL_96: {len(missing)}", flush=True)
    by_sub = {}
    for r, s in missing:
        by_sub.setdefault(s, []).append(r)
    print(f"  subsystems contributing: {len(by_sub)}", flush=True)

    # Add-back via the v7 base at hpf=96 with a reaction mask: zero out everything
    # not in BL's reaction set ∪ this subsystem's LD-only extras. This way the S
    # matrix is consistent across all add-backs and matches the parent topology.
    base96 = load_base(96)
    base_rxns = base96["rxns"]
    Snnz_base = np.diff(base96["S"].tocsc().indptr); is_ex_base = Snnz_base == 1

    def mask_biomax(allowed_rxns):
        lb = base96["lb"].copy(); ub = base96["ub"].copy()
        lb_, ub_, idx_b, _ = apply_context_bounds(lb, ub, base_rxns, is_ex_base, "BL", 96)
        # Zero out reactions not in allowed set
        allowed = set(allowed_rxns)
        for i, r in enumerate(base_rxns):
            if r not in allowed:
                lb_[i] = 0.0; ub_[i] = 0.0
        return fba_biomass(base96, lb_, ub_, idx_b, get_duals=False)[0]

    bl_baseline_in_base = mask_biomax(bl_rxns)
    full_ceiling      = mask_biomax(bl_rxns.union(r for r, _ in missing))
    print(f"  baseline (BL_rxns via base, BL env): {bl_baseline_in_base:.4f} "
          f"(reference BL_96 child FBA was 8.535)", flush=True)
    print(f"  ceiling  (BL_rxns ∪ all LD-only extras): {full_ceiling:.4f}", flush=True)

    rows_bk = []
    for sub, extras in sorted(by_sub.items(), key=lambda x: -len(x[1])):
        bm = mask_biomax(bl_rxns.union(extras))
        gain = bm - bl_baseline_in_base
        rows_bk.append(dict(
            subsystem=sub, n_extras=len(extras),
            biomax_bl_plus_subsystem=bm, baseline=bl_baseline_in_base,
            ceiling_full_addback=full_ceiling,
            gain_vs_baseline=gain,
            frac_of_full_recovery=(bm - bl_baseline_in_base) / max(full_ceiling - bl_baseline_in_base, 1e-9),
        ))
        print(f"  +{sub!r}: n_extras={len(extras):3d}  biomax={bm:.4f}  gain={gain:+.4f}",
              flush=True)
    bk_df = pd.DataFrame(rows_bk).sort_values("gain_vs_baseline", ascending=False)
    bk_df.to_csv(OUT / "tier0_bl96_addback.csv", index=False)
    print(f"[write] tier0_bl96_addback.csv  rows={len(bk_df)}", flush=True)
else:
    # Report the binding exchange and its fc ratio
    print("  envelope-dominated at BL_96; naming top binding exchange from Part 1.", flush=True)
    bl96_pc = pc_df[pc_df.cell == "BL_96"].iloc[0]
    print(f"  BL_96 top binding (child): {bl96_pc.child_top_rxn} contrib={bl96_pc.child_top_contrib:.3f}",
          flush=True)
    rxn_id = bl96_pc.child_top_rxn
    sel = MB[(MB.ex_rxn == rxn_id) & (MB.hpf == 96)]
    if len(sel):
        fcs = sel.set_index("condition")["fc"].to_dict()
        if "BL" in fcs and "LD" in fcs:
            print(f"  {rxn_id} fc BL/LD = {fcs['BL']:.4f} / {fcs['LD']:.4f}  "
                  f"= {fcs['BL']/fcs['LD']:.3f}×", flush=True)

print("\n=== done ===", flush=True)
