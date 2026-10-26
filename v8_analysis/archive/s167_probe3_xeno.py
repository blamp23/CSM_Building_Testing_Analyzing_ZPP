"""
s167 -- Probe 3: generic xenobiotic burden.

Scratch model only. Add to each cell identically:
  X_ex:    ∅ -> X_c                     (exchange, lb=-u, ub=0; UPTAKE forced)
  Phase1:  X_c + NADPH + O2 + H+ -> X-OH + NADP+ + H2O
  Phase2a: X-OH + UDP-GlcA -> X-gluc + UDP          [glucuronidation variant]
  Phase2b: X-OH + GSH -> X-SG                       [GSH-conjugation variant]
  Exp_a:   X-gluc -> ∅
  Exp_b:   X-SG -> ∅

Variant 3a uses Phase2a + Exp_a. Variant 3b uses Phase2b + Exp_b. Each gives
its own set of scratch reactions appended to the base S.

For each cell × variant × u in {0, 0.5, 1, 2, 5}: biomass_max(u). Slope =
(biomax(u=0) - biomax(u=5)) / 5  (cost per unit load).

Output: probe3_xeno_burden.csv
"""
import sys, time
from pathlib import Path
import h5py, numpy as np, pandas as pd
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds, compute_biomax

V7_EXTRACTED = REPO / "reviewer_packet_v4_slim/v7_models"
OUT = REPO / "grace_v4/results/v7_susceptibility"
OUT.mkdir(parents=True, exist_ok=True)

BIOMASS_RXN = "MAR00021"
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = [f"{c}_{h}" for c in CONDS for h in HPFS]

MET = {
    "NADPH":    "MAM02555c",
    "NADP":     "MAM02554c",
    "O2":       "MAM02630c",
    "H":        "MAM02039c",
    "H2O":      "MAM02040c",
    "UDP_GlcA": "MAM03109c",
    "UDP":      "MAM03106c",
    "GSH":      "MAM02026c",
}

U_LIST = (0.0, 0.5, 1.0, 2.0, 5.0)


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
        try: mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        except Exception: mets = []
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, S=S, lb=lb, ub=ub)


def _add_scratch(mdl, variant):
    """Return (S_aug_csr, lb_aug, ub_aug, scratch_names, scratch_js) where
    scratch reactions are appended as columns. Variant 'a' uses glucuronide,
    'b' uses GSH-conjugate."""
    S = mdl["S"].tocsc()
    n_met, n_rxn = S.shape
    mets = list(mdl["mets"])

    # Add scratch metabolites: X_c, X-OH_c, X-gluc_c or X-SG_c
    scratch_mets = ["Xparent_c", "XOH_c"]
    if variant == "a":
        scratch_mets.append("Xgluc_c")
    else:
        scratch_mets.append("XSG_c")
    met_idx = {m: i for i, m in enumerate(mets)}
    for sm in scratch_mets:
        met_idx[sm] = len(mets)
        mets.append(sm)
    n_met_new = len(mets)

    # Prepare dense new rows for existing columns = zero
    pad_rows = n_met_new - n_met
    S_padded = sp.vstack([S, sp.csc_matrix((pad_rows, n_rxn))], format="csc")

    # Build scratch columns
    def _col(stoich_dict):
        """stoich_dict: {met_name_or_id: coef}; met_name_or_id uses either
        the met_idx keys above or raw MAM IDs."""
        rows = []; vals = []
        for k, c in stoich_dict.items():
            # Resolve: scratch name or MAM id
            if k in met_idx:
                rows.append(met_idx[k])
            elif k in MET and MET[k] in met_idx:
                rows.append(met_idx[MET[k]])
            else:
                raise KeyError(k)
            vals.append(c)
        return sp.csc_matrix((vals, (rows, [0] * len(rows))), shape=(n_met_new, 1))

    scratch_cols = []
    scratch_names = []
    # 1. X_ex: ∅ -> X_c
    scratch_cols.append(_col({"Xparent_c": +1.0}))
    scratch_names.append("Xparent_ex")
    # 2. Phase1: X + NADPH + O2 + H -> X-OH + NADP + H2O
    scratch_cols.append(_col({
        "Xparent_c": -1.0, "NADPH": -1.0, "O2": -1.0, "H": -1.0,
        "XOH_c": +1.0, "NADP": +1.0, "H2O": +1.0,
    }))
    scratch_names.append("Phase1")
    if variant == "a":
        # 3a. Phase2a: X-OH + UDP-GlcA -> X-gluc + UDP
        scratch_cols.append(_col({
            "XOH_c": -1.0, "UDP_GlcA": -1.0,
            "Xgluc_c": +1.0, "UDP": +1.0,
        }))
        scratch_names.append("Phase2a")
        # 4a. Export X-gluc -> ∅
        scratch_cols.append(_col({"Xgluc_c": -1.0}))
        scratch_names.append("Export_a")
    else:
        # 3b. Phase2b: X-OH + GSH -> X-SG
        scratch_cols.append(_col({"XOH_c": -1.0, "GSH": -1.0, "XSG_c": +1.0}))
        scratch_names.append("Phase2b")
        # 4b. Export X-SG -> ∅
        scratch_cols.append(_col({"XSG_c": -1.0}))
        scratch_names.append("Export_b")

    S_aug = sp.hstack([S_padded] + scratch_cols, format="csr")
    k = len(scratch_cols)
    scratch_js = list(range(n_rxn, n_rxn + k))
    # Default bounds for scratch rxns: internal reactions ±1000; X_ex later pinned to -u
    lb_aug = np.concatenate([mdl["lb"], np.array([-1000.0] + [0.0] * (k - 1))])
    ub_aug = np.concatenate([mdl["ub"], np.array([0.0] + [1000.0] * (k - 1))])
    # Phase1 forward only (irreversible CYP)
    # already 0..1000 by default; good
    return S_aug, lb_aug, ub_aug, scratch_names, scratch_js


def _biomass_at(mdl, S_aug, lb, ub, idx_bio, u, scratch_js, scratch_names):
    """Pin X_ex to -u (uptake of u units), maximise biomass. Return biomass."""
    lb2 = lb.copy(); ub2 = ub.copy()
    xex_j = scratch_js[scratch_names.index("Xparent_ex")]
    # X_ex: lb=-u, ub=-u → forced uptake = u (since +1 in S means production of X_c,
    # positive flux produces X_c; to force uptake we pin flux = u)
    lb2[xex_j] = u
    ub2[xex_j] = u
    n_rxn = S_aug.shape[1]
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(n_rxn, lb=lb2, ub=ub2)
    gm.addMConstr(S_aug, v, "=", np.zeros(S_aug.shape[0]))
    obj = np.zeros(n_rxn); obj[idx_bio] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    bm = gm.ObjVal if gm.Status == GRB.OPTIMAL else np.nan
    gm.dispose(); env.dispose()
    return bm


rows = []
t0 = time.time()
for cell in CELLS:
    cond, h = cell.split("_"); hpf = int(h)
    mdl = _load(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex, cond, hpf,
        S=mdl["S"],
    )
    non_ex = ~is_ex
    lb[(lb <= -999.9) & non_ex] = -1e6
    ub[(ub >= 999.9) & non_ex] = +1e6
    for variant, label in [("a", "Phase2_glucuronide"), ("b", "Phase2_GSH")]:
        S_aug, lb_aug, ub_aug, sc_names, sc_js = _add_scratch(
            dict(mdl, lb=lb, ub=ub), variant)
        bm_by_u = {}
        for u in U_LIST:
            bm = _biomass_at(mdl, S_aug, lb_aug, ub_aug, idx[BIOMASS_RXN],
                             u, sc_js, sc_names)
            bm_by_u[u] = bm
            rows.append(dict(cell=cell, condition=cond, hpf=hpf,
                             variant=label, u=u, biomass_max=bm))
        # Slope = (bm(0) - bm(5))/5
        slope = (bm_by_u[0.0] - bm_by_u[5.0]) / 5.0 if np.isfinite(bm_by_u[0.0]) and np.isfinite(bm_by_u[5.0]) else np.nan
        rows.append(dict(cell=cell, condition=cond, hpf=hpf,
                         variant=label, u=np.nan, biomass_max=np.nan,
                         slope=slope))
        dt = time.time() - t0
        print(f"  {cell:7s}  {label:20s}  bm[0]={bm_by_u[0.0]:.3f}  bm[5]={bm_by_u[5.0]:.3f}  slope={slope:+.4f}  ({dt:.1f}s)",
              flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "probe3_xeno_burden.csv", index=False)
print(f"\n[write] {OUT/'probe3_xeno_burden.csv'}  ({len(df)} rows)")

# verdict: slope (growth cost per unit X) per condition x hpf x variant
slope_df = df[df.u.isna()][["cell", "condition", "hpf", "variant", "slope"]]
piv_a = slope_df[slope_df.variant == "Phase2_glucuronide"].pivot_table(
    index="hpf", columns="condition", values="slope", aggfunc="mean")
piv_b = slope_df[slope_df.variant == "Phase2_GSH"].pivot_table(
    index="hpf", columns="condition", values="slope", aggfunc="mean")
print("\n--- slope (biomass cost per unit X, Phase2a glucuronide) ---")
print(piv_a.round(4).to_string())
print("\n--- slope (biomass cost per unit X, Phase2b GSH) ---")
print(piv_b.round(4).to_string())
