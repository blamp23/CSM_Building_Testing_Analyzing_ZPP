"""
s275 -- 2x2 "metabolic capacity across development" figure.

Panels A, B: ATP ceilings from existing s268 / s269 outputs.
Panels C, D: Max xenobiotic clearance rate for glucuronide and GSH routes,
             subject to biomass >= 0.9 x biomax (free biomass ceiling).

For C and D we build a scratch Phase I + Phase II + export pathway on each
cell's extracted model and maximize xenobiotic uptake.
"""
import sys, h5py
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FUNCTIONS, V7_MODELS as V7, FIGURES as OUT, TABLES, SUSCEPTIBILITY
sys.path.insert(0, str(FUNCTIONS))
from apply_context_bounds import apply_context_bounds, compute_biomax
OUT.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

ATP_RXN = "MAR03964"
BIOMASS_RXN = "MAR00021"
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)

# Required metabolites for Phase I and Phase II
# NADPH cytosolic, O2 cytosolic, GSH cytosolic, UDP-glucuronate cytosolic
NADPH = "MAM02555c"
O2    = "MAM02630c"
H     = "MAM02039c"
H2O   = "MAM02040c"
NADP  = "MAM02554c"
GSH   = "MAM02026c"
UDP_GLCA = "MAM03109c"
UDP   = "MAM03106c"


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
        mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = int(ir.max()) + 1 if len(ir) else 0
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, S=S, lb=lb, ub=ub)


def _max_x_uptake(mdl, lb, ub, idx, biomass_floor, route):
    """Build scratch pathway and maximize X uptake while biomass >= floor.

    Scratch columns appended to S:
        X_ex:  ∅ -> X_c                     (lb=-1000, ub=0 by convention: v<=0 means uptake)
        PhaseI: X + NADPH + O2 + H -> X-OH + NADP + H2O
        PhaseIIa (glucuronide): X-OH + UDP-GlcA -> X-gluc + UDP
        PhaseIIb (GSH):         X-OH + GSH -> X-SG
        Export: X-gluc or X-SG -> ∅

    Scratch mets: X_c, X_OH_c, X_gluc_c or X_SG_c
    """
    S = mdl["S"].tolil()
    n_met_orig, n_rxn_orig = S.shape
    mets_new = mdl["mets"] + ["X_c", "X_OH_c",
                              ("X_gluc_c" if route == "glucuronide" else "X_SG_c")]
    n_met_new = len(mets_new)
    # Pad S with 3 new metabolite rows
    for _ in range(3):
        S.resize((S.shape[0] + 1, S.shape[1]))
    X_i = n_met_orig
    X_OH_i = n_met_orig + 1
    X_conj_i = n_met_orig + 2
    # Required cofactor met indices
    def _mi(m):
        return mdl["mets"].index(m) if m in mdl["mets"] else None
    nadph_i = _mi(NADPH); o2_i = _mi(O2); h_i = _mi(H); h2o_i = _mi(H2O); nadp_i = _mi(NADP)
    gsh_i = _mi(GSH); udp_glca_i = _mi(UDP_GLCA); udp_i = _mi(UDP)

    # Add 4 scratch reactions
    # Each reaction = one new column in S
    def _add_rxn(stoich_dict):
        """stoich_dict maps met_idx -> coeff."""
        nonlocal S
        col = sp.lil_matrix((n_met_new, 1))
        for mi, coeff in stoich_dict.items():
            if mi is not None:
                col[mi, 0] = coeff
        S = sp.hstack([S.tocsr(), col.tocsr()]).tolil()

    # X_ex: ∅ -> X_c (produces X_c when v > 0; we'll constrain v in [0, 1000])
    _add_rxn({X_i: +1.0})
    j_X_ex = n_rxn_orig

    # Phase I: X + NADPH + O2 + H -> X-OH + NADP + H2O
    _add_rxn({X_i: -1.0, nadph_i: -1.0, o2_i: -1.0, h_i: -1.0,
              X_OH_i: +1.0, nadp_i: +1.0, h2o_i: +1.0})
    j_phI = n_rxn_orig + 1

    # Phase II
    if route == "glucuronide":
        _add_rxn({X_OH_i: -1.0, udp_glca_i: -1.0, X_conj_i: +1.0, udp_i: +1.0})
    else:
        _add_rxn({X_OH_i: -1.0, gsh_i: -1.0, X_conj_i: +1.0})
    j_phII = n_rxn_orig + 2

    # Export: X_conj -> ∅
    _add_rxn({X_conj_i: -1.0})
    j_export = n_rxn_orig + 3

    # Build Gurobi model
    S_csr = S.tocsr()
    n_rxn_new = S_csr.shape[1]
    lb_new = np.concatenate([lb, [0.0, 0.0, 0.0, 0.0]])
    ub_new = np.concatenate([ub, [1000.0, 1000.0, 1000.0, 1000.0]])
    # Apply biomass floor
    bio_j = idx[BIOMASS_RXN]
    lb_new[bio_j] = biomass_floor

    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model(env=env)
    v = gm.addMVar(n_rxn_new, lb=lb_new, ub=ub_new)
    gm.addMConstr(S_csr, v, "=", np.zeros(S_csr.shape[0]))
    obj = np.zeros(n_rxn_new); obj[j_X_ex] = 1.0
    gm.setMObjective(None, obj, 0.0, sense=GRB.MAXIMIZE)
    gm.optimize()
    if gm.Status != GRB.OPTIMAL:
        gm.dispose(); env.dispose()
        return np.nan
    val = gm.ObjVal
    gm.dispose(); env.dispose()
    return float(val)


# ---------- Compute Panels C and D ----------
print("[prep] loading cells + computing biomax + max X clearance")
results = []
for cond in CONDS:
    for h in HPFS:
        cell = f"{cond}_{h}"
        mdl = _load(V7 / f"trans_rfastcormics_{cell}.mat")
        is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
        lb, ub, idx, _ = apply_context_bounds(
            mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex, cond, h,
            S=mdl["S"])
        non_ex = ~is_ex
        lb[(lb <= -999.9) & non_ex] = -1e6
        ub[(ub >= 999.9) & non_ex] = +1e6
        bm, _ = compute_biomax(mdl, lb, ub, idx)
        floor = 0.9 * bm
        v_gluc = _max_x_uptake(mdl, lb, ub, idx, floor, "glucuronide")
        v_gsh  = _max_x_uptake(mdl, lb, ub, idx, floor, "GSH")
        results.append(dict(cell=cell, condition=cond, hpf=h,
                            biomax=bm, biomass_floor=floor,
                            max_X_glucuronide=v_gluc,
                            max_X_GSH=v_gsh))
        print(f"  {cell:7s}  biomax={bm:6.2f}  floor={floor:6.2f}  "
              f"max_X_gluc={v_gluc:7.2f}  max_X_GSH={v_gsh:7.2f}", flush=True)

df_xeno = pd.DataFrame(results)
df_xeno.to_csv(TABLES / "xenobiotic_capacity.csv", index=False)
print(f"\n[write] {TABLES/'xenobiotic_capacity.csv'}")

# ---------- Pull ATP data for Panels A and B ----------
atp_full = pd.read_csv(SUSCEPTIBILITY / "atp_ceiling_physio_ocr.csv")
atp_glc  = pd.read_csv(SUSCEPTIBILITY / "atp_ceiling_glucose_only.csv")
atp_full_free = atp_full[atp_full.regime == "free"].copy()
atp_glc_free = atp_glc[atp_glc.regime == "free"].copy()
atp_full_free["hpf"] = atp_full_free.cell.str.split("_").str[1].astype(int)
atp_full_free["cond"] = atp_full_free.cell.str.split("_").str[0]
atp_glc_free["hpf"]  = atp_glc_free.cell.str.split("_").str[1].astype(int)
atp_glc_free["cond"] = atp_glc_free.cell.str.split("_").str[0]

# ---------- Build the 2x2 figure ----------
cond_color  = {"BL": "#2E7D32", "D": "#333333", "LD": "#E65100"}
cond_marker = {"BL": "o", "D": "s", "LD": "^"}

fig, axes = plt.subplots(2, 2, figsize=(11, 9))

# Panel A: ATP full envelope
ax = axes[0, 0]
piv = atp_full_free.pivot_table(index="hpf", columns="cond", values="atp_max")
for c in CONDS:
    ax.plot(piv.index, piv[c].values, marker=cond_marker[c], color=cond_color[c],
            linewidth=2, markersize=10, label=c, markeredgecolor="black", markeredgewidth=0.7)
ax.set_xlabel("hpf", fontsize=12)
ax.set_ylabel("ATP max flux (mmol / gDW / hr)", fontsize=12)
ax.set_xticks(HPFS)
ax.grid(alpha=0.3)
ax.legend(fontsize=10, loc="upper left")

# Panel B: ATP glucose only
ax = axes[0, 1]
piv = atp_glc_free.pivot_table(index="hpf", columns="cond", values="atp_max")
for c in CONDS:
    ax.plot(piv.index, piv[c].values, marker=cond_marker[c], color=cond_color[c],
            linewidth=2, markersize=10, label=c, markeredgecolor="black", markeredgewidth=0.7)
ax.set_xlabel("hpf", fontsize=12)
ax.set_ylabel("ATP max flux (mmol / gDW / hr)", fontsize=12)
ax.set_xticks(HPFS)
ax.grid(alpha=0.3)

# Panel C: xenobiotic clearance glucuronide
ax = axes[1, 0]
piv = df_xeno.pivot_table(index="hpf", columns="condition", values="max_X_glucuronide")
for c in CONDS:
    ax.plot(piv.index, piv[c].values, marker=cond_marker[c], color=cond_color[c],
            linewidth=2, markersize=10, label=c, markeredgecolor="black", markeredgewidth=0.7)
ax.set_xlabel("hpf", fontsize=12)
ax.set_ylabel("max X clearance (mmol / gDW / hr)", fontsize=12)
ax.set_xticks(HPFS)
ax.grid(alpha=0.3)

# Panel D: xenobiotic clearance GSH
ax = axes[1, 1]
piv = df_xeno.pivot_table(index="hpf", columns="condition", values="max_X_GSH")
for c in CONDS:
    ax.plot(piv.index, piv[c].values, marker=cond_marker[c], color=cond_color[c],
            linewidth=2, markersize=10, label=c, markeredgecolor="black", markeredgewidth=0.7)
ax.set_xlabel("hpf", fontsize=12)
ax.set_ylabel("max X clearance (mmol / gDW / hr)", fontsize=12)
ax.set_xticks(HPFS)
ax.grid(alpha=0.3)

# Panel labels (A, B, C, D only — no in-figure descriptions)
for ax, label in zip(axes.flat, ("A", "B", "C", "D")):
    ax.text(0.02, 0.98, label, transform=ax.transAxes,
            fontsize=16, fontweight="bold", va="top")

fig.tight_layout()
out_pdf = OUT / "metabolic_capacity_developmental.pdf"
out_png = OUT / "metabolic_capacity_developmental.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\n[write] {out_pdf}")
print(f"[write] {out_png}")

# ---------- Print table for chapter ----------
print("\n=== xenobiotic clearance capacity ===")
print(df_xeno[["cell","condition","hpf","biomax","max_X_glucuronide","max_X_GSH"]].to_string(index=False))
