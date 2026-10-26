"""
s189 -- Table 5b: 2×2 topology × envelope decomposition for all 12 F2
adjacent-stage pairs (within-regime), same structure as Table 5.

Columns: pair, gap_raw, topology, envelope, interaction, verdict.
Gates (locked): |interaction|/|gap_raw| > 0.30 → not separable;
                |env|/(|topo|+|env|) ≥ 0.70 → envelope-dominated;
                |topo|/(|topo|+|env|) ≥ 0.70 → topology-dominated;
                else → mixed.
Frozen ±1000 bounds. Free biomass (no matched floor; Table 5 convention).

Seven of the 12 pairs are already captured in s178_lipid_relay_baseline.csv
(same (A, env_A), (A, env_B), (B, env_A), (B, env_B) biomax quadruple).
This script computes the five missing pairs and merges all 12 into
grace_v4/results/v7_final_addenda/table5b_developmental_2x2.csv.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import h5py

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds, compute_biomax

OUT = REPO / "grace_v4/results/v7_final_addenda"
V7 = REPO / "reviewer_packet_v4_slim/v7_models"

# All 12 F2 pairs (A = earlier, B = later)
PAIRS = [
    ("BL_24", "BL_48"),  ("BL_48", "BL_72"), ("BL_72", "BL_96"), ("BL_96", "BL_120"),
    ("D_24",  "D_48"),   ("D_48",  "D_72"),  ("D_72",  "D_96"),  ("D_96",  "D_120"),
    ("LD_24", "LD_48"),  ("LD_48", "LD_72"), ("LD_72", "LD_96"), ("LD_96", "LD_120"),
]


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
        mets = [_resolve(f, r) for r in np.asarray(m["mets"]).ravel()]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir = np.asarray(S_grp["ir"]); jc = np.asarray(S_grp["jc"])
        n_met = max(int(ir.max()) + 1 if len(ir) else 0, len(mets))
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, len(jc) - 1))
        lb = np.asarray(m["lb"]).ravel().astype(float)
        ub = np.asarray(m["ub"]).ravel().astype(float)
    return dict(rxns=rxns, mets=mets, S=S, lb=lb, ub=ub)


def _biomax_corner(model_cell, env_cell):
    """Biomax with model_cell's network under env_cell's envelope.
    Frozen ±1000 defaults (no widening); free biomass."""
    mdl = _load_model(V7 / f"trans_rfastcormics_{model_cell}.mat")
    env_cond, env_hpf = env_cell.split("_"); env_hpf = int(env_hpf)
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex,
        env_cond, env_hpf, S=mdl["S"])
    bm, _ = compute_biomax(mdl, lb, ub, idx)
    return float(bm)


def _classify(topo, env, inter, total):
    if abs(total) < 1e-6: return "no material gap"
    if abs(inter) / abs(total) > 0.30: return "not separable"
    denom = abs(topo) + abs(env)
    if denom < 1e-9: return "mixed"
    env_sh = abs(env) / denom
    topo_sh = abs(topo) / denom
    if env_sh >= 0.70:  return "envelope"
    if topo_sh >= 0.70: return "topology"
    return "mixed"


# ---- Load the 7 pairs we already have (from s178 baseline) ----
baseline = pd.read_csv(OUT / "s178_lipid_relay_baseline.csv")
# collapse to unique (A, B) pairs since baseline has subsystem-level dup rows
baseline = baseline.drop_duplicates(subset=["A", "B"])[
    ["A", "B", "biomax_AA", "biomax_AB", "biomax_BA", "biomax_BB"]]
have = {(r.A, r.B): r for _, r in baseline.iterrows()}

rows = []
for A, B in PAIRS:
    if (A, B) in have:
        r = have[(A, B)]
        AA, AB, BA, BB = r.biomax_AA, r.biomax_AB, r.biomax_BA, r.biomax_BB
        source = "s178 baseline (lipid relay)"
    else:
        print(f"  computing new corners for {A} → {B}", flush=True)
        AA = _biomax_corner(A, A)
        AB = _biomax_corner(A, B)
        BA = _biomax_corner(B, A)
        BB = _biomax_corner(B, B)
        source = "s189 fresh"
    topo = AA - BA
    env  = AA - AB
    total = AA - BB
    inter = total - topo - env
    verdict = _classify(topo, env, inter, total)
    rows.append(dict(
        pair=f"{A}→{B}", A=A, B=B,
        biomax_AA=round(AA, 3), biomax_AB=round(AB, 3),
        biomax_BA=round(BA, 3), biomax_BB=round(BB, 3),
        gap_raw=round(total, 3),
        topology=round(topo, 3),
        envelope=round(env, 3),
        interaction=round(inter, 3),
        verdict=verdict,
        source=source,
    ))

df = pd.DataFrame(rows)
df.to_csv(OUT / "table5b_developmental_2x2.csv", index=False)
print(f"\n[write] {OUT/'table5b_developmental_2x2.csv'}  ({len(df)} rows)")
print()
print(df[["pair","gap_raw","topology","envelope","interaction","verdict"]].to_string(index=False))
print()
print("--- verdict counts ---")
print(df.verdict.value_counts().to_string())
