"""
s185 -- lipid-relay baseline from existing pooled data (AA, BB corners).

Before dispatching the 14 swap-corner sampling runs, compute everything the
preregistration (s178_lipid_relay_prereg.md) needs that can be derived from
already-sampled data:

  - s_AA and s_BB per pair (from existing grace_v4/results/pooled_v7/*)
  - total = s_AA − s_BB  (comparable to the Tier A direction × W1 sign)
  - biomax for all four corners (AA known; AB, BA, BB computed fresh)
  - matched floor = 0.9 × min(biomax_AA, biomax_AB, biomax_BA, biomax_BB)

Writes grace_v4/results/v7_final_addenda/s178_lipid_relay_baseline.csv

The swap-corner samples (AB, BA, s_AB, s_BA, decomposition, verdict) are
left for s184 → s186. This script lets us set expectations and budgets
before any sampling runs.
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import h5py
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
from apply_context_bounds import apply_context_bounds, compute_biomax

R = REPO / "grace_v4/results"
OUT = R / "v7_final_addenda"
V7 = REPO / "reviewer_packet_v4_slim/v7_models"

PAIRS = [
    # (A, B, subsystem)
    ("D_24",  "D_48",  "Fatty acid activation (cytosolic)"),
    ("LD_24", "LD_48", "Fatty acid activation (cytosolic)"),
    ("D_48",  "D_72",  "Fatty acid activation (cytosolic)"),
    ("LD_48", "LD_72", "Fatty acid activation (cytosolic)"),
    ("BL_48", "BL_72", "Fatty acid oxidation"),
    ("LD_48", "LD_72", "Fatty acid oxidation"),
    ("BL_48", "BL_72", "Arachidonic acid metabolism"),
    ("LD_48", "LD_72", "Arachidonic acid metabolism"),
    ("D_72",  "D_96",  "Arachidonic acid metabolism"),
    ("LD_72", "LD_96", "Arachidonic acid metabolism"),
    ("D_48",  "D_72",  "Leukotriene metabolism"),
    ("LD_48", "LD_72", "Leukotriene metabolism"),
    ("D_72",  "D_96",  "Leukotriene metabolism"),
    ("LD_72", "LD_96", "Leukotriene metabolism"),
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


# ---------- Pre-load all needed models once ----------
cells = sorted(set(a for a, b, s in PAIRS) | set(b for a, b, s in PAIRS))
models = {c: _load_model(V7 / f"trans_rfastcormics_{c}.mat") for c in cells}
subsets = {c: pd.read_csv(R / f"subsets_v7/{c}_subset_membership.csv") for c in cells}
pooled  = {c: pd.read_csv(R / f"pooled_v7/{c}_pooled.csv").set_index("rxn") for c in cells}
print(f"[prep] loaded {len(models)} models, subsets, and pooled samples", flush=True)


def _biomax_corner(model_cell, env_cell):
    """Biomax for (model_cell network, env_cell bounds). Wang internal
    ±1000 defaults INTACT to match frozen pooled_v7 (user flag 2026-10-14)."""
    mdl = models[model_cell]
    env_cond, env_hpf = env_cell.split("_"); env_hpf = int(env_hpf)
    is_ex = (np.diff(mdl["S"].tocsc().indptr) == 1)
    lb, ub, idx, _ = apply_context_bounds(
        mdl["lb"].copy(), mdl["ub"].copy(), mdl["rxns"], is_ex,
        env_cond, env_hpf, S=mdl["S"])
    bm, _ = compute_biomax(mdl, lb, ub, idx)
    return bm


def _s_corner(cell, subsystem, mu_pair):
    """signed median log2(|v|/mu_pair) for pooled samples on representatives
    of `subsystem` in `cell`."""
    sub = subsets[cell]
    reps = sub.loc[sub.subsystem == subsystem, "representative"].unique()
    pool = pooled[cell]
    rep_in = [r for r in reps if r in pool.index]
    if not rep_in:
        return np.nan, 0.0, 0
    flat = pool.loc[rep_in].values.flatten()
    nonzero = flat[flat != 0]
    zero_frac = 1 - nonzero.size / flat.size if flat.size else 0.0
    if nonzero.size == 0:
        return np.nan, zero_frac, 0
    s = np.median(np.log2(np.abs(nonzero) / mu_pair))
    return s, zero_frac, nonzero.size


# ---------- Walk pairs ----------
rows = []
for a, b, subsystem in PAIRS:
    print(f"\n[pair] {a} → {b}  ({subsystem})", flush=True)

    # Biomax for all four corners
    bm_AA = _biomax_corner(a, a)
    bm_BB = _biomax_corner(b, b)
    bm_AB = _biomax_corner(a, b)
    bm_BA = _biomax_corner(b, a)
    matched_floor = 0.9 * min(bm_AA, bm_AB, bm_BA, bm_BB)

    print(f"  biomax  AA={bm_AA:.3f}  AB={bm_AB:.3f}  BA={bm_BA:.3f}  BB={bm_BB:.3f}", flush=True)
    print(f"  matched floor (0.9 × min)  = {matched_floor:.3f}", flush=True)

    # mu_pair reference: use pool-mean |v| across AA + BB subsystem flux
    # (will extend to all 4 corners once AB, BA sampled).
    sub_a = subsets[a]; sub_b = subsets[b]
    reps_a = sub_a.loc[sub_a.subsystem == subsystem, "representative"].unique()
    reps_b = sub_b.loc[sub_b.subsystem == subsystem, "representative"].unique()
    rep_a_in = [r for r in reps_a if r in pooled[a].index]
    rep_b_in = [r for r in reps_b if r in pooled[b].index]
    all_vals = np.concatenate([
        pooled[a].loc[rep_a_in].values.flatten() if rep_a_in else np.array([]),
        pooled[b].loc[rep_b_in].values.flatten() if rep_b_in else np.array([]),
    ])
    nz_all = all_vals[all_vals != 0]
    mu_pair_AA_BB = np.abs(nz_all).mean() if nz_all.size else 1.0

    s_AA, zf_AA, n_AA = _s_corner(a, subsystem, mu_pair_AA_BB)
    s_BB, zf_BB, n_BB = _s_corner(b, subsystem, mu_pair_AA_BB)
    total = s_AA - s_BB

    print(f"  s_AA={s_AA:+.3f}  s_BB={s_BB:+.3f}  total={total:+.3f}  "
          f"zf_A={zf_AA:.2f}  zf_B={zf_BB:.2f}", flush=True)

    rows.append(dict(
        A=a, B=b, subsystem=subsystem,
        biomax_AA=bm_AA, biomax_AB=bm_AB, biomax_BA=bm_BA, biomax_BB=bm_BB,
        matched_floor=matched_floor,
        mu_pair_AA_BB=mu_pair_AA_BB,
        s_AA=s_AA, s_BB=s_BB,
        total_AA_minus_BB=total,
        zero_frac_AA=zf_AA, zero_frac_BB=zf_BB,
        n_nonzero_AA=n_AA, n_nonzero_BB=n_BB,
        subsystem_rep_count_A=len(rep_a_in),
        subsystem_rep_count_B=len(rep_b_in),
    ))

df = pd.DataFrame(rows)
df.to_csv(OUT / "s178_lipid_relay_baseline.csv", index=False)
print(f"\n[write] {OUT/'s178_lipid_relay_baseline.csv'}  ({len(df)} rows)")

# Quick summary
print("\n--- compute budget for s184 swap sampling ---")
unique_corners = set()
for a, b, _ in PAIRS:
    unique_corners.add((a, b)); unique_corners.add((b, a))
print(f"  unique swap corners (model, env):  {len(unique_corners)}")
print(f"  estimated per-corner cost:          8–15 min (polytope + 4-chain hopsy + pool)")
print(f"  estimated total compute (local):    {len(unique_corners)*10}–{len(unique_corners)*15} min")
print()
print("--- totals (|s_AA − s_BB|; sign = AA > BB) ---")
for _, r in df.iterrows():
    arrow = "↑" if r.total_AA_minus_BB > 0 else "↓"
    print(f"  {r.A:7s}→{r.B:7s}  {r.subsystem[:40]:40s}  total {arrow} {abs(r.total_AA_minus_BB):.3f}  "
          f"floor {r.matched_floor:.2f}")
