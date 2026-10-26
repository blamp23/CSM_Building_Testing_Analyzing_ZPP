"""
s186 -- aggregate the 4 swap corners per pair, decompose, apply the locked
decision rule from s178_lipid_relay_prereg.md, write s178_lipid_relay.{csv,md}.

Reads:
  - grace_v4/results/pooled_v7/{cell}_pooled.csv            (AA, BB)
  - grace_v4/results/v7_final_addenda/lipid_relay_samples/
      {A}_mdl__{B}_env_pooled.csv                           (AB)
      {B}_mdl__{A}_env_pooled.csv                           (BA)
  - grace_v4/results/subsets_v7/{cell}_subset_membership.csv for subsystem mapping

Statistic (per preregistration):
  s_corner = median_{rxn, sample} log2(|v| / mu_pair)
  where mu_pair = pool-mean |v| across the 4 corners' subsystem flux.

Decomposition:
  topology    = s_AA − s_BA
  envelope    = s_AA − s_AB
  total       = s_AA − s_BB
  interaction = total − topology − envelope

Decision gates (locked):
  1. |interaction| / |total| > 0.30            → not separable
  2. |envelope| / (|topo|+|env|) ≥ 0.70         → relayed
  3. |topology| / (|topo|+|env|) ≥ 0.70         → generated (needs GPR audit)
  4. else                                       → mixed

Writes s178_lipid_relay.csv + s178_lipid_relay.md.
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
OUT = R / "v7_final_addenda"
POOL_STD = R / "pooled_v7"
POOL_SWAP = OUT / "lipid_relay_samples"
SUBSETS = R / "subsets_v7"

PAIRS = [
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


def _reps(cell, subsystem):
    sub = pd.read_csv(SUBSETS / f"{cell}_subset_membership.csv")
    return sub.loc[sub.subsystem == subsystem, "representative"].unique().tolist()


def _flux_values(cell, subsystem, pool_fp_or_chain_dir, tag=None):
    """Flattened flux values across (reps × samples) for the subsystem.
    Accepts either a single pooled CSV path (pool_v7 standard runs) or a
    directory of chain CSVs (swap-corner runs); when a dir, `tag` is
    required and chains are concatenated column-wise."""
    reps = _reps(cell, subsystem)
    if pool_fp_or_chain_dir.is_dir():
        import glob
        chain_fps = sorted(pool_fp_or_chain_dir.glob(f"{tag}_ch*.csv"))
        if not chain_fps: return np.array([])
        frames = [pd.read_csv(fp).set_index("rxn") for fp in chain_fps]
        pool = pd.concat(frames, axis=1)
    else:
        pool = pd.read_csv(pool_fp_or_chain_dir).set_index("rxn")
    rep_in = [r for r in reps if r in pool.index]
    if not rep_in: return np.array([])
    return pool.loc[rep_in].values.flatten()


def _s(flux_vals, mu_pair):
    nz = flux_vals[flux_vals != 0]
    if nz.size == 0: return (np.nan, 1.0, 0)
    s = np.median(np.log2(np.abs(nz) / mu_pair))
    return (float(s), 1 - nz.size / flux_vals.size, int(nz.size))


def _classify(topo, env, inter, total):
    """Apply locked decision gates in order."""
    if abs(total) < 1e-9:
        return "no material gap"
    if abs(inter) / abs(total) > 0.30:
        return "not separable"
    denom = abs(topo) + abs(env)
    if denom < 1e-9:
        return "mixed"
    env_sh = abs(env) / denom
    topo_sh = abs(topo) / denom
    if env_sh >= 0.70:  return "relayed"
    if topo_sh >= 0.70: return "generated (needs GPR audit)"
    return "mixed"


CHAIN_DIR = POOL_SWAP / "chains"
rows = []
for A, B, subsystem in PAIRS:
    AA_fp = POOL_STD / f"{A}_pooled.csv"
    BB_fp = POOL_STD / f"{B}_pooled.csv"
    AB_tag = f"{A}_mdl__{B}_env"
    BA_tag = f"{B}_mdl__{A}_env"
    # Count completed chains for AB and BA
    import glob
    AB_chains = sorted(CHAIN_DIR.glob(f"{AB_tag}_ch*.csv")) if CHAIN_DIR.exists() else []
    BA_chains = sorted(CHAIN_DIR.glob(f"{BA_tag}_ch*.csv")) if CHAIN_DIR.exists() else []
    missing = []
    if not AA_fp.exists(): missing.append(str(AA_fp))
    if not BB_fp.exists(): missing.append(str(BB_fp))
    if len(AB_chains) < 4: missing.append(f"{AB_tag}_ch*.csv ({len(AB_chains)}/4)")
    if len(BA_chains) < 4: missing.append(f"{BA_tag}_ch*.csv ({len(BA_chains)}/4)")
    if missing:
        rows.append(dict(A=A, B=B, subsystem=subsystem, verdict="pending",
                         missing=";".join(missing)))
        continue

    flux_AA = _flux_values(A, subsystem, AA_fp)
    flux_AB = _flux_values(A, subsystem, CHAIN_DIR, AB_tag)
    flux_BA = _flux_values(B, subsystem, CHAIN_DIR, BA_tag)
    flux_BB = _flux_values(B, subsystem, BB_fp)
    nz_all = np.concatenate([flux_AA, flux_AB, flux_BA, flux_BB])
    nz_all = nz_all[nz_all != 0]
    mu_pair = float(np.abs(nz_all).mean()) if nz_all.size else 1.0
    s_AA, zf_AA, n_AA = _s(flux_AA, mu_pair)
    s_AB, zf_AB, n_AB = _s(flux_AB, mu_pair)
    s_BA, zf_BA, n_BA = _s(flux_BA, mu_pair)
    s_BB, zf_BB, n_BB = _s(flux_BB, mu_pair)
    topo  = s_AA - s_BA
    env   = s_AA - s_AB
    total = s_AA - s_BB
    inter = total - topo - env
    verdict = _classify(topo, env, inter, total)
    rows.append(dict(
        A=A, B=B, subsystem=subsystem,
        mu_pair=mu_pair,
        s_AA=s_AA, s_AB=s_AB, s_BA=s_BA, s_BB=s_BB,
        topology=topo, envelope=env, interaction=inter, total=total,
        verdict=verdict,
        zero_frac_AA=zf_AA, zero_frac_AB=zf_AB,
        zero_frac_BA=zf_BA, zero_frac_BB=zf_BB,
        n_nz_AA=n_AA, n_nz_AB=n_AB, n_nz_BA=n_BA, n_nz_BB=n_BB,
    ))

df = pd.DataFrame(rows)
df.to_csv(OUT / "s178_lipid_relay.csv", index=False)
print(f"[write] {OUT/'s178_lipid_relay.csv'}  ({len(df)} rows)")

# ---------- Markdown summary ----------
md = ["# s178 — Tier A lipid relay test (per s178_lipid_relay_prereg.md)",
      "",
      "Decomposition of the 14 Tier A F2 lipid transitions into topology, envelope,",
      "and interaction components via the §3.5 four-corner swap on **sampled flux**.",
      "The statistic is a signed median of log₂(|v|/μ) with μ pooled across the 4",
      "corners (common reference). Thresholds are locked in the preregistration.",
      "",
      "| A | B | subsystem | s_AA | s_AB | s_BA | s_BB | topology | envelope | interaction | total | verdict |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|",]
for _, r in df.iterrows():
    if r.verdict == "pending":
        md.append(f"| {r.A} | {r.B} | {r.subsystem} | — | — | — | — | — | — | — | — | pending |")
        continue
    md.append(f"| {r.A} | {r.B} | {r.subsystem} | "
              f"{r.s_AA:+.3f} | {r.s_AB:+.3f} | {r.s_BA:+.3f} | {r.s_BB:+.3f} | "
              f"{r.topology:+.3f} | {r.envelope:+.3f} | {r.interaction:+.3f} | {r.total:+.3f} | "
              f"**{r.verdict}** |")
md += ["",
       "## Verdict counts",
       ""]
counts = df.verdict.value_counts(dropna=False).to_dict()
for v, n in sorted(counts.items(), key=lambda x: -x[1]):
    md.append(f"- {v}: **{n}**")
md += ["",
       "## Decision rule (locked, per preregistration)",
       "",
       "Applied in order to the signed decomposition on sampled flux:",
       "",
       "1. `|interaction| / |total| > 0.30` → **not separable** (BL_96 precedent)",
       "2. `|envelope| / (|topo| + |env|) ≥ 0.70` → **relayed**",
       "3. `|topology| / (|topo| + |env|) ≥ 0.70` → **generated** (requires GPR audit)",
       "4. else → **mixed**",
       "",
       "## Mandatory follow-up on any 'generated' verdict",
       "",
       "Run the GPR + zFPKM audit from s170 / s177 on reactions in the "
       "A \\ B subsystem set. A 'generated' verdict is biological only if the "
       "carrying reactions are TRANSCRIPTIONAL_PARTIAL. If UNKNOWN_ZONE or "
       "GENELESS, restate as inseparable from extraction choice.",
       ""]
(OUT / "s178_lipid_relay.md").write_text("\n".join(md))
print(f"[write] {OUT/'s178_lipid_relay.md'}")
print()
print(df[["A","B","subsystem","topology","envelope","interaction","total","verdict"]]
      .to_string(index=False))
