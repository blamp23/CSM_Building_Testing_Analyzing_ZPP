"""
s181 -- F2 late relay decompose (D 96→120, LD 96→120).

Same scheme as s186_lipid_relay_decompose.py; different PAIRS table.
F2 convention: A = earlier stage. Covered subsystems come from Tier A at
these transitions:
    D   96→120 : Nucleotide, Pyrimidine, Phenylalanine, tyrosine and tryptophan biosynthesis
    LD  96→120 : Arginine and proline metabolism, Glycolysis / Gluconeogenesis,
                 Phenylalanine, tyrosine and tryptophan biosynthesis,
                 Valine, leucine, and isoleucine metabolism

Writes grace_v4/results/v7_final_addenda/s181_f2_late_relay.{csv,md}.
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
CHAIN_DIR = POOL_SWAP / "chains"

PAIRS = [
    # (A = earlier stage, B = later stage, subsystem)
    ("D_96",  "D_120",  "Nucleotide metabolism"),
    ("D_96",  "D_120",  "Pyrimidine metabolism"),
    ("D_96",  "D_120",  "Phenylalanine, tyrosine and tryptophan biosynthesis"),
    ("LD_96", "LD_120", "Arginine and proline metabolism"),
    ("LD_96", "LD_120", "Glycolysis / Gluconeogenesis"),
    ("LD_96", "LD_120", "Phenylalanine, tyrosine and tryptophan biosynthesis"),
    ("LD_96", "LD_120", "Valine, leucine, and isoleucine metabolism"),
]


def _reps(cell, subsystem):
    sub = pd.read_csv(SUBSETS / f"{cell}_subset_membership.csv")
    return sub.loc[sub.subsystem == subsystem, "representative"].unique().tolist()


def _flux_values(cell, subsystem, pool_fp_or_chain_dir, tag=None):
    reps = _reps(cell, subsystem)
    if pool_fp_or_chain_dir.is_dir():
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
    if abs(total) < 1e-9: return "no material gap"
    if abs(inter) / abs(total) > 0.30: return "not separable"
    denom = abs(topo) + abs(env)
    if denom < 1e-9: return "mixed"
    if abs(env) / denom >= 0.70: return "relayed"
    if abs(topo) / denom >= 0.70: return "generated (needs GPR audit)"
    return "mixed"


rows = []
for A, B, subsystem in PAIRS:
    AA_fp = POOL_STD / f"{A}_pooled.csv"
    BB_fp = POOL_STD / f"{B}_pooled.csv"
    AB_tag = f"{A}_mdl__{B}_env"
    BA_tag = f"{B}_mdl__{A}_env"
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
df.to_csv(OUT / "s181_f2_late_relay.csv", index=False)
print(f"[write] {OUT/'s181_f2_late_relay.csv'}  ({len(df)} rows)")

md = ["# s181 — F2 late relay test  (D / LD at 96→120 hpf)",
      "",
      "Signed median of log₂(|v|/μ) across subsystem representatives,",
      "pooled across all four corners of each pair. F2 convention:",
      "A = earlier stage. Thresholds locked (|interaction|/|total| > 0.30 →",
      "not separable; else 70 % envelope → relayed; 70 % topology → generated).",
      "",
      "| A | B | subsystem | s_AA | s_AB | s_BA | s_BB | topology | envelope | interaction | total | verdict |",
      "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for _, r in df.iterrows():
    if r.verdict == "pending":
        md.append(f"| {r.A} | {r.B} | {r.subsystem} | — | — | — | — | — | — | — | — | pending |")
        continue
    md.append(f"| {r.A} | {r.B} | {r.subsystem} | "
              f"{r.s_AA:+.3f} | {r.s_AB:+.3f} | {r.s_BA:+.3f} | {r.s_BB:+.3f} | "
              f"{r.topology:+.3f} | {r.envelope:+.3f} | {r.interaction:+.3f} | {r.total:+.3f} | "
              f"**{r.verdict}** |")
md += ["", "## Verdict counts", ""]
for v, n in sorted(df.verdict.value_counts(dropna=False).to_dict().items(), key=lambda x: -x[1]):
    md.append(f"- {v}: **{n}**")
md += ["",
       "## Convention notes",
       "- F2 cross-stage: envelope swap moves BOTH 105-metabolite fc AND OCR cap",
       "  (OCR[24]=97 → OCR[120]=313). Report envelope column as “measured environment",
       "  (metabolome fc + OCR cap)” in the chapter.",
       "- A \\ B topology audit (§8.6 GPR class) is MANDATORY on any 'generated' verdict.",
       ""]
(OUT / "s181_f2_late_relay.md").write_text("\n".join(md))
print(f"[write] {OUT/'s181_f2_late_relay.md'}")
print()
if len(df):
    print(df[["A","B","subsystem","topology","envelope","interaction","total","verdict"]]
          .to_string(index=False))
