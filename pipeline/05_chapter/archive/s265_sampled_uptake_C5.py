"""
C5: Sampled uptake of inositol, β-carotene, every row in the B3 offending
set, and glucose for scale. From existing v7 pooled samples.

Negative-signed flux = uptake (model convention). Median and 95th
percentile of signed flux per cell. Also summary stats per row across
cells.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from s250_atp_sweep import CELLS

OUT = REPO / "grace_v4/results/atp_sweep_v7"
POOL_DIR = REPO / "grace_v4/results/pooled_v7"
SUBS_DIR = REPO / "grace_v4/results/subsets_v7"

TARGETS = {
    "MAR09361": "inositol",
    "MAR09276": "beta-carotene",
    "MAR09152": "alpha-tocotrienol",
    "MAR09154": "gamma-tocotrienol",
    "MAR09034": "glucose (measured scale)",
}

rows = []
for cell in CELLS:
    sub = pd.read_csv(SUBS_DIR / f"{cell}_subset_membership.csv")
    pool = pd.read_csv(POOL_DIR / f"{cell}_pooled.csv")
    pool_idx = pool.set_index(pool.columns[0])
    for mar, name in TARGETS.items():
        r = sub[sub.rxn == mar]
        if r.empty:
            rows.append(dict(cell=cell, rxn=mar, metabolite=name,
                             representative="", expansion_coef=np.nan,
                             median=np.nan, p05=np.nan, p95=np.nan,
                             fraction_negative=np.nan, note="not in cell"))
            continue
        rep = r.iloc[0].representative
        coef = r.iloc[0].expansion_coef
        if rep not in pool_idx.index:
            rows.append(dict(cell=cell, rxn=mar, metabolite=name,
                             representative=rep, expansion_coef=coef,
                             median=np.nan, p05=np.nan, p95=np.nan,
                             fraction_negative=np.nan,
                             note="representative not in pool"))
            continue
        signed = pool_idx.loc[rep].values * coef
        rows.append(dict(cell=cell, rxn=mar, metabolite=name,
                         representative=rep, expansion_coef=coef,
                         median=float(np.median(signed)),
                         p05=float(np.percentile(signed, 5)),
                         p95=float(np.percentile(signed, 95)),
                         min=float(signed.min()), max=float(signed.max()),
                         fraction_negative=float((signed < -1e-6).mean()),
                         note=""))

df = pd.DataFrame(rows)
df.to_csv(OUT / "sampled_uptake_offending.csv", index=False)
print(f"[write] sampled_uptake_offending.csv  rows={len(df)}", flush=True)

# Pivot: cell x metabolite, median
piv_med = df.pivot(index="cell", columns="metabolite", values="median")
piv_p95 = df.pivot(index="cell", columns="metabolite", values="p95")
piv_med.to_csv(OUT / "sampled_uptake_median_pivot.csv")
piv_p95.to_csv(OUT / "sampled_uptake_p95_pivot.csv")

print("\n=== Median signed flux per cell (negative = uptake) ===")
print(piv_med.to_string())
print("\n=== 95th-pct signed flux per cell ===")
print(piv_p95.to_string())
print("\n=== Fraction of samples with uptake (flux < 0) ===")
pf = df.pivot(index="cell", columns="metabolite", values="fraction_negative")
print(pf.to_string())
