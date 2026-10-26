"""
s20b_rank_rhat -- rank-normalized split-Rhat + bulk/tail ESS (Vehtari 2021)
on EXISTING chain CSVs.  No resampling.

Flexible layout: pass --outdir + --suffix so this reads
    {outdir}/flux_samples_{cell}{suffix}_ch{1..n_chains}.csv

Examples (from grace_v4/):
    python scripts/s20b_rank_rhat.py --outdir results/hopsy_v7_local --cells all
    python scripts/s20b_rank_rhat.py --outdir ../ocr_anchored_extraction/results/chrr_pilot --cells LD_72 --suffix _optA
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd, arviz as az

ap = argparse.ArgumentParser()
ap.add_argument("--outdir", required=True)
ap.add_argument("--cells", required=True, help="comma-sep list or 'all'")
ap.add_argument("--n_chains", type=int, default=4)
ap.add_argument("--suffix", default="")
ap.add_argument("--label", default="")
a = ap.parse_args()

ALL_CELLS = [f"{c}_{h}" for c in ("BL","D","LD") for h in (24,48,72,96,120)]
cells = ALL_CELLS if a.cells == "all" else [c.strip() for c in a.cells.split(",")]
outdir = Path(a.outdir).resolve()

label = a.label or outdir.name
rows = []
for cell in cells:
    chain_paths = [outdir / f"flux_samples_{cell}{a.suffix}_ch{ch}.csv"
                   for ch in range(1, a.n_chains + 1)]
    missing = [p for p in chain_paths if not p.exists()]
    if missing:
        rows.append(dict(label=label, cell=cell, status=f"missing {len(missing)}/{a.n_chains}"))
        continue
    chains = [pd.read_csv(p, index_col=0) for p in chain_paths]
    rxns = chains[0].index
    # align shape (chain, draw, rxn)
    x = np.stack([c.reindex(rxns).values.T for c in chains])
    var_ok = x.reshape(-1, x.shape[2]).std(axis=0) > 1e-12
    x_ok, rxns_ok = x[:, :, var_ok], rxns[var_ok]
    ds = az.convert_to_dataset({"v": x_ok})
    rhat  = az.rhat(ds, method="rank")["v"].values
    ess_b = az.ess(ds, method="bulk")["v"].values
    ess_t = az.ess(ds, method="tail")["v"].values
    pd.DataFrame({"rxn": rxns_ok, "rhat_rank": rhat,
                  "ess_bulk": ess_b, "ess_tail": ess_t}) \
      .to_csv(outdir / f"{cell}_rhat_rank.csv", index=False)
    rows.append(dict(
        label=label, cell=cell,
        n_rxns=len(rxns), n_nonfixed=len(rxns_ok),
        n_samples_per_chain=x.shape[1],
        rhat_median=float(np.median(rhat)),
        rhat_max=float(np.max(rhat)),
        pct_rhat_lt_1_01=float(100 * np.mean(rhat < 1.01)),
        pct_rhat_lt_1_05=float(100 * np.mean(rhat < 1.05)),
        pct_rhat_gt_1_10=float(100 * np.mean(rhat > 1.10)),
        pct_rhat_gt_1_20=float(100 * np.mean(rhat > 1.20)),
        ess_bulk_median=float(np.median(ess_b)),
        ess_bulk_p10=float(np.quantile(ess_b, .1)),
        ess_tail_median=float(np.median(ess_t)),
        pct_ess_bulk_gt_400=float(100 * np.mean(ess_b > 400)),
        pct_ess_bulk_lt_100=float(100 * np.mean(ess_b < 100)),
    ))
    r = rows[-1]
    print(f"[{label}:{cell}] "
          f"n_nonfixed={r['n_nonfixed']}  n_samp/chain={r['n_samples_per_chain']}  "
          f"rhat_rank: median={r['rhat_median']:.3f} max={r['rhat_max']:.2f}  "
          f"<1.01={r['pct_rhat_lt_1_01']:.1f}%  <1.05={r['pct_rhat_lt_1_05']:.1f}%  "
          f">1.10={r['pct_rhat_gt_1_10']:.1f}%  >1.20={r['pct_rhat_gt_1_20']:.1f}%  "
          f"| ESS_bulk med={r['ess_bulk_median']:.0f}  p10={r['ess_bulk_p10']:.0f}  "
          f">400={r['pct_ess_bulk_gt_400']:.1f}%  <100={r['pct_ess_bulk_lt_100']:.1f}%",
          flush=True)

df = pd.DataFrame(rows)
out_csv = outdir / f"rank_rhat_summary_{label}.csv"
df.round(3).to_csv(out_csv, index=False)
print(f"\nwrote {out_csv}")
