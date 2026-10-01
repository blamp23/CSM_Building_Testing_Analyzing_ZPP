"""
Split-chain check for reviewer v4 (task #105).

One SLURM array task per contrast (three total):
  1: LD_72_halfA vs LD_72_halfB     (negative control -> every z should be ~0)
  2: D_72_halfA  vs LD_72_halfA     (headline computed from chains {1,2})
  3: D_72_halfB  vs LD_72_halfB     (headline computed from chains {3,4})

Each task pools its own inputs from results/hopsy_loopless/ with skip-if-exists,
then runs scripts/two_tier_subsystem_test.py at ndraw=2000. Concurrent execution
across three nodes; wall time per task ~= one contrast in s52b_twotier.

Report step (scripts/s60_split_chain_report.py) runs after the array finishes
and produces cross_half_headline.csv + neg_control_summary.txt.

Usage:
    python s60_split_chain.py --task {1,2,3}
"""
import argparse, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

RESULTS = Path("results")
LOOPLESS = RESULTS / "hopsy_loopless"
SUBSETS = RESULTS / "subsets_v4"
OUT = RESULTS / "split_chain_v4"
OUT.mkdir(parents=True, exist_ok=True)

MEMB = SUBSETS / "LD_72_subset_membership.csv"

TASKS = {
    1: ("LD72_halfA_vs_halfB", ("LD_72", "halfA"), ("LD_72", "halfB")),
    2: ("D72_vs_LD72_halfA",   ("D_72",  "halfA"), ("LD_72", "halfA")),
    3: ("D72_vs_LD72_halfB",   ("D_72",  "halfB"), ("LD_72", "halfB")),
}
HALF_CHAINS = {"halfA": [1, 2], "halfB": [3, 4]}


def pool_half(cell, tag):
    out_csv = OUT / f"reps_pooled_{cell}_{tag}.csv"
    if out_csv.exists():
        print(f"  skip pool {cell}_{tag} (exists)")
        return out_csv
    memb = pd.read_csv(SUBSETS / f"{cell}_subset_membership.csv")
    rep_ids = memb.drop_duplicates("representative").representative.tolist()
    parts = []
    for ch in HALF_CHAINS[tag]:
        fp = LOOPLESS / f"flux_samples_{cell}_ch{ch}.csv"
        if not fp.exists():
            sys.exit(f"missing {fp}")
        df = pd.read_csv(fp, index_col=0)
        parts.append(df.reindex(rep_ids).values)
    pooled = np.concatenate(parts, axis=1)
    dfo = pd.DataFrame(pooled, index=rep_ids,
                       columns=[f"s{i+1}" for i in range(pooled.shape[1])])
    dfo.index.name = "rxn"
    dfo.to_csv(out_csv, float_format="%.5g")
    print(f"  pooled {cell}_{tag}: {dfo.shape} -> {out_csv.name}")
    return out_csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", type=int, required=True, choices=[1, 2, 3])
    ap.add_argument("--ndraw", type=int, default=2000)
    args = ap.parse_args()

    name, (cell_a, tag_a), (cell_b, tag_b) = TASKS[args.task]
    print(f"==== task {args.task}: {name} ====")

    a_csv = pool_half(cell_a, tag_a)
    b_csv = pool_half(cell_b, tag_b)

    out_csv = OUT / f"two_tier_{name}.csv"
    print(f"  running two-tier -> {out_csv.name}")
    subprocess.run([sys.executable, "-u", "scripts/two_tier_subsystem_test.py",
                    str(a_csv), str(b_csv), str(MEMB), str(out_csv),
                    "--ndraw", str(args.ndraw)], check=True)
    print(f"==== task {args.task} DONE ====")


if __name__ == "__main__":
    main()
