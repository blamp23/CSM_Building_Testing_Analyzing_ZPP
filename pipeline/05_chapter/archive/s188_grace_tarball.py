"""
s188 -- package a lipid-relay (or F1/F2-late) batch for Grace dispatch.

Produces a single tarball:
  lipid_relay_<batch>_<date>.tar.gz
    s184_lipid_relay_swap_sample.py
    run_lipid_relay.slurm                   (SLURM array file for this batch)
    manifest.csv                            (one row per sampling task)
    polytopes/*.pkl                         (polytope caches for this batch)
    README_grace.md                         (how to dispatch)

Usage:
  python s188_grace_tarball.py --batch lipid14
  python s188_grace_tarball.py --batch F1_24_48
  python s188_grace_tarball.py --batch F2_late
  python s188_grace_tarball.py --batch F1_72_96
  python s188_grace_tarball.py --batch F1_120
  python s188_grace_tarball.py --batch all          (one giant tarball)
"""
import argparse, shutil, subprocess, datetime, csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
POLY_DIR = REPO / "grace_v4/results/v7_final_addenda/lipid_relay_samples/polytopes"
S184 = REPO / "grace_v4/scripts/s184_lipid_relay_swap_sample.py"

BATCH_CORNERS = {
    "lipid14": [
        ("D_24",  "D_48",  8.48), ("D_48",  "D_24",  8.48),
        ("LD_24", "LD_48", 8.48), ("LD_48", "LD_24", 8.48),
        ("D_48",  "D_72",  8.36), ("D_72",  "D_48",  8.36),
        ("LD_48", "LD_72", 8.43), ("LD_72", "LD_48", 8.43),
        ("BL_48", "BL_72", 8.23), ("BL_72", "BL_48", 8.23),
        ("D_72",  "D_96",  10.11), ("D_96", "D_72",  10.11),
        ("LD_72", "LD_96", 11.41), ("LD_96","LD_72", 11.41),
    ],
    "F1_24_48": [
        ("BL_24", "D_24",  7.5), ("D_24",  "BL_24", 7.5),
        ("D_24",  "LD_24", 7.5), ("LD_24", "D_24",  7.5),
        ("BL_48", "LD_48", 6.5), ("LD_48", "BL_48", 6.5),
        ("D_48",  "LD_48", 6.5), ("LD_48", "D_48",  6.5),
    ],
    "F2_late": [
        ("D_96",  "D_120", 10.0), ("D_120", "D_96",  10.0),
        ("LD_96", "LD_120",11.0), ("LD_120","LD_96", 11.0),
    ],
    "F1_72_96": [
        ("BL_72", "D_72",  7.5), ("D_72",  "BL_72", 7.5),
        ("BL_72", "LD_72", 7.5), ("LD_72", "BL_72", 7.5),
        ("D_72",  "LD_72", 7.5), ("LD_72", "D_72",  7.5),
        ("BL_96", "LD_96", 7.5), ("LD_96", "BL_96", 7.5),
        ("D_96",  "LD_96", 10.0), ("LD_96", "D_96", 10.0),
    ],
    "F1_120": [
        ("BL_120","D_120", 10.0), ("D_120", "BL_120",10.0),
        ("BL_120","LD_120",10.0), ("LD_120","BL_120",10.0),
        ("D_120", "LD_120",10.0), ("LD_120","D_120", 10.0),
    ],
}


SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --job-name={batch}
#SBATCH --array=1-{n_tasks}%12
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/{batch}_%A_%a.out
#SBATCH --error=logs/{batch}_%A_%a.err

source "$HOME/.bashrc"
export V7_REPO=$SCRATCH/{batch}
cd "$V7_REPO"

# Manifest: one row per task
MANIFEST=$(cat manifest.csv | tail -n +2)
IDX=$((SLURM_ARRAY_TASK_ID))
LINE=$(echo "$MANIFEST" | sed -n "${{IDX}}p")
IFS=',' read -r MDL ENV FLOOR CHAIN <<<"$LINE"
echo "[task $SLURM_ARRAY_TASK_ID] model=$MDL env=$ENV floor=$FLOOR chain=$CHAIN"

python s184_lipid_relay_swap_sample.py \\
    --model_cell "$MDL" --env_cell "$ENV" --floor "$FLOOR" \\
    --out_dir "." --chain_only "$CHAIN" \\
    --n_chains 4 --n_samples 500 --thinning 300
"""


README = """# Grace dispatch — lipid-relay / F1 / F2-late swap sampling

## Contents
- `s184_lipid_relay_swap_sample.py`  : chain sampler (`--chain_only N` mode).
- `polytopes/*.pkl`                  : pre-rounded polytope caches (one per model×env combo).
- `manifest.csv`                     : one row per sampling task (model, env, floor, chain).
- `run_{batch}.slurm`                : SLURM array file; dispatches 4 chains per corner.

## Dispatch
```
mkdir -p $SCRATCH/{batch} && cd $SCRATCH/{batch}
scp ... .  # copy this tarball
tar xzf lipid_relay_{batch}_*.tar.gz
mkdir -p logs
sbatch run_{batch}.slurm
```

Each chain writes `chains/{{model_cell}}_mdl__{{env_cell}}_env_ch{{N}}.csv`.
Download the `chains/` directory to local
`grace_v4/results/v7_final_addenda/lipid_relay_samples/chains/`, then run
`python grace_v4/scripts/s186_lipid_relay_decompose.py` (updated to pool
the four chains per corner).

## Convention notes
- Internal reaction bounds: Wang-default ±1000 (match frozen pooled_v7).
- Matched floor: pre-computed per pair from `s178_lipid_relay_baseline.csv`;
  stored in manifest.csv column 3.
- Thinning 300; 500 samples/chain; 4 chains per corner.
- hopsy billiard-walk; CycleFreeFlux NOT applied per sample (chapter method
  of record §11.14 uses pooled mean + drop zeros, robust to per-sample loops).
"""


def pack(batch_name):
    corners = BATCH_CORNERS[batch_name]
    out_dir = REPO / "grace_v4/results/v7_final_addenda" / f"grace_tarball_{batch_name}"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "polytopes").mkdir(exist_ok=True)

    # --- Copy polytopes
    missing = []
    for mdl, env, _ in corners:
        tag = f"{mdl}_mdl__{env}_env"
        src = POLY_DIR / f"{tag}_polytope.pkl"
        if not src.exists():
            missing.append(tag)
            continue
        shutil.copy(src, out_dir / "polytopes" / src.name)
        # also copy meta if present
        meta_src = POLY_DIR / f"{tag}_meta.json"
        if meta_src.exists():
            shutil.copy(meta_src, out_dir / "polytopes" / meta_src.name)
    if missing:
        print(f"  [warn] {len(missing)} polytopes missing — not yet built:")
        for t in missing: print(f"    {t}")
        print(f"  Partial tarball will be created.  Rerun after missing polytopes build.")

    # --- Copy sampler script
    shutil.copy(S184, out_dir / "s184_lipid_relay_swap_sample.py")

    # --- Write manifest
    tasks = []
    for mdl, env, floor in corners:
        tag = f"{mdl}_mdl__{env}_env"
        if (out_dir / "polytopes" / f"{tag}_polytope.pkl").exists():
            for chain in (1, 2, 3, 4):
                tasks.append((mdl, env, floor, chain))
    with open(out_dir / "manifest.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model_cell", "env_cell", "floor", "chain"])
        for t in tasks:
            w.writerow(t)

    # --- Write SLURM
    slurm = SLURM_TEMPLATE.format(batch=batch_name, n_tasks=len(tasks))
    (out_dir / f"run_{batch_name}.slurm").write_text(slurm)

    # --- README
    (out_dir / "README_grace.md").write_text(README.replace("{batch}", batch_name))

    # --- Tar
    date = datetime.date.today().isoformat().replace("-", "")
    tarball = REPO / "grace_v4/results/v7_final_addenda" / f"lipid_relay_{batch_name}_{date}.tar.gz"
    subprocess.run(["tar", "czf", str(tarball), "-C", str(out_dir.parent),
                    out_dir.name], check=True)
    print(f"\n[pack] {tarball}")
    print(f"       {len(tasks)} tasks  ({len(corners)} polytopes × 4 chains)")
    if missing:
        print(f"       (!) {len(missing)} polytopes still pending")
    return tarball, missing


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", choices=list(BATCH_CORNERS) + ["all"], required=True)
    args = ap.parse_args()
    if args.batch == "all":
        for b in BATCH_CORNERS:
            print(f"\n============ {b} ============")
            pack(b)
    else:
        pack(args.batch)
