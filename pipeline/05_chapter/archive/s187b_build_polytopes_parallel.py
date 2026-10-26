"""
s187b -- parallel polytope builder, 3-wide with Gurobi thread limits.

Each worker runs one s184 preprocess with GRB_THREADS=4. On a 14-core /
36 GB Mac this gives ~2× throughput over the sequential s187 driver with
~12 GB RAM headroom.

Usage:
  python s187b_build_polytopes_parallel.py --batch all
  python s187b_build_polytopes_parallel.py --batch F1_24_48 --jobs 3
  python s187b_build_polytopes_parallel.py --batch F2_late     # D/LD 96→120
  python s187b_build_polytopes_parallel.py --batch F1_72_96
  python s187b_build_polytopes_parallel.py --batch F1_120

Smallest-model-first ordering within each batch so a partial run leaves
the most useful corners cached.
"""
import argparse, subprocess, time, os
from pathlib import Path
from multiprocessing import Pool

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "grace_v4/results/v7_final_addenda/lipid_relay_samples"

# ---------- All 28 additional swap corners, grouped by batch -----------
# Each entry: (model_cell, env_cell, floor, batch_tag, approx_rxn_count)
# Floors are 0.9 × min(biomax_AA, biomax_AB, biomax_BA, biomax_BB) under
# frozen ±1000. Computed in s185 for the pairs where we have it; new
# pairs below use a safe lower bound (will be refined after biomax check
# at build time if needed — s184 enforces feasibility).

CORNERS = [
    # ======= Batch F1_24_48 — 8 polytopes =======
    # BL_24↔D_24
    ("BL_24", "D_24",  7.5, "F1_24_48"),  # placeholder floor; s184 verifies
    ("D_24",  "BL_24", 7.5, "F1_24_48"),
    # D_24↔LD_24
    ("D_24",  "LD_24", 7.5, "F1_24_48"),
    ("LD_24", "D_24",  7.5, "F1_24_48"),
    # BL_48↔LD_48
    ("BL_48", "LD_48", 6.5, "F1_24_48"),
    ("LD_48", "BL_48", 6.5, "F1_24_48"),
    # D_48↔LD_48
    ("D_48",  "LD_48", 6.5, "F1_24_48"),
    ("LD_48", "D_48",  6.5, "F1_24_48"),

    # ======= Batch F2_late — 4 polytopes =======
    ("D_96",  "D_120", 10.0, "F2_late"),
    ("D_120", "D_96",  10.0, "F2_late"),
    ("LD_96", "LD_120",11.0, "F2_late"),
    ("LD_120","LD_96", 11.0, "F2_late"),

    # ======= Batch F1_72_96 — 10 polytopes =======
    # BL_72↔D_72
    ("BL_72", "D_72",  7.5, "F1_72_96"),
    ("D_72",  "BL_72", 7.5, "F1_72_96"),
    # BL_72↔LD_72
    ("BL_72", "LD_72", 7.5, "F1_72_96"),
    ("LD_72", "BL_72", 7.5, "F1_72_96"),
    # D_72↔LD_72
    ("D_72",  "LD_72", 7.5, "F1_72_96"),
    ("LD_72", "D_72",  7.5, "F1_72_96"),
    # BL_96↔LD_96
    ("BL_96", "LD_96", 7.5, "F1_72_96"),
    ("LD_96", "BL_96", 7.5, "F1_72_96"),
    # D_96↔LD_96
    ("D_96",  "LD_96", 10.0, "F1_72_96"),
    ("LD_96", "D_96",  10.0, "F1_72_96"),

    # ======= Batch F1_120 — 6 polytopes =======
    # BL_120↔D_120
    ("BL_120","D_120", 10.0, "F1_120"),
    ("D_120", "BL_120",10.0, "F1_120"),
    # BL_120↔LD_120
    ("BL_120","LD_120",10.0, "F1_120"),
    ("LD_120","BL_120",10.0, "F1_120"),
    # D_120↔LD_120
    ("D_120", "LD_120",10.0, "F1_120"),
    ("LD_120","D_120", 10.0, "F1_120"),
]

BATCHES = ["F1_24_48", "F2_late", "F1_72_96", "F1_120"]


def _worker(args):
    """Build one polytope corner. Returns (tag, status, elapsed_min)."""
    model_cell, env_cell, floor, batch = args
    tag = f"{model_cell}_mdl__{env_cell}_env"
    cache = OUT / "polytopes" / f"{tag}_polytope.pkl"
    if cache.exists():
        return (tag, "cached", 0.0, batch)
    env = os.environ.copy()
    env["GRB_THREADS"] = "4"   # limit per-process Gurobi parallelism
    t0 = time.time()
    cmd = [
        "python",
        str(REPO / "grace_v4/scripts/s184_lipid_relay_swap_sample.py"),
        "--model_cell", model_cell,
        "--env_cell", env_cell,
        "--floor", str(floor),
        "--out_dir", str(OUT),
        "--preprocess_only",
    ]
    try:
        subprocess.run(cmd, env=env, check=True, cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       timeout=3600)
        return (tag, "OK" if cache.exists() else "FAIL_no_cache",
                (time.time() - t0) / 60, batch)
    except subprocess.CalledProcessError as e:
        return (tag, f"FAIL_rc={e.returncode}", (time.time() - t0) / 60, batch)
    except subprocess.TimeoutExpired:
        return (tag, "TIMEOUT_1h", 60.0, batch)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", choices=BATCHES + ["all"], default="all",
                    help="which batch to run")
    ap.add_argument("--jobs", type=int, default=3,
                    help="parallel workers (default 3 = ~12 GB RAM, ~12 cores)")
    args = ap.parse_args()

    tasks = CORNERS if args.batch == "all" else [c for c in CORNERS if c[3] == args.batch]
    # smallest-first proxy: sort by env hpf, then model hpf (early hpf = smaller)
    def _hpf(c): return int(c.split("_")[1])
    tasks.sort(key=lambda t: (_hpf(t[0]), _hpf(t[1])))

    print(f"[s187b] batch={args.batch}  tasks={len(tasks)}  jobs={args.jobs}")
    for t in tasks:
        print(f"  queue  {t[0]:7s} mdl + {t[1]:7s} env  floor={t[2]}  batch={t[3]}")
    print()
    t0 = time.time()
    with Pool(args.jobs) as pool:
        for i, (tag, status, mins, batch) in enumerate(pool.imap_unordered(_worker, tasks), 1):
            elapsed_total = (time.time() - t0) / 60
            print(f"  [{i:2d}/{len(tasks)}]  {tag}  {status}  ({mins:.1f} min, "
                  f"total {elapsed_total:.1f} min, batch={batch})", flush=True)
    print(f"\n[s187b] DONE  total {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
