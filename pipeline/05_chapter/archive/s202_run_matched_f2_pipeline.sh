#!/usr/bin/env bash
# s202 -- run the matched-floor F2 pipeline end to end:
#   1. CFF (Gurobi) on 9 cells × 4 chains in hopsy_v7_matched_f2 → hopsy_v7_matched_f2_loopless
#   2. Pool each cell's 4 chains → pooled_v7_matched_f2/{cell}_pooled.csv
#   3. Run 6 F2 contrasts via two_tier → stats_v7_matched_f2/{A}_vs_{B}.csv
# After this: run `python scripts/s85_export_tier_tables.py` to update floor_resolved.

set -euo pipefail
cd "$(dirname "$0")/.."

CELLS="BL_24 BL_48 BL_72 D_24 D_48 D_72 LD_24 LD_48 LD_72"
IN_DIR="results/hopsy_v7_matched_f2"
LL_DIR="results/hopsy_v7_matched_f2_loopless"
POOL_DIR="results/pooled_v7_matched_f2"
STATS_DIR="results/stats_v7_matched_f2"

mkdir -p "$LL_DIR" "$POOL_DIR" "$STATS_DIR" logs/matched_f2

# ----- 1. CFF per cell (parallel across cells; chains serial within cell) -----
echo "=== STAGE 1: CFF (Gurobi) ==="
echo "$CELLS" | tr ' ' '\n' | xargs -I{} -P 4 bash -c '
  cell="$1"
  log="logs/matched_f2/cff_${cell}.log"
  if ls results/hopsy_v7_matched_f2_loopless/flux_samples_${cell}_ch4.csv >/dev/null 2>&1; then
    echo "[skip CFF] $cell"
    exit 0
  fi
  echo "[cff start] $cell"
  python -u scripts/s51b_cyclefreeflux_gurobi.py \
    --cell "$cell" \
    --outdir_in hopsy_v7_matched_f2 \
    --outdir_out hopsy_v7_matched_f2_loopless \
    > "$log" 2>&1 && echo "[cff done]  $cell" || { echo "[cff FAIL] $cell see $log"; exit 1; }
' _ {}

# ----- 2. Pool chains per cell -----
echo "=== STAGE 2: pool chains ==="
for cell in $CELLS; do
  outfp="$POOL_DIR/${cell}_pooled.csv"
  if [ -f "$outfp" ]; then
    echo "[skip pool] $cell"
    continue
  fi
  python -u - <<PY
import pandas as pd
from pathlib import Path
cell = "$cell"
indir = Path("$LL_DIR")
dfs = []
for ch in (1,2,3,4):
    d = pd.read_csv(indir / f"flux_samples_{cell}_ch{ch}.csv", index_col=0)
    d.columns = [f"ch{ch}_{c}" for c in d.columns]
    dfs.append(d)
pooled = pd.concat(dfs, axis=1)
pooled.index.name = "rxn"
pooled.to_csv("$outfp")
print(f"[pool done] {cell}  shape={pooled.shape}")
PY
done

# ----- 3. Run 6 F2 contrasts -----
echo "=== STAGE 3: six F2 matched-floor contrasts ==="
F2_CONTRASTS=(
  "BL_24 BL_48"  "BL_48 BL_72"
  "D_24  D_48"   "D_48  D_72"
  "LD_24 LD_48"  "LD_48 LD_72"
)
for pair in "${F2_CONTRASTS[@]}"; do
  A=$(echo $pair | awk '{print $1}')
  B=$(echo $pair | awk '{print $2}')
  OUT="$STATS_DIR/${A}_vs_${B}.csv"
  LOG="logs/matched_f2/stats_${A}_vs_${B}.log"
  if [ -f "$OUT" ]; then
    echo "[skip stats] ${A} vs ${B}"
    continue
  fi
  echo "[stats start] ${A} vs ${B}"
  python -u scripts/two_tier_subsystem_test_fast_mean_dropzero.py \
    "$POOL_DIR/${A}_pooled.csv" \
    "$POOL_DIR/${B}_pooled.csv" \
    "results/subsets_v7/${A}_subset_membership.csv" \
    "$OUT" > "$LOG" 2>&1 && echo "[stats done] ${A} vs ${B}" \
      || { echo "[stats FAIL] ${A} vs ${B} see $LOG"; exit 1; }
done

echo "=== DONE. Next: python scripts/s85_export_tier_tables.py (needs patch to read stats_v7_matched_f2) ==="
