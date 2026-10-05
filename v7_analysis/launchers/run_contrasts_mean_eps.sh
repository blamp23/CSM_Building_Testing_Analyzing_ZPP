#!/usr/bin/env bash
# Launch all 42 contrasts (27 primary + 15 neg-control split-chain) via xargs -P 12.
# Inputs:
#   results/pooled_v7/{cell}_pooled.csv       (primary contrasts)
#   results/pooled_v7/{cell}_halfA.csv        (neg control)
#   results/pooled_v7/{cell}_halfB.csv
#   results/subsets_v7/{cell}_subset_membership.csv
# Outputs:
#   results/stats_v7/{A}_vs_{B}.csv
#   results/stats_v7/neg_{cell}.csv
#   logs/stats_v7/{A}_vs_{B}.log

set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results/stats_v7 logs/stats_v7

HPFS="24 48 72 96 120"
CONDS="BL D LD"

# Build contrast list: lines of "A B membership out_prefix"
TMP=$(mktemp)
# A) three-way per hpf: BL-D, BL-LD, D-LD at each hpf
for h in $HPFS; do
  for pair in "BL D" "BL LD" "D LD"; do
    A="$(echo $pair | cut -d' ' -f1)_${h}"
    B="$(echo $pair | cut -d' ' -f2)_${h}"
    echo "$A $B $A ${A}_vs_${B}" >> $TMP
  done
done
# B) within-condition consecutive timewise
for c in $CONDS; do
  prev=""
  for h in $HPFS; do
    if [ -n "$prev" ]; then
      A="${c}_${prev}"; B="${c}_${h}"
      echo "$A $B $A ${A}_vs_${B}" >> $TMP
    fi
    prev=$h
  done
done
# C) neg-control split-chain (halfA vs halfB per cell), using cell's membership
NEG=$(mktemp)
for c in $CONDS; do
  for h in $HPFS; do
    cell="${c}_${h}"
    echo "$cell $cell neg_${cell}" >> $NEG
  done
done

echo "Primary contrasts: $(wc -l < $TMP)"
echo "Neg controls:      $(wc -l < $NEG)"

# Primary launcher: positional args A B MEM OUT
cat $TMP | xargs -n4 -P 12 bash -c '
  A=$1; B=$2; MEM=$3; OUT=$4
  OUTCSV=results/stats_v7/${OUT}.csv
  LOG=logs/stats_v7/${OUT}.log
  if [ -f "$OUTCSV" ]; then echo "[skip exists] $OUT"; exit 0; fi
  python -u scripts/two_tier_subsystem_test_fast.py \
    results/pooled_v7/${A}_pooled.csv \
    results/pooled_v7/${B}_pooled.csv \
    results/subsets_v7/${MEM}_subset_membership.csv \
    $OUTCSV > $LOG 2>&1 && echo "[done] $OUT" || echo "[FAIL] $OUT ($LOG)"
' _

# Neg control launcher: positional args cell cell outname
cat $NEG | xargs -n3 -P 12 bash -c '
  CELL=$1; OUT=$3
  OUTCSV=results/stats_v7/${OUT}.csv
  LOG=logs/stats_v7/${OUT}.log
  if [ -f "$OUTCSV" ]; then echo "[skip exists] $OUT"; exit 0; fi
  python -u scripts/two_tier_subsystem_test_fast.py \
    results/pooled_v7/${CELL}_halfA.csv \
    results/pooled_v7/${CELL}_halfB.csv \
    results/subsets_v7/${CELL}_subset_membership.csv \
    $OUTCSV > $LOG 2>&1 && echo "[done] $OUT" || echo "[FAIL] $OUT ($LOG)"
' _

echo "ALL CONTRASTS COMPLETE"
rm -f $TMP $NEG
