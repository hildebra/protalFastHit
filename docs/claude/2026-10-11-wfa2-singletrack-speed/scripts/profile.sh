#!/usr/bin/env bash
# profile.sh VARIANT: callgrind of one wfareplay round of VARIANT per sample; WFA2's functions, self and inclusive, as
# shares of wavefront_align's inclusive instructions (plus the aligner's setup, outside wavefront_align).
set -uo pipefail
W=$HOME/wfasimd; C=$W/calls; R=$W/replay; P=$W/prof; V=${1:-singletrack}; BIN=${BIN:-$R/wfareplay}; TAG=${TAG:-base}
mkdir -p $P
fn() { sed -nE 's/^ *([0-9,]+) \([^)]*\) +.*:([A-Za-z_][A-Za-z0-9_]*) \[.*/\1 \2/p' "$1" | tr -d ','; }
for s in ${SETS:-pe pb ont}; do
  o=$P/$s.$TAG.$(echo $V | tr / h)
  taskset -c 0-3 valgrind --tool=callgrind --callgrind-out-file=$o.out $BIN $C/$s.tsv 1 $P/res $V > /dev/null 2>&1
  callgrind_annotate --inclusive=yes --threshold=100 $o.out 2>/dev/null > $o.incl
  callgrind_annotate --threshold=100 $o.out 2>/dev/null > $o.self
  total=$(fn $o.incl | awk '$2=="wavefront_align" {print $1; exit}')
  echo "-- $s $V ($TAG): wavefront_align $total instructions (inclusive shares / self shares)"
  fn $o.incl | awk -v t=$total '$2 ~ /^(wavefront_(extend|compute|heuristic|termination|backtrace|aligner_reset|aligner_init|sequences_init|unialign_terminate|slab)|wf_heuristic)/ && $1/t > 0.02 && !seen[$2]++ {printf "   incl %5.1f%%  %s\n", 100*$1/t, $2}' | head -14
  fn $o.self | awk -v t=$total '$1/t > 0.02 {printf "   self %5.1f%%  %s\n", 100*$1/t, $2}' | head -14
done
