#!/usr/bin/env bash
# run_variants.sh: library variants built (variants.sh), then each sample's calls replayed by every variant binary in
# turn, ROUNDS times (binaries alternated, each running classic and singletrack); median thread CPU seconds per variant
# and mode; per-call results compared with the base build.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfasimd
W=$HOME/wfasimd; C=$W/calls; O=$W/varout; ROUNDS=${ROUNDS:-5}
mkdir -p $O
bash $S/variants.sh base "" x86-64
bash $S/variants.sh bounds bounds x86-64
bash $S/variants.sh bounds_term bounds,term x86-64
bash $S/variants.sh bounds_term_v2 bounds,term x86-64-v2
bash $S/variants.sh bounds_term_v3 bounds,term x86-64-v3
bash $S/variants.sh base_v3 "" x86-64-v3
NAMES="base bounds bounds_term bounds_term_v2 bounds_term_v3 base_v3"
for s in pe pb ont; do
  rm -f $O/$s.*.times
  for r in $(seq 1 $ROUNDS); do
    for n in $NAMES; do
      taskset -c 0-3 $W/libs/$n/wfareplay $C/$s.tsv 1 $O/$s.$n classic singletrack singletrack/2 2>/dev/null \
        | awk -v n=$n '/median/ {print n, $1, $3}' >> $O/$s.times
    done
  done
  echo "-- $s (median of $ROUNDS rounds, thread CPU seconds)"
  for n in $NAMES; do
    for m in classic singletrack singletrack/2; do
      awk -v n=$n -v m=$m '$1==n && $2==m {print $3}' $O/$s.times | sort -n | awk -v n=$n -v m=$m '{a[NR]=$1} END {printf "   %-16s %-14s %.4f\n", n, m, a[int((NR+1)/2)]}'
    done
  done
  for n in $NAMES; do
    for m in classic singletrack singletrack_half; do
      cmp -s $O/$s.base.$m $O/$s.$n.$m || echo "   $s $n $m: results DIFFER from base"
    done
  done
done
echo "VARDONE $(date +%T)"
