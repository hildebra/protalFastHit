#!/usr/bin/env bash
# timing.sh: base, bounds_term_v2 and bt_v3x (v3 without WFA2's AVX2 extend kernels) replayed in turn, ROUNDS rounds,
# each run timing classic, singletrack and singletrack/2 in the same process; medians of thread CPU seconds.
set -uo pipefail
W=$HOME/wfasimd; C=$W/calls; O=$W/varout; ROUNDS=${ROUNDS:-9}
NAMES="base bounds_term bounds_term_v2 bt_v3x"
for s in pe pb ont; do
  rm -f $O/$s.t2
  for r in $(seq 1 $ROUNDS); do
    for n in $NAMES; do
      taskset -c 0-3 $W/libs/$n/wfareplay $C/$s.tsv 1 $O/t2.$s.$n classic singletrack singletrack/2 classic/2 2>/dev/null \
        | awk -v n=$n '/median/ {print n, $1, $3}' >> $O/$s.t2
    done
  done
  echo "-- $s: median thread CPU seconds of $ROUNDS rounds (vs base classic | vs base singletrack)"
  bc=$(awk '$1=="base" && $2=="classic" {print $3}' $O/$s.t2 | sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')
  bs=$(awk '$1=="base" && $2=="singletrack" {print $3}' $O/$s.t2 | sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')
  for n in $NAMES; do
    for m in classic classic/2 singletrack singletrack/2; do
      awk -v n=$n -v m=$m '$1==n && $2==m {print $3}' $O/$s.t2 | sort -n | awk -v n=$n -v m=$m -v bc=$bc -v bs=$bs \
        '{a[NR]=$1} END {x=a[int((NR+1)/2)]; printf "   %-16s %-14s %.4f  %+5.0f%% | %+5.0f%%\n", n, m, x, 100*(x-bc)/bc, 100*(x-bs)/bs}'
    done
  done
done
uptime
echo "TIMINGDONE $(date +%T)"
