#!/usr/bin/env bash
# massif.sh PROTAL TAG N [EXTRA...] - heap profile (valgrind massif) of a --profile_only run over N copies of the w900
# SAM (memory.sh's cohort folders), to see what a sample keeps (docs/claude/2026-10-07-sam-combine). One thread: valgrind
# runs threads one at a time anyway. Writes $M/massif_TAG_N.out and the last snapshot's top allocation sites.
set -u
PROTAL=$1; TAG=$2; N=$3; shift 3
DB=${DB:-$HOME/protal-perf/db900n}
M=${M:-$HOME/samcombine/mem}
o=$M/massif_${TAG}_$N
rm -rf "$o"
taskset -c 0-3 nice -n 5 valgrind --tool=massif --massif-out-file="$o.out" --threshold=0.5 --detailed-freq=1 \
  "$PROTAL" --db "$DB" --profile_only "$M/cohort$N/*.sam.zst" -o "$o" -t 1 --no_qcmsa "$@" > "$o.log" 2>&1
echo "exit $?"
ms_print --threshold=0.5 "$o.out" > "$o.txt"
# The snapshot with the most heap after the profiling (the last detailed one before teardown is the end of the run).
grep -n "^ *[0-9]* *[0-9,]* *[0-9,]* *[0-9,]* *[0-9,]* *[0-9,]*$" "$o.txt" | tail -3
