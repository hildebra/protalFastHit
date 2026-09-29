#!/bin/bash
# Transparent huge pages for the malloc'd index, without a code change (glibc >= 2.35 tunable),
# alternated with runs without: thp_test.sh DB PLAIN_READS_DIR LABEL [THREADS] [REPS]
# Prints each run's AnonHugePages (from /proc) to confirm the tunable took effect.
source "$(dirname "$0")/env.sh"
db=$1; src=$2; ds=$3; t=${4:-1}; reps=${5:-2}
check_thp() {
  for i in $(seq 1 120); do
    pid=$(pgrep -n -f "$BIN --db $db"); [ -n "$pid" ] && grep -q "Start parallel" $PERF_DIR/runs/$1/stdout.log 2>/dev/null && break; sleep 0.5
  done
  [ -n "$pid" ] && grep AnonHugePages /proc/$pid/smaps_rollup 2>/dev/null | sed "s/^/    $1 /"
}
for rep in $(seq 1 $reps); do
  for mode in off on; do
    label=thp_${mode}_${ds}_t${t}_$rep
    check_thp $label &
    if [ $mode = on ]; then export GLIBC_TUNABLES=glibc.malloc.hugetlb=1; else unset GLIBC_TUNABLES; fi
    bash $here/run_one.sh $label $db $src $t --no_profile --verbose | grep -E "rc=|Seeding|Load Index took|Aligning reads"
    wait
  done
done
