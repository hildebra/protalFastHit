#!/usr/bin/env bash
# The scheduling change (protal-par2): the multi-sample outputs against 0.7.1's (compare_profiles.sh's), and the
# eight full-database SAMs profiled together at -t 6, alternated with 0.7.1, three times (rows into speed/runs.tsv).
set -uo pipefail
cp $HOME/mt-work/sched/src/build/protal $HOME/mt-work/bin/protal-par2
C=$HOME/mt-work/compare
B=$HOME/bench071/runs
# The SAM lists as compare_profiles.sh built them (its associative arrays, in their order).
eval "$(sed -n '/^declare -A FULL=/,/^DBMISSING=/p' "$(dirname "$0")/compare_profiles.sh")"
for set in FULL MISSING; do
  declare -n M=$set
  db=$([ $set = FULL ] && echo $DBFULL || echo $DBMISSING)
  all=""; prefixes=""
  for sample in "${!M[@]}"; do all+="${all:+,}${M[$sample]}"; prefixes+="${prefixes:+,}$sample"; done
  for t in 6 16 3; do
    d=$C/$set.all.par2_$t; rm -rf $d; mkdir -p $d
    nice -n 5 $HOME/mt-work/bin/protal-par2 --db $db --profile_only $all --prefix $prefixes -o $d/out -t $t --no_qcmsa > $d/log 2>&1
    if diff -r -q $C/$set.all.base/out $d/out > /dev/null; then echo "SAME $set.all t$t"; else echo "DIFFER $set.all t$t"; fi
  done
done
