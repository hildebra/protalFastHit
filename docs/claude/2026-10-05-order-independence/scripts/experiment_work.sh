#!/usr/bin/env bash
# After the change: the unit tests; the change four times and with the scrambled maps (scramble_maps.pl) twice; 48e8cd0
# with the scrambled maps twice (what map order did before); each compared with run 1 of its build, and the change with
# 48e8cd0 (what changes once); the profiling stage alone on the change's run-1 SAMs with the reads shuffled.
HERE=$(cd "$(dirname "$0")" && pwd); W=$HOME/det-order; R=$W/runs
echo "== unit tests"; $W/work/tree/build/tests/protal_tests 2>&1 | tail -3
$HERE/runs.sh work 4
$HERE/runs.sh work-reorder 2
$HERE/runs.sh base-reorder 2
grep -h '^TEST' $R/work-reorder/pe1.log $R/base-reorder/pe1.log
for t in pe pb; do
  for i in 2 3 4; do $HERE/compare.sh $R/work/${t}1 $R/work/$t$i "work $t run 1 vs $i"; done
  for i in 1 2; do $HERE/compare.sh $R/work/${t}1 $R/work-reorder/$t$i "work $t run 1 vs scrambled maps run $i"; done
  for i in 1 2; do $HERE/compare.sh $R/base/${t}1 $R/base-reorder/$t$i "base $t run 1 vs scrambled maps run $i"; done
  for i in 1 2 3 4; do $HERE/compare.sh $R/base/${t}$i $R/work/${t}1 "base $t run $i vs work run 1"; done
  LINES_SHOWN=12 $HERE/profile_shuffled.sh work $R/work/${t}1 $t $R/work/shuffled_$t
done
