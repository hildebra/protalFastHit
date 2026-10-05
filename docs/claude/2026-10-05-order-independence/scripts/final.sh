#!/usr/bin/env bash
# The final state of the change: rebuilt, unit tests, four runs (each compared with run 1 and with 48e8cd0's runs), the
# shuffled-SAM profiles, the profiling A/B against 48e8cd0.
HERE=$(cd "$(dirname "$0")" && pwd); R=$HOME/det-order/runs
$HERE/build.sh work tree || exit 1
$HERE/check_tests.sh work
$HERE/runs.sh work 4
for t in pe pb; do
  for i in 2 3 4; do $HERE/compare.sh $R/work/${t}1 $R/work/$t$i "work $t run 1 vs $i"; done
  $HERE/compare.sh $R/work/${t}1 $R/work-reorder/${t}1 "work $t run 1 vs scrambled maps run 1 (built before the single lookup)"
  for i in 1 2 3 4; do $HERE/compare.sh $R/base/${t}$i $R/work/${t}1 "base $t run $i vs work run 1"; done
  LINES_SHOWN=12 $HERE/profile_shuffled.sh work $R/work/${t}1 $t $R/work/shuffled_$t
done
$HERE/profile_ab.sh
