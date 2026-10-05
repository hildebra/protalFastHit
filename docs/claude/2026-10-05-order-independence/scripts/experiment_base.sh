#!/usr/bin/env bash
# Before the change: 48e8cd0 four times, and with the genome map in another order (reorder_genomes.pl) twice; each
# compared with run 1, and the profiling stage alone on run 1's SAMs with the reads shuffled.
HERE=$(cd "$(dirname "$0")" && pwd); R=$HOME/det-order/runs
$HERE/runs.sh base 4
$HERE/runs.sh base-reorder 2
for t in pe pb; do
  for i in 2 3 4; do $HERE/compare.sh $R/base/${t}1 $R/base/$t$i "base $t run 1 vs $i"; done
  for i in 1 2; do $HERE/compare.sh $R/base/${t}1 $R/base-reorder/$t$i "base $t run 1 vs reordered genome map run $i"; done
  $HERE/profile_shuffled.sh base $R/base/${t}1 $t $R/base/shuffled_$t
done
