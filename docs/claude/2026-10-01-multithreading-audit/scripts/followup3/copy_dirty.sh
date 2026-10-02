#!/usr/bin/env bash
# As in a run: the records file freshly written (its pages dirty, not synced) before it is appended behind a
# 2 MB header; copy_file_range on 1, 2 and 4 threads, alternated, five passes, niced.
set -uo pipefail
W=$HOME/mt-work/copy
cat /proc/loadavg
for rep in 1 2 3 4 5; do
  for t in 1 2 4; do
    rm -f $W/fresh.partial; cp $HOME/mt-audit/runs/b6/s1.sam.zst $W/fresh.partial
    nice $W/copybench $W/fresh.partial $W/out.sam.zst cfr $t 2097152
  done
done
rm -f $W/fresh.partial
echo DIRTY DONE
