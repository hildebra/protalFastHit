#!/usr/bin/env bash
# The committed src (a212559; was run as f359d49 + its diff): clean build and unit tests, f359d49's own build, and byte comparisons
# of --profile_only outputs for three samples at 1, 6 and 8 threads.
set -uo pipefail
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
bash $SP/check.sh final2 "" a212559 || exit 1
bash $SP/build_wt.sh fbase2 /dev/null f359d49 protal || exit 1
W=$HOME/mt-work/fcmp; rm -rf $W; mkdir -p $W
B=$HOME/bench071/runs; DB=$HOME/bench071/V071/protal_db
for s in v071.full.pe.rl150_p500000_s_1/rl150_p500000_s_1 v071.full.ont.ont_b90000000_s_1/ont_b90000000_s_1 v071.full.se.rl150_p500000_s_1/rl150_p500000_s_1; do
  n=$(basename $s)
  nice $HOME/mt-work/fbase2/src/build/protal --db $DB --profile_only $B/$s.sam.zst --prefix $n -o $W/base -t 6 --no_qcmsa > /dev/null 2>&1
  for t in 1 6 8; do
    rm -rf $W/new; nice $HOME/mt-work/final2/src/build/protal --db $DB --profile_only $B/$s.sam.zst --prefix $n -o $W/new -t $t --no_qcmsa > /dev/null 2>&1
    diff -r -q $W/base $W/new > /dev/null && echo "SAME $n t$t" || echo "DIFFER $n t$t"
  done
  rm -rf $W/base
done
echo FINAL DONE
