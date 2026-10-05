#!/usr/bin/env bash
# The profiling stage's cost, 48e8cd0 against the change: --profile_only on run 1's SAMs at 6 threads, 5 alternated rounds;
# the "Profiling took" line of each.
W=$HOME/det-order; R=$W/runs; DB=$HOME/bench071/V073/protal_db/database.protal
for t in pe pb; do
  for round in 1 2 3 4 5; do
    for b in base work; do
      rm -rf $W/ab/$b$t
      $W/$b/tree/build/protal --db $DB --profile_only $R/base/${t}1/s.sam.zst --read_type $t --prefix s -o $W/ab/$b$t -t 6 --no_qcmsa > $W/ab.log 2>&1
      echo "$t $b $round $(grep -h '^Profiling took' $W/ab.log | sed 's/Profiling took //')"
    done
  done
done
