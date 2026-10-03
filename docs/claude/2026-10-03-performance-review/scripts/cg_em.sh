#!/usr/bin/env bash
# callgrind of --profile_only on the 500k SAM with the EM on indices (em build), one thread.
set -uo pipefail
W=$HOME/mt-work/perf4; B=$W/em/build/protal; DB=$HOME/bench071/V073/protal_db
rm -rf $W/cgp.em
valgrind --tool=callgrind --callgrind-out-file=$W/cg.em.out --compress-strings=no --compress-pos=no \
  $B --db $DB --profile_only $W/o.pe500k_t1/s.sam.zst --prefix s -o $W/cgp.em -t 1 --no_qcmsa > $W/cg.em.log 2>&1
callgrind_annotate --inclusive=yes $W/cg.em.out 2>/dev/null | head -60 > $W/cg.em.incl.txt
grep -E 'PROGRAM TOTALS|AbundanceWeightedShares|ApplySampleContext|pow@|CongenerDistances::Compared' $W/cg.em.incl.txt | cut -c1-160
