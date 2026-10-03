#!/usr/bin/env bash
W=$HOME/perf-pass2; B=$W/work/build/protal; DB=$HOME/bench071/V073/protal_db; R=$W/reads
rm -rf $W/cga2; valgrind --tool=callgrind --callgrind-out-file=$W/cg2.align.out --compress-strings=no --compress-pos=no \
  $B --db $DB -1 $R/pe100k_R1.fq.gz -2 $R/pe100k_R2.fq.gz --read_type pe --no_profile --prefix s -o $W/cga2 -t 1 --no_qcmsa > $W/cg2.align.log 2>&1
callgrind_annotate --inclusive=yes $W/cg2.align.out 2>/dev/null | head -60 > $W/cg2.align.incl.txt
rm -rf $W/cgp2; valgrind --tool=callgrind --callgrind-out-file=$W/cg2.prof.out --compress-strings=no --compress-pos=no \
  $B --db $DB --profile_only $W/a.head.1/s.sam.zst --prefix s -o $W/cgp2 -t 1 --no_qcmsa > $W/cg2.prof.log 2>&1
callgrind_annotate --inclusive=yes $W/cg2.prof.out 2>/dev/null | head -60 > $W/cg2.prof.incl.txt
grep -m1 "PROGRAM TOTALS" $W/cg.align.incl.txt $W/cg2.align.incl.txt $W/cg.prof.incl.txt $W/cg2.prof.incl.txt
