#!/usr/bin/env bash
# callgrind of HEAD: aligning 100k pairs (1 thread, no profiling), and profiling the 500k-pair SAM (1 thread).
set -uo pipefail
W=$HOME/perf-pass2; B=$W/head/build/protal; DB=$HOME/bench071/V073/protal_db; R=$W/reads
rm -rf $W/cga; valgrind --tool=callgrind --callgrind-out-file=$W/cg.align.out --compress-strings=no --compress-pos=no \
  $B --db $DB -1 $R/pe100k_R1.fq.gz -2 $R/pe100k_R2.fq.gz --read_type pe --no_profile --prefix s -o $W/cga -t 1 --no_qcmsa > $W/cg.align.log 2>&1
callgrind_annotate --inclusive=yes $W/cg.align.out 2>/dev/null | head -120 > $W/cg.align.incl.txt
callgrind_annotate $W/cg.align.out 2>/dev/null | head -80 > $W/cg.align.self.txt
rm -rf $W/cgp; valgrind --tool=callgrind --callgrind-out-file=$W/cg.prof.out --compress-strings=no --compress-pos=no \
  $B --db $DB --profile_only $W/o.pe500k_t1/s.sam.zst --prefix s -o $W/cgp -t 1 --no_qcmsa > $W/cg.prof.log 2>&1
callgrind_annotate --inclusive=yes $W/cg.prof.out 2>/dev/null | head -120 > $W/cg.prof.incl.txt
callgrind_annotate $W/cg.prof.out 2>/dev/null | head -80 > $W/cg.prof.self.txt
echo done
