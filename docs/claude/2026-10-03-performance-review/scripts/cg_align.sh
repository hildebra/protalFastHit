#!/usr/bin/env bash
# callgrind of aligning 100k pairs (first 100k of rl150_p500000_s_1), --no_profile, one thread.
set -uo pipefail
W=$HOME/mt-work/perf4; B=$W/ref/build/protal; DB=$HOME/bench071/V073/protal_db
P=$HOME/bench071/samples/points/rl150_p500000/sim/reads
mkdir -p $W/reads
[ -s $W/reads/pe100k_R1.fq.gz ] || { zcat $P/rl150_p500000_s_1_R1.fq.gz | head -400000 | gzip -1 > $W/reads/pe100k_R1.fq.gz; zcat $P/rl150_p500000_s_1_R2.fq.gz | head -400000 | gzip -1 > $W/reads/pe100k_R2.fq.gz; }
rm -rf $W/cga.pe100k
valgrind --tool=callgrind --callgrind-out-file=$W/cg.align.out --compress-strings=no --compress-pos=no \
  $B --db $DB -1 $W/reads/pe100k_R1.fq.gz -2 $W/reads/pe100k_R2.fq.gz --read_type pe --no_profile --prefix s -o $W/cga.pe100k -t 1 --no_qcmsa > $W/cg.align.log 2>&1
callgrind_annotate --inclusive=yes $W/cg.align.out 2>/dev/null | head -150 > $W/cg.align.incl.txt
callgrind_annotate $W/cg.align.out 2>/dev/null | head -150 > $W/cg.align.self.txt
echo "== callgrind align done"
head -30 $W/cg.align.self.txt
