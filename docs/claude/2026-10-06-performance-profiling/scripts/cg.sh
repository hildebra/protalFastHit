#!/usr/bin/env bash
# callgrind of HEAD, one thread: aligning 100k pairs (no profiling), profiling the 500k-pair SAM,
# PacBio 3 Mb and ONT 3 Mb whole runs, single-end 100k reads (no profiling).
set -uo pipefail
W=$HOME/perf6; B=$W/head/build/protal; DB=$HOME/bench071/V075/protal_db
P=$HOME/bench071/samples/points
mkdir -p $W/reads $W/cg
PE1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; PE2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
[ -s $W/reads/pe100k_R1.fq.gz ] || { zcat $PE1 | head -400000 | gzip -1 > $W/reads/pe100k_R1.fq.gz; zcat $PE2 | head -400000 | gzip -1 > $W/reads/pe100k_R2.fq.gz; }
cgrun() { # name args...
  local n=$1; shift
  rm -rf $W/cg/o.$n
  valgrind --tool=callgrind --callgrind-out-file=$W/cg/$n.out --compress-strings=no --compress-pos=no \
    $B --db $DB "$@" --prefix s -o $W/cg/o.$n -t 1 --no_qcmsa --verbose > $W/cg/$n.log 2>&1
  callgrind_annotate --inclusive=yes $W/cg/$n.out 2>/dev/null | head -200 > $W/cg/$n.incl.txt
  callgrind_annotate $W/cg/$n.out 2>/dev/null | head -150 > $W/cg/$n.self.txt
  echo "$n: $(grep -m1 'PROGRAM TOTALS' $W/cg/$n.self.txt) $(date +%T)"
}
cgrun align -1 $W/reads/pe100k_R1.fq.gz -2 $W/reads/pe100k_R2.fq.gz --read_type pe --no_profile
cgrun prof --profile_only $W/o.pe_t6.r1/s.sam.zst
cgrun pb -1 $P/pb_b3000000/sim/reads/pb_b3000000_s_1.fq.gz --read_type pb
cgrun ont -1 $P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz --read_type ont
cgrun se -1 $W/reads/pe100k_R1.fq.gz --read_type se --no_profile
echo CGDONE
