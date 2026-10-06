#!/usr/bin/env bash
# callgrind with the -g build (line tables): aligning 100k pairs, PacBio 3 Mb, profiling the 500k SAM; then the
# 1- and 6-thread timed runs (round 2) of the plain build.
set -uo pipefail
W=$HOME/perf6; B=$W/headg/build/protal; DB=$HOME/bench071/V075/protal_db
P=$HOME/bench071/samples/points
mkdir -p $W/cgg
cgrun() { # name args...
  local n=$1; shift
  rm -rf $W/cgg/o.$n
  valgrind --tool=callgrind --callgrind-out-file=$W/cgg/$n.out \
    $B --db $DB "$@" --prefix s -o $W/cgg/o.$n -t 1 --no_qcmsa --verbose > $W/cgg/$n.log 2>&1
  echo "$n: done $(date +%T)"
}
cgrun align -1 $W/reads/pe100k_R1.fq.gz -2 $W/reads/pe100k_R2.fq.gz --read_type pe --no_profile
cgrun pb -1 $P/pb_b3000000/sim/reads/pb_b3000000_s_1.fq.gz --read_type pb --no_profile
cgrun prof --profile_only $W/o.pe_t6.r1/s.sam.zst
echo CGGDONE
bash $W/runs.sh 2
THREADS=6 bash $W/runs.sh 3
echo RUNSDONE
