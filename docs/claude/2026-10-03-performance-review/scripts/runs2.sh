#!/usr/bin/env bash
# Output-handler check (plain SAM, 3 threads) and the profiling split of the 5M-pair SAM at 6 and 1 threads.
set -uo pipefail
W=$HOME/mt-work/perf4; B=$W/ref/build/protal; DB=$HOME/bench071/V073/protal_db
P=$HOME/bench071/samples/points
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
run() { local n=$1 t=$2; shift 2; rm -rf $W/o.$n
  /usr/bin/time -v $B --db $DB "$@" --prefix s -o $W/o.$n -t $t --no_qcmsa --verbose > $W/o.$n.log 2> $W/o.$n.time
  echo "== $n t=$t: $(grep -E 'Elapsed|User time' $W/o.$n.time | tr -s ' ' | tr '\n' ';')"
  grep -E "^(Output handler|Alignment handler|Aligning reads|Thread 0 Profile sample|Profiling|Run protal) took" $W/o.$n.log; }
run pe500k_t6_plain 6 -1 $R1 -2 $R2 --read_type pe --sam_format sam --no_profile
run pe500k_t6_np 6 -1 $R1 -2 $R2 --read_type pe --no_profile
run pe500k_t3_np 3 -1 $R1 -2 $R2 --read_type pe --no_profile
run pe500k_t1_np 1 -1 $R1 -2 $R2 --read_type pe --no_profile
run prof5M_t6 6 --profile_only $W/o.pe5M_t6/s.sam.zst
run prof5M_t1 1 --profile_only $W/o.pe5M_t6/s.sam.zst
run prof500k_t6 6 --profile_only $W/o.pe500k_t1/s.sam.zst
echo "SAM records pe500k: $(zstd -dc $W/o.pe500k_t1/s.sam.zst | grep -vc '^@')  pe5M: $(zstd -dc $W/o.pe5M_t6/s.sam.zst | grep -vc '^@')"
