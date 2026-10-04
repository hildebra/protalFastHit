#!/usr/bin/env bash
set -uo pipefail
W=$HOME/perf-gtdb; H=$W/head/build/protal; N=$W/work/build/protal; X=$W/e2e
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
run() { local bin=$1 o=$2; shift 2; rm -rf $o; $bin --db $X/db/database.protal -1 $R1 -2 $R2 --read_type pe --prefix s -o $o -t 6 --no_qcmsa "$@" > $o.log 2>&1; }
cmpout() { diff -rq -x '*_runtime.tsv' -x '*.sam.zst' -x '*.statistics.tsv' $1 $2 > /dev/null && echo same || echo DIFF; }
for i in 1 2 3 4; do run $N $X/rs$i --sequential_load; run $N $X/rc$i; done
for i in 2 3 4; do echo "sequential 1 vs $i: $(cmpout $X/rs1 $X/rs$i); concurrent 1 vs $i: $(cmpout $X/rc1 $X/rc$i); sequential 1 vs concurrent $i: $(cmpout $X/rs1 $X/rc$i)"; done
echo "== profile-only on one SAM, 4 times at 6 threads (the profiling stage alone)"
for i in 1 2 3 4; do rm -rf $X/po$i; $N --db $X/db/database.protal --profile_only $X/rs1/s.sam.zst --read_type pe --prefix s -o $X/po$i -t 6 --no_qcmsa --sequential_load > $X/po$i.log 2>&1; done
for i in 2 3 4; do echo "profile-only 1 vs $i: $(cmpout $X/po1 $X/po$i)"; done
