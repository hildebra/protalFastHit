#!/usr/bin/env bash
set -uo pipefail
W=$HOME/perf-gtdb; H=$W/head/build/protal; N=$W/work/build/protal; X=$W/e2e; ORIG=$HOME/bench071/V073/protal_db/database.protal
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
run() { local bin=$1 db=$2 o=$3; shift 3; rm -rf $o; $bin --db $db --prefix s -o $o -t 6 --no_qcmsa "$@" > $o.log 2>&1; }
same() { a=$(cmp -s <(zstd -dc $1/s.sam.zst | sort) <(zstd -dc $2/s.sam.zst | sort) && echo SAM-same || echo SAM-DIFF)
  b=$(diff -rq -x '*_runtime.tsv' -x '*.sam.zst' -x '*.statistics.tsv' $1 $2 > /dev/null && echo out-same || echo OUT-DIFF); echo "$3: $a $b"; }
run $N $ORIG $X/po -1 $PB --read_type pb --sequential_load
run $N $X/folder $X/pf -1 $PB --read_type pb --sequential_load
run $H $ORIG $X/hp2 -1 $PB --read_type pb
same $X/hp $X/hp2 "head vs head again"
same $X/hp $X/po "head vs work on the original file (text tables)"
same $X/hp $X/pf "head vs work on the folder (text tables)"
same $X/po $X/ps1 "work original (text) vs work with gene_table.bin"
echo "== pe determinism"
run $N $X/db/database.protal $X/pec1 -1 $R1 -2 $R2 --read_type pe
run $N $X/db/database.protal $X/pec2 -1 $R1 -2 $R2 --read_type pe
run $N $X/db/database.protal $X/pes1 -1 $R1 -2 $R2 --read_type pe --sequential_load
run $H $ORIG $X/peh2 -1 $R1 -2 $R2 --read_type pe
same $X/pec1 $X/pec2 "pe concurrent1 vs concurrent2"; same $X/pes1 $X/pec1 "pe sequential vs concurrent"; same $X/h $X/peh2 "pe head vs head again"; same $X/h $X/pes1 "pe head vs sequential"
