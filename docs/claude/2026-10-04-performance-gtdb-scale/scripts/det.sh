#!/usr/bin/env bash
set -uo pipefail
W=$HOME/perf-gtdb; H=$W/head/build/protal; N=$W/work/build/protal; X=$W/e2e; DB=$X/db/database.protal
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
run() { local bin=$1 o=$2; shift 2; rm -rf $o; $bin --db $DB -1 $PB --read_type pb --prefix s -o $o -t 6 --no_qcmsa "$@" > $o.log 2>&1; }
same() { a=$(cmp -s <(zstd -dc $1/s.sam.zst | sort) <(zstd -dc $2/s.sam.zst | sort) && echo SAM-same || echo SAM-DIFF)
  b=$(diff -rq -x '*_runtime.tsv' -x '*.sam.zst' -x '*.statistics.tsv' $1 $2 > /dev/null && echo out-same || echo OUT-DIFF); echo "$3: $a $b"; }
run $N $X/ps1 --sequential_load; run $N $X/ps2 --sequential_load
run $N $X/pc1; run $N $X/pc2; run $N $X/pc3
same $X/hp $X/ps1 "head vs sequential1"; same $X/ps1 $X/ps2 "sequential1 vs sequential2"
same $X/pc1 $X/pc2 "concurrent1 vs concurrent2"; same $X/pc1 $X/pc3 "concurrent1 vs concurrent3"; same $X/ps1 $X/pc1 "sequential1 vs concurrent1"
