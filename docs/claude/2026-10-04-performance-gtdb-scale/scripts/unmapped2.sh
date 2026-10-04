#!/usr/bin/env bash
# As unmapped.sh, but the 45M unmapped records' ZF tags name 1-8 taxids drawn from 1..140000 (as on GTDB r226, where the
# failed taxa of a chunk are many thousands of distinct species): does the per-chunk evidence merge become the cost?
set -uo pipefail
W=$HOME/perf-gtdb; N=$W/work/build/protal; DB=$HOME/bench071/V073/protal_db; X=$W/unmapped
if [ ! -e $X/big2.sam ]; then
  cp $X/deep.sam $X/big2.sam
  awk 'BEGIN { srand(11); for (i = 0; i < 45000000; i++) { n = 1 + int(rand() * 8); zf = 1 + int(rand() * 140000)
       for (j = 1; j < n; j++) zf = zf "," (1 + int(rand() * 140000))
       printf "A01234:567:HXXXXXXXX:1:%d:%d:%d\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:%s\n", 1101 + int(i / 40000), 1000 + int(rand() * 30000), 1000 + int(rand() * 40000), zf } }' >> $X/big2.sam
fi
ls -la $X/big2.sam | awk '{ print $5, $9 }'
prof() { local name=$1 t=$2 f=$3; rm -rf $X/o; /usr/bin/time -f "%e s wall, %U s user" $N --db $DB --profile_only $f --read_type pe --prefix s -o $X/o -t $t --no_qcmsa > $X/o.log 2> $X/o.time
  echo "$name (t=$t): $(tail -1 $X/o.time); $(grep -E '^Profiling sample' $X/o.log | cut -c1-120)"; }
prof "big2 plain .sam, 45M unmapped, ZF over 140k taxids" 6 $X/big2.sam
prof "big plain .sam, 45M unmapped, ZF over 760 taxids" 6 $X/big.sam
prof "big2 plain .sam" 1 $X/big2.sam
