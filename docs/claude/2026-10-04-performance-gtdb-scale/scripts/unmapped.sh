#!/usr/bin/env bash
# Where the profiling stage's "reading the SAM" goes with 45M unmapped records (as the r226 run's pe1 SAM): the deep
# sample's SAM (4.66M records, 393k unmapped) plain and .sam.zst, and the same with 45M synthetic unmapped records
# appended (ZF tags of 1-4 taxids), plain and zstd -T6; --profile_only at 6 and 1 threads.
set -uo pipefail
W=$HOME/perf-gtdb; N=$W/work/build/protal; DB=$HOME/bench071/V073/protal_db; X=$W/unmapped; mkdir -p $X
SAM=$W/d.work.1/s.sam.zst
[ -e $X/deep.sam ] || zstd -dc $SAM > $X/deep.sam
if [ ! -e $X/big.sam ]; then
  cp $X/deep.sam $X/big.sam
  awk 'BEGIN { srand(7); for (i = 0; i < 45000000; i++) { n = 1 + int(rand() * 4); zf = 1 + int(rand() * 760)
       for (j = 1; j < n; j++) zf = zf "," (1 + int(rand() * 760))
       printf "A01234:567:HXXXXXXXX:1:%d:%d:%d\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:%s\n", 1101 + int(i / 40000), 1000 + int(rand() * 30000), 1000 + int(rand() * 40000), zf } }' >> $X/big.sam
  ls -la $X/big.sam
fi
[ -e $X/big.sam.zst ] || zstd -q -T6 -3 $X/big.sam -o $X/big.sam.zst
ls -la $X/*.sam* | awk '{ print $5, $9 }'
prof() { local name=$1 t=$2 f=$3; rm -rf $X/o; /usr/bin/time -f "%e s wall, %U s user" $N --db $DB --profile_only $f --read_type pe --prefix s -o $X/o -t $t --no_qcmsa > $X/o.log 2> $X/o.time
  echo "$name (t=$t): $(tail -1 $X/o.time); $(grep -E '^Profiling sample' $X/o.log | cut -c1-200)"; }
prof "deep .sam.zst (protal's)" 6 $SAM
prof "deep plain .sam" 6 $X/deep.sam
prof "big plain .sam, 45M unmapped" 6 $X/big.sam
prof "big .sam.zst (zstd -T6, one stream)" 6 $X/big.sam.zst
prof "big plain .sam" 1 $X/big.sam
