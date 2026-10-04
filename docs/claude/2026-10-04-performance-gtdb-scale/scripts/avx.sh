#!/usr/bin/env bash
# The working tree + the WFA2 AVX2 extend-kernel patch (1fced40) built in ~/perf-gtdb/avx; PacBio and short-read runs against ~/perf-gtdb/work.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
W=$HOME/perf-gtdb; E=$W/avx
rm -rf $E; mkdir -p $E
rsync -a --exclude '/build/' $W/work/ $E/
cd $E && patch -p1 < $S/wfa_avx2.patch > $W/avx.patch.log 2>&1 || { echo "PATCH FAILED"; cat $W/avx.patch.log; exit 1; }
nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF > $W/avx.configure.log 2>&1 && \
  nice cmake --build build --target protal -j 6 > $W/avx.build.log 2>&1 && echo "OK build" || { echo FAIL; grep -E 'error' -A3 $W/avx.build.log | head -30; exit 1; }
grep 'AVX2 chosen' $W/avx.configure.log
DB=$HOME/bench071/V073/protal_db
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
run() { local name=$1 bin=$2 t=$3 o=$4; shift 4; rm -rf $o
  /usr/bin/time -f "%e s wall, %U s user" $bin --db $DB "$@" --prefix s -o $o -t $t --no_qcmsa --no_profile > $o.log 2> $o.time
  echo "$name: $(tail -1 $o.time); $(grep -E '^Aligning reads took' $o.log)"; }
for round in 1 2 3; do
  run "pb work $round" $W/work/build/protal 6 $W/x.work.$round -1 $PB --read_type pb
  run "pb avx  $round" $E/build/protal 6 $W/x.avx.$round -1 $PB --read_type pb
done
for round in 1 2; do
  run "pe work $round" $W/work/build/protal 1 $W/y.work.$round -1 $R1 -2 $R2 --read_type pe
  run "pe avx  $round" $E/build/protal 1 $W/y.avx.$round -1 $R1 -2 $R2 --read_type pe
done
cmp -s <(zstd -dc $W/x.work.1/s.sam.zst | sort) <(zstd -dc $W/x.avx.1/s.sam.zst | sort) && echo "pb SAM identical as sets" || echo "pb SAM DIFFERS"
cmp -s <(zstd -dc $W/y.work.1/s.sam.zst) <(zstd -dc $W/y.avx.1/s.sam.zst) && echo "pe SAM identical" || echo "pe SAM DIFFERS"
