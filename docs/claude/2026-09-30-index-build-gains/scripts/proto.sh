#!/bin/bash
# The prototype (proto.py on perf) against perf on the 765-species world, alternated: phase times,
# the index byte for byte and unique_kmers.tsv as sorted rows.
set -e
S=$(cd "$(dirname "$0")" && pwd)  # this folder; head56.tar, syncmer.patch, huge.patch: see README
B=$HOME/protal-scale; W=/tmp/scale
if [ ! -x $B/proto/build/protal ]; then
  rm -rf $B/proto; mkdir -p $B/proto/src; tar -xf $S/head56.tar -C $B/proto/src
  python3 $S/instrument.py $B/proto/src
  (cd $B/proto/src && patch -p1 -s < $S/syncmer.patch && patch -p1 -s < $S/huge.patch)
  python3 $S/proto.py $B/proto/src
  cmake -S $B/proto/src -B $B/proto/build -DCMAKE_BUILD_TYPE=Release > $B/proto/configure.log 2>&1
  cmake --build $B/proto/build -j 8 --target protal > $B/proto/make.log 2>&1 || { tail -30 $B/proto/make.log; exit 1; }
fi
for rep in 1 2; do
  for t in 1 8; do
    for v in perf proto; do
      D=$W/db_$v; rm -rf $D; cp -r $W/in100 $D
      $B/$v/build/protal --build --no_profile --no_compress -t $t --db $D --reference $D/reference.fna \
          --full_reference $D/full_reference.fna 2>&1 | python3 $S/stamp.py > $W/plog_${v}_${t}_$rep
      python3 $S/phases.py "$v p100 t$t" $W/plog_${v}_${t}_$rep
      echo "    index $(md5sum < $D/index.prx | cut -c1-8), unique_kmers sorted $(sort $D/unique_kmers.tsv | md5sum | cut -c1-8), as written $(md5sum < $D/unique_kmers.tsv | cut -c1-8)"
    done
    rm -rf $W/db_perf $W/db_proto
  done
done
# corrected counters (the prototype counts only lookups that found entries) at each size, 1 thread
for p in 25 50; do
  D=$W/db_proto; rm -rf $D; cp -r $W/in$p $D
  $B/proto/build/protal --build --no_profile --no_compress -t 1 --db $D --reference $D/reference.fna \
      --full_reference $D/full_reference.fna 2>&1 | python3 $S/stamp.py > $W/plog_proto_p${p}
  python3 $S/phases.py "proto p$p t1" $W/plog_proto_p${p}
  rm -rf $D
done
