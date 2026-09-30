#!/bin/bash
# zstd level against time and size for the index and the reference of the 765-species tuning
# world: protal --compress_db --no_bundle on a copy of an uncompressed build, 8 threads.
set -e
S=$(cd "$(dirname "$0")" && pwd)  # this folder; head56.tar, syncmer.patch, huge.patch: see README
P=$HOME/protal-scale/base/build/protal
W=/tmp/scale
if [ ! -f $W/raw/index.prx ]; then
  rm -rf $W/raw; cp -r $W/in100 $W/raw
  $P --build --no_profile --no_compress -t 8 --db $W/raw --reference $W/raw/reference.fna \
     --full_reference $W/raw/full_reference.fna > $W/raw.log 2>&1
  rm -f $W/raw/full_reference.fna
fi
echo "raw: index.prx $(du -m $W/raw/index.prx | cut -f1) MB, reference.fna $(du -m $W/raw/reference.fna | cut -f1) MB"
for L in 3 9 12 15 17 19; do
  D=$W/lvl; rm -rf $D; cp -r $W/raw $D
  $P --compress_db --no_bundle --compress_level $L -t 8 --db $D 2>&1 | python3 $S/stamp.py > $W/lvl$L.log
  ti=$(grep -m1 "Compress index.prx took" $W/lvl$L.log | awk '{print $1}')
  t0=$(grep -m1 "Write .*index.prx.zst from" $W/lvl$L.log | awk '{print $1}')
  tr=$(grep -m1 "Compress reference took\|Reference .*reference.fna.zst:" $W/lvl$L.log | awk '{print $1}')
  echo "level $L: index $(du -k $D/index.prx.zst | cut -f1) KB, reference $(du -k $D/reference.fna.zst | cut -f1) KB;" \
       "log marks: index start $t0 end $ti, reference end $tr, total $(tail -n 1 $W/lvl$L.log | awk '{print $1}') s"
done
rm -rf $W/lvl
