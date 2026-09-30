#!/bin/bash
# Unpacks the mini DB into ~/audit6/io/db and simulates reads.
set -eu
S=$(dirname "$0")
P=~/strain-build/bin/protal
W=~/audit6/io
mkdir -p $W/db $W/reads
if [ ! -e $W/db/reference.map ]; then
  $P --unpack_db --db ~/strain-build/mini_db/protal_db/database.protal --unpack_dir $W/db -t 2 > $W/unpack.log 2>&1
fi
ls -la $W/db
cp ~/strain-build/mini_db/protal_db/database.protal $W/bundle.protal
cd $W/db
REF=reference.fna
[ -e reference.fna ] || { zstd -dcq reference.fna.zst > $W/reference.plain.fna; REF=$W/reference.plain.fna; }
[ -e $W/reference.plain.fna ] || cp reference.fna $W/reference.plain.fna
python3 $S/mkreads.py reads $W/reference.plain.fna $W/reads/sa 12 1
python3 $S/mkreads.py reads $W/reference.plain.fna $W/reads/sb 6 2
wc -l $W/reads/*.fq
head -4 $W/reads/sa_R1.fq
grep -c '>' $W/reference.plain.fna
wc -l $W/db/reference.map
head -3 $W/db/reference.map
ls $W/db
