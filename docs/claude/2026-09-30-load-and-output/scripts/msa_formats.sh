#!/bin/bash
# Strain MSAs compressed as gzip, BGZF (libdeflate), zstd -3 and -19: a read set split into samples
# (the same species in each), protal with qcmsa, then every raw and filtered MSA compressed.
#   msa_formats.sh DB R1.fq.gz R2.fq.gz SAMPLES PAIRS_PER_SAMPLE
# Needs the to_bgzf helper (built from to_bgzf.cpp) in $OUT, and $PROTAL_SRC/scripts/qcmsa.py.
source "$(dirname "$0")/env.sh"
db=$1; w1=$2; w2=$3; samples=$4; pairs=$5
O=$OUT/msa_formats; rm -rf $O; mkdir -p $O/reads
for m in 1 2; do
  f=$w1; [ $m = 2 ] && f=$w2
  zcat $f | awk -v o=$O/reads -v m=$m -v lines=$((4 * pairs)) -v n=$samples \
    '{ s = int((NR - 1) / lines) + 1; if (s <= n) print > (o "/s" s "_R" m ".fq") }'
done
printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\n' $O/out > $O/samples.map
for s in $(seq 1 $samples); do printf 's%d\ts%d\t%s\t%s\n' $s $s $O/reads/s${s}_R1.fq $O/reads/s${s}_R2.fq >> $O/samples.map; done
PROTAL_QCMSA_SCRIPT=$PROTAL_SRC/scripts/qcmsa.py $BIN --db $db --map $O/samples.map -t 8 > $O/run.log 2>&1
echo "protal rc=$?; $(grep -oE 'Strain-level MSAs took [0-9a-z ]+' $O/run.log)"
S=$O/out/strains
echo "raw MSAs: $(ls $S/*.raw.msa.fna 2>/dev/null | wc -l), rows of the largest: $(grep -c '>' $(ls -S $S/*.raw.msa.fna | head -1))"
cat $S/*.raw.msa.fna > $O/raw.fna
ls $S/*.msa.fna | grep -v '\.raw\.' | xargs -r cat > $O/filtered.fna
now() { date +%s.%N; }
el() { echo "$(now) - $1" | bc; }
for set in raw filtered; do
  f=$O/$set.fna; size=$(stat -c %s $f)
  [ $size -gt 0 ] || continue
  echo "== $set MSAs: $size bytes"
  t=$(now); gzip -6 -c $f > $f.gz; c=$(el $t)
  echo "  gzip -6            ratio $(echo "scale=1; $size / $(stat -c %s $f.gz)" | bc)  compress ${c}s"
  t=$(now); $OUT/to_bgzf $f $f.bgzf.gz 1; c=$(el $t)
  echo "  BGZF libdeflate 6  ratio $(echo "scale=1; $size / $(stat -c %s $f.bgzf.gz)" | bc)  compress ${c}s"
  for lvl in 3 19; do
    t=$(now); zstd -q -f -$lvl -T1 $f -o $f.$lvl.zst; c=$(el $t)
    echo "  zstd -$lvl           ratio $(echo "scale=1; $size / $(stat -c %s $f.$lvl.zst)" | bc)  compress ${c}s"
  done
done
big=$(ls -S $S/*.raw.msa.fna | head -1); name=$(basename $big .raw.msa.fna)
gzip -6 -c $big > $O/big.raw.msa.fna.gz; cp $big $O/big.raw.msa.fna
for v in big.raw.msa.fna big.raw.msa.fna.gz; do
  t=$(now); python3 $PROTAL_SRC/scripts/qcmsa.py $O/$v $S/$name.raw.partition.txt $(ls $S/$name*meta* | head -1) --prefix $O/q > /dev/null 2>&1
  echo "qcmsa on $v: rc=$? $(el $t) s"
done
