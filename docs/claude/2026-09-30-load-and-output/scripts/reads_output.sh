#!/bin/bash
# Reads in and SAM out, 8 threads: the alignment stage with gzipped and plain reads (alternated),
# inflating alone, pigz on the SAM, and profiling that SAM plain and gzipped.
#   reads_output.sh DB NAME GZ_DIR PLAIN_DIR [NAME GZ_DIR PLAIN_DIR ...]
# Each GZ_DIR holds *_R1.fq.gz/*_R2.fq.gz, each PLAIN_DIR the same reads as *_R1.fq/*_R2.fq. The SAM
# of the last NAME is compressed and profiled.
source "$(dirname "$0")/env.sh"
db=$1; shift
now() { date +%s.%N; }
el() { echo "$(now) - $1" | bc; }
sets=()
while [ $# -ge 3 ]; do sets+=("$1:$2:$3"); shift 3; done
first=${sets[0]}; gzdir=$(echo $first | cut -d: -f2); plaindir=$(echo $first | cut -d: -f3)
echo "== inflating R1 of $(echo $first | cut -d: -f1) alone"
for rep in 1 2; do
  t=$(now); zcat $(ls $gzdir/*_R1.fq.gz | head -1) > /dev/null; a=$(el $t)
  t=$(now); pigz -dc -p 4 $(ls $gzdir/*_R1.fq.gz | head -1) > /dev/null; b=$(el $t)
  echo "rep $rep: zcat $a s, pigz -dc $b s"
done
echo "== alignment stage, 8 threads, --no_profile"
for rep in 1 2; do
  for s in "${sets[@]}"; do
    name=$(echo $s | cut -d: -f1)
    for kind in gz plain; do
      dir=$(echo $s | cut -d: -f2); ext=.fq.gz; [ $kind = plain ] && { dir=$(echo $s | cut -d: -f3); ext=.fq; }
      v=${name}_$kind; rm -rf $OUT/o_$v
      /usr/bin/time -f "%e s wall, %U s user, %S s sys" $BIN --db $db -1 $(ls $dir/*_R1$ext | head -1) -2 $(ls $dir/*_R2$ext | head -1) \
        -o $OUT/o_$v -t 8 --no_qcmsa --no_profile --prefix $v > $OUT/$v.log 2> $OUT/$v.time
      rt=$OUT/o_$v/misc/${v}_runtime.tsv
      echo -e "rep $rep $v\t$(grep 'Aligning reads took' $OUT/$v.log | sed 's/Aligning reads took //')\treader/thread $(awk -F'\t' '$1 == "Sequence reader" { print $4 }' $rt)\toutput/thread $(awk -F'\t' '$1 == "Output handler" { print $4 }' $rt)\t$(cat $OUT/$v.time)\tload $(cut -d' ' -f1 /proc/loadavg)"
    done
  done
done
last=${sets[${#sets[@]}-1]}; name=$(echo $last | cut -d: -f1)
S=$OUT/o_${name}_gz/${name}_gz.sam
echo "== the SAM of ${name}_gz: $(stat -c %s $S) bytes, $(grep -c '^@' $S) header lines, $(grep -vc '^@' $S) records"
for p in 1 8; do t=$(now); pigz -p $p -c $S > $OUT/${name}_gz.sam.gz; echo "pigz -p $p: $(el $t) s, $(stat -c %s $OUT/${name}_gz.sam.gz) bytes"; done
echo "== profiling it (--profile_only, 1 thread; one sample is profiled by one thread)"
for rep in 1 2; do
  for f in $S $OUT/${name}_gz.sam.gz; do
    rm -rf $OUT/po; mkdir -p $OUT/po
    /usr/bin/time -f "%e s wall, %U s user" $BIN --db $db --profile_only $f --prefix $name -o $OUT/po -t 1 --no_qcmsa > $OUT/po.log 2> $OUT/po.time
    echo "rep $rep $(basename $f): $(grep -E 'Profiling took' $OUT/po.log) | $(tail -1 $OUT/po.time)"
  done
done
