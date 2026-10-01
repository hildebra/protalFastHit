#!/bin/bash
# Do two builds write the same alignments? same_output.sh BIN_A BIN_B DB READS_DIR PAIRS
# Runs both on the first PAIRS pairs (one thread, alignment only) and compares the decompressed SAMs.
source "$(dirname "$0")/env.sh"
A=$1; B=$2; db=$3; reads=$4; pairs=$5
d=$PERF_DIR/runs2/same; rm -rf $d; mkdir -p $d
head -n $((pairs*4)) $(ls $reads/*_R1.fq* | head -1) > $d/r1.fq; head -n $((pairs*4)) $(ls $reads/*_R2.fq* | head -1) > $d/r2.fq
for x in a b; do
  bin=$A; [ $x = b ] && bin=$B
  $bin --db $db -1 $d/r1.fq -2 $d/r2.fq -o $d/$x -t 1 --no_qcmsa --no_profile > $d/$x.log 2>&1
done
f() { ls $d/$1/*.sam* 2>/dev/null | head -1; }
fa=$(f a); fb=$(f b)
dec() { case $1 in *.zst) zstd -dc $1 ;; *.gz) zcat $1 ;; *) cat $1 ;; esac; }
dec $fa | grep -v '^@' | sort > $d/a.sam; dec $fb | grep -v '^@' | sort > $d/b.sam
echo "records: $(wc -l < $d/a.sam) / $(wc -l < $d/b.sam); differing lines: $(diff $d/a.sam $d/b.sam | grep -c '^[<>]')"
