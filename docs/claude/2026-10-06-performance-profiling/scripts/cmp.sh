#!/usr/bin/env bash
# cmp.sh <name>: callgrind of ~/perf6/<name>/build/protal against HEAD's on aligning 100k pairs and PacBio 3 Mb
# (no profiling), and the SAM records compared (decompressed, the @PG/@CO header lines left out).
set -uo pipefail
N=$1; W=$HOME/perf6; DB=$HOME/bench071/V075/protal_db; P=$HOME/bench071/samples/points
mkdir -p $W/cmp_$N
one() { # binary tag set args...
  local b=$1 tag=$2 set=$3; shift 3
  local o=$W/cmp_$N/$set.$tag
  rm -rf $o
  valgrind --tool=callgrind --callgrind-out-file=$o.cg $b --db $DB "$@" --prefix s -o $o -t 1 --no_qcmsa --no_profile > $o.log 2>&1
  callgrind_annotate --inclusive=yes $o.cg 2>/dev/null > $o.incl.txt
  echo "$set $tag: $(grep -m1 'PROGRAM TOTALS' $o.incl.txt)"
}
for set in pe pb; do
  if [ $set = pe ]; then args="-1 $W/reads/pe100k_R1.fq.gz -2 $W/reads/pe100k_R2.fq.gz --read_type pe"; else args="-1 $P/pb_b3000000/sim/reads/pb_b3000000_s_1.fq.gz --read_type pb"; fi
  [ -s $W/cmp_head_$set.incl.txt ] || { one $W/head/build/protal head $set $args; cp $W/cmp_$N/$set.head.incl.txt $W/cmp_head_$set.incl.txt; cp -r $W/cmp_$N/$set.head $W/cmp_head_$set.out; }
  one $W/$N/build/protal $N $set $args
  a=$W/cmp_head_$set.out/s.sam.zst; b=$W/cmp_$N/$set.$N/s.sam.zst
  if cmp -s <(zstd -dc $a | grep -v '^@PG\|^@CO') <(zstd -dc $b | grep -v '^@PG\|^@CO'); then echo "$set: SAM records identical ($(zstd -dc $a | grep -vc '^@'))"; else echo "$set: SAM DIFFERS"; fi
  for f in RunPairedEnd RunLongReads MayAlign AlignAnchor ChainAnchorFinder\<protal::KmerLookupSM\>::operator; do
    h=$(grep -m1 "$f" $W/cmp_head_$set.incl.txt | awk '{print $1}'); n=$(grep -m1 "$f" $W/cmp_$N/$set.$N.incl.txt | awk '{print $1}')
    [ -n "$h$n" ] && echo "  $f: head $h, $N $n"
  done
done
