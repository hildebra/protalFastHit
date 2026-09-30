#!/bin/bash
# callgrind, 1 thread: reading pairs as single-member gzip, BGZF and plain; writing a SAM as .sam.gz,
# .sam.zst and .sam.  gzip_cg.sh DB PLAIN_R1.fq PLAIN_R2.fq R1.fq.gz R2.fq.gz
# (the plain pair for reading, e.g. 50k pairs of mix; the gzipped pair for writing, e.g. 20k of w900)
source "$(dirname "$0")/env.sh"
db=$1; p1=$2; p2=$3; w1=$4; w2=$5
O=$OUT/gzip_cg; mkdir -p $O
g++ -std=c++20 -O2 -I$PROTAL_SRC/src/IO $here/to_bgzf.cpp -o $OUT/to_bgzf -ldeflate -pthread || exit 1
gzip -c $p1 > $O/r1.fq.gz; gzip -c $p2 > $O/r2.fq.gz
$OUT/to_bgzf $p1 $O/r1.bgzf.fq.gz 4; $OUT/to_bgzf $p2 $O/r2.bgzf.fq.gz 4
cg() {
  local name=$1; shift
  rm -rf $O/o; valgrind --tool=callgrind --callgrind-out-file=$O/cg_$name.out $BIN --db $db "$@" -o $O/o -t 1 --no_qcmsa --no_profile > /dev/null 2>&1
  callgrind_annotate --inclusive=yes --threshold=99.9 $O/cg_$name.out 2>/dev/null > $O/cg_$name.txt
  echo -e "$name\ttotal $(grep 'PROGRAM TOTALS' $O/cg_$name.txt | awk '{print $1}')\tinflate $(grep -E 'ThreadedGzStreambuf::Inflate' $O/cg_$name.txt | grep -v "'2" | head -1 | awk '{print $1}')\tSamOutput::Write $(grep -E 'SamOutput::Write' $O/cg_$name.txt | grep -v "'2" | head -1 | awk '{print $1}')"
}
cg in_gzip -1 $O/r1.fq.gz -2 $O/r2.fq.gz --sam_format sam
cg in_bgzf -1 $O/r1.bgzf.fq.gz -2 $O/r2.bgzf.fq.gz --sam_format sam
cg in_plain -1 $p1 -2 $p2 --sam_format sam
for f in gz zst sam; do cg out_$f -1 $w1 -2 $w2 --sam_format $f; done
