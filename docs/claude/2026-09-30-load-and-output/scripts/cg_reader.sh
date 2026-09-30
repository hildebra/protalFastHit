#!/bin/bash
# callgrind, 50k pairs, plain and gzipped, 1 thread: the FASTQ reader under its lock (LoadBatch),
# outside it (NextSequence) and the inflating thread.  cg_reader.sh DB PLAIN_READS_DIR
# PLAIN_READS_DIR holds *_R1.fq and *_R2.fq.
source "$(dirname "$0")/env.sh"
db=$1; src=$2
head -n 200000 $(ls $src/*_R1.fq | head -1) > $OUT/m50_R1.fq; head -n 200000 $(ls $src/*_R2.fq | head -1) > $OUT/m50_R2.fq
gzip -c $OUT/m50_R1.fq > $OUT/m50_R1.fq.gz; gzip -c $OUT/m50_R2.fq > $OUT/m50_R2.fq.gz
for v in plain gz; do
  s=""; [ $v = gz ] && s=.gz
  rm -rf $OUT/cgr
  valgrind --tool=callgrind --callgrind-out-file=$OUT/cg_reader_$v.out $BIN --db $db -1 $OUT/m50_R1.fq$s -2 $OUT/m50_R2.fq$s \
    -o $OUT/cgr -t 1 --no_qcmsa --no_profile > /dev/null 2>&1
  callgrind_annotate --inclusive=yes --threshold=99.9 $OUT/cg_reader_$v.out 2>/dev/null > $OUT/cg_reader_$v.txt
  echo "== $v"; grep "PROGRAM TOTALS" $OUT/cg_reader_$v.txt
  grep -E "RunPairedEnd<|LoadBatch|ReadNextSequence|ThreadedGzStreambuf::Inflate|crc32_z|SimpleKmerHandler.*operator|SimpleAlignmentHandler::operator|OutputHandler.*operator" \
    $OUT/cg_reader_$v.txt | cg_lines | awk '!seen[$2]++'
done
