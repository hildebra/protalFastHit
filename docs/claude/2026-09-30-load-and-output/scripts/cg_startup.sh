#!/bin/bash
# callgrind of a 1000-pair run (1 thread, with profiling; startup.sh makes the reads): where the
# start-up and shut-down instructions go.  cg_startup.sh DB [DB ...]
source "$(dirname "$0")/env.sh"
for db in "$@"; do
  name=$(basename $db)
  rm -rf $OUT/cgo
  valgrind --tool=callgrind --callgrind-out-file=$OUT/cg_start_$name.out $BIN --db $db -1 $OUT/tiny_R1.fq.gz \
    -2 $OUT/tiny_R2.fq.gz -o $OUT/cgo -t 1 --no_qcmsa > $OUT/cg_start_$name.log 2>&1
  callgrind_annotate --inclusive=yes --threshold=99.9 $OUT/cg_start_$name.out 2>/dev/null > $OUT/cg_start_$name.txt
  echo "== $name"; grep "PROGRAM TOTALS" $OUT/cg_start_$name.txt
  grep -E "Seedmap::Load |DecodeChunk|ZSTD_decompressDCtx|LoadModel|LoadAllGenomes|ParallelReadFrames|ProtalDB::ProtalDB|LoadUniqueKmers|LoadPositionMap|WriteSamHeader|IntTaxonomy::Load|CheckReference|ProfileWrapper|RunPairedEnd<" \
    $OUT/cg_start_$name.txt | cg_lines | awk '!seen[$2]++'
done
