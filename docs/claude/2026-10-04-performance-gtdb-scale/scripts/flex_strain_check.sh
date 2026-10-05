#!/usr/bin/env bash
# The AVX2 flex scan and the parallel strain stage: build + all tests + the scan's bench; then runs on the v0.7.3 world
# (pe 500k, PacBio 90 Mb, 6 threads) with the AVX2 scan and with PROTAL_FLEX_SCAN=scalar (identical outputs? seeding
# time, 3 alternated rounds); then a cohort of both from their SAMs at 6 threads and at 1 (identical strain outputs? time).
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh '*' 8 || exit 1
W=$HOME/perf-gtdb; N=$W/work/build/protal; DB=$W/e2e/db/database.protal; X=$W/avx; rm -rf $X; mkdir -p $X
echo "== bench"; PROTAL_FLEX_BENCH=1 $W/work/build/tests/protal_tests --gtest_filter='FlexScan.*' 2>&1 | grep -E 'blocks|PASSED|FAILED'
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
same() { a=$(cmp -s <(zstd -dc $1/$3.sam.zst | sort) <(zstd -dc $2/$3.sam.zst | sort) && echo SAM-same || echo SAM-DIFF)
  b=$(diff -rq -x '*_runtime.tsv' -x '*.sam.zst' -x '*.statistics.tsv' $1 $2 > /dev/null && echo out-same || echo OUT-DIFF); echo "$4: $a $b"; }
for round in 1 2 3; do
  for v in avx2 scalar; do
    for s in pe pb; do
      o=$X/$s.$v.$round; rm -rf $o
      if [ $s = pe ]; then reads=(-1 $R1 -2 $R2); else reads=(-1 $PB); fi
      if [ $v = scalar ]; then PROTAL_FLEX_SCAN=scalar $N --db $DB "${reads[@]}" --read_type $s --prefix $s -o $o -t 6 --no_qcmsa > $o.log 2>&1
      else $N --db $DB "${reads[@]}" --read_type $s --prefix $s -o $o -t 6 --no_qcmsa > $o.log 2>&1; fi
      seeding=$(awk -F'\t' '$1 == "Seeding" { printf "%.2f", $4 }' $o/misc/${s}_runtime.tsv)
      echo "$s $v $round: seeding $seeding s per thread; $(grep -hE '^Aligning reads took|^Run protal took' $o.log | tr '\n' ' ')"
    done
  done
done
same $X/pe.avx2.1 $X/pe.scalar.1 pe "pe AVX2 vs scalar"; same $X/pb.avx2.1 $X/pb.scalar.1 pb "pb AVX2 vs scalar"
same $X/pe.avx2.1 $X/pe.avx2.2 pe "pe AVX2 run 1 vs 2"
grep -h 'seeding:' $X/pe.avx2.1.log $X/pe.scalar.1.log | cut -c1-120
echo "== cohort (strain stage) at 6 and 1 threads"
mkdir -p $X/sams; cp $X/pe.avx2.1/pe.sam.zst $X/pb.avx2.1/pb.sam.zst $X/sams/
for t in 6 1; do
  o=$X/cohort.t$t; rm -rf $o; map=$X/cohort.t$t.map
  { printf '#OUTPUT_DIR\t%s\n' "$o"; printf '#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tREAD_TYPE\n'
    printf 'pe\t%s\t%s\t%s\tpe\tpe\n' "$R1" "$R2" "$X/sams/pe.sam.zst"; printf 'pb\t%s\t-\t%s\tpb\tpb\n' "$PB" "$X/sams/pb.sam.zst"; } > $map
  /usr/bin/time -f "%e s wall" $N --db $DB --map $map -t $t > $o.log 2> $o.time
  echo "threads $t: $(tail -1 $o.time); $(grep -hE '^Strain-level MSAs took|^Strain-level MSAs of|^Building the strain MSAs took|^qcMSA took|All alignments' $o.log | tr '\n' ' ')"
done
diff -rq -x '*_runtime.tsv' -x '*.statistics.tsv' $X/cohort.t6/strains $X/cohort.t1/strains > $X/strains.diff && echo "strain outputs identical at 6 and 1 threads ($(ls $X/cohort.t6/strains | wc -l) files)" || { echo "STRAIN OUTPUTS DIFFER"; head $X/strains.diff; }
diff -rq -x '*_runtime.tsv' -x '*.statistics.tsv' -x strains $X/cohort.t6 $X/cohort.t1 > /dev/null && echo "other cohort outputs identical" || echo "OTHER COHORT OUTPUTS DIFFER"
echo "== log order: species names in the 6-thread log follow the species list"
diff <(grep -E '^s__' $X/cohort.t6.log | head -20) <(grep -E '^s__' $X/cohort.t1.log | head -20) > /dev/null && echo "same order" || echo "ORDER DIFFERS"
