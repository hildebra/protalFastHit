#!/usr/bin/env bash
# The unaligned reads' counts in the SAM header: all tests; pe and pb runs by default (header counts) and with --write_unmapped_reads,
# compared with the earlier runs (avx/*.avx2.1, unmapped records): the profiles must be the same, the SAM only lose its unmapped records.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh "${TESTS:-*}" 14 || exit 1
W=$HOME/perf-gtdb; N=$W/work/build/protal; DB=$W/e2e/db/database.protal; X=$W/hdr; rm -rf $X; mkdir -p $X
P=$HOME/bench071/samples/points; R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
outs() { diff -rq -x '*_runtime.tsv' -x '*.sam.zst' -x '*.statistics.tsv' $1 $2 > /dev/null && echo out-same || { echo OUT-DIFF; diff -rq -x '*_runtime.tsv' -x '*.sam.zst' -x '*.statistics.tsv' $1 $2 | head -5; }; }
recs() { zstd -dc $1 | grep -v '^@' | sort; }
for s in pe pb; do
  if [ $s = pe ]; then reads=(-1 $R1 -2 $R2); else reads=(-1 $PB); fi
  for v in counts records; do
    o=$X/$s.$v; extra=(); [ $v = records ] && extra=(--write_unmapped_reads)
    $N --db $DB "${reads[@]}" --read_type $s --prefix $s -o $o -t 6 --no_qcmsa "${extra[@]}" > $o.log 2>&1
    echo "$s $v: $(outs $o $W/avx/$s.avx2.1) $(stat -c %s $o/$s.sam.zst) bytes (before $(stat -c %s $W/avx/$s.avx2.1/$s.sam.zst)); $(grep -h '^SAM header:' $o.log | cut -c1-230)"
  done
  cmp -s <(recs $X/$s.records/$s.sam.zst) <(recs $W/avx/$s.avx2.1/$s.sam.zst) && echo "$s --write_unmapped_reads: records as before" || echo "$s --write_unmapped_reads: RECORDS DIFFER"
  cmp -s <(recs $X/$s.counts/$s.sam.zst) <(recs $W/avx/$s.avx2.1/$s.sam.zst | awk -F'\t' '!and($2,4) || $3 != "*"') && echo "$s counts: records as before less the unmapped ones" || echo "$s counts: RECORDS DIFFER"
  zstd -dc $X/$s.counts/$s.sam.zst | grep -c '^@CO	protal failed candidates' | sed "s/^/$s header lines with counts: /"
done
# Profiles from the SAMs alone: the new one (header) and the old one (records), at 6 and 1 threads.
for t in 6 1; do for v in counts old; do
  o=$X/cohort.$v.t$t; sam=$X/pe.counts/pe.sam.zst; [ $v = old ] && sam=$W/avx/pe.avx2.1/pe.sam.zst
  mkdir -p $o/s; cp $sam $o/s/pe.sam.zst
  { printf '#OUTPUT_DIR\t%s\n' "$o"; printf '#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tREAD_TYPE\n'; printf 'pe\t%s\t%s\t%s\tpe\tpe\n' "$R1" "$R2" "$o/s/pe.sam.zst"; } > $o.map
  $N --db $DB --map $o.map -t $t --no_strains > $o.log 2>&1
  cmp -s $o/pe.profile $W/avx/pe.avx2.1/pe.profile && cmp -s $o/pe.profile.log $W/avx/pe.avx2.1/pe.profile.log && echo "from the $v SAM at $t threads: profile same" || echo "from the $v SAM at $t threads: PROFILE DIFFERS"
done; done
