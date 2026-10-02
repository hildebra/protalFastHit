#!/usr/bin/env bash
# The clones build again with the final patch (NoteRecord no longer cloned), then callgrind of it on all three
# workloads (cl2), as oneb_cg.sh, and the comparison of its outputs with base's.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/build_wt.sh oneb $S/oneb_full.patch d11381f "protal protal_tests simulate_metagenomes" || exit 1
W=$HOME/mt-work/onebcg; I=$HOME/mt-work/isacg
DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
B=$HOME/mt-work/oneb/src/build/protal
declare -A ARGS=(
  [pe100k]="-1 $I/r1.fq -2 $I/r2.fq --no_profile"
  [ont3M]="-1 $P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz --read_type ont --no_profile"
  [prof500k]="--profile_only $HOME/bench071/runs/v071.full.pe.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst" )
one() {
  local w=$1 b=cl2; rm -rf $W/out.$w.$b
  nice valgrind --tool=callgrind --callgrind-out-file=$W/cg.$w.$b $B --db $DB ${ARGS[$w]} --prefix s -o $W/out.$w.$b -t 1 --no_qcmsa > $W/log.$w.$b 2>&1
  callgrind_annotate --threshold=100 $W/cg.$w.$b 2>/dev/null | c++filt > $W/$w.$b.txt
  echo "done $w $b: $(grep 'PROGRAM TOTALS' $W/$w.$b.txt)"
  r=$(diff -r -q -x '*_runtime.tsv' -x '*.sam.zst' $W/out.$w.base $W/out.$w.$b > /dev/null && echo same || echo DIFFER)
  s=""; [ -f $W/out.$w.base/s.sam.zst ] && { cmp -s <(zstdcat $W/out.$w.base/s.sam.zst) <(zstdcat $W/out.$w.$b/s.sam.zst) && s="SAM text same" || s="SAM text DIFFER"; }
  echo "$w $b vs base: outputs $r $s"
}
for w in pe100k ont3M prof500k; do one $w & done; wait
echo REBUILD DONE
