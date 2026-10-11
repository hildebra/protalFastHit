#!/usr/bin/env bash
# wfadiff.sh: wfadiff.cpp against the reference tree's WFA2-lib (v2.3.6) and the work tree's (the fork, patched; its
# library rebuilt from the current sources first), classic in both and Singletrack in the work tree: the three outputs
# compared line by line (status, score, validity, operations).
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad
R=$HOME/stk/repro
(cd $HOME/stk/work && taskset -c 0-3 nice cmake --build build -j4 --target wfa_lib > $R/wfa_lib.log 2>&1) || { echo "FAIL wfa_lib"; tail $R/wfa_lib.log; }
for t in ref work; do
  H=$HOME/stk/$t; def=""; [ $t = work ] && def="-DSINGLETRACK"
  g++ -O2 -std=c++17 $def -I$H/lib -I$H/lib/wfa2-lib $S/wfadiff.cpp $H/build/src/libwfa_lib.a -o $R/wfadiff_$t 2>&1 | grep error | head -5
done
taskset -c 0-3 $R/wfadiff_ref > $R/diff_ref.tsv
taskset -c 0-3 $R/wfadiff_work > $R/diff_work.tsv
taskset -c 0-3 $R/wfadiff_work singletrack > $R/diff_single.tsv
for f in ref work single; do
  echo "$f: $(wc -l < $R/diff_$f.tsv) alignments, $(awk -F'\t' '$3==0' $R/diff_$f.tsv | wc -l) completed, $(awk -F'\t' '$5==0' $R/diff_$f.tsv | wc -l) with invalid operations"
done
echo "ref vs work classic: $(paste $R/diff_ref.tsv $R/diff_work.tsv | awk -F'\t' '$3!=$9 || $4!=$10 || $6!=$12' | wc -l) lines differ"
paste $R/diff_ref.tsv $R/diff_work.tsv | awk -F'\t' '$3!=$9 || $4!=$10 || $6!=$12 {print "  kind " $1 " case " $2 ": status " $3 "/" $9 " score " $4 "/" $10 " valid " $5 "/" $11}' | head -10
echo "work classic vs singletrack: outcome or score differ in $(paste $R/diff_work.tsv $R/diff_single.tsv | awk -F'\t' '$3!=$9 || $4!=$10' | wc -l), operations in $(paste $R/diff_work.tsv $R/diff_single.tsv | awk -F'\t' '$6!=$12' | wc -l)"
paste $R/diff_work.tsv $R/diff_single.tsv | awk -F'\t' '$3!=$9 || $4!=$10 {print "  kind " $1 " case " $2 ": status " $3 "/" $9 " score " $4 "/" $10 " valid " $5 "/" $11}' | head -10
for k in 0 1 2 3; do
  echo "  kind $k: $(paste $R/diff_work.tsv $R/diff_single.tsv | awk -F'\t' -v k=$k '$1==k && $3==0' | wc -l) completed, other operations $(paste $R/diff_work.tsv $R/diff_single.tsv | awk -F'\t' -v k=$k '$1==k && $6!=$12' | wc -l)"
done
