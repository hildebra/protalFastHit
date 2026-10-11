#!/usr/bin/env bash
# nonull.sh: a scratch copy of the work tree's WFA2-lib with wf_heuristic_nullify_boundaries (the fork's nulling of the
# elements next to a heuristic-pruned wavefront) disabled, through wfadiff: classic against v2.3.6 (diff_ref.tsv),
# Singletrack against classic.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad
R=$HOME/stk/repro; rm -rf $R/inc2; mkdir -p $R/inc2; cp -r $HOME/stk/work/lib/wfa2-lib $R/inc2/wfa2-lib
L=$R/inc2/wfa2-lib
perl -0pi -e 's/(  if \(wavefront == NULL \|\| wavefront->null\) return;)/  return;  \/\/ disabled (nonull.sh)\n$1/' $L/wavefront/wavefront_heuristic.c
grep -n "disabled (nonull.sh)" $L/wavefront/wavefront_heuristic.c
cd $L && rm -f *.o
for f in $(sed -n '/set(WFA2_SOURCES/,/)/p' $HOME/stk/work/lib/wfa2-lib.cmake | grep -oE '[a-z_0-9]+/[a-z_0-9]+\.(c|cpp)') bindings/cpp/WFAligner.cpp; do
  case $f in *.cpp) g++ -O2 -std=c++17 -I. -c $f -o $(basename $f).o ;; *) gcc -O2 -I. -c $f -o $(basename $f).o ;; esac
done
ar rcs libwfa.a *.o
g++ -O2 -std=c++17 -DSINGLETRACK -I$R/inc2 -I$L $S/wfadiff.cpp $L/libwfa.a -o $R/wfadiff_nonull 2>&1 | grep error | head
taskset -c 0-3 $R/wfadiff_nonull > $R/diff_nonull.tsv
taskset -c 0-3 $R/wfadiff_nonull singletrack > $R/diff_nonull_single.tsv
for f in nonull nonull_single; do
  echo "$f: $(awk -F'\t' '$3==0' $R/diff_$f.tsv | wc -l) completed, $(awk -F'\t' '$5==0' $R/diff_$f.tsv | wc -l) invalid, $(awk -F'\t' '$3==-300' $R/diff_$f.tsv | wc -l) unattainable"
done
echo "ref vs nonull classic: $(paste $R/diff_ref.tsv $R/diff_nonull.tsv | awk -F'\t' '$3!=$9 || $4!=$10 || $6!=$12' | wc -l) lines differ"
echo "nonull classic vs singletrack: outcome or score differ in $(paste $R/diff_nonull.tsv $R/diff_nonull_single.tsv | awk -F'\t' '$3!=$9 || $4!=$10' | wc -l), operations in $(paste $R/diff_nonull.tsv $R/diff_nonull_single.tsv | awk -F'\t' '$6!=$12' | wc -l)"
paste $R/diff_nonull.tsv $R/diff_nonull_single.tsv | awk -F'\t' '$3!=$9 || $4!=$10 {print "  kind " $1 " case " $2 ": status " $3 "/" $9 " score " $4 "/" $10}' | head
