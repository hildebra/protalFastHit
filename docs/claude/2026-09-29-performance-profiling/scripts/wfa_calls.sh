#!/bin/bash
# Every WFA call of a run, with the anchor chain it starts from, and the anchored-extension
# prototype next to it: wfa_calls.sh DB READS_DIR LABEL
# Builds $PERF_DIR/instr from $PERF_DIR/src with wfa_calls.patch (logging only; outputs unchanged),
# runs 1 thread on READS_DIR/*_R{1,2}.fq.gz and tabulates the log.
# Log columns: 1 read length 2 links 3 exact bases in links 4 bases left of the chain 5 right of it
# 6 indel the seeds imply 7-8 dovetail left/right 9 window length 10-11 free read ends 12 aligned
# 13 score 14 mismatches 15-16 internal I/D operations 17 ns 18 mismatches on the anchor diagonal
# 19 extension aligned (-1: not a single-link anchor) 20 extension score 21 extension ns
set -u
source "$(dirname "$0")/env.sh"
db=$1; reads=$2; label=$3
I=$PERF_DIR/instr
if [ ! -x $I/build/protal_avx2 ]; then
  rm -rf $I; mkdir -p $I; cp -a $PERF_DIR/src $I/src
  (cd $I/src && patch -p1 < $here/wfa_calls.patch) || exit 1
  cmake -S $I/src -B $I/build -G Ninja -DCMAKE_BUILD_TYPE=Release > $I/cmake.log 2>&1
  cmake --build $I/build --target protal_avx2 -j "$(nproc)" > $I/make.log 2>&1 || exit 1
fi
PROTAL_WFA_LOG=$I/wfa_$label.tsv $I/build/protal_avx2 --db $db -1 $(ls $reads/*_R1.fq.gz) -2 $(ls $reads/*_R2.fq.gz) \
  -o $I/out_$label -t 1 --no_qcmsa --no_profile > $I/run_$label.log 2>&1
f=$I/wfa_$label.tsv.0
echo "$label: $(wc -l < $f) WFA calls"
awk -F'\t' '
{ n++; T += $17
  if (!$12) c = "1 failed"; else if ($15 + $16 > 0) c = "2 aligned, internal indel"; else if ($14 == 0) c = "3 aligned, exact"
  else if ($14 == 1) c = "4 aligned, 1 mismatch"; else if ($14 <= 5) c = "5 aligned, 2-5 mismatches"; else c = "6 aligned, >5 mismatches"
  cn[c]++; ct[c] += $17
  if ($2 == 1) { s++; full += $17; ext += $21
    if ($12 && $19 == 1) { both++; same += $20 == $13; better += $20 > $13; worse += $20 < $13 }
    else if ($12) lost++; else if ($19 == 1) gained++ } }
END {
  printf "WFA %.3f s, %.0f ns per call\n", T / 1e9, T / n
  for (c in cn) printf "  %-28s %7d calls %5.1f%%  %6.0f ns/call  %5.1f%% of WFA time\n", substr(c, 3), cn[c], 100*cn[c]/n, ct[c]/cn[c], 100*ct[c]/T
  printf "single-link anchors %d (%.1f%%): full WFA %.3f s, anchored extension %.3f s (%.2fx); both aligned %d: same score %d, extension better %d, worse %d; lost %d, gained %d\n",
    s, 100*s/n, full/1e9, ext/1e9, full/ext, both, same, better, worse, lost, gained }' $f | sort
