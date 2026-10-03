#!/usr/bin/env bash
# After ab_packed.sh: what differs in the two profile logs, the read file names of the other worlds,
# a one-thread A/B on mix (all files byte for byte), and the RSS per stage (mem_trace.sh of the
# 2026-09-30 memory report) for the baseline and the packed build.
WORK=${WORK:-$HOME/protal-pack}
PERF=${PERF:-$HOME/protal-perf}
SRC=${SRC:-/mnt/c/Users/hildebra/Documents/locDev/protal}
MEM=$SRC/docs/claude/2026-09-30-memory-profiling/scripts
OUT=$WORK/ab
echo "== differing files (mix, 6 threads) =="
for f in mix_R.profile.gene.log mix_R.profile.genes.log; do
  echo "-- $f: $(diff "$OUT/mix.base/$f" "$OUT/mix.new/$f" | grep -c '^<') lines differ"
  diff "$OUT/mix.base/$f" "$OUT/mix.new/$f" | head -6 | cut -c1-200
done
echo "== read files =="
ls "$PERF/reads/w900" "$PERF/reads/dense_w" | head -12
echo "== mix, 1 thread, every file byte for byte =="
r1=$(ls "$PERF"/reads/mix/*_R1.fq* | head -1); r2=$(ls "$PERF"/reads/mix/*_R2.fq* | head -1)
for v in base new; do
  bin=$WORK/base/build/protal; [ $v = new ] && bin=$WORK/build/protal
  rm -rf "$OUT/mix1.$v"
  "$bin" --db "$PERF/db900n" -1 "$r1" -2 "$r2" -o "$OUT/mix1.$v" -t 1 --no_qcmsa > "$OUT/mix1.$v.stdout" 2>&1
done
cd "$OUT/mix1.base" && for f in $(find . -type f | sort); do
  case "$f" in *_runtime.tsv) continue ;; esac
  cmp -s "$f" "$OUT/mix1.new/$f" || echo "  differs: $f"
done
echo "  compared $(find . -type f | wc -l) files"
cd "$WORK" || exit 1
echo "== RSS per stage (mix, 6 threads) =="
for v in base new; do
  bin=$WORK/base/build/protal; [ $v = new ] && bin=$WORK/build/protal
  OUT="$WORK/trace" BIN="$bin" bash "$MEM/mem_trace.sh" "mix_$v" "$PERF/db900n" "$PERF/reads/mix" 6 > /dev/null 2>&1
  echo "-- $v"; bash "$MEM/stage_peaks.sh" "$WORK/trace/mix_$v" 2>&1 | head -12
done
