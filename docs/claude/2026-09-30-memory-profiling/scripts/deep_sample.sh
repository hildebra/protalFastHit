export OUT=~/protal-mem/runs
S=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-09-30-memory-profiling/scripts
r=~/protal-perf/reads/w900/reads
mkdir -p ~/protal-mem/reads/w900x8
for m in 1 2; do rm -f ~/protal-mem/reads/w900x8/w900x8_R$m.fq.gz; for i in 1 2 3 4 5 6 7 8; do cat $r/w900_1_R$m.fq.gz >> ~/protal-mem/reads/w900x8/w900x8_R$m.fq.gz; done; done
bash $S/mem_trace.sh w900x8_t8 ~/protal-perf/db900 ~/protal-mem/reads/w900x8 8 >/dev/null 2>&1
grep -E "took|Output alignments|Processed" $OUT/w900x8_t8/stdout.ts | cut -c1-120
ls -la $OUT/w900x8_t8/out/alignments 2>/dev/null
bash $S/stage_peaks.sh $OUT/w900x8_t8
