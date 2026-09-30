export OUT=~/protal-mem/runs
S=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-09-30-memory-profiling/scripts
r=~/protal-perf/reads/w900/reads
for t in 1 8; do echo "== default t$t"; bash $S/stage_peaks.sh $OUT/db900_w900_t$t; done
bash $S/mem_trace.sh w900_nostrain_t1 ~/protal-perf/db900 $r 1 --no_strains >/dev/null 2>&1; echo "== --no_strains t1"; bash $S/stage_peaks.sh $OUT/w900_nostrain_t1
for t in 2 4 16 32; do bash $S/mem_trace.sh w900_t$t ~/protal-perf/db900 $r $t >/dev/null 2>&1; echo "== default t$t"; bash $S/stage_peaks.sh $OUT/w900_t$t; done
