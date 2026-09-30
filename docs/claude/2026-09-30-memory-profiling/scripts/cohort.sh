export OUT=~/protal-mem/runs
S=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-09-30-memory-profiling/scripts
n=$1; t=$2
r=~/protal-perf/reads/w900/reads
m=~/protal-mem/cohort$n.map
{ printf '#OUTPUT_DIR\t%s/cohort%s_out\n#INPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\n' $OUT $n $r
  for i in $(seq 1 $n); do printf 's%s\tw900_1_R1.fq.gz\tw900_1_R2.fq.gz\ts%s.sam.zst\ts%s\ts%s.profile\n' $i $i $i $i; done; } > $m
MAP=$m bash $S/mem_trace.sh cohort${n}_t$t ~/protal-perf/db900 x $t $EXTRA >/dev/null 2>&1
grep -E "took|MSA" $OUT/cohort${n}_t$t/stdout.ts | cut -c1-110 | tail -6
bash $S/stage_peaks.sh $OUT/cohort${n}_t$t
ls $OUT/cohort${n}_out/strains 2>/dev/null | head -3; du -sh $OUT/cohort${n}_out/alignments
