#!/bin/bash
# Checks of the final binary (the vote's floor rule built in, default --mate_guidance consensus):
#  1. its default run gives the records of the evaluation's consensus_floor runs (the same rule switched on by an
#     environment variable in the evaluation build): two points of the finished database;
#  2. thread independence: one sample on 1 and on 4 threads, records and profile;
#  3. guide reproduces the build's own runs (protal 0.8.4 at 5369c99, the training database): the rl150 test points'
#     original SAMs against this binary's --mate_guidance guide runs (run_all.sh's tdb/guide).
# usage: check_final.sh <protal> <runs folder> <out folder>     (WSL; 4 cores)
set -u
PROTAL=$1; RUNS=$2; OUT=$3
HERE=$(cd "$(dirname "$0")" && pwd)
PY=$HOME/micromamba/envs/protal-db-build/bin/python
FIN=$HOME/c4test/db/protal_db
WORK=$HOME/c4test/db/work
P=$WORK/test/points
mkdir -p $OUT
POINTS="rl150_p100000 rl250_p100000" bash $HERE/run_sims.sh $PROTAL consensus $OUT/final > /dev/null
specs=()
for point in rl150_p100000 rl250_p100000; do
  for sam in $OUT/final/consensus/$point/alignments/*.sam.zst; do
    s=$(basename $sam .sam.zst)
    specs+=("$s:$RUNS/fin/consensus_floor/$point/alignments/$s.sam.zst:$sam")
  done
done
$PY -I $HERE/compare_records.py $OUT/final_vs_floor.tsv floor final $FIN "${specs[@]}" 2> /dev/null
echo "1. final against consensus_floor: $(tail -1 $OUT/final_vs_floor.tsv | cut -f1-6)"
R=$HOME/mc/sims/rl150_p100000/sim/reads
for t in 1 4; do
  taskset -c 0-3 nice -n 5 $PROTAL --db $FIN -1 $R/rl150_p100000_s_1_R1.fq.zst -2 $R/rl150_p100000_s_1_R2.fq.zst --prefix t$t -o $OUT/threads \
      -t $t --no_strains --no_qcmsa --force > $OUT/threads_t$t.log 2>&1
done
zstdcat $OUT/threads/t1.sam.zst | grep -v '^@' | sort > $OUT/t1.sam
zstdcat $OUT/threads/t4.sam.zst | grep -v '^@' | sort > $OUT/t4.sam
echo "2. records on 1 and 4 threads: $(wc -l < $OUT/t1.sam) and $(wc -l < $OUT/t4.sam), differing lines $(diff $OUT/t1.sam $OUT/t4.sam | grep -c '^[<>]'); profiles $(cmp -s $OUT/threads/t1.profile $OUT/threads/t4.profile && echo identical || echo DIFFER)"
rm -f $OUT/t1.sam $OUT/t4.sam
specs=()
for point in rl150_p10000 rl150_p100000 rl150_p1000000; do
  for sam in $P/$point/protal/alignments/*.sam.zst; do
    s=$(basename $sam .sam.zst)
    specs+=("$s:$sam:$RUNS/tdb/guide/$point/alignments/$s.sam.zst")
  done
done
$PY -I $HERE/compare_records.py $OUT/guide_vs_084.tsv original guide $WORK/training_db --work $WORK "${specs[@]}" 2> /dev/null
echo "3. guide against the build's 0.8.4 runs: $(tail -1 $OUT/guide_vs_084.tsv | cut -f1-6)"
