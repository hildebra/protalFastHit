#!/usr/bin/env bash
# Build + all unit tests, then scripts/measure_performance.sh on the v0.7.3 world (pe 500k + PacBio 90 Mb), one repeat,
# with the cohort run and qcMSA.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
bash $S/build_work.sh '*' 8 || exit 1
W=$HOME/perf-gtdb; cd $W/work || exit 1
P=$HOME/bench071/samples/points
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
rm -rf $W/mp
eval "$($HOME/micromamba/bin/micromamba shell hook -s bash 2> /dev/null)" 2> /dev/null; micromamba activate protal 2> /dev/null || true
command -v python3; python3 -c 'import numpy; print("numpy", numpy.__version__)' 2>&1 | tail -1
THREADS=6 REPEATS=1 PROTAL=$W/work/build/protal bash scripts/measure_performance.sh $W/mp $W/e2e/db/database.protal pe:$R1:$R2 pb:$PB 2>&1 | tail -30
echo "== runs.tsv"; cat $W/mp/runs.tsv | cut -f1,3,5,12,13,26-33
echo "== cohort.tsv"; cat $W/mp/cohort.tsv
echo "== cohort log"; grep -E 'All alignments|took|Strain-level MSAs of|SAM header|qcmsa\]|seeding' $W/mp/logs/cohort.1.log | grep -v 'Profiling sample' | head -20 | cut -c1-200
echo "== pe log seeding/header"; grep -E 'seeding:|SAM header' $W/mp/logs/pe1.1.log | cut -c1-400
