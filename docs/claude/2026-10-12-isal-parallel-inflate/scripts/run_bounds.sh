#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
rm -rf ~/isalpar_build
bash $D/build.sh
mkdir -p $D/logs
cat ~/gzpar/data/R1.member.gz > /dev/null
cd $D
taskset -c 0-3 nice -n 5 $D/bin/isalpar bounds ~/gzpar/data/R1.member.gz $D/logs/R1.bounds.tsv 2>&1 | tee $D/logs/bounds.log
head -3 $D/logs/R1.bounds.tsv; wc -l $D/logs/R1.bounds.tsv
taskset -c 0-3 nice -n 5 $D/bin/isalpar seq ~/gzpar/data/R1.member.gz 2>&1 | tee $D/logs/seq_first.log
taskset -c 0-3 nice -n 5 $D/bin/isalpar seq ~/gzpar/data/R1.member.gz 128 2>&1 | tee -a $D/logs/seq_first.log
