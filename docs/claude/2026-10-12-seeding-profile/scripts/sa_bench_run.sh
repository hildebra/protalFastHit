#!/usr/bin/env bash
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
bash "$A/sa_bench_compile.sh" | tail -2
W=$HOME/seed; S=$W/samples; R=$W/cg
cat /proc/loadavg
taskset -c 0-3 nice -n 5 $W/head/build/seedbench $W/db/database.protal "$S/pe1m/reads/pe1m_1_R1.fq.gz" "$S/pe1m/reads/pe1m_1_R2.fq.gz" 500000 4 3 > "$R/seedbench2.tsv" 2> "$R/seedbench2.log"
echo "exit $?"; cat "$R/seedbench2.log" "$R/seedbench2.tsv"; cat /proc/loadavg
