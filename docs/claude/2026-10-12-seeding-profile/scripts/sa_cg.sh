#!/usr/bin/env bash
# Exact: the instructions of the seed and anchor finding, by function and by source line, for 100,000 pairs aligned on
# one thread (--no_profile) by HEAD built with -g; and seedbench on 4 pinned cores. After sa_setup.sh.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
until grep -q "SETUP DONE" "$A/sa_setup.out" 2> /dev/null; do sleep 15; done
W=$HOME/seed; B=$W/head/build; S=$W/samples; R=$W/cg; mkdir -p "$R"
DB=$W/db/database.protal
echo "[$(date +%T)] the run's seeding line, 1M pairs, 4 threads"
rm -rf "$R/run1m"
taskset -c 0-3 nice -n 5 "$B/protal" --db "$DB" -1 "$S/pe1m/reads/pe1m_1_R1.fq.gz" -2 "$S/pe1m/reads/pe1m_1_R2.fq.gz" -o "$R/run1m" -t 4 --no_profile \
    > "$R/run1m.log" 2>&1 || echo "FAILED run1m"
grep -E "seeding:|Aligning reads took|Load Index took" "$R/run1m.log" | cut -c1-700
for f in "$R/run1m"/misc/*_runtime.tsv; do [ -e "$f" ] && sort -t$'\t' -k3,3nr "$f" | head -14; done
echo "[$(date +%T)] callgrind, 100k pairs, 1 thread"
rm -rf "$R/out"
taskset -c 1 nice -n 5 valgrind --tool=callgrind --dump-instr=yes --callgrind-out-file="$R/callgrind.out" \
    "$B/protal" --db "$DB" -1 "$S/pe100k_R1.fq.gz" -2 "$S/pe100k_R2.fq.gz" -o "$R/out" -t 1 --no_profile > "$R/cg.log" 2>&1 || echo "FAILED callgrind"
callgrind_annotate --inclusive=yes --threshold=99.9 "$R/callgrind.out" 2> /dev/null > "$R/inclusive.txt"
callgrind_annotate --inclusive=no --threshold=99 "$R/callgrind.out" 2> /dev/null > "$R/self.txt"
SRC=$W/head/src/src
callgrind_annotate --auto=yes --show-percs=yes --context=0 --include="$SRC" --include="$SRC/Core" --include="$SRC/Hash" "$R/callgrind.out" \
    "$SRC/Core/ChainAnchorFinder.h" "$SRC/Hash/KmerLookup.h" "$SRC/Hash/FlexScan.h" "$SRC/Hash/Seedmap.h" 2> /dev/null > "$R/lines.txt"
head -60 "$R/inclusive.txt" | cut -c1-200
grep -E "seeding:" "$R/cg.log" | cut -c1-500
echo "[$(date +%T)] seedbench, 500k pairs, 4 threads"
taskset -c 0-3 nice -n 5 "$B/seedbench" "$DB" "$S/pe1m/reads/pe1m_1_R1.fq.gz" "$S/pe1m/reads/pe1m_1_R2.fq.gz" 500000 4 3 > "$R/seedbench.tsv" 2> "$R/seedbench.log"
echo "exit $?"; cat "$R/seedbench.log" "$R/seedbench.tsv"
echo "[$(date +%T)] CG DONE"
