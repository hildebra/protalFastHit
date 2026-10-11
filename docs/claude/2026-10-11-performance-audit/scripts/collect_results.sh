#!/usr/bin/env bash
# Copies the audit's local results into the report's results/local (small text files only).
D=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-11-performance-audit/results/local
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
R=$HOME/perf9/runs
mkdir -p "$D"
cp "$R/align/runs.tsv" "$D/align_runs.tsv"
cp "$R/tables/tables.tsv" "$R/tables/startup.tsv" "$R/tables/table_counts.txt" "$D/" 2> /dev/null
cp "$R/tables/table_files.txt" "$D/table_files.txt" 2> /dev/null
cp "$R/cohort/cohort.tsv" "$D/cohort.tsv"
cp "$R/cohort_evidence/evidence.tsv" "$D/cohort_evidence.tsv"
[ -s "$R/profiling/profiling.tsv" ] && cp "$R/profiling/profiling.tsv" "$R/profiling/sam_finish.tsv" "$D/" 2> /dev/null
bash "$A/cg_funcs.sh" > "$D/callgrind_functions.txt" 2>&1
grep -h "PROGRAM TOTALS" "$R/callgrind/"*.top.txt > /dev/null 2>&1
bash "$A/summarize.sh" > "$D/summary.txt" 2>&1
# the massif profile: the heap per snapshot and the late snapshot's tree (two levels)
bash "$A/massif_last.sh" > "$D/massif_cohort24.txt" 2>&1
grep -E "^(Memory after|Profiling took|Run protal took)" "$R/massif/massif.log" >> "$D/massif_cohort24.txt"
# the local database and builds
{ cat "$HOME/perf9/arcdb.time"; grep -E "took|Memory after|^Congener gaps:|^Column weights( codes)?" "$HOME/perf9/arcdb.build.log"; } > "$D/arc_build.txt"
{ cat "$HOME/perf9/db078.time"; tail -5 "$HOME/perf9/db078.build.log"; } > "$D/db078_build.txt"
for v in v078 v080 v081; do echo "$v $(cat "$HOME/perf9/$v/commit") $("$HOME/perf9/$v/build/protal" --version 2>&1 | head -1)"; done > "$D/binaries.txt"
ls -la "$D"
