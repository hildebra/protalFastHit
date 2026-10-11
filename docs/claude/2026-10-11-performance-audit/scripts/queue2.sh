#!/usr/bin/env bash
# After queue1 (its m5 failed on the SAM paths): m5 again, on cores 0-3.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
W=$HOME/perf9
until grep -q "ALL DONE" "$A/queue1.out" 2> /dev/null; do sleep 20; done
echo "[$(date +%T)] m5 profiling (again)"
bash "$A/m5_profiling.sh" > "$W/m5b.out" 2>&1
echo "[$(date +%T)] $(tail -1 "$W/m5b.out")"
echo "[$(date +%T)] ALL DONE"
