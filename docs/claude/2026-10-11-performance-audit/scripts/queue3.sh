#!/usr/bin/env bash
# The evidence experiment: rebuild tim (with the evidence drop), then m4b; on cores 0-3 beside the callgrind runs.
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
echo "[$(date +%T)] build tim"; bash "$A/build_tim.sh"
echo "[$(date +%T)] m4b"; bash "$A/m4b_evidence.sh" > "$HOME/perf9/m4b.out" 2>&1; tail -1 "$HOME/perf9/m4b.out"
cat "$HOME/perf9/runs/cohort_evidence/evidence.tsv" "$HOME/perf9/runs/cohort_evidence/profiles.md5"
echo "[$(date +%T)] ALL DONE"
