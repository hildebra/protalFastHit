#!/usr/bin/env bash
# Alternated runs of binaries given as NAME=PATH on mapped and mixed, at 1 and 6 threads; prints the read timer and the parse sums.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/37e6327e-2743-4622-91c7-90dfd0315ee2/scratchpad
cd ~/perf-gtdb/samread
for round in 1 2 3; do for b in "$@"; do name=${b%%=*}; bin=${b#*=}; for f in mapped mixed; do for t in ${THREADS:-1 6}; do
  PROTAL=$bin TAG=$name.$round FILES=$f THREADS=$t bash $S/samread/run.sh > /dev/null
  o=out.$name.$round.$f.t$t.log
  echo "$name r$round $f t$t: $(grep -oE 'reading the SAM [0-9.]+s' $o) | $(grep -oE 'the reading thread: [0-9.]+ s reading, [0-9.]+ s cutting|parsing [0-9.]+ s' $o | tr '\n' ' ')"
done; done; done; done
