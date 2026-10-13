#!/bin/bash
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
bash $D/build.sh >/dev/null
cat ~/gzpar/data/R1.member.gz > /dev/null
for r in 1 2; do for impl in 0 1; do
  echo "MR_IMPL=$impl round $r"
  MR_IMPL=$impl taskset -c 0-3 nice -n 5 $D/bin/isalpar cost ~/gzpar/data/R1.member.gz $D/logs/R1.bounds.tsv 4 - 2>&1 | grep -v markers
done; done | tee $D/logs/mr_impl.log
