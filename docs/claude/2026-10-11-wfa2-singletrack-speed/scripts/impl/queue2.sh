#!/usr/bin/env bash
# queue2.sh: after queue1 (whole runs, e2e): lto.sh, then asan.sh.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfaspeed
until grep -q QUEUE1DONE $S/queue1.log; do sleep 10; done
echo "== lto ($(date +%T))"
bash $S/lto.sh 9 5
echo "== asan ($(date +%T))"
bash $S/asan.sh
echo "QUEUE2DONE $(date +%T)"
