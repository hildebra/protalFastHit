#!/usr/bin/env bash
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
for v in ref lr reseed; do echo "== $v"; bash $S/perf3/calls.sh $HOME/mt-work/perf3/cg.ont_b3000000.$v | head -16; done
