#!/bin/bash
for t in test test_congeners; do echo "=== $t"; ~/protal-train/bin/python -W ignore /mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/fbe611b7-59a8-4a3b-b159-58842b2426b7/scratchpad/tune/tune_eval.py ~/tune ~/audit4/bin/protal $t "$@"; done
