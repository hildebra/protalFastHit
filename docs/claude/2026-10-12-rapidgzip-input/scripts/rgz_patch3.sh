#!/bin/bash
# Clean 2.32.1 + the fork's stopping points (0007, fuzzy); write the resulting diff for review.
set -uo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
cd ~/rgz/isal
git reset -q --hard HEAD~1 2>/dev/null || true
git log --oneline | head -3
git reset -q --hard
git clean -qfd
patch -p1 -l --fuzz=3 --no-backup-if-mismatch < ~/rgz/split/0007
git diff > $SP/isal_0007_applied.diff
wc -l $SP/isal_0007_applied.diff
