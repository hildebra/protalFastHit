#!/usr/bin/env bash
# The clones build (protal, protal_tests) and an x86-64-v2 baseline build of d11381f, one after the other.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/build_wt.sh oneb $S/clones1.patch d11381f "protal protal_tests"
bash $S/build_wt.sh oneb-v2 $S/v2.patch d11381f "protal"
echo BUILDS DONE
