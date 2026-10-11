#!/usr/bin/env bash
# repro3.sh REF READ: repro3.cpp built against the reference tree (WFA2-lib v2.3.6) and the work tree (the fork, patched).
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad
for t in ref work; do
  H=$HOME/stk/$t
  g++ -O2 -std=c++17 -I$H/lib -I$H/lib/wfa2-lib $S/repro3.cpp $H/build/src/libwfa_lib.a -o $HOME/stk/repro/repro3_$t 2>&1 | grep error | head -5
  echo "-- $t"; $HOME/stk/repro/repro3_$t "$1" "$2"
done
