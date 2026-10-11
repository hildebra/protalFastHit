#!/usr/bin/env bash
# repro.sh: compiles repro.cpp against ~/stk/work and runs it on the failing pair.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad
H=$HOME/stk/work; P=${1:-$HOME/stk/repro/pattern.txt}; T=${2:-$HOME/stk/repro/text.txt}
g++ -O2 -std=c++20 -I$H/src -I$H/src/IO -I$H/src/Hash -I$H/src/Core -I$H/src/SequenceUtils -I$H/src/Alignment -I$H/src/Utilities -I$H/src/Taxonomy -I$H/lib -I$H/lib/robin -I$H/lib/tsl -I$H/lib/wfa2-lib -I$H/lib/wfa2-lib/wavefront -I$H/lib/wfa2-lib/utils -I$H/lib/gzstream -I$H/build/zlib-ng -I$H/lib/zlib-ng -I$H/build/generated \
  $S/repro.cpp $H/build/src/libwfa_lib.a -o $HOME/stk/repro/repro 2>&1 | grep -E "error" | head
taskset -c 0-3 $HOME/stk/repro/repro $P $T
