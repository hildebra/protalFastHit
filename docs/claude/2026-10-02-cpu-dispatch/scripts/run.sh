#!/usr/bin/env bash
# run.sh BINARY: the instruction-set scan of a protal build, its lzcnt/tzcnt, its dynamic libraries, and the
# launcher's choice with fake /proc/cpuinfo contents. Results into ../results.
set -uo pipefail
D=$(cd "$(dirname "$0")" && pwd); O=$D/../results; mkdir -p $O
bash $D/isa_check.sh "$1" | c++filt > $O/isa.txt
bash $D/lzcnt_check.sh "$1" > $O/lzcnt_tzcnt.txt
echo "lzcnt instructions in the whole binary: $(objdump -d --no-show-raw-insn "$1" | grep -c -P '\tlzcnt')" >> $O/lzcnt_tzcnt.txt
ldd "$1" > $O/ldd.txt
bash $D/launcher_test.sh > $O/launcher.txt
