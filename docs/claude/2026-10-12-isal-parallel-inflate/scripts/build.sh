#!/bin/bash
set -e
D=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad/isal_parallel
B=$D/bin
mkdir -p $B
cp $D/isalpar.cpp $B/
cd $B
taskset -c 0-3 nice -n 5 g++ -O2 -march=native -std=c++17 -pthread -Wall -Wextra -o isalpar isalpar.cpp \
  -I ~/rgzbuild/build/isa-l-install/include ~/rgzbuild/build/isa-l-install/lib/libisal.a
ls -la $B/isalpar
