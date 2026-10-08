#!/usr/bin/env bash
G=$HOME/gzpar; D=$G/data; T="taskset -c 0-3 nice -n 5"; W=$G/work
SC=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/cdfd1fc3-105f-41f2-9438-012dde2a0acc/scratchpad/gz
cd $G/rb && cp $SC/onefile.cpp . && $T g++ -std=c++20 -O3 -march=x86-64-v2 -I$W/src -I$W/src/IO -I$W/src/Utilities -I$W/lib onefile.cpp -lisal -lzstd -lpthread -o onefile || exit 1
cat $D/R1.bgzf.gz $D/R1.member.gz > /dev/null
uptime
for rep in 1 2 3; do
  echo -e "member\t$($T ./onefile $D/R1.member.gz)"
  for n in 1 2 3; do echo -e "bgzf\t$(PROTAL_INFLATE_THREADS=$n $T ./onefile $D/R1.bgzf.gz)"; done
done
uptime
