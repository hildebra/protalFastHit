#!/usr/bin/env bash
# pfbench against the head07c build's libraries (07c371b; the seeding code is unchanged since), pinned, niced.
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad/perf4
T=$HOME/mt-work/head07c/src; B=$T/build; W=$HOME/mt-work/pfbench; mkdir -p $W
INC=""; for d in src src/IO src/Hash src/Core src/Alignment src/Utilities src/Profiling src/SequenceUtils src/Taxonomy src/SNPs lib lib/robin lib/tsl lib/wfa2-lib lib/wfa2-lib/wavefront lib/wfa2-lib/utils lib/gzstream; do INC="$INC -I$T/$d"; done
g++ -std=c++20 -O3 -march=x86-64 -mtune=generic -ffp-contract=off -fopenmp $INC -I$B/generated -I$B/zlib-ng -I$T/lib/zlib-ng $SP/pfbench.cpp \
  $B/src/libprotal_lib.a $B/src/gzstream/libgzstream_lib.a $B/zlib-ng/libz-ng.a -lzstd -ldeflate -lpthread -o $W/pfbench 2> $W/build.log || { echo "FAIL build"; grep -m5 -A3 error $W/build.log; exit 1; }
cat /proc/loadavg
DB=$HOME/bench071/V071/protal_db/database.protal
R=$HOME/bench071/samples/points/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz
taskset -c ${CPU:-2} nice $W/pfbench $DB $R ${READS:-400000} ${ROUNDS:-3}
cat /proc/loadavg
