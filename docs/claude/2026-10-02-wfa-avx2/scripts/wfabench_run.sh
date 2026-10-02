#!/usr/bin/env bash
# wfabench against the wfa build's libraries (f07acb3 + wfa.patch), niced, after the timed whole runs.
SP=$(cd "$(dirname "$0")" && pwd)  # this folder

T=$HOME/mt-work/wfa/src; B=$T/build; W=$HOME/mt-work/wfabench; mkdir -p $W
WFA=$(find $B -name 'libwfa_lib.a' | head -1)
INC=""; for d in src src/IO src/Hash src/Core src/Alignment src/Utilities src/Profiling src/SequenceUtils src/Taxonomy src/SNPs lib lib/robin lib/tsl lib/wfa2-lib lib/wfa2-lib/wavefront lib/wfa2-lib/utils lib/gzstream; do INC="$INC -I$T/$d"; done
g++ -std=c++20 -O3 -march=x86-64 -fopenmp $INC -I$B/generated -I$B/zlib-ng -I$T/lib/zlib-ng $SP/wfabench.cpp $WFA -o $W/wfabench || { echo "FAIL build"; exit 1; }
cat /proc/loadavg
nice $W/wfabench 7
for d in 0.01 0.03 0.05 0.10 0.15 0.20 0.30; do echo "divergence $d"; nice $W/wfabench 5 $d; done
echo WFABENCH DONE
