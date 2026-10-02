#!/usr/bin/env bash
# The records of the 5M-pair run (its .sam.zst, 400 MB, as a stand-in for the records file) appended behind a
# 2 MB header: copy_file_range and pread/pwrite on 1, 2, 4 and 6 threads, three passes, niced.
set -uo pipefail
SP=$(cd "$(dirname "$0")" && pwd)  # this folder
W=$HOME/mt-work/copy; mkdir -p $W
g++ -std=c++20 -O2 -pthread $SP/copybench.cpp -o $W/copybench || exit 1
cp $HOME/mt-audit/runs/b6/s1.sam.zst $W/records.partial
sync; cat $W/records.partial > /dev/null
cat /proc/loadavg
for rep in 1 2 3; do
  for mode in cfr rw; do for t in 1 2 4 6; do
    nice $W/copybench $W/records.partial $W/out.sam.zst $mode $t 2097152
  done; done
done
echo COPY DONE
