#!/bin/bash
# Each unit of the split prototype and the original main.cpp compiled alone, one after the other, on one core:
# per-unit times free of the other units' contention.
W=~/ctime
D=$W/split
R=$W/res
cd $D
BASE=$(cat base_cmd.txt)
uptime
one() {
    local name=$1 src=$2
    /usr/bin/time -f "$name: %e s wall, %U s user, %M KB" taskset -c 2 nice -n 5 $BASE -ftime-report -c $src -o $D/seq_$name.o 2> $D/seqft_$name.txt
    tail -1 $D/seqft_$name.txt
}
one original $W/src/src/main.cpp
for u in main split_align split_buildindex split_options split_profile split_strains split_buildmaint; do
    one $u $u.cpp
done
uptime
for u in original main split_align split_buildindex split_options split_profile split_strains split_buildmaint; do
    f=seqft_$u.txt
    echo "$u $(grep -E '^ TOTAL' $f) parse=$(grep 'phase parsing' $f | awk -F: '{print $2}' | awk '{print $1}') deferred=$(grep 'phase lang. deferred' $f | awk -F: '{print $2}' | awk '{print $1}') opt=$(grep 'phase opt and generate' $f | awk -F: '{print $2}' | awk '{print $1}') text=$(size seq_$u.o | awk 'NR==2 {print $1}')"
done
