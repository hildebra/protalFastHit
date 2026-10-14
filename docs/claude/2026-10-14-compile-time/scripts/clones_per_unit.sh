#!/bin/bash
# split_strains.cpp and split_buildmaint.cpp with and without the x86-64-v3 clones: are the profiler's SAM parsing
# (PrepareMAPQ), the index decoder (DecodeChunk) and their callees still compiled into them?
D=~/ctime/split
cd $D
BASE=$(cat base_cmd.txt)
for u in split_strains split_buildmaint; do
    /usr/bin/time -f "$u noclones: %e s wall, %U s user" taskset -c 3 nice -n 5 $BASE -DPROTAL_NO_CLONES -c $u.cpp -o nc_$u.o
    for o in $u.o nc_$u.o; do
        echo "$o: text $(size $o | awk 'NR==2 {print $1}'), PrepareMAPQ $(nm -C $o | grep -c 'PrepareMAPQ'), DecodeChunk $(nm -C $o | grep -c 'DecodeChunk'), NoteRecord $(nm -C $o | grep -c 'NoteRecord'), PackAvx2 $(nm -C $o | grep -c 'PackAvx2'), ancestry::Decode $(nm -C $o | grep -c 'ancestry::Decode('), regex _Compiler $(nm -C $o | grep -c '_Compiler<std::__cxx11::regex_traits')"
    done
done
