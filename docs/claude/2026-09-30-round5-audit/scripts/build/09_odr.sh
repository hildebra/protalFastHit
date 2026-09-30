#!/bin/bash
# protal_avx2 = main.cpp at -march=x86-64-v3 + libprotal_lib.a at -march=x86-64: how many inline
# (weak/COMDAT) functions are defined in both, i.e. compiled under two ISAs with one copy kept?
set -u
A=~/audit6/build
B=$A/src/build
M=$B/CMakeFiles/protal_avx2.dir/src/main.cpp.o
L=$B/src/libprotal_lib.a
nm --defined-only $M | awk '$2=="W"||$2=="V"||$2=="u"{print $3}' | sort -u > $A/odr_main_weak.txt
nm --defined-only $L 2>/dev/null | awk '$2=="W"||$2=="V"||$2=="u"{print $3}' | sort -u > $A/odr_lib_weak.txt
echo "weak symbols: main.cpp.o $(wc -l < $A/odr_main_weak.txt), libprotal_lib.a $(wc -l < $A/odr_lib_weak.txt), in both $(comm -12 $A/odr_main_weak.txt $A/odr_lib_weak.txt | wc -l)"
echo "-- sample of shared protal:: symbols"
comm -12 $A/odr_main_weak.txt $A/odr_lib_weak.txt | c++filt | grep 'protal::' | head -8
echo "-- link order (main.cpp.o before the archive: its copies win)"
tr ' ' '\n' < $B/CMakeFiles/protal_avx2.dir/link.txt | grep -nE 'main.cpp.o|libprotal_lib.a'
echo "-- objects in libprotal_lib.a"; ar t $L
echo "-- ymm instructions: protal vs protal_avx2"
for b in protal protal_avx2; do echo "$b: $(objdump -d --no-show-raw-insn $B/$b | grep -c ymm) ymm instructions"; done
echo DONE
