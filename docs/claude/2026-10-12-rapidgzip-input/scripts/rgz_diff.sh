#!/bin/bash
# Diff rapidgzip 0.16.0's ISA-L fork against upstream ISA-L 2.32.1 (igzip, crc, include).
set -euo pipefail
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/0ff4a1d7-3599-4615-a178-e4ca253109af/scratchpad
W=~/rgz
mkdir -p $W
cd $W
rm -rf up fork
mkdir up
tar -xzf ~/isal/isa-l-2.32.1.tar.gz -C up
ls up
cp -r $SP/rapidgzip-dl/x/rapidgzip-0.16.0/src/external/isa-l fork
U=$(ls -d up/isa-l-*)
echo "upstream dir: $U"
for d in igzip crc include; do
  diff -ruN $U/$d fork/$d > diff_$d.txt || true
  echo "== $d: $(grep -c '^diff ' diff_$d.txt) files differ, $(wc -l < diff_$d.txt) lines"
  grep '^diff \|^Only' diff_$d.txt | sed 's/^diff -ruN //' | head -60
done
