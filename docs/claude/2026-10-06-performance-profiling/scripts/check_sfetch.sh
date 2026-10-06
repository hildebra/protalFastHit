#!/usr/bin/env bash
# check_sfetch.sh: lib/static-deps.cmake's own nasm and ISA-L, from the local tarballs (no nasm on this PATH, no ISA-L on
# the system: every binary links the fetched one), then protal from that build tree against the ISA-L build's SAM.
# (check_isal.sh's block of the same, again after nasm kept its tarball's file times.)
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/perf8; DB=$HOME/bench071/V075/protal_db; RD=$HOME/perf6/reads
cp $REPO/lib/static-deps.cmake $W/isal/lib/static-deps.cmake
S=$W/isal/build-sfetch; rm -rf $S
( cd $W/isal && nice cmake -S . -B $S -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DPROTAL_STATIC_FETCH_DEPS=ON -DPROTAL_ISAL_URL=$HOME/isal/isa-l-2.32.1.tar.gz -DPROTAL_NASM_URL=$HOME/isal/nasm-2.16.03.tar.xz \
    > $W/sfetch.configure.log 2>&1 ); echo "configure exit $?"
grep -E "nasm for ISA-L|ISA-L:|Static binaries" $W/sfetch.configure.log | sed 's/^/    /'
nice cmake --build $S --target protal -j 5 > $W/sfetch.build.log 2>&1; echo "build exit $?"
grep -E "error|Error" $W/sfetch.build.log | head -5
ls -l $S/static-deps/nasm-install/bin/nasm $S/static-deps/isal-install/lib/libisal.a 2>&1 | awk '{print "    " $5 " " $NF}'
$S/protal --version 2>&1 | head -1 | sed 's/^/    /'
rm -rf $W/cmpi/sfetch
$S/protal --db $DB -1 $RD/pe100k_R1.fq.gz -2 $RD/pe100k_R2.fq.gz --read_type pe --prefix s -o $W/cmpi/sfetch -t 6 --no_qcmsa > $W/cmpi/sfetch.log 2>&1
echo "    run exit $?"
s1=$(zstd -dc $W/cmpi/sfetch/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
s2=$(zstd -dc $W/cmpi/pe100k.isal/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
echo "    SAM records as the ISA-L build's: $([ "$s1" = "$s2" ] && echo yes || echo NO)"
echo "SFETCHDONE $(date +%T)"
