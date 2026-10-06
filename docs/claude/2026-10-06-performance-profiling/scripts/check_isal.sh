#!/usr/bin/env bash
# check_isal.sh: the ISA-L build (~/perf8/isal: items 1 and 3 plus ISA-L for all gzip) against the work build (items 1
# and 3, zlib-ng and libdeflate): unit tests; whole runs on single-member gzip (pe100k, the PacBio sample) and BGZF
# input (the simulator's), SAM records and every other file identical; a .sam.gz written with ISA-L, read by zcat and
# profiled as the .sam.zst; the end-to-end tests; then the build paths: ISA-L and nasm built by lib/static-deps.cmake
# from local tarballs (no nasm on this PATH, no system ISA-L: protal links the fetched one) and a fully static protal
# against ~/isal/prefix's libisal.a and the system's libzstd.a.
set -uo pipefail
W=$HOME/perf8; K=$W/work/build; I=$W/isal/build
DB=$HOME/bench071/V075/protal_db; P=$HOME/bench071/samples/points; RD=$HOME/perf6/reads
mkdir -p $W/cmpi
echo "== unit tests, ISA-L build ($(date +%T))"
$I/tests/protal_tests --gtest_brief=1 > $W/unit.isal.log 2>&1; echo "exit $?"
grep -E "^\[  (PASSED|FAILED|SKIPPED)|tests ran|FAILED  \]" $W/unit.isal.log | head -30
$I/tests/protal_tests --gtest_filter='ThreadedGzStream.*:SamFile.*:Compressor.*:Bgzf*' > $W/unit.isal.gz.log 2>&1
grep -E "^\[       OK|FAILED" $W/unit.isal.gz.log

echo "== whole runs, 6 threads ($(date +%T))"
run() { local b=$1 o=$2; shift 2; rm -rf $o; $b --db $DB "$@" --prefix s -o $o -t 6 --no_qcmsa --verbose > $o.log 2>&1; echo "    $(basename $o): exit $?"; }
for set in pe100k pe se pb; do
  case $set in
    pe100k) args="-1 $RD/pe100k_R1.fq.gz -2 $RD/pe100k_R2.fq.gz --read_type pe" ;;   # single-member gzip
    pe) args="-1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz -2 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz --read_type pe" ;;  # BGZF
    se) args="-1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz --read_type se" ;;
    pb) args="-1 $P/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz --read_type pb" ;;  # single-member gzip with a name
  esac
  echo "-- $set"
  run $K/protal $W/cmpi/$set.work $args
  run $I/protal $W/cmpi/$set.isal $args
  a=$W/cmpi/$set.work; b=$W/cmpi/$set.isal
  sa=$(zstd -dc $a/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
  sb=$(zstd -dc $b/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
  n=$(zstd -dc $a/s.sam.zst | grep -vc '^@')
  diffs=$(cd $a && find . -type f ! -name '*.sam.zst' ! -name '*.err' ! -name '*_runtime.tsv' | while read -r f; do cmp -s "$f" "$b/$f" || echo "$f"; done | tr '\n' ' ')
  echo "    work vs isal: SAM records $([ "$sa" = "$sb" ] && echo "identical ($n)" || echo "DIFFER ($sa vs $sb)"); other files differing: ${diffs:-none}"
  grep -E "reads? (pairs )?read|Sequence reader took" $b.log | head -2 | cut -c1-120 | sed 's/^/    /'
done

echo "== .sam.gz written with ISA-L ($(date +%T))"
run $I/protal $W/cmpi/gz.isal -1 $RD/pe100k_R1.fq.gz -2 $RD/pe100k_R2.fq.gz --read_type pe --sam_format gz
g=$W/cmpi/gz.isal/s.sam.gz; z=$W/cmpi/pe100k.isal/s.sam.zst
echo "    $(ls -l $g | awk '{print $5}') bytes (the .sam.zst: $(ls -l $z | awk '{print $5}')); gzip -t: $(gzip -t $g && echo ok)"
sg=$(zcat $g | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12); sz=$(zstd -dc $z | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
echo "    records as the .sam.zst's: $([ "$sg" = "$sz" ] && echo yes || echo "NO ($sg vs $sz)")"
diffs=$(cd $W/cmpi/pe100k.isal && find . -type f ! -name '*.sam.*' ! -name '*.err' ! -name '*_runtime.tsv' | while read -r f; do cmp -s "$f" "$W/cmpi/gz.isal/$f" || echo "$f"; done | tr '\n' ' ')
echo "    profile files differing from the .sam.zst run: ${diffs:-none}"

echo "== end-to-end tests on the mini database ($(date +%T))"
cd $W/isal && PROTAL_TEST_DB=$HOME/buildmem/new/data/mini_db/protal_db PROTAL=$I/protal SIMULATE=$I/simulate_metagenomes \
  timeout 3000 python3 -m unittest tests/e2e/test_protal_e2e.py > $W/e2e.isal.log 2>&1; echo "exit $?"
tail -4 $W/e2e.isal.log

echo "== lib/static-deps.cmake: nasm and ISA-L built from the local tarballs, protal linked to that ISA-L ($(date +%T))"
S=$W/isal/build-sfetch; rm -rf $S
( cd $W/isal && PATH=$(echo "$PATH" | tr ':' '\n' | grep -v isal | paste -sd:) nice cmake -S . -B $S -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DPROTAL_STATIC_FETCH_DEPS=ON -DPROTAL_ISAL_URL=$HOME/isal/isa-l-2.32.1.tar.gz -DPROTAL_NASM_URL=$HOME/isal/nasm-2.16.03.tar.xz \
    > $W/sfetch.configure.log 2>&1 ); echo "configure exit $?"
grep -E "nasm for ISA-L|ISA-L:|Static binaries" $W/sfetch.configure.log | sed 's/^/    /'
nice cmake --build $S --target protal -j 5 > $W/sfetch.build.log 2>&1; echo "build exit $?"
ls -l $S/static-deps/nasm-install/bin/nasm $S/static-deps/isal-install/lib/libisal.a 2>&1 | awk '{print "    " $5 " " $NF}'
$S/protal --version 2>&1 | head -1 | sed 's/^/    /'
run $S/protal $W/cmpi/sfetch -1 $RD/pe100k_R1.fq.gz -2 $RD/pe100k_R2.fq.gz --read_type pe
s1=$(zstd -dc $W/cmpi/sfetch/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
s2=$(zstd -dc $W/cmpi/pe100k.isal/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
echo "    SAM records as the ISA-L build's: $([ "$s1" = "$s2" ] && echo yes || echo NO)"

echo "== a fully static protal: ~/isal/prefix's libisal.a, the system's libzstd.a ($(date +%T))"
T=$W/isal/build-static; rm -rf $T
( cd $W/isal && nice cmake -S . -B $T -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$HOME/isal/prefix > $W/static.configure.log 2>&1 ); echo "configure exit $?"
nice cmake --build $T --target protal_static -j 5 > $W/static.build.log 2>&1; echo "build exit $?"
sb=$(ls $T/protal_*_static | head -1)
echo "    $(basename $sb): $(file -b $sb | cut -c1-80)"
run $sb $W/cmpi/static -1 $RD/pe100k_R1.fq.gz -2 $RD/pe100k_R2.fq.gz --read_type pe
s3=$(zstd -dc $W/cmpi/static/s.sam.zst | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-12)
echo "    SAM records as the ISA-L build's: $([ "$s3" = "$s2" ] && echo yes || echo NO)"
echo "CHECKISALDONE $(date +%T)"
