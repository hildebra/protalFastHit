#!/bin/bash
# Experiment 9: the ASan/UBSan build on (a) long reads given as single-end (e4: 20, 70 and 140 kb
# reads through the short-read path), (b) the 200 kb HiFi and ONT reads of e2 (chunked), and
# (c) the long-read unit tests.
S=$(dirname "$0")
W=~/audit6/longreads
PA=$W/build-asan/protal
D=$W/mini_db
OUT=$S/out_e9.txt
export ASAN_OPTIONS=detect_leaks=0
export UBSAN_OPTIONS=print_stacktrace=1
mkdir -p $W/e9 && cd $W/e9
{
rm -rf out
$PA --db $D -1 $W/e4/mixed_se.fq --model_se $D/model_pe.xml -o out -t 2 --sam_format sam --no_qcmsa --no_strains --prefix mixed_se > mixed_se.log 2>&1
echo "(a) se with long reads rc=$?"
grep -E "runtime error|AddressSanitizer|SUMMARY" mixed_se.log | sort | uniq -c | head -10
grep -A12 "runtime error" mixed_se.log | grep -E "^\s+#[0-9]" | head -12
$PA --db $D -1 $W/e2/h200000.fq,$W/e2/o200000.fq --prefix h200000,o200000 --read_type pb --model_pb $D/model_pe.xml -o out -t 2 --sam_format sam --no_qcmsa --no_strains > long.log 2>&1
echo "(b) 200 kb reads (as pb) rc=$?"
grep -E "runtime error|AddressSanitizer|SUMMARY" long.log | sort | uniq -c | head -10
$W/build-asan/tests/protal_tests --gtest_filter='ChunkRead.*:LongRead*:ReadConsensus.*:SeqReader.*:Options.*:SampleMap.*' > tests.log 2>&1
echo "(c) unit tests rc=$?"
grep -E "^\[  (PASSED|FAILED)|tests ran|FAILED" tests.log | head
grep -E "runtime error|AddressSanitizer" tests.log | sort | uniq -c | head
} > $OUT 2>&1
echo done
