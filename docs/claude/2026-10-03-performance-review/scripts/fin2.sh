#!/usr/bin/env bash
# The fin tree with the corrected test (and the comment in SampleContext.h): tests again; then alternated timings
# on a quiet machine: --profile_only of the 5M SAM (EM build vs final), the four-sample map (reference vs final),
# and a single 500k-pair run at one thread (reference vs final).
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/mt-work/perf4; E=$W/fin; DB=$HOME/bench071/V073/protal_db; P=$HOME/bench071/samples/points
cp $REPO/tests/test_SampleContext.cpp $E/tests/; cp $REPO/src/Profiling/SampleContext.h $E/src/Profiling/
cd $E && nice cmake --build build --target protal protal_tests -j 5 > $W/fin2.build.log 2>&1 && echo "OK fin2 build" || { echo "FAIL fin2 build"; grep -E 'error' -A3 $W/fin2.build.log | head -40; exit 1; }
$E/build/tests/protal_tests --gtest_filter='SampleContext.*' 2>&1 | tail -2
ctest --test-dir build -j 4 > $W/fin2.ctest.log 2>&1; echo "ctest: $(grep -E 'tests passed|tests failed' $W/fin2.ctest.log | tail -1)"
B=$E/build/protal; EM=$W/em/build/protal; R=$W/ref/build/protal
echo "--- alternated: --profile_only 5M SAM, 6 threads (EM build | final), 3 rounds"
for r in 1 2 3; do for v in em fin; do bin=$EM; [ $v = fin ] && bin=$B; rm -rf $W/t.$v.5M
  /usr/bin/time -f "%e" $bin --db $DB --profile_only $W/o.pe5M_t6/s.sam.zst --prefix s -o $W/t.$v.5M -t 6 --no_qcmsa > $W/t.$v.5M.log 2> $W/t.$v.5M.time
  printf "   round %d %s: %s s wall; %s\n" $r $v "$(tail -1 $W/t.$v.5M.time)" "$(grep -E '^Profiling took' $W/t.$v.5M.log)"; done; done
echo "--- alternated: four-sample map, 6 threads (reference 27423c6 | final), 2 rounds"
for r in 1 2; do for v in ref fin; do bin=$R; [ $v = fin ] && bin=$B; rm -rf $W/t.$v.map4; mkdir -p $W/t.$v.map4
  { printf '#OUTPUT_DIR\t%s\n#INPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\n' $W/t.$v.map4 $P
    for i in 1 2 3 4; do printf 's%d\trl150_p500000/sim/reads/rl150_p500000_s_%d_R1.fq.gz\trl150_p500000/sim/reads/rl150_p500000_s_%d_R2.fq.gz\ts%d\n' $i $i $i $i; done; } > $W/t.$v.map4.tsv
  /usr/bin/time -f "%e" $bin --db $DB --map $W/t.$v.map4.tsv --read_type pe -t 6 --no_qcmsa > $W/t.$v.map4.log 2> $W/t.$v.map4.time
  printf "   round %d %s: %s s wall; %s\n" $r $v "$(tail -1 $W/t.$v.map4.time)" "$(grep -E '^(Processing all samples|Profiling|Strain-level MSAs) took' $W/t.$v.map4.log | sed 's/ took / /' | tr '\n' ';')"; done; done
echo "--- alternated: 500k pairs, whole run, 1 thread (reference | final), 2 rounds"
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
for r in 1 2; do for v in ref fin; do bin=$R; [ $v = fin ] && bin=$B; rm -rf $W/t.$v.pe1
  /usr/bin/time -f "%e" $bin --db $DB -1 $R1 -2 $R2 --read_type pe --prefix s -o $W/t.$v.pe1 -t 1 --no_qcmsa > $W/t.$v.pe1.log 2> $W/t.$v.pe1.time
  printf "   round %d %s: %s s wall; %s\n" $r $v "$(tail -1 $W/t.$v.pe1.time)" "$(grep -E '^(Aligning reads|Profiling) took' $W/t.$v.pe1.log | sed 's/ took / /' | tr '\n' ';')"; done; done
