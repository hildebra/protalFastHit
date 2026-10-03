#!/usr/bin/env bash
# The final tree (HEAD + final.patch: items 2-5) built in ~/mt-work/perf4/fin: tests; SAM text against the reference
# build (pe, ONT, 1 thread); --profile_only outputs against the EM build's (c.em.*: item 2 changes them, the
# differences are summarised with compare_profiles.sh); the four-sample map against the reference build's; ctest; e2e.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
SCR=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/85c3e3d5-6df2-4cf2-a5bc-63aafb61f818/scratchpad
W=$HOME/mt-work/perf4; E=$W/fin; DB=$HOME/bench071/V073/protal_db; P=$HOME/bench071/samples/points
rm -rf $E; mkdir -p $E
git -C $REPO archive HEAD | tar -x -C $E
(cd $E && patch -p1 -s < $SCR/final.patch) || { echo "FAIL patch"; exit 1; }
cd $E && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/fin.configure.log 2>&1 && \
  nice cmake --build build --target protal protal_tests simulate_metagenomes -j 5 > $W/fin.build.log 2>&1 && echo "OK fin build" || { echo "FAIL fin build"; grep -E 'error' -A3 $W/fin.build.log | head -60; exit 1; }
B=$E/build/protal
$E/build/tests/protal_tests --gtest_filter='SampleContext.*' 2>&1 | tail -2
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
ONT=$P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz
for spec in "pe -1 $R1 -2 $R2 --read_type pe" "ont -1 $ONT --read_type ont"; do
  set -- $spec; n=$1; shift; rm -rf $W/a.fin_${n}_t1
  $B --db $DB "$@" --prefix s -o $W/a.fin_${n}_t1 -t 1 --no_profile --no_qcmsa > $W/a.fin_${n}_t1.log 2>&1
  if cmp -s <(zstd -dc $W/a.ref_${n}_t1/s.sam.zst) <(zstd -dc $W/a.fin_${n}_t1/s.sam.zst); then echo "== fin $n t1: SAM text identical to the reference build's"; else echo "== fin $n t1: SAM TEXT DIFFERS"; fi
done
for spec in "pe500k_t1 1 $W/o.pe500k_t1/s.sam.zst" "pe500k_t6 6 $W/o.pe500k_t1/s.sam.zst" "pe5M_t6 6 $W/o.pe5M_t6/s.sam.zst" "ont90M_t6 6 $W/o.ont90M_t6/s.sam.zst" "pb90M_t6 6 $W/o.pb90M_t6/s.sam.zst"; do
  set -- $spec; n=$1; t=$2; sam=$3; rm -rf $W/c.fin.$n
  /usr/bin/time -f "%e s wall, %U s user" $B --db $DB --profile_only $sam --prefix s -o $W/c.fin.$n -t $t --no_qcmsa --verbose > $W/c.fin.$n.log 2> $W/c.fin.$n.time
  echo "== $n t=$t: $(tail -1 $W/c.fin.$n.time); $(grep -E '^Profiling took' $W/c.fin.$n.log) (EM build: $(grep -E '^Profiling took' $W/c.em.$n.log | sed 's/Profiling took //'))"
  if diff -r $W/c.em.$n $W/c.fin.$n > /dev/null; then echo "   outputs identical to the EM build's"; else echo "   outputs differ from the EM build's (item 2):"; bash $SCR/compare_profiles.sh $W/c.em.$n $W/c.fin.$n; fi
done
rm -rf $W/m.fin4; mkdir -p $W/m.fin4
{ printf '#OUTPUT_DIR\t%s\n#INPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\n' $W/m.fin4 $P
  for i in 1 2 3 4; do printf 's%d\trl150_p500000/sim/reads/rl150_p500000_s_%d_R1.fq.gz\trl150_p500000/sim/reads/rl150_p500000_s_%d_R2.fq.gz\ts%d\n' $i $i $i $i; done; } > $W/m.fin4.tsv
/usr/bin/time -f "%e s wall, %U s user" $B --db $DB --map $W/m.fin4.tsv --read_type pe -t 6 --no_qcmsa --verbose > $W/m.fin4.log 2> $W/m.fin4.time
echo "== map4 (final): $(tail -1 $W/m.fin4.time); $(grep -E '^(Processing all samples|Profiling|Strain-level MSAs|Run protal) took' $W/m.fin4.log | sed 's/ took / /' | tr '\n' ';') (reference run: 15.2 s)"
for i in 1 2 3 4; do echo "   sample s$i vs the reference build's map4:"; bash $SCR/compare_profiles.sh $W/o.map4 $W/m.fin4 s$i; done
if diff -r $W/o.map4/strains $W/m.fin4/strains > /dev/null; then echo "   strain MSAs identical to the reference build's"; else echo "   strain MSAs differ: $(diff -rq $W/o.map4/strains $W/m.fin4/strains | wc -l) files"; fi
ctest --test-dir build --output-on-failure -j 4 > $W/fin.ctest.log 2>&1; echo "ctest: $(grep -E 'tests passed|tests failed' $W/fin.ctest.log | tail -1)"
PROTAL_TEST_DB=$HOME/mt-work/head07c/mini_db/protal_db PROTAL=$B SIMULATE=$E/build/simulate_metagenomes \
  python3 -m unittest tests/e2e/test_protal_e2e.py > $W/fin.e2e.log 2>&1; echo "e2e (mini database): $(grep -E '^Ran |^OK|^FAILED' $W/fin.e2e.log | tr '\n' ' ')"
grep -E '^(FAIL|ERROR):' $W/fin.e2e.log | head -5
