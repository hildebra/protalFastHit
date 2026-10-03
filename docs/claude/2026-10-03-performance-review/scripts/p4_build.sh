#!/usr/bin/env bash
# Item 4 (profiling ahead) on top of 3 and 5: HEAD + p4.patch built in ~/mt-work/perf4/p4; outputs of map runs
# (4 and 8 samples) and single-sample runs compared between the default and --profile_after_alignment, and with the
# reference build's map4 outputs; times; callgrind of --profile_only (sketches) on the 500k SAM.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
SCR=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/85c3e3d5-6df2-4cf2-a5bc-63aafb61f818/scratchpad
W=$HOME/mt-work/perf4; E=$W/p4; DB=$HOME/bench071/V073/protal_db; P=$HOME/bench071/samples/points
rm -rf $E; mkdir -p $E
git -C $REPO archive HEAD | tar -x -C $E
(cd $E && patch -p1 -s < $SCR/p4.patch) || { echo "FAIL patch"; exit 1; }
cd $E && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/p4.configure.log 2>&1 && \
  nice cmake --build build --target protal protal_tests simulate_metagenomes -j 5 > $W/p4.build.log 2>&1 && echo "OK p4 build" || { echo "FAIL p4 build"; grep -E 'error' -A3 $W/p4.build.log | head -60; exit 1; }
B=$E/build/protal
# callgrind: --profile_only on the 500k SAM (the sketches of item 3)
rm -rf $W/cgp.p4
valgrind --tool=callgrind --callgrind-out-file=$W/cg.p4.prof.out --compress-strings=no --compress-pos=no \
  $B --db $DB --profile_only $W/o.pe500k_t1/s.sam.zst --prefix s -o $W/cgp.p4 -t 1 --no_qcmsa > $W/cg.p4.prof.log 2>&1
callgrind_annotate --inclusive=yes $W/cg.p4.prof.out 2>/dev/null | head -80 > $W/cg.p4.prof.incl.txt
echo "--- callgrind --profile_only 500k (em build | p4):"
for f in "PROGRAM TOTALS" "ApplySampleContext" "AbundanceWeightedShares" "CongenerDistances::Compared" "GeneSketch" "MakeSketch"; do
  printf "   %-40s %s | %s\n" "$f" "$(grep -m1 -F "$f" $W/cg.em.incl.txt | awk '{print $1}')" "$(grep -m1 -F "$f" $W/cg.p4.prof.incl.txt | awk '{print $1}')"; done
# maps
mk_map() { # out n
  local out=$1 n=$2
  { printf '#OUTPUT_DIR\t%s\n#INPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\n' $out $P
    for i in 1 2 3 4; do printf 'a%d\trl150_p500000/sim/reads/rl150_p500000_s_%d_R1.fq.gz\trl150_p500000/sim/reads/rl150_p500000_s_%d_R2.fq.gz\ta%d\n' $i $i $i $i; done
    if [ $n = 8 ]; then for i in 1 2 3 4; do printf 'b%d\trl100_p500000/sim/reads/rl100_p500000_s_%d_R1.fq.gz\trl100_p500000/sim/reads/rl100_p500000_s_%d_R2.fq.gz\tb%d\n' $i $i $i $i; done; fi; }
}
run_map() { # name n args...
  local n=$1 k=$2; shift 2; rm -rf $W/m.$n; mkdir -p $W/m.$n; mk_map $W/m.$n $k > $W/m.$n.tsv
  /usr/bin/time -f "%e s wall, %U s user" $B --db $DB --map $W/m.$n.tsv --read_type pe -t 6 --no_qcmsa --verbose "$@" > $W/m.$n.log 2> $W/m.$n.time
  echo "   $n: $(tail -1 $W/m.$n.time); $(grep -E '^(Processing all samples|Profiling|Strain-level MSAs|Run protal) took' $W/m.$n.log | sed 's/ took / /' | tr '\n' ';')"
  grep -E '^Profiled while' $W/m.$n.log
}
cmp_dirs() { # a b label
  if diff -r -x '*_runtime.tsv' -x '*.log' -x '*.err' $1 $2 > $W/m.$3.diff; then echo "== $3: outputs identical ($(find $2 -type f | wc -l) files)"; else echo "== $3: OUTPUTS DIFFER"; head -8 $W/m.$3.diff; fi
}
for round in 1 2; do
  run_map map4_ahead_$round 4
  run_map map4_after_$round 4 --profile_after_alignment
done
cmp_dirs $W/m.map4_after_1 $W/m.map4_ahead_1 map4_ahead_vs_after
cmp_dirs $W/m.map4_after_1 $W/m.map4_ahead_2 map4_ahead2_vs_after
# against the reference build's four-sample map (profiles and strains; the sample names differ: s1.. vs a1..)
for i in 1 2 3 4; do
  for f in profile profile.log profile.gene.log profile.genes.log; do
    cmp -s $W/o.map4/profiles/s$i.$f $W/m.map4_ahead_1/profiles/a$i.$f 2>/dev/null || cmp -s $W/o.map4/s$i.$f $W/m.map4_ahead_1/a$i.$f 2>/dev/null || echo "   differs or missing: sample $i $f"
  done
done
echo "   (compared the reference map4's profiles with the ahead run's, sample by sample; any difference is listed above)"
ls $W/o.map4 | head -3; ls $W/m.map4_ahead_1 | head -3
run_map map8_ahead 8
run_map map8_after 8 --profile_after_alignment
cmp_dirs $W/m.map8_after $W/m.map8_ahead map8_ahead_vs_after
# single sample: the only sample is not profiled ahead (the stage takes it on all threads)
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
for v in ahead after; do rm -rf $W/s1.$v; args=""; [ $v = after ] && args="--profile_after_alignment"
  /usr/bin/time -f "%e s wall" $B --db $DB -1 $R1 -2 $R2 --read_type pe --prefix s -o $W/s1.$v -t 6 --no_qcmsa --verbose $args > $W/s1.$v.log 2> $W/s1.$v.time
  echo "   single $v: $(tail -1 $W/s1.$v.time); $(grep -E '^(Aligning reads|Profiling|Run protal) took' $W/s1.$v.log | sed 's/ took / /' | tr '\n' ';')"; grep -E '^Profiled while' $W/s1.$v.log; done
cmp_dirs $W/s1.after $W/s1.ahead single_ahead_vs_after
ctest --test-dir build --output-on-failure -j 4 > $W/p4.ctest.log 2>&1; echo "ctest: $(grep -E 'tests passed|tests failed' $W/p4.ctest.log | tail -1)"
# SAM text of the final alignment code (items 5a, 5b, 5d kept; 5c and 5e reverted) against the reference build, 1 thread.
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
ONT=$P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz
for spec in "pe -1 $R1 -2 $R2 --read_type pe" "ont -1 $ONT --read_type ont"; do
  set -- $spec; n=$1; shift; rm -rf $W/a.p4_${n}_t1
  /usr/bin/time -f "%e s wall" $B --db $DB "$@" --prefix s -o $W/a.p4_${n}_t1 -t 1 --no_profile --no_qcmsa --verbose > $W/a.p4_${n}_t1.log 2> $W/a.p4_${n}_t1.time
  if cmp -s <(zstd -dc $W/a.ref_${n}_t1/s.sam.zst) <(zstd -dc $W/a.p4_${n}_t1/s.sam.zst); then echo "== p4 $n t1: SAM text identical to the reference build's ($(tail -1 $W/a.p4_${n}_t1.time); $(grep -E '^Aligning reads took' $W/a.p4_${n}_t1.log))"; else echo "== p4 $n t1: SAM TEXT DIFFERS"; fi
done
rm -rf $W/cga.p4
valgrind --tool=callgrind --callgrind-out-file=$W/cg.p4.align.out --compress-strings=no --compress-pos=no \
  $B --db $DB -1 $W/reads/pe100k_R1.fq.gz -2 $W/reads/pe100k_R2.fq.gz --read_type pe --no_profile --prefix s -o $W/cga.p4 -t 1 --no_qcmsa > $W/cg.p4.align.log 2>&1
callgrind_annotate --inclusive=yes $W/cg.p4.align.out 2>/dev/null | head -150 > $W/cg.p4.align.incl.txt
echo "--- callgrind aligning 100k pairs (ref | p4 final):"
for f in "PROGRAM TOTALS" "RunPairedEnd" "SimpleAlignmentHandler::operator" "ChainAnchorFinder<protal::KmerLookupSM>::operator" "IsAlignmentValid" "gene_neighbours::Table::Assess" "AlignmentEdits" "NextCompressedCigar"; do
  printf "   %-50s %s | %s\n" "$f" "$(grep -m1 -F "$f" $W/cg.align.incl.txt | awk '{print $1}')" "$(grep -m1 -F "$f" $W/cg.p4.align.incl.txt | awk '{print $1}')"; done
