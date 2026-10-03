#!/usr/bin/env bash
# Items 3 and 5: HEAD (4cae2f4) + p5.patch built in ~/mt-work/perf4/p5; SampleContext tests; SAM text compared with
# the reference build's at 1 thread (pe, se, ONT, PacBio) and sorted at 6 threads (pe); --profile_only outputs compared
# with the EM build's (c.em.*); stage timers; callgrind of aligning 100k pairs.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
SCR=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/85c3e3d5-6df2-4cf2-a5bc-63aafb61f818/scratchpad
W=$HOME/mt-work/perf4; E=$W/p5; R=$W/ref/build/protal; DB=$HOME/bench071/V073/protal_db
P=$HOME/bench071/samples/points
rm -rf $E; mkdir -p $E
git -C $REPO archive HEAD | tar -x -C $E
(cd $E && patch -p1 -s < $SCR/p5.patch) || { echo "FAIL patch"; exit 1; }
cd $E && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/p5.configure.log 2>&1 && \
  nice cmake --build build --target protal protal_tests -j 5 > $W/p5.build.log 2>&1 && echo "OK p5 build" || { echo "FAIL p5 build"; grep -E 'error' -A3 $W/p5.build.log | head -60; exit 1; }
B=$E/build/protal
$E/build/tests/protal_tests --gtest_filter='SampleContext.*:GeneNeighbours.*:AnchoredAlignment.*:PairJoin.*' 2>&1 | tail -3
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
ONT=$P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz; PB=$P/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
align() { # build name threads args...
  local bin=$1 n=$2 t=$3; shift 3; rm -rf $W/a.$n
  /usr/bin/time -f "%e s wall, %U s user" $bin --db $DB "$@" --prefix s -o $W/a.$n -t $t --no_profile --no_qcmsa --verbose > $W/a.$n.log 2> $W/a.$n.time
  echo "   $n: $(tail -1 $W/a.$n.time); $(grep -E '^Aligning reads took' $W/a.$n.log)"
}
samcmp() { # name_ref name_new sorted?
  local a=$W/a.$1/s.sam.zst b=$W/a.$2/s.sam.zst
  if [ "${3:-}" = sorted ]; then
    if cmp -s <(zstd -dc $a | grep -v '^@' | sort) <(zstd -dc $b | grep -v '^@' | sort) && cmp -s <(zstd -dc $a | grep '^@' | sort) <(zstd -dc $b | grep '^@' | sort); then echo "== $2: sorted SAM identical to $1"; else echo "== $2: SORTED SAM DIFFERS from $1"; fi
  else
    if cmp -s <(zstd -dc $a) <(zstd -dc $b); then echo "== $2: SAM text identical to $1 ($(zstd -dc $b | grep -vc '^@') records)"; else echo "== $2: SAM TEXT DIFFERS from $1"; diff <(zstd -dc $a) <(zstd -dc $b) | head -6; fi
  fi
}
align $R ref_pe_t1 1 -1 $R1 -2 $R2 --read_type pe
align $B p5_pe_t1 1 -1 $R1 -2 $R2 --read_type pe
samcmp ref_pe_t1 p5_pe_t1
align $R ref_se_t1 1 -1 $R1 --read_type se
align $B p5_se_t1 1 -1 $R1 --read_type se
samcmp ref_se_t1 p5_se_t1
align $R ref_ont_t1 1 -1 $ONT --read_type ont
align $B p5_ont_t1 1 -1 $ONT --read_type ont
samcmp ref_ont_t1 p5_ont_t1
align $R ref_pb_t1 1 -1 $PB --read_type pb
align $B p5_pb_t1 1 -1 $PB --read_type pb
samcmp ref_pb_t1 p5_pb_t1
align $R ref_pe_t6 6 -1 $R1 -2 $R2 --read_type pe
align $B p5_pe_t6 6 -1 $R1 -2 $R2 --read_type pe
samcmp ref_pe_t6 p5_pe_t6 sorted
echo "--- stage timers, 500k pairs, 1 thread (ref | p5):"
for s in "Seeding" "Sorting Seeds" "Alignment handler" "Output handler" "Joining alignment pairs and sorting" "Aligning reads"; do
  printf "   %-36s %s | %s\n" "$s" "$(grep -m1 "^\s*$s took" $W/a.ref_pe_t1.log | sed 's/.*took //')" "$(grep -m1 "^\s*$s took" $W/a.p5_pe_t1.log | sed 's/.*took //')"; done
prof() { # name threads sam
  local n=$1 t=$2 sam=$3; rm -rf $W/c.p5.$n
  /usr/bin/time -f "%e s wall, %U s user" $B --db $DB --profile_only $sam --prefix s -o $W/c.p5.$n -t $t --no_qcmsa --verbose > $W/c.p5.$n.log 2> $W/c.p5.$n.time
  echo "   p5 $n t=$t: $(tail -1 $W/c.p5.$n.time); $(grep -E '^Profiling took' $W/c.p5.$n.log) (em build: $(grep -E '^Profiling took' $W/c.em.$n.log))"
  if diff -r $W/c.em.$n $W/c.p5.$n > $W/c.p5.$n.diff; then echo "== $n: profile outputs identical to the EM build's"; else echo "== $n: PROFILE OUTPUTS DIFFER"; head -10 $W/c.p5.$n.diff; fi
}
prof pe500k_t1 1 $W/o.pe500k_t1/s.sam.zst
prof pe500k_t6 6 $W/o.pe500k_t1/s.sam.zst
prof pe5M_t6 6 $W/o.pe5M_t6/s.sam.zst
prof ont90M_t6 6 $W/o.ont90M_t6/s.sam.zst
prof pb90M_t6 6 $W/o.pb90M_t6/s.sam.zst
rm -rf $W/cga.p5
valgrind --tool=callgrind --callgrind-out-file=$W/cg.p5.align.out --compress-strings=no --compress-pos=no \
  $B --db $DB -1 $W/reads/pe100k_R1.fq.gz -2 $W/reads/pe100k_R2.fq.gz --read_type pe --no_profile --prefix s -o $W/cga.p5 -t 1 --no_qcmsa > $W/cg.p5.align.log 2>&1
callgrind_annotate --inclusive=yes $W/cg.p5.align.out 2>/dev/null | head -150 > $W/cg.p5.align.incl.txt
callgrind_annotate $W/cg.p5.align.out 2>/dev/null | head -150 > $W/cg.p5.align.self.txt
echo "--- callgrind aligning 100k pairs (ref | p5):"
for f in "PROGRAM TOTALS" "RunPairedEnd" "SimpleAlignmentHandler::operator" "AnchoredAligner::Align" "wavefront_align" "ChainAnchorFinder<protal::KmerLookupSM>::operator" "ProtalPairedOutputHandler" "IsAlignmentValid" "gene_neighbours::Table::Assess" "wavefront_slab_reap_repurpose" "wavefront_slab_allocate" "__introsort_loop" "AlignmentEdits" "GeneSketch"; do
  printf "   %-50s %s | %s\n" "$f" "$(grep -m1 -F "$f" $W/cg.align.incl.txt | awk '{print $1}')" "$(grep -m1 -F "$f" $W/cg.p5.align.incl.txt | awk '{print $1}')"; done
