#!/usr/bin/env bash
# The output-identical changes (id.patch on HEAD): unit tests (the two new ones shown), then outputs against the
# reference: 100k pairs -t 1 (SAM text), 500k pairs -t 6 (sorted SAM, profile from the same SAM), single-end, the
# long reads (SAM text, -t 1); callgrind of the 100k pairs.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/build_wt.sh perf3-id $S/perf3/id.patch $(cat $HOME/mt-work/perf3/ref/COMMIT) "protal protal_tests" || exit 1
B=$HOME/mt-work/perf3-id/src/build
(cd $B && ./tests/protal_tests --gtest_filter='AnchoredAlignment.UngappedFlanks*:IndexCodec.LongRuns*' 2>&1 | grep -E "aligned,|OK|FAILED|PASSED|Failure" | head -10)
(cd $B && nice ctest -j 4 > $HOME/mt-work/perf3-id/ctest.log 2>&1; echo "ctest: $(grep -E 'tests passed|tests failed' $HOME/mt-work/perf3-id/ctest.log)")
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points; I=$HOME/mt-work/isacg
REF=$W/ref/build/protal; ID=$B/protal; R5=$P/rl150_p500000/sim/reads/rl150_p500000_s_1; S5=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz
run() { local name=$1; shift; for b in ref id; do bin=$REF; [ $b = id ] && bin=$ID; rm -rf $W/i.$name.$b; nice $bin --db $DB "$@" --prefix s -o $W/i.$name.$b --no_qcmsa > $W/i.$name.$b.log 2>&1 || echo "FAIL $name $b"; done; }
run pe100k -1 $I/r1.fq -2 $I/r2.fq --no_profile -t 1
cmp -s <(zstdcat $W/i.pe100k.ref/s.sam.zst) <(zstdcat $W/i.pe100k.id/s.sam.zst) && echo "100k pairs -t 1: SAM text same" || echo "100k pairs -t 1: SAM DIFFER"
run pe500k -1 ${R5}_R1.fq.gz -2 ${R5}_R2.fq.gz -t 6
cmp -s <(zstdcat $W/i.pe500k.ref/s.sam.zst | sort) <(zstdcat $W/i.pe500k.id/s.sam.zst | sort) && echo "500k pairs -t 6: sorted SAM same" || echo "500k pairs -t 6: SAM DIFFER"
rm -rf $W/i.pe500k.reprof; nice $REF --db $DB --profile_only $W/i.pe500k.id/s.sam.zst --prefix s -o $W/i.pe500k.reprof -t 6 --no_qcmsa > /dev/null 2>&1
diff -r -q -x '*_runtime.tsv' -x '*.sam.zst' -x '*.log' $W/i.pe500k.id $W/i.pe500k.reprof > /dev/null && echo "500k pairs: its SAM profiled by the reference gives its profile" || echo "500k pairs: profile DIFFERS from the reference's on its SAM"
run se500k -1 $S5 --read_type se --no_profile -t 1
cmp -s <(zstdcat $W/i.se500k.ref/s.sam.zst) <(zstdcat $W/i.se500k.id/s.sam.zst) && echo "single-end 500k -t 1: SAM text same" || echo "single-end: SAM DIFFER"
for s in ont_b3000000 pb_b3000000; do t=ont; [ $s = pb_b3000000 ] && t=pb
  run $s -1 $P/$s/sim/reads/${s}_s_1.fq.gz --read_type $t --no_profile -t 1
  cmp -s <(zstdcat $W/i.$s.ref/s.sam.zst) <(zstdcat $W/i.$s.id/s.sam.zst) && echo "$s -t 1: SAM text same" || echo "$s: SAM DIFFER"
done
run hifi3M -1 $W/hifi/hifi_3M.fq.gz --read_type pb --no_profile -t 1
cmp -s <(zstdcat $W/i.hifi3M.ref/s.sam.zst) <(zstdcat $W/i.hifi3M.id/s.sam.zst) && echo "HiFi 3 Mb -t 1: SAM text same" || echo "HiFi: SAM DIFFER"
nice valgrind --tool=callgrind --callgrind-out-file=$W/cg.pe100k.id $ID --db $DB -1 $I/r1.fq -2 $I/r2.fq --no_profile --prefix s -o $W/cg.o.pe100k.id -t 1 --no_qcmsa > /dev/null 2>&1
callgrind_annotate --inclusive=yes --threshold=99.5 $W/cg.pe100k.id 2>/dev/null | c++filt > $W/cg.pe100k.id.incl.txt
echo "== pe100k.id: $(grep 'PROGRAM TOTALS' $W/cg.pe100k.id.incl.txt | awk '{print $1}')"
grep -E 'RunPairedEnd|AnchoredAligner::(Align|Flank)|alignEndsFree|DecodeChunk|Seedmap::Load ' $W/cg.pe100k.id.incl.txt | sed -E 's/\[\/home[^]]*\]//; s/\(([^()]|\([^()]*\))*\)//g' | awk '{ $1 = $1; print }' | cut -c1-110 | awk '!seen[$2]++' | head -8
echo ID DONE
