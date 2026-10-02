#!/usr/bin/env bash
# The combined prototype (long reads through their chain with re-seeding and a piece aligner; the ungapped one-mismatch
# flank; the empty-block fill of the index decoder): short reads byte-identical to the reference (100k pairs at -t 1,
# 500k pairs at -t 6 sorted); long-read records against the reference (Nanopore, PacBio, HiFi 3 Mb); callgrind.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/perf3/build_exp.sh lr || exit 1
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points; I=$HOME/mt-work/isacg
REF=$W/ref/build/protal; LR=$HOME/mt-work/perf3-lr/src/build/protal
R5=$P/rl150_p500000/sim/reads/rl150_p500000_s_1
for b in ref lr; do bin=$REF; [ $b = lr ] && bin=$LR
  rm -rf $W/pe.$b $W/pe5.$b
  nice $bin --db $DB -1 $I/r1.fq -2 $I/r2.fq --no_profile --prefix s -o $W/pe.$b -t 1 --no_qcmsa > $W/pe.$b.log 2>&1
  nice $bin --db $DB -1 ${R5}_R1.fq.gz -2 ${R5}_R2.fq.gz --prefix s -o $W/pe5.$b -t 6 --no_qcmsa > $W/pe5.$b.log 2>&1
done
cmp -s <(zstdcat $W/pe.ref/s.sam.zst) <(zstdcat $W/pe.lr/s.sam.zst) && echo "100k pairs -t 1: SAM text same" || echo "100k pairs -t 1: SAM DIFFER"
cmp -s <(zstdcat $W/pe5.ref/s.sam.zst | sort) <(zstdcat $W/pe5.lr/s.sam.zst | sort) && echo "500k pairs -t 6: sorted SAM same" || echo "500k pairs -t 6: SAM DIFFER"
diff -r -q -x '*_runtime.tsv' -x '*.sam.zst' $W/pe5.ref $W/pe5.lr > /dev/null && echo "500k pairs -t 6: profile and other outputs same" || echo "500k pairs -t 6: outputs DIFFER"
declare -A RD=( [ont_b3000000]="$P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz ont" [pb_b3000000]="$P/pb_b3000000/sim/reads/pb_b3000000_s_1.fq.gz pb" [hifi3M]="$W/hifi/hifi_3M.fq.gz pb" )
for s in ont_b3000000 pb_b3000000 hifi3M; do set -- ${RD[$s]}
  for b in ref lr4; do bin=$REF; [ $b = lr4 ] && bin=$LR
    rm -rf $W/o.$s.$b; nice $bin --db $DB -1 $1 --read_type $2 --no_profile --prefix s -o $W/o.$s.$b -t 1 --no_qcmsa > /dev/null 2>&1
  done
  echo "== $s"; bash $S/perf3/lrcmp.sh $W/o.$s.ref/s.sam.zst $W/o.$s.lr4/s.sam.zst
done
cg() { local name=$1 bin=$2; shift 2; nice valgrind --tool=callgrind --callgrind-out-file=$W/cg.$name $bin --db $DB "$@" --prefix s -o $W/cg.o.$name -t 1 --no_qcmsa > /dev/null 2>&1; }
cg pe100k.lr4 $LR -1 $I/r1.fq -2 $I/r2.fq --no_profile &
cg ont_b3000000.lr4 $LR -1 $P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz --read_type ont --no_profile &
cg pb_b3000000.lr4 $LR -1 $P/pb_b3000000/sim/reads/pb_b3000000_s_1.fq.gz --read_type pb --no_profile &
cg hifi3M.lr4 $LR -1 $W/hifi/hifi_3M.fq.gz --read_type pb --no_profile &
cg hifi3M.ref $REF -1 $W/hifi/hifi_3M.fq.gz --read_type pb --no_profile &
wait
for n in pe100k.lr4 ont_b3000000.ref ont_b3000000.lr4 pb_b3000000.ref pb_b3000000.lr4 hifi3M.ref hifi3M.lr4; do
  callgrind_annotate --inclusive=yes --threshold=99.5 $W/cg.$n 2>/dev/null | c++filt > $W/cg.$n.incl.txt
  echo "== $n: $(grep 'PROGRAM TOTALS' $W/cg.$n.incl.txt | awk '{print $1}')"
  grep -E 'RunLongReads|RunPairedEnd|AnchoredAligner::(Align|Flank|Gap|Reseed|Piece)|Between|alignEndsFree|alignEnd2End|DecodeChunk|Seedmap::Load ' $W/cg.$n.incl.txt | sed -E 's/\[\/home[^]]*\]//; s/\(([^()]|\([^()]*\))*\)//g' | awk '{ $1 = $1; print }' | cut -c1-110 | awk '!seen[$2]++' | head -12
done
echo LR4 DONE
