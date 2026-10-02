#!/usr/bin/env bash
# callgrind of the reference (HEAD) and the long-read prototype on the Nanopore and PacBio 3 Mb samples, one
# thread, SAM only; inclusive counts of the alignment parts.
W=$HOME/mt-work/perf3; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
B=${B:-lr}
declare -A BIN=( [ref]=$W/ref/build/protal [$B]=$HOME/mt-work/perf3-$B/src/build/protal )
one() { local s=$1 b=$2 t=ont; [ $s = pb_b3000000 ] && t=pb
  rm -rf $W/cg.o.$s.$b
  nice valgrind --tool=callgrind --callgrind-out-file=$W/cg.$s.$b ${BIN[$b]} --db $DB -1 $P/$s/sim/reads/${s}_s_1.fq.gz --read_type $t --no_profile --prefix s -o $W/cg.o.$s.$b -t 1 --no_qcmsa > /dev/null 2>&1
  callgrind_annotate --inclusive=yes --threshold=99.5 $W/cg.$s.$b 2>/dev/null | c++filt > $W/cg.$s.$b.incl.txt
  callgrind_annotate --threshold=100 $W/cg.$s.$b 2>/dev/null | c++filt > $W/cg.$s.$b.self.txt
}
for s in ont_b3000000 pb_b3000000; do for b in ${BUILDS:-ref $B}; do one $s $b & done; done; wait
for s in ont_b3000000 pb_b3000000; do for b in ${BUILDS:-ref $B}; do
  echo "== $s $b: $(grep 'PROGRAM TOTALS' $W/cg.$s.$b.incl.txt | awk '{print $1}')"
  grep -E 'RunLongReads|AlignAnchor|AnchoredAligner::(Align|Flank|Gap)|Between|alignEndsFree|wavefront_align |LongReadAligner.*::Align|Seedmap::Load ' $W/cg.$s.$b.incl.txt | sed -E 's/\[\/home[^]]*\]//; s/\(([^()]|\([^()]*\))*\)//g' | awk '{ $1 = $1; print }' | cut -c1-140 | head -12
done; done
