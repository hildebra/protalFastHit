#!/usr/bin/env bash
# collect.sh <dest>: the run summaries, logs, callgrind summaries and bench outputs of this round into <dest>/results.
set -uo pipefail
D=$1; W=$HOME/perf6
mkdir -p $D/results/logs $D/results/callgrind $D/results/bench $D/results/prototypes
cp $W/head.commit $D/results/
bash $W/summary.sh > $D/results/runs.txt
for f in $W/o.*.log $W/o.*.time; do cp $f $D/results/logs/; done
for n in align prof pb ont se; do
  head -120 $W/cg/$n.incl.txt > $D/results/callgrind/$n.incl.txt
  head -90 $W/cg/$n.self.txt > $D/results/callgrind/$n.self.txt
  grep -E "took|seeding:|candidate|SAM read|read EM" $W/cg/$n.log > $D/results/callgrind/$n.log.txt
done
for f in Alignment/AlignmentScreen.h SequenceUtils/KmerIterator.h Hash/KmerLookup.h Hash/FlexScan.h Core/ChainAnchorFinder.h SequenceUtils/GenomeLoader.h; do
  n=$(basename $f .h)
  bash $W/annot.sh $W/cgg/align.out $f 30 > $D/results/callgrind/lines_align_$n.txt
done
bash $W/annot.sh $W/cgg/pb.out Alignment/AlignmentScreen.h 30 > $D/results/callgrind/lines_pb_AlignmentScreen.txt
cp $W/bench/*.txt $D/results/bench/ 2>/dev/null
for p in screen seed; do
  for s in pe pb; do
    [ -f $W/cmp_$p/$s.$p.incl.txt ] && head -150 $W/cmp_$p/$s.$p.incl.txt > $D/results/prototypes/${s}_$p.incl.txt
  done
done
head -150 $W/cmp_head_pe.incl.txt > $D/results/prototypes/pe_head.incl.txt
head -150 $W/cmp_head_pb.incl.txt > $D/results/prototypes/pb_head.incl.txt
[ -d $W/drift ] && for t in pass2 head; do head -150 $W/drift/$t.incl.txt > $D/results/callgrind/drift_$t.incl.txt 2>/dev/null; done
ls -R $D/results | head -80
