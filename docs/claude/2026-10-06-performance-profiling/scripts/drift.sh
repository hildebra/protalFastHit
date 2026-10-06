#!/usr/bin/env bash
# The drift since performance pass 2: its build (0e4c5a4, ~/perf-pass2/head) and HEAD on the same 100k pairs and the
# v0.7.3 database it used, callgrind, aligning only; then both profiling the 500k-pair SAM HEAD wrote on V073.
set -uo pipefail
W=$HOME/perf6; DB=$HOME/bench071/V073/protal_db; R=$W/reads
mkdir -p $W/drift
for tag in pass2 head; do
  b=$HOME/perf-pass2/head/build/protal; [ $tag = head ] && b=$W/head/build/protal
  rm -rf $W/drift/o.$tag
  valgrind --tool=callgrind --callgrind-out-file=$W/drift/$tag.cg $b --db $DB -1 $R/pe100k_R1.fq.gz -2 $R/pe100k_R2.fq.gz \
    --read_type pe --no_profile --prefix s -o $W/drift/o.$tag -t 1 --no_qcmsa > $W/drift/$tag.log 2>&1
  callgrind_annotate --inclusive=yes $W/drift/$tag.cg 2>/dev/null > $W/drift/$tag.incl.txt
  echo "$tag: $(grep -m1 'PROGRAM TOTALS' $W/drift/$tag.incl.txt)"
  for f in 'RunPairedEnd' 'Seedmap::Load' 'SimpleAlignmentHandler::operator' 'ChainAnchorFinder<protal::KmerLookupSM>::operator' 'OutputHandler<false>::operator' 'AlignmentScreen::MayAlign' 'GuideMate' 'JoinAlignmentPairs' 'SimpleKmerHandler<protal::ClosedSyncmer>::operator' 'ReadGeneNeighbours' 'LoadAllGenomes' 'ThreadedGzStreambuf::Inflate'; do
    echo "   $f: $(grep -m1 -F "$f" $W/drift/$tag.incl.txt | awk '{print $1}')"
  done
done
echo DRIFTDONE
