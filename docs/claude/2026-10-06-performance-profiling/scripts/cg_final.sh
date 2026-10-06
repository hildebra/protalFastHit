#!/usr/bin/env bash
# cg_final.sh: callgrind (one thread) of HEAD (~/perf8/ref) and the final build (~/perf8/isal: shared seeds, packed read
# strands, ISA-L) aligning 100k pairs (single-member gzip input) and PacBio 3 Mb; the functions this round changed.
set -uo pipefail
W=$HOME/perf8; DB=$HOME/bench071/V075/protal_db; P=$HOME/bench071/samples/points; RD=$HOME/perf6/reads
mkdir -p $W/cgf
cg() { local b=$1 tag=$2; shift 2; rm -rf $W/cgf/$tag
  valgrind --tool=callgrind --callgrind-out-file=$W/cgf/$tag.out $b --db $DB "$@" --prefix s -o $W/cgf/$tag -t 1 --no_qcmsa --no_profile > $W/cgf/$tag.log 2>&1
  callgrind_annotate --inclusive=yes --threshold=100 $W/cgf/$tag.out 2>/dev/null > $W/cgf/$tag.incl.txt
  echo "    $tag: $(grep -m1 'PROGRAM TOTALS' $W/cgf/$tag.incl.txt)"; }
for build in ref isal; do
  cg $W/$build/build/protal pe.$build -1 $RD/pe100k_R1.fq.gz -2 $RD/pe100k_R2.fq.gz --read_type pe
  cg $W/$build/build/protal pb.$build -1 $P/pb_b3000000/sim/reads/pb_b3000000_s_1.fq.gz --read_type pb
done
for set in pe pb; do
  for f in 'RunPairedEnd' 'RunLongReads' 'ChainAnchorFinder<protal::KmerLookupSM>::operator()' 'ChainAnchorFinder<protal::KmerLookupSM>::Sort' \
           'SharedSeeds::Sort' 'SimpleAlignmentHandler::AlignAnchor' 'AlignmentScreen::MayAlignPacked' 'AlignmentScreen::SharesEnough' \
           'ThreadedGzStreambuf::Inflate' 'isal_inflate' 'zng_inflate' 'libdeflate_deflate_decompress' 'isal_inflate_stateless'; do
    h=$(grep -m1 -F "$f" $W/cgf/$set.ref.incl.txt | awk '{print $1}'); w=$(grep -m1 -F "$f" $W/cgf/$set.isal.incl.txt | awk '{print $1}')
    [ -n "$h$w" ] && echo "    $set $f: ref ${h:--}, final ${w:--}"
  done
done
for b in ref isal; do grep -h -E "refused by the k-mer screen" $W/cgf/pe.$b.log | cut -c1-200 | sed "s/^/    $b: /"; done
echo "CGFINALDONE $(date +%T)"
