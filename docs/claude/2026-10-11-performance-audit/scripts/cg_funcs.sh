#!/usr/bin/env bash
# Inclusive instruction counts of chosen functions in each callgrind profile of the audit (names matched loosely).
R=$HOME/perf9/runs/callgrind
names=(
  "PROGRAM TOTALS"
  "protal::LoadIndex("
  "LoadAllGenomes("
  "protal::LoadStrainAlleles("
  "ReadGeneNeighbours("
  "BindGeneNeighbours("
  "classify::RunPairedEnd<"
  "classify::RunSingleEnd<"
  "classify::RunLongReads<"
  "ChainAnchorFinder<protal::KmerLookupSM>::operator()"
  "KmerLookupSM::GetFromLookup"
  "SimpleAlignmentHandler::operator()"
  "SimpleAlignmentHandler::ScoreAlleles"
  "SimpleAlignmentHandler::TryCongeners"
  "SimpleAlignmentHandler::SettleBySites"
  "LongReadAligner::Align"
  "IndelBound::MayAlign"
  "SimpleAlignmentHandler::AlignAnchor"
  "AlignmentScreen::MayAlign"
  "IndelBound::LcsAbove"
  "WFA2Wrapper2::Alignment"
  "ProfileSample("
  "ProfileSam"
  "RecordEvidenceCollector::NoteRecord"
  "AncestrySitesOf"
  "ColumnWeightsOf"
  "ApplySampleContext"
  "ReleaseReadData"
  "malloc_trim"
  "LoadColumnWeights("
  "LoadAncestrySites("
  "LoadCongenerGaps("
)
for f in "$R"/callgrind.out.*; do
  name=${f##*callgrind.out.}
  [ -s "$R/$name.all.txt" ] || callgrind_annotate --inclusive=yes --threshold=100 "$f" 2> /dev/null > "$R/$name.all.txt"
  echo "== $name"
  for n in "${names[@]}"; do
    line=$(grep -F -m1 "$n" "$R/$name.all.txt")
    [ -n "$line" ] && printf '%s\t%s\n' "$(echo "$line" | awk '{print $1}')" "$n"
  done
done
