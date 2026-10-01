#!/bin/bash
# Stage breakdown (inclusive instructions) of a callgrind.sh run, and call counts of the hot
# functions: cg_stage.sh LABEL PAIRS
source "$(dirname "$0")/env.sh"
out=$PERF_DIR/cg/$1; pairs=$2
f=$out/inclusive_all.txt; tree=$out/caller_tree.txt
[ -s $f ] || callgrind_annotate --inclusive=yes --threshold=100 $out/callgrind.out 2>/dev/null > $f
[ -s $tree ] || callgrind_annotate --tree=caller --inclusive=yes --threshold=100 $out/callgrind.out 2>/dev/null > $tree
tot=$(grep "PROGRAM TOTALS" $f | awk '{gsub(",","",$1); print $1}')
get() { grep -F "$1" $f | grep -v "'2" | head -1 | awk '{v=$1; gsub(",","",v); print v}'; }
rp=$(get "protal::Statistics protal::classify::RunPairedEnd<")
row() { v=$(get "$2"); [ -z "$v" ] && v=0; awk -v n="$1" -v v=$v -v t=$tot -v p=$pairs -v rp=$rp 'BEGIN{printf "%-36s %8.2f G  %6.2f%% of run  %6.1f%% of loop  %8.0f Ir/pair\n", n, v/1e9, 100*v/t, 100*v/rp, v/p}'; }
calls() { awk -v fn="$1" '
  /^ *[0-9,]+ .*< / { if (match($0, /\(([0-9,]+)x\)/)) { s = substr($0, RSTART+1, RLENGTH-3); gsub(",", "", s); pend += s } next }
  /\*  / { if (index($0, fn) > 0) total += pend; pend = 0; next }
  { pend = 0 }
  END { print total + 0 }' $tree; }
echo "$1: total $(awk -v t=$tot 'BEGIN{printf "%.2f", t/1e9}') G instructions, $pairs pairs"
row "index load (Seedmap::Load)" "protal::Seedmap::Load("
row "  zstd decompress" "ZSTD_decompressDCtx"
row "  column decode" "protal::index_codec::detail::DecodeChunk"
row "model parse (cPMML from_string)" "cpmml::Model::from_string"
row "  exceptions (__cxa_throw)" "???:__cxa_throw"
row "genomes (LoadAllGenomes)" "protal::GenomeLoader::LoadAllGenomes(int) ["
row "alignment loop (RunPairedEnd)" "protal::Statistics protal::classify::RunPairedEnd<"
row "  FASTQ reader (SeqReaderPE)" "protal::SeqReaderPE::operator()"
row "    gzread" "???:gzread"
row "  syncmers (SimpleKmerHandler)" "protal::SimpleKmerHandler<protal::ClosedSyncmer>::operator()(std::basic_string_view<char, std::char_traits<char> > const&&"
row "  seeds+anchors (ChainAnchorFinder)" "protal::ChainAnchorFinder<protal::KmerLookupSM>::operator()("
row "  alignment (SimpleAlignmentHandler)" "protal::SimpleAlignmentHandler::operator()"
row "    WFA (alignEndsFree)" "wfa::WFAligner::alignEndsFree"
row "    WF-adaptive cut-off" "???:wavefront_heuristic_cufoff"
row "  SAM output handler" "protal::ProtalPairedOutputHandler<false>::operator()"
row "profiling (ProfileWrapper)" "protal::ProfileWrapper(protal::Options&, protal::ProtalDB&, protal::profiler::TaxonFilterForest const&) ["
row "  SAM parse (ReadSamGroups)" "ReadSamGroups<"
for fn in "wfa::WFAligner::alignEndsFree" "protal::Seedmap::Get(" "protal::KmerLookupSM::GetFromLookup" "???:__cxa_throw"; do
  printf "  calls %-40s %10d  (%.2f per read)\n" "$fn" "$(calls "$fn")" "$(awk -v c=$(calls "$fn") -v p=$pairs 'BEGIN{print c/(2*p)}')"
done
