#!/bin/bash
# Turns the gdb backtraces of sample.sh into two tables: the innermost frame, and the stage of the first
# frame (from the innermost out) that belongs to one. sample_report.sh RAW OUTDIR
raw=$1; out=$2
awk '
function stage(s) {
  if (s ~ /ReverseComplement|Complement/) return "reverse complement"
  if (s ~ /ScanWindows|FillSmers|SimpleKmerHandler/) return "syncmers"
  if (s ~ /steady_clock|clock_gettime|Benchmark/) return "stage timers"
  if (s ~ /Seedmap::Get|Seedmap::|KmerLookupSM|FindSeeds|LookupPointer|LookupResult/) return "seeding (index lookups)"
  if (s ~ /ChainAnchorFinder|ChainAlignmentAnchor|ChainingStrategy|RecoverAnchors/) return "anchors"
  if (s ~ /wavefront|WFAligner|WFA2Wrapper|wfa::/) return "WFA"
  if (s ~ /AnchoredAligner/) return "anchored aligner (non-WFA)"
  if (s ~ /GeneSequence|Gene::Sequence|packed::|Unpack/) return "gene decode"
  if (s ~ /PostProcess|AlignmentInfo|GetInstructionCounts|IsAlignmentValid|SimpleAlignmentHandler/) return "alignment handler (rest)"
  if (s ~ /OutputHandler|SamEntry|ArtoSAM|ExtractSNPs|VariantHandler|SamOutput|BufferedString/) return "SAM output"
  if (s ~ /FastxReader|SeqReader|NextFastq|Uppercase|ThreadedGz|memchr/) return "FASTQ reader"
  if (s ~ /JoinAlignmentPairs|SortAlignmentPairs/) return "pairing/sorting"
  if (s ~ /RunPairedEnd/) return "read loop (inlined rest)"
  return ""
}
/^===/ { if (n) flush(); n = 0; skip = 0; next }
/^\[Switching/ { skip = 1; next }
/^#[0-9]+ / { if (skip && $1 == "#0") { skip = 0; next } skip = 0
  line = $0; sub(/^#[0-9]+ +(0x[0-9a-f]+ in )?/, "", line); sub(/ \(.*/, "", line); fr[n++] = line; next }
function flush(   i, st) { top[fr[0]]++; total++
  st = ""; for (i = 0; i < n && st == ""; i++) st = stage(fr[i]); if (st == "") st = "other: " substr(fr[0], 1, 40); cat[st]++ }
END { if (n) flush()
  print "samples: " total; print "\nby stage:"; for (c in cat) printf "%5d %5.1f%%  %s\n", cat[c], 100*cat[c]/total, c | "sort -rn | head -25"; close("sort -rn | head -25")
  print "\ninnermost frame:"; for (f in top) printf "%5d %5.1f%%  %s\n", top[f], 100*top[f]/total, substr(f, 1, 110) | "sort -rn | head -25" }' $raw > $out/report.txt
cat $out/report.txt
