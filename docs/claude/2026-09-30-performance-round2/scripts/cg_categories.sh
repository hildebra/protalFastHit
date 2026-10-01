#!/bin/bash
# Per-pair instructions by category, from cg_diff.sh's perpair.tsv: cg_categories.sh LABEL
# A function is classified by its name (the part before its parameter list), not its signature: the first
# version matched the signature, and every function taking a std::string landed in "strings".
source "$(dirname "$0")/env.sh"
f=$PERF_DIR/cg/$1/perpair.tsv
awk -F'\t' '
function cat(full,   n, p) {
  n = full; sub(/^[^:]*:/, "", n)            # drop the source file
  p = index(n, "("); if (p > 0) n = substr(n, 1, p - 1)
  if (n ~ /wavefront_|wfa::|WFAligner/) return "WFA (alignment core)"
  if (n ~ /inflate|zng_|crc32|adler|ThreadedGz/) return "gzip inflate (zlib-ng)"
  if (n ~ /SimpleKmerHandler|KmerIterator/) return "syncmer extraction"
  if (n ~ /ReverseComplement|Complement/) return "reverse complement"
  if (n ~ /basic_string|char_traits|to_string|__to_chars|charconv/) return "std::string operations"
  if (n ~ /malloc|operator new|operator delete|_int_free|^free$|tcache/) return "malloc/free"
  if (n ~ /memmove|memcpy|memset|memchr|memcmp/) return "memcpy/memset/memchr"
  if (n ~ /FastxReader|Uppercase|SeqReader/) return "FASTQ parsing"
  if (n ~ /ChainAnchorFinder|KmerLookup|Seedmap|ChainingStrategy|robin|SeedingStrategy|LookupResult|LookupPointer/) return "seeding + anchors"
  if (n ~ /AnchoredAligner|SimpleAlignmentHandler|AlignmentInfo|PostProcess|IsAlignmentValid|GeneSequence|UnpackRange|Unpack/) return "alignment handler (not WFA)"
  if (n ~ /OutputHandler|ArtoSAM|ExtractSNPs|SamEntry|SamOutput|ReferenceOf|NextCompressedCigar|JoinAlignmentPairs|SortAlignmentPairs|Bitscore/) return "pairing, SAM output"
  if (n ~ /SequenceRange|VariantHandler|Profiler|profiler|strtol|strtod|Strain|LineSplitter|ReadSamGroups/) return "profiling stage"
  if (n ~ /steady_clock|clock_gettime|chrono|Benchmark/) return "stage timers"
  if (n ~ /RunPairedEnd|RunSingleEnd/) return "read loop (inlined)"
  return "other"
}
{ c = cat($2); s[c] += $1; t += $1 }
END { for (c in s) printf "%8d %5.1f%%  %s\n", s[c], 100*s[c]/t, c; printf "%8d        total\n", t }' $f | sort -rn
