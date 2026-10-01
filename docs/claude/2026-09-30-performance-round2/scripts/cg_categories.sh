#!/bin/bash
# Per-pair instructions by category, from cg_diff.sh's perpair.tsv: cg_categories.sh LABEL
source "$(dirname "$0")/env.sh"
f=$PERF_DIR/cg/$1/perpair.tsv
awk -F'\t' '
function cat(n) {
  if (n ~ /wavefront_|wfa::|WFAligner/) return "WFA (alignment core)"
  if (n ~ /inflate|zng_|crc32|adler|ThreadedGz/) return "gzip inflate (zlib-ng)"
  if (n ~ /KmerIterator|avx2intrin|avxintrin|emmintrin|ScanWindows|stl_vector.h:protal::SimpleKmerHandler|vector.tcc:protal::SimpleKmerHandler/) return "syncmer extraction"
  if (n ~ /ReverseComplement|KmerUtils\.h|basic_string|char_traits|string_view/) return "strings: reverse complement, copies"
  if (n ~ /malloc|operator new|operator delete|_int_free|free$/) return "malloc/free"
  if (n ~ /memmove|memcpy|memset|memchr/) return "memcpy/memset/memchr"
  if (n ~ /FastxReader|Uppercase/) return "FASTQ parsing"
  if (n ~ /ChainAnchorFinder|KmerLookup|Seedmap|ChainingStrategy|robin|SeedingStrategy/) return "seeding + anchors"
  if (n ~ /GetInstructionCountsAndCompress|AlignmentUtils|SNPUtils|AnchoredAlignment|AlignmentStrategy|PackedSequence|GeneSequence|AlignmentOutputHandler|SamHandler|SamFile|ScoreAlign/) return "alignment post-processing, SAM"
  if (n ~ /SequenceRange|VariantHandler|Profiler|strtol|strtod|charconv|Strain|Profiling/) return "profiling stage"
  if (n ~ /steady_clock|clock_gettime|chrono|Benchmark/) return "stage timers"
  return "other"
}
{ c = cat($2); s[c] += $1; t += $1 }
END { for (c in s) printf "%8d %5.1f%%  %s\n", s[c], 100*s[c]/t, c; printf "%8d        total\n", t }' $f | sort -rn
