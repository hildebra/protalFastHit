// wfa_issue_repro.c: two ends-free gap-affine alignments with WFA2-lib's default attributes (match 0, mismatch 4,
// gap opening 6, gap extension 2, wf-adaptive 10/50/1), the text's (read's) end free, everything else not.
//   ./repro 1 singletrack   case 1 with wavefront_memory_singletrack (only in builds with Singletrack)
//   ./repro 2 high          case 2 with wavefront_memory_high
// Prints the status, the score, and whether the CIGAR is an alignment of the two sequences (every M on equal
// characters, every X on different ones, both sequences consumed) and what it costs.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wavefront/wavefront_align.h"

static const char* const kPattern1 =
    "AACCTACTAGAGCACTTACCCGACCTAGACGACTGCGTCAAGTCCGCCGTACCTTCCTGCCAACCGGAACATGAGGCCCATTTCGAGCAGCGCTTTCC"
    "AGAACGGGAACCCCAGCTTCTGCTACCGGTCCCTCGTCGGCCCTATCAGGAACCGGTCGAACGCGTCGTACCGGTGGTGCGGGGAGTGCCGTGGCCAC"
    "CCTTTCGGCAGGCACCACCAGTCGCCGCCGTCGCGTACCTTCCGCGCTTTATCCATCCTACTCGGCAGGTGGTCTTATTACCCAAACAGCAAGCGGCG"
    "GCGCAGACGGTGGTA";
static const char* const kText1 =
    "CCGCTACCCCCGACGGCCGTCCCGCATCCGCTTGTTGACTGGCGCTGAGTCTACGCCCGCCAGGCCGCCGCGCCTAGCGCAGTACGCAGTCAAGGAGGA"
    "CCAGCCGATGTTGCCGTCAGTCTGGCCGGAAACCGGGCGGACCGACGGGTGGCCGAGCTCGTGGTGGTGCCCTTAGCAGGCGGGCCGGCGGGCGACGAC"
    "TATGAGCAGTATCCCCCCTATCTAGCGCCTGTAGACAAAGGTGTCTGCGCGAGCGAGCGACGGAAGCCGCTGCGCCCGAAAGCGGAAGTACACCCACCG"
    "CGGGTGGTAGTCAAGCGTCAGGCCTAGCCGAAGCTCCAACTCCATGGTGTTGCCTACCAGCCCGTCCGGCGGGAACCTACGCTAGTGGTGCATCACTGG"
    "CGGCTTGAACGGCATCTAGTCGCCGCGCCCGGCGCCGAGCGGCCAGGAAGCTGGGGCCCAGCACTCGTCGTCGCAGCCGTCTTACGGGCCGGGCGGGTG"
    "GCGTCACCTCAGGGTGTACACGCGGAAAAGGCCTAGGGCGACTATCTCA";
static const char* const kPattern2 =
    "AGATAAACGCCCTTGCCGCTTTCACGCATCGTAATCCCTGAGCCCTTTCTTAATTACATGAATGTCACCCCTCGTGAATCCGTATTCTTAGTACCTGTA"
    "TCCGAATGGAATCTACACCGGGCATCGGTACCGATCAGTTTCTCGGCGTTCCT";
static const char* const kText2 =
    "CCTGACGCGATCTAGGTAGAACATGGGATTTGGGATGAATTGGCCCCGACAGCTAACCCGGGGCCAGTTCTTTATGGCGTCTAACGAATTGAGTTTGCAA"
    "TAAAATCAGATTAGACGGGCAAAAGCTCGTCTAAAAGAAAAAGGACCTTATGTTGCCTCGGTCCCAAGCATTCTCTTCTGTCGGCCTCGGGGAAAATCTC"
    "GCGCGATAATATGCTGCC";

int main(int argc, char** argv) {
  if (argc < 3) { fprintf(stderr, "usage: %s 1|2 high|singletrack\n", argv[0]); return 2; }
  const char* const pattern = argv[1][0] == '1' ? kPattern1 : kPattern2;
  const char* const text = argv[1][0] == '1' ? kText1 : kText2;
  const int plen = (int)strlen(pattern), tlen = (int)strlen(text);
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;  // gap-affine 0/4/6/2, wf-adaptive 10/50/1
  attributes.alignment_scope = compute_alignment;
  attributes.memory_mode = wavefront_memory_high;
  if (strcmp(argv[2], "singletrack") == 0) {
#ifdef HAVE_SINGLETRACK
    attributes.memory_mode = wavefront_memory_singletrack;
#else
    fprintf(stderr, "this WFA2-lib has no Singletrack\n");
    return 2;
#endif
  }
  attributes.alignment_form.span = alignment_endsfree;
  attributes.alignment_form.pattern_begin_free = 0;
  attributes.alignment_form.pattern_end_free = 0;
  attributes.alignment_form.text_begin_free = 0;
  attributes.alignment_form.text_end_free = tlen;
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  const int status = wavefront_align(wf_aligner, pattern, plen, text, tlen);
  const cigar_t* const cigar = wf_aligner->cigar;
  printf("case %s, %s: pattern %d, text %d bases; status %d, score %d\n", argv[1], argv[2], plen, tlen, status, cigar->score);
  // Walk the CIGAR.
  int v = 0, h = 0, i, bad = -1, bad_v = 0, bad_h = 0, cost = 0, in_gap = 0;
  for (i = cigar->begin_offset; i < cigar->end_offset; ++i) {
    const char op = cigar->operations[i];
    if (op == 'M' || op == 'X') {
      if (bad < 0 && (v >= plen || h >= tlen || (pattern[v] == text[h]) != (op == 'M'))) {
        bad = i - cigar->begin_offset; bad_v = v; bad_h = h;
      }
      if (op == 'X') cost += 4;
      ++v; ++h; in_gap = 0;
    } else if (op == 'D' || op == 'I') {
      const int trailing_free = (op == 'I' && v == plen);  // the text's free end
      if (!trailing_free) { cost += (in_gap == op ? 2 : 8); }
      in_gap = op;
      if (op == 'D') ++v; else ++h;
    }
  }
  if (bad >= 0) {
    printf("CIGAR is no alignment of the sequences: operation %d is '%c' on pattern[%d] = '%c', text[%d] = '%c'\n", bad,
           cigar->operations[cigar->begin_offset + bad], bad_v, pattern[bad_v], bad_h, text[bad_h]);
  }
  printf("CIGAR %s: consumed pattern %d/%d, text %d/%d; its cost %d\n", bad >= 0 ? "INVALID" : "valid", v, plen, h, tlen, cost);
  wavefront_aligner_delete(wf_aligner);
  return 0;
}
