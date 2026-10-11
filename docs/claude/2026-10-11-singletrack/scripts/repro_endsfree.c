// repro_endsfree.c: an ends-free gap-affine alignment with the default attributes (match 0, mismatch 4, gap opening 6,
// gap extension 2, wf-adaptive 10/50/1), only the text's end free. Prints the status and score, and checks that the
// CIGAR is an alignment of the two sequences.
//   cc -O2 -I. repro_endsfree.c lib/libwfa.a -lm -o repro_endsfree && ./repro_endsfree high && ./repro_endsfree singletrack
#include <stdio.h>
#include <string.h>
#include "wavefront/wavefront_align.h"

int main(int argc, char** argv) {
  const char* const pattern =
      "AGATAAACGCCCTTGCCGCTTTCACGCATCGTAATCCCTGAGCCCTTTCTTAATTACATGAATGTCACCCCTCGTGAATCCGTATTCTTAGTACCTGTA"
      "TCCGAATGGAATCTACACCGGGCATCGGTACCGATCAGTTTCTCGGCGTTCCT";
  const char* const text =
      "CCTGACGCGATCTAGGTAGAACATGGGATTTGGGATGAATTGGCCCCGACAGCTAACCCGGGGCCAGTTCTTTATGGCGTCTAACGAATTGAGTTTGCAA"
      "TAAAATCAGATTAGACGGGCAAAAGCTCGTCTAAAAGAAAAAGGACCTTATGTTGCCTCGGTCCCAAGCATTCTCTTCTGTCGGCCTCGGGGAAAATCTC"
      "GCGCGATAATATGCTGCC";
  const int plen = (int)strlen(pattern), tlen = (int)strlen(text);
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.alignment_scope = compute_alignment;
  attributes.memory_mode = (argc > 1 && strcmp(argv[1], "singletrack") == 0) ?
      wavefront_memory_singletrack : wavefront_memory_high;
  attributes.alignment_form.span = alignment_endsfree;
  attributes.alignment_form.pattern_begin_free = 0;
  attributes.alignment_form.pattern_end_free = 0;
  attributes.alignment_form.text_begin_free = 0;
  attributes.alignment_form.text_end_free = tlen;
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  const int status = wavefront_align(wf_aligner, pattern, plen, text, tlen);
  const cigar_t* const cigar = wf_aligner->cigar;
  printf("status %d, score %d\n", status, cigar->score);
  int v = 0, h = 0, i;
  for (i = cigar->begin_offset; i < cigar->end_offset; ++i) {
    const char op = cigar->operations[i];
    if ((op == 'M' || op == 'X') && (pattern[v] == text[h]) != (op == 'M')) {
      printf("invalid CIGAR: operation %d is '%c', but pattern[%d] = '%c' and text[%d] = '%c'\n",
             i - cigar->begin_offset, op, v, pattern[v], h, text[h]);
      break;
    }
    if (op != 'I') ++v;
    if (op != 'D') ++h;
  }
  if (i == cigar->end_offset) printf("CIGAR valid (pattern %d/%d, text %d/%d)\n", v, plen, h, tlen);
  wavefront_aligner_delete(wf_aligner);
  return 0;
}
