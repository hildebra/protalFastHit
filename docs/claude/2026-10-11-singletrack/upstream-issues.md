# Draft: two issues for tmokveld/WFA2-lib

Drafted 2026-10-11 for the user to file at <https://github.com/tmokveld/WFA2-lib/issues>; not sent.

**How the drafts were checked**:
- The reproducer is `scripts/repro_endsfree.c`, run against the fork at `85d063f` built with its own `make clean all`
  (gcc 13.3, Ubuntu 24.04, x86-64).
- It was also run against smarco/WFA2-lib v2.3.6 and `main` (`3282f83`), and against protal's patched copy
  (`lib/wfa2-lib`).
- Background is in [README.md](README.md), section 2.

Both issues use the same program and the same pair of sequences. Each issue is written to stand alone.

---

## Issue 1

**Title:** Singletrack backtrace exits the process (`exit(-1)`) on an ends-free alignment pruned by wf-adaptive

**Body:**

With `wavefront_memory_singletrack`, an ends-free gap-affine alignment can end the whole process inside the
backtrace. The cause is that the M-only backtrace finds no path back to the start. It prints
`[WFA::Backtrace] I?/D?-Beginning backtrace error` and the two sequences, then calls `exit(-1)`
(`wavefront_backtrace_affine_m_only` in `wavefront/wavefront_backtrace.c`).

It needs the default wf-adaptive heuristic. With `wf_heuristic_none` the same alignment completes, with the optimal
score.

We vendor the fork in [protal](https://protal.earlham.ac.uk) (a metagenomic profiler). It ended a PacBio run on a read's
flank: the read's end free, the reference's start not.

**To reproduce** (fork at 85d063f, `make clean all`, gcc 13.3, x86-64):

```c
// repro_endsfree.c: cc -O2 -I. repro_endsfree.c lib/libwfa.a -lm -o repro_endsfree
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
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;  // gap-affine 0/4/6/2, wf-adaptive 10/50/1
  attributes.alignment_scope = compute_alignment;
  attributes.memory_mode = (argc > 1 && strcmp(argv[1], "singletrack") == 0) ?
      wavefront_memory_singletrack : wavefront_memory_high;
  attributes.alignment_form.span = alignment_endsfree;
  attributes.alignment_form.pattern_begin_free = 0;
  attributes.alignment_form.pattern_end_free = 0;
  attributes.alignment_form.text_begin_free = 0;
  attributes.alignment_form.text_end_free = tlen;  // only the text's end free
  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
  const int status = wavefront_align(wf_aligner, pattern, plen, text, tlen);
  const cigar_t* const cigar = wf_aligner->cigar;
  printf("status %d, score %d\n", status, cigar->score);
  int v = 0, h = 0, i;
  for (i = cigar->begin_offset; i < cigar->end_offset; ++i) {  // is the CIGAR an alignment of the sequences?
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
```

```
$ ./repro_endsfree singletrack
[WFA::Backtrace] I?/D?-Beginning backtrace error
>AGATAAACGCCCTTGCCGCTTTCACGCATCG...
<CCTGACGCGATCTAGGTAGAACATGGGATTT...
$ echo $?
255
```

With `attributes.heuristic.strategy = wf_heuristic_none` it completes: `status 0, score -298`, and the CIGAR is valid.
A second pair, 309 × 544 bp, behaves the same; we can send it.

**What we think happens:**
1. wf-adaptive (`wavefront_heuristic_cufoff`) narrows a wavefront's `lo`/`hi` after it is computed, but leaves the
   pruned offsets in the array. `wf_heuristic_nullify_boundaries` nulls only `lo-1` and `hi+1`.
2. An ends-free output wavefront can be wider than its inputs (`wavefront_compute_endsfree_limits`). Its kernels then
   read the inputs over that wider range.
3. `wavefront_compute_init_ends_wf_higher` and `_lower` null only beyond `max(wf_elements_init_max, hi)` (and
   `min(wf_elements_init_min, lo)`). So the kernels read the pruned offsets between `hi+2` and `wf_elements_init_max`,
   and a state can descend from a pruned M state.
4. `wavefront_backtrace_affine_m_only` looks a predecessor up with `wavefront_backtrace_m_offset`, which accepts only
   `k` within the active `lo`/`hi`. For such a state it finds none, searches gap lengths until the score runs out, and
   ends in the error branch.

We traced this in the 309 × 544 bp pair. There the path ends `...XXXDDDMMM` before the free end of the text, and the
3-base deletion opened from an M state outside the pruned bounds of its wavefront.

**A fix that works for us** (we would gladly open a PR):
- In `wavefront_backtrace_affine_m_only`, when no predecessor is found within `lo`/`hi`, try the same cell again with
  the lookups over `[wf_elements_init_min, wf_elements_init_max]`, i.e. what the kernels read. Results that work now
  stay as they are: the wider lookup is only tried after the active bounds gave no path.
- Not exiting from the library: return `false` from the backtrace and fail the alignment, e.g. as
  `WF_STATUS_UNATTAINABLE`, so that a caller can recover. We realign such a pair with `wavefront_memory_high`.

With this, Singletrack gave the classic traceback's outcome and score, and an alignment of that score, in all of
33,500 synthetic short-read, long-read, ends-free-flank and end-to-end alignments, and in all our whole-run tests. It
never needed the fallback.

---

## Issue 2

**Title:** High-memory backtrace returns an invalid CIGAR after wf-adaptive pruning (since caecb60, active M-wavefront
bounds)

**Body:**

With `wavefront_memory_high` (the classic backtrace) and the default wf-adaptive heuristic, an ends-free gap-affine
alignment returns the right score but a CIGAR that is no alignment of the two sequences: an `M` on two different bases.

| Version | Score | CIGAR |
|---|---|---|
| fork at 85d063f | 322 | invalid |
| smarco/WFA2-lib v2.3.6 | 322 | valid |
| smarco/WFA2-lib main (3282f83) | 322 | valid |

**To reproduce:** the program below (fork at 85d063f, `make clean all`, gcc 13.3, x86-64):

```c
[the same repro_endsfree.c as in issue 1]
```

```
$ ./repro_endsfree high
status 0, score -322
invalid CIGAR: operation 151 is 'M', but pattern[149] = 'C' and text[101] = 'A'
```

Built against upstream v2.3.6 or main (3282f83):

```
status 0, score -322
CIGAR valid (pattern 152/152, text 218/218)
```

With `wf_heuristic_none` the fork's CIGAR is valid too (`score -298`). A second pair, 309 × 544 bp, gives an invalid
`M` at operation 358; we can send it.

**What we think happens:**
- `caecb60` ("Use active M-wavefront bounds during traceback") made `wavefront_backtrace_m_offset` accept only `k`
  within an M-wavefront's active `lo`/`hi`.
- The classic backtrace uses it for its mismatch and gap-opening lookups (`wavefront_backtrace_misms`,
  `wavefront_backtrace_del1_open`, `_ins1_open`, ...). The gap-extension lookups still use
  `wf_elements_init_min/max`, as upstream does for all of them.
- After wf-adaptive pruning, the compute kernels of a wider ends-free wavefront still read pruned offsets beyond `hi+1`
  (`wf_heuristic_nullify_boundaries` nulls only `lo-1`/`hi+1`; `wavefront_compute_init_ends_*` null only beyond the
  initialized range). So a state can descend from a pruned one.
- When its true predecessor is such a state, the active-bounds lookup misses it. The backtrace then takes another
  branch and writes operations that do not match the sequences.

**A fix that works for us:** the classic backtrace's M lookups over `[wf_elements_init_min, wf_elements_init_max]`
again, as in upstream (an element that is `WAVEFRONT_OFFSET_NULL` matches nothing).
- With this alone, no invalid CIGAR in 33,500 synthetic alignments.
- With `wf_heuristic_nullify_boundaries` also off (next paragraph), the fork's high-memory results were identical to
  upstream v2.3.6's: status, score and CIGAR in all 33,500.

**Related observation**, perhaps intended: `wf_heuristic_nullify_boundaries` (new in the fork) changes heuristic
results relative to upstream.
- In 180 of the 33,500 synthetic alignments the result differed from v2.3.6 (in 71 the score), always worse: 482
  against 486, 762 against 778. All were ends-free flanks pruned by wf-adaptive and given a large step budget.
- If the aim was to stop gaps opening from pruned states, nulling only `lo-1`/`hi+1` does not do that for wider
  ends-free outputs (see above). If it was not needed for something else, dropping it would match upstream again.
