# WFA2-lib, built for protal as the static library `wfa_lib`.
#
# lib/wfa2-lib is WFA2-lib v2.3.6 with the Singletrack traceback (wavefront_memory_singletrack;
# López-Villellas et al., bioRxiv 10.1101/2025.10.31.685625), from T. Mokveld's fork
#   https://github.com/tmokveld/WFA2-lib at commit 85d063f0c4121d8331b91d7f953b4caa7024d2f8 (2026-10-05),
# which merges upstream's main (smarco/WFA2-lib 3282f83, after v2.3.6) and adds Singletrack for
# end-to-end, ends-free and extension alignments with the adaptive, drop and band heuristics, besides
# fixes; without its img/, tests/, tools/, examples/ and .github/ directories. Until 2026-10-11 it was
# upstream's v2.3.6 tarball. protal's classic traceback (--wfa_traceback classic) is WFA2's high-memory
# mode, as before; its outputs were the same as with v2.3.6 (docs/claude/2026-10-11-singletrack).
#
# protal's changes to it (comments starting "protal"; docs/claude/2026-10-11-singletrack). wf-adaptive
# narrows a wavefront's lo/hi but leaves the pruned offsets, and the kernels of a later, wider ends-free
# wavefront still read them, so a state can descend from a pruned one:
#  1. wavefront_heuristic.c: the fork's nulling of the elements next to a pruned wavefront is off, as in
#     v2.3.6 (with it, heuristic-pruned alignments past protal's budgets scored worse than with v2.3.6).
#  2. wavefront_backtrace.c: the classic backtrace's M lookups read any initialized element, as in v2.3.6
#     (the fork's active bounds gave operations that are no alignment of the sequences).
#  3. wavefront_backtrace.c, .h, wavefront_unialign.c: Singletrack's M-only backtrace looks in the
#     initialized elements too when lo/hi give no path (it ended the process, exit(-1), on a PacBio
#     flank with the read's end free); one that still finds none fails the alignment
#     (WF_STATUS_UNATTAINABLE) instead of exiting, and WFA2Wrapper2 aligns it with the classic traceback.
# With them, 33,500 protal-like alignments gave v2.3.6's results byte for byte (classic) and its scores
# (Singletrack). Keep them when updating the library, unless upstream has the fixes (tests
# WfaTraceback.*, scripts in docs/claude/2026-10-11-singletrack).
#
# Speed changes (2026-10-11, docs/claude/2026-10-11-wfa2-singletrack-speed), none changing a result:
#  4. wavefront_heuristic.c: the loops read lo/hi into locals, so GCC can compute their trip count and
#     vectorize them (wf-adaptive's distances were 9-20% of WFA2's instructions, all scalar).
#  5. wavefront_extend_kernels.c: the ends-free extend kernel calls wavefront_termination_endsfree only
#     for a diagonal at the end of the text or the pattern (it was called for every diagonal).
#  6. utils/commons.h WFA_CLONE_V3 on wf-adaptive's distance loops and wavefront_compute_affine_idm: a
#     second copy for x86-64-v3, chosen at run time (not -march: WFA2's AVX2 kernels change ends-free ties).
#  7. wavefront_slab.c, .h: the reset at each alignment frees only the wavefronts handed out since the
#     last one, when the slab holds nothing to discard (it walked every wavefront ever allocated), and
#     rewrites only the part of the free stack the alignment changed.
# WFA2Wrapper2 also aligns with the penalties divided by their greatest common divisor (4/6/2 -> 2/3/1).
# The per-step bookkeeping (init_ends, trim_ends, allocation: ~45% of short reads' WFA2 instructions after
# these) was profiled too; compiling the library with -flto did not help (+1-2%).
#
# Upstream's own CMakeLists.txt is not used: in Release builds it adds -march=native (the baseline
# protal binary must run on any x86-64), and it needs pkg-config and builds shared libraries, a
# benchmark tool and a test. The source list and include directories below are upstream's
# (wfa2lib_SOURCE + the C++ binding); keep them in sync when updating the library. Compiler flags
# come from protal's build type. WFA_PARALLEL stays undefined: protal parallelises over reads.

set(WFA2_DIR ${CMAKE_CURRENT_LIST_DIR}/wfa2-lib)

set(WFA2_SOURCES
        wavefront/wavefront_align.c
        wavefront/wavefront_aligner.c
        wavefront/wavefront_attributes.c
        wavefront/wavefront_backtrace_buffer.c
        wavefront/wavefront_backtrace.c
        wavefront/wavefront_backtrace_offload.c
        wavefront/wavefront_bialign.c
        wavefront/wavefront_bialigner.c
        wavefront/wavefront.c
        wavefront/wavefront_components.c
        wavefront/wavefront_compute_affine2p.c
        wavefront/wavefront_compute_affine.c
        wavefront/wavefront_compute.c
        wavefront/wavefront_compute_edit.c
        wavefront/wavefront_compute_linear.c
        wavefront/wavefront_debug.c
        wavefront/wavefront_display.c
        wavefront/wavefront_extend.c
        wavefront/wavefront_heuristic.c
        wavefront/wavefront_pcigar.c
        wavefront/wavefront_penalties.c
        wavefront/wavefront_plot.c
        wavefront/wavefront_sequences.c
        wavefront/wavefront_slab.c
        wavefront/wavefront_unialign.c
        wavefront/wavefront_termination.c
        wavefront/wavefront_extend_kernels_avx.c
        wavefront/wavefront_extend_kernels.c
        system/mm_stack.c
        system/mm_allocator.c
        system/profiler_counter.c
        system/profiler_timer.c
        utils/bitmap.c
        utils/dna_text.c
        utils/sequence_buffer.c
        utils/vector.c
        utils/commons.c
        utils/heatmap.c
        alignment/affine2p_penalties.c
        alignment/affine_penalties.c
        alignment/cigar.c
        alignment/cigar_utils.c
        alignment/score_matrix.c
        bindings/cpp/WFAligner.cpp)
list(TRANSFORM WFA2_SOURCES PREPEND ${WFA2_DIR}/)

add_library(wfa_lib STATIC ${WFA2_SOURCES})
target_include_directories(wfa_lib PUBLIC ${WFA2_DIR} ${WFA2_DIR}/wavefront ${WFA2_DIR}/utils)
target_compile_definitions(wfa_lib PRIVATE _FILE_OFFSET_BITS=64)
# -DPROTAL_NO_CLONES (given in CMAKE_CXX_FLAGS, as the sanitizer builds do) also compiles WFA2's cloned
# loops once (WFA_CLONE_V3, patch 6 above).
if(CMAKE_CXX_FLAGS MATCHES "PROTAL_NO_CLONES" OR CMAKE_C_FLAGS MATCHES "PROTAL_NO_CLONES")
    target_compile_definitions(wfa_lib PRIVATE PROTAL_NO_CLONES)
endif()
# UBSan is off for WFA2 (ASan stays on): it compares sequences in unaligned 8-byte blocks and shifts
# negative offsets left, idioms UBSan reports and that would stop a sanitizer run with
# halt_on_error=1 as soon as a test aligns. They are not protal's code; the flag does nothing
# outside sanitizer builds.
target_compile_options(wfa_lib PRIVATE -fno-sanitize=undefined)
set_target_properties(wfa_lib PROPERTIES LINKER_LANGUAGE CXX)
