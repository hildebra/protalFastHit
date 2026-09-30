# WFA2-lib, built for protal as the static library `wfa_lib`.
#
# lib/wfa2-lib is WFA2-lib v2.3.6, unmodified, from
#   https://github.com/smarco/WFA2-lib/archive/refs/tags/v2.3.6.tar.gz
#   (sha256 fd74c4bfdd5764ae8668cdeee7a80bfa35583b1980e261daa9dbf425bf12bb0b)
# without its img/, tests/, tools/, examples/ and .github/ directories.
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
# UBSan is off for WFA2 (ASan stays on): it compares sequences in unaligned 8-byte blocks and shifts
# negative offsets left, idioms UBSan reports and that would stop a sanitizer run with
# halt_on_error=1 as soon as a test aligns. They are not protal's code; the flag does nothing
# outside sanitizer builds.
target_compile_options(wfa_lib PRIVATE -fno-sanitize=undefined)
set_target_properties(wfa_lib PROPERTIES LINKER_LANGUAGE CXX)
