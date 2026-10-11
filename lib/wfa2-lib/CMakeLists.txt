# WFA2-lib CMake build
# For Debian currently with
#
#   cd build
#   cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DOPENMP=1 ..
#   cmake --build .
#   ctest .
#   make install
# See below option statements and the README for build information

cmake_minimum_required(VERSION 3.16)
project(wfa2lib)

set(CMAKE_CXX_STANDARD 17)

include(FeatureSummary)
include(GNUInstallDirs)

find_package(PkgConfig REQUIRED)

feature_summary(
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  WHAT REQUIRED_PACKAGES_NOT_FOUND)

# ---- Options

option(OPENMP "Enable OpenMP" OFF) # enables WFA_PARALLEL
option(PROFILING "Enable profiling" OFF)
option(ASAN "Use address sanitiser" OFF)
option(EXTRA_FLAGS "Add optimization flags for C/C++ compiler" OFF)

# include(CheckIPOSupported) # adds lto
# check_ipo_supported(RESULT ipo_supported OUTPUT output)

option(BUILD_SHARED_LIB "Build shared WFA2 library" OFF)
option(BUILD_CXX_BINDING_LIB "Build C++ binding library" OFF)
option(WFA_BUILD_PERF_TESTS "Enable timing-sensitive performance CTest tests" OFF)

# ---- Dependencies

if(OPENMP)
  include(FindOpenMP)
  set(OPTIMIZE_FLAGS "-DWFA_PARALLEL")
endif(OPENMP)

if(EXTRA_FLAGS)
  set(OPTIMIZE_FLAGS "${OPTIMIZE_FLAGS} ${EXTRA_FLAGS}")
endif(EXTRA_FLAGS)

find_package(Threads)
set_package_properties(Threads PROPERTIES TYPE REQUIRED)

# ---- Build switches
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${ipo_supported})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
          "Choose the type of build, options are: Release|Debug|RelWithDebInfo (for distros)." FORCE)
endif()

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  set(OPTIMIZE_FLAGS "${OPTIMIZE_FLAGS} -D_FILE_OFFSET_BITS=64")
endif()

if ((${CMAKE_BUILD_TYPE} MATCHES Release) OR (${CMAKE_BUILD_TYPE} MATCHES RelWithDebInfo))
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPTIMIZE_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAGS}")
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPTIMIZE_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAGS}")
  add_definitions(-Wfatal-errors)
endif ()

if (ASAN)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
endif(ASAN)

if(PROFILING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif(PROFILING)

if(GPROF)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif(GPROF)

# ---- Include files

file(GLOB INCLUDES
  wavefront/*.h*
  )
file(GLOB UTILS_INCLUDES
  utils/*.h*
  )
file(GLOB ALIGNMENT_INCLUDES
  alignment/*.h*
  )
file(GLOB SYSTEM_INCLUDES
  system/*.h*
  )

set(wfa2lib_SOURCE
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
)

add_library(wfa2_static STATIC ${wfa2lib_SOURCE})
set_target_properties(wfa2_static PROPERTIES OUTPUT_NAME wfa2)
target_include_directories(wfa2_static PUBLIC . wavefront utils)
add_library(wfa2::wfa2_static ALIAS wfa2_static)

if (BUILD_SHARED_LIB)
    add_library(wfa2 SHARED ${wfa2lib_SOURCE})
    set_target_properties(wfa2 PROPERTIES SOVERSION 0)
    target_include_directories(wfa2 PUBLIC . wavefront utils)
    add_library(wfa2::wfa2 ALIAS wfa2)
endif ()

if(OPENMP)
  target_link_libraries(wfa2_static PRIVATE OpenMP::OpenMP_C)
  if (BUILD_SHARED_LIB)
    target_link_libraries(wfa2 PRIVATE OpenMP::OpenMP_C)
  endif ()
endif(OPENMP)

# ---- C++ binding library

if (BUILD_CXX_BINDING_LIB)
    set(wfa2cpp_SOURCE
      bindings/cpp/WFAligner.cpp
    )
    file(GLOB CPP_INCLUDES
      bindings/cpp/*.h*
      )
    add_library(wfa2cpp_static STATIC ${wfa2cpp_SOURCE})
    set_target_properties(wfa2cpp_static PROPERTIES OUTPUT_NAME wfa2cpp)
    target_link_libraries(wfa2cpp_static PUBLIC wfa2_static)
    add_library(wfa2::wfa2cpp_static ALIAS wfa2cpp_static)
    if (BUILD_SHARED_LIB)
      add_library(wfa2cpp SHARED ${wfa2cpp_SOURCE})
      set_target_properties(wfa2cpp PROPERTIES SOVERSION 0)
      target_link_libraries(wfa2cpp PUBLIC wfa2)
      add_library(wfa2::wfa2cpp ALIAS wfa2cpp)
    endif ()

    if(OPENMP)
      target_link_libraries(wfa2cpp_static PRIVATE OpenMP::OpenMP_CXX)
      if (BUILD_SHARED_LIB)
        target_link_libraries(wfa2cpp PRIVATE OpenMP::OpenMP_CXX)
      endif ()
    endif(OPENMP)
endif ()

# ---- Get version

file (STRINGS "VERSION.txt" BUILD_NUMBER)
add_definitions(-DWFA2LIB_VERSION="${BUILD_NUMBER}")
add_definitions(-DVERSION="${BUILD_NUMBER}")

set(wfa2lib_LIBS
)

# add_dependencies(wfa2lib ${wfa2lib_DEPS})

# ---- Build all
# add_subdirectory(tools/align_benchmark)


add_executable(align_benchmark
  tools/align_benchmark/align_benchmark_params.c
  tools/align_benchmark/align_benchmark.c
  tools/align_benchmark/benchmark/benchmark_check.c
  tools/align_benchmark/benchmark/benchmark_edit.c
  tools/align_benchmark/benchmark/benchmark_gap_affine2p.c
  tools/align_benchmark/benchmark/benchmark_gap_affine.c
  tools/align_benchmark/benchmark/benchmark_gap_linear.c
  tools/align_benchmark/benchmark/benchmark_indel.c
  tools/align_benchmark/benchmark/benchmark_utils.c
  system/mm_allocator.c
  system/mm_stack.c
  system/profiler_counter.c
  system/profiler_timer.c
  utils/bitmap.c
  utils/commons.c
  utils/dna_text.c
  utils/heatmap.c
  utils/sequence_buffer.c
  utils/vector.c
  alignment/affine2p_penalties.c
  alignment/affine_penalties.c
  alignment/cigar.c
  alignment/score_matrix.c
  tools/align_benchmark/gap_affine/affine_matrix.c
  tools/align_benchmark/gap_affine/swg.c
  tools/align_benchmark/gap_linear/nw.c
  tools/align_benchmark/gap_affine2p/affine2p_dp.c
  tools/align_benchmark/gap_affine2p/affine2p_matrix.c
  tools/align_benchmark/edit/edit_bpm.c
  tools/align_benchmark/edit/edit_dp.c
  tools/align_benchmark/indel/indel_dp.c
)
target_include_directories(align_benchmark PUBLIC . wavefront system utils tools/align_benchmark)
target_link_libraries(align_benchmark m)
if (BUILD_SHARED_LIB)
  add_dependencies(align_benchmark wfa2)
  target_link_libraries(align_benchmark wfa2)
else ()
  add_dependencies(align_benchmark wfa2_static)
  target_link_libraries(align_benchmark wfa2_static)
endif ()
if(OPENMP)
  target_link_libraries(align_benchmark OpenMP::OpenMP_C)
endif(OPENMP)

# ---- Test

enable_testing()

function(add_wfa_test)
  add_test(
    NAME wfa2lib
    COMMAND ./tests/wfa.utest.sh ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endfunction()

add_wfa_test()

add_test(
  NAME wfa2lib_singletrack
  COMMAND ./tests/wfa.singletrack.utest.sh ${CMAKE_CURRENT_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
  NAME wfa2lib_biwfa_endsfree
  COMMAND ./tests/wfa.biwfa.endsfree.proof.py ${CMAKE_CURRENT_BINARY_DIR}/align_benchmark
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_executable(wfa_packed2bits_regression
  tests/wfa_packed2bits_regression.c
)
target_link_libraries(wfa_packed2bits_regression wfa2_static m)
add_test(
  NAME wfa_packed2bits_regression
  COMMAND wfa_packed2bits_regression
)

add_executable(wfa_backtrace_active_bounds_regression
  tests/wfa_backtrace_active_bounds_regression.c
)
target_link_libraries(wfa_backtrace_active_bounds_regression wfa2_static m)
add_test(
  NAME wfa_backtrace_active_bounds_regression
  COMMAND wfa_backtrace_active_bounds_regression
)

add_executable(wfa_extend_kernel_regression
  tests/wfa_extend_kernel_regression.c
)
target_link_libraries(wfa_extend_kernel_regression wfa2_static m)
add_test(
  NAME wfa_extend_kernel_regression
  COMMAND wfa_extend_kernel_regression
)

add_executable(wfa_zdrop_piggyback_regression
  tests/wfa_zdrop_piggyback_regression.c
)
target_link_libraries(wfa_zdrop_piggyback_regression wfa2_static m)
add_test(
  NAME wfa_zdrop_piggyback_regression
  COMMAND wfa_zdrop_piggyback_regression
)

add_executable(wfa_pcigar_extend_regression
  tests/wfa_pcigar_extend_regression.c
)
target_link_libraries(wfa_pcigar_extend_regression wfa2_static m)
add_test(
  NAME wfa_pcigar_extend_regression
  COMMAND wfa_pcigar_extend_regression
)

add_executable(wfa_cigar_sam_regression
  tests/wfa_cigar_sam_regression.c
)
target_link_libraries(wfa_cigar_sam_regression wfa2_static m)
add_test(
  NAME wfa_cigar_sam_regression
  COMMAND wfa_cigar_sam_regression
)

add_executable(wfa_vector_regression
  tests/wfa_vector_regression.c
)
target_link_libraries(wfa_vector_regression wfa2_static m)
add_test(
  NAME wfa_vector_regression
  COMMAND wfa_vector_regression
)

add_executable(wfa_vector_debug_regression
  tests/wfa_vector_debug_regression.c
  utils/vector.c
)
target_compile_definitions(wfa_vector_debug_regression PRIVATE VECTOR_DEBUG)
target_include_directories(wfa_vector_debug_regression PUBLIC . utils)
target_link_libraries(wfa_vector_debug_regression m)
add_test(
  NAME wfa_vector_debug_out_of_range
  COMMAND wfa_vector_debug_regression out_of_range
)
set_tests_properties(wfa_vector_debug_out_of_range PROPERTIES WILL_FAIL TRUE)
add_test(
  NAME wfa_vector_debug_type_mismatch
  COMMAND wfa_vector_debug_regression type_mismatch
)
set_tests_properties(wfa_vector_debug_type_mismatch PROPERTIES WILL_FAIL TRUE)
add_test(
  NAME wfa_vector_debug_empty_last
  COMMAND wfa_vector_debug_regression empty_last
)
set_tests_properties(wfa_vector_debug_empty_last PROPERTIES WILL_FAIL TRUE)

add_executable(wfa_vector_append_microbenchmark
  tests/wfa_vector_append_microbenchmark.c
)
target_link_libraries(wfa_vector_append_microbenchmark wfa2_static m)

add_executable(wfa_swg_endsfree_mapping_regression
  tests/wfa_swg_endsfree_mapping_regression.c
  tools/align_benchmark/gap_affine/swg.c
  tools/align_benchmark/gap_affine/affine_matrix.c
  alignment/cigar.c
  system/mm_allocator.c
  utils/commons.c
  utils/vector.c
)
target_include_directories(wfa_swg_endsfree_mapping_regression PUBLIC
  . wavefront system utils tools/align_benchmark
)
target_link_libraries(wfa_swg_endsfree_mapping_regression m)
add_test(
  NAME wfa_swg_endsfree_mapping_regression
  COMMAND wfa_swg_endsfree_mapping_regression
)

add_executable(wfa_biwfa_endsfree_rejection_regression
  tests/wfa_biwfa_endsfree_rejection_regression.c
)
target_link_libraries(wfa_biwfa_endsfree_rejection_regression wfa2_static m)
add_test(
  NAME wfa_biwfa_endsfree_rejection_regression
  COMMAND wfa_biwfa_endsfree_rejection_regression
)

add_executable(wfa_biwfa_endsfree_endpoint_regression
  tests/wfa_biwfa_endsfree_endpoint_regression.c
)
target_link_libraries(wfa_biwfa_endsfree_endpoint_regression wfa2_static m)
add_test(
  NAME wfa_biwfa_endsfree_endpoint_regression
  COMMAND wfa_biwfa_endsfree_endpoint_regression
)

add_executable(wfa_biwfa_endsfree_score_perf_regression
  tests/wfa_biwfa_endsfree_score_perf_regression.c
)
target_link_libraries(wfa_biwfa_endsfree_score_perf_regression wfa2_static m)
add_test(
  NAME wfa_biwfa_endsfree_score_regression
  COMMAND wfa_biwfa_endsfree_score_perf_regression
)
if(WFA_BUILD_PERF_TESTS)
  add_test(
    NAME wfa_biwfa_endsfree_score_perf_regression
    COMMAND ${CMAKE_COMMAND} -E env WFA_RUN_PERF_REGRESSION=1 $<TARGET_FILE:wfa_biwfa_endsfree_score_perf_regression>
  )
  set_tests_properties(wfa_biwfa_endsfree_score_perf_regression
    PROPERTIES LABELS perf)
endif()

add_executable(wfa_singletrack_heuristic_api_regression
  tests/wfa_singletrack_heuristic_api_regression.c
)
target_link_libraries(wfa_singletrack_heuristic_api_regression wfa2_static m)
add_test(
  NAME wfa_singletrack_heuristic_api_regression
  COMMAND wfa_singletrack_heuristic_api_regression
)

add_executable(wfa_cpp_singletrack_regression
  tests/wfa_cpp_singletrack_regression.cpp
  bindings/cpp/WFAligner.cpp
)
target_link_libraries(wfa_cpp_singletrack_regression wfa2_static m)
add_test(
  NAME wfa_cpp_singletrack_regression
  COMMAND wfa_cpp_singletrack_regression
)

add_executable(wfa_cpp_upstream_regression
  tests/wfa_cpp_upstream_regression.cpp
  bindings/cpp/WFAligner.cpp
)
target_link_libraries(wfa_cpp_upstream_regression wfa2_static m)
add_test(
  NAME wfa_cpp_upstream_regression
  COMMAND wfa_cpp_upstream_regression
)

add_executable(wfa_banded_biwfa_regression
  tests/wfa_banded_biwfa_regression.c
)
target_link_libraries(wfa_banded_biwfa_regression wfa2_static m)
add_test(
  NAME wfa_banded_biwfa_regression
  COMMAND wfa_banded_biwfa_regression
)

# ---- Install

# Do not install anything when used with FetchContent
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  install(TARGETS wfa2_static ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )
  if (BUILD_SHARED_LIB)
      install(TARGETS wfa2 ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )
  endif ()

  install(FILES ${INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/wavefront)
  install(FILES ${UTILS_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/utils)
  install(FILES ${ALIGNMENT_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/alignment)
  install(FILES ${SYSTEM_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/system)

  if (BUILD_CXX_BINDING_LIB)
      install(TARGETS wfa2cpp_static ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
      if (BUILD_SHARED_LIB)
          install(TARGETS wfa2cpp ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
      endif ()
      install(FILES ${CPP_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/bindings/cpp)
  endif ()

  configure_file(${CMAKE_SOURCE_DIR}/libwfa2.pc.in ${CMAKE_BINARY_DIR}/libwfa2.pc @ONLY)
  install(FILES ${CMAKE_BINARY_DIR}/libwfa2.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
endif()
