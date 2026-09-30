# zlib-ng, protal's gzip library (gzip input that is not BGZF, gzstream), built as the static
# library target `zlib-ng` (libz-ng.a) with its native API: zlib-ng.h, zng_ functions.
#
# lib/zlib-ng is zlib-ng 2.3.3, unmodified, from
#   https://github.com/zlib-ng/zlib-ng/archive/refs/tags/2.3.3.tar.gz
#   (sha256 f9c65aa9c852eb8255b636fd9f07ce1c406f061ec19a2e7d508b318ca0c907d1, the release
#   conda-forge packages) without its test/ and doc/ directories.
# To update: replace lib/zlib-ng with the new release, less test/ and doc/, and update the
# version and sha256 here.
#
# Upstream's CMakeLists.txt is used as it is: it detects the compiler's SIMD support and builds
# every variant (SSE2 to AVX-512 on x86-64, NEON on ARM, ...), chosen at run time, so the baseline
# protal binary still runs on any x86-64. Being static and linked in, it is no run-time dependency,
# and the static build needs nothing more. The native API (not ZLIB_COMPAT) cannot be mixed up
# with a system zlib: a leftover #include <zlib.h> or gzopen() fails to compile or link.
#
# The settings below apply to this subdirectory only (they are restored after it): a static
# library only, no tests (they would fetch googletest), no install rules, no zlib-named aliases.

set(PROTAL_ZNG_SETTINGS
        BUILD_SHARED_LIBS OFF
        BUILD_TESTING OFF
        CMAKE_POSITION_INDEPENDENT_CODE ON
        ZLIB_COMPAT OFF
        ZLIB_ALIASES OFF
        WITH_GZFILEOP ON
        WITH_NATIVE_INSTRUCTIONS OFF
        WITH_RUNTIME_CPU_DETECTION ON
        SKIP_INSTALL_ALL ON)
set(PROTAL_ZNG_VARIABLES "")
while (PROTAL_ZNG_SETTINGS)
    list(POP_FRONT PROTAL_ZNG_SETTINGS variable value)
    list(APPEND PROTAL_ZNG_VARIABLES ${variable})
    if (DEFINED ${variable})
        set(PROTAL_ZNG_SAVED_${variable} "${${variable}}")
    endif()
    set(${variable} ${value})
endwhile()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/zlib-ng ${CMAKE_BINARY_DIR}/zlib-ng EXCLUDE_FROM_ALL)
# zlib-ng.h declares gzFile and the zng_gz* functions only with WITH_GZFILEOP, which upstream's
# pkg-config file passes on to users but its CMake target does not.
target_compile_definitions(zlib-ng INTERFACE WITH_GZFILEOP)

foreach (variable ${PROTAL_ZNG_VARIABLES})
    if (DEFINED PROTAL_ZNG_SAVED_${variable})
        set(${variable} "${PROTAL_ZNG_SAVED_${variable}}")
        unset(PROTAL_ZNG_SAVED_${variable})
    else()
        unset(${variable})
    endif()
endforeach()
