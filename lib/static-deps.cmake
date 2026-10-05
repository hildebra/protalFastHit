# zstd and libdeflate for the static binaries (protal_static, simulate_metagenomes_static), built
# in the build tree, with -DPROTAL_STATIC_FETCH_DEPS=ON (as `just static` sets): the system then
# needs no libzstd.a or libdeflate.a. Defines the imported static libraries `protal_static_zstd`
# and `protal_static_deflate` and PROTAL_STATIC_DEPS_INCLUDE_DIR, their headers.
#
# The release tarballs are downloaded when the build first needs them (the sha256 is checked):
#   https://github.com/facebook/zstd/releases/download/v1.5.7/zstd-1.5.7.tar.gz
#   https://github.com/ebiggers/libdeflate/releases/download/v1.26/libdeflate-1.26.tar.gz
# Without network access, download them elsewhere and pass their paths as -DPROTAL_ZSTD_URL=... and
# -DPROTAL_LIBDEFLATE_URL=... . To update: change the versions and sha256 here.
#
# Both are built with their own CMakeLists.txt, as static libraries only, with protal's compilers
# and no -march: they choose their SSE/AVX2/BMI2 code at run time, so the static binaries still run
# on any x86-64. zstd keeps its multithreading support (protal compresses with several workers).

include(ExternalProject)
# Extracted files get the time of the extraction, so a new tarball rebuilds them.
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(PROTAL_ZSTD_URL "https://github.com/facebook/zstd/releases/download/v1.5.7/zstd-1.5.7.tar.gz"
        CACHE STRING "zstd 1.5.7 release tarball (URL or local path) for the static binaries")
set(PROTAL_ZSTD_SHA256 eb33e51f49a15e023950cd7825ca74a4a2b43db8354825ac24fc1b7ee09e6fa3)
set(PROTAL_LIBDEFLATE_URL "https://github.com/ebiggers/libdeflate/releases/download/v1.26/libdeflate-1.26.tar.gz"
        CACHE STRING "libdeflate 1.26 release tarball (URL or local path) for the static binaries")
set(PROTAL_LIBDEFLATE_SHA256 125856d4656e0feab660f94842f835923410c9281fedbcee64598c918da42b5a)

set(PROTAL_STATIC_DEPS_DIR ${CMAKE_BINARY_DIR}/static-deps)
set(PROTAL_STATIC_DEPS_INCLUDE_DIR ${PROTAL_STATIC_DEPS_DIR}/include)
# Release whatever protal's build type: the libraries are only linked, never debugged here.
# lib, not lib64 (GNUInstallDirs picks lib64 on Red Hat-like systems).
set(PROTAL_STATIC_DEPS_CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_INSTALL_PREFIX=${PROTAL_STATIC_DEPS_DIR}
        -DCMAKE_INSTALL_LIBDIR=lib
        -DCMAKE_INSTALL_INCLUDEDIR=include)

ExternalProject_Add(protal_static_zstd_build
        URL ${PROTAL_ZSTD_URL}
        URL_HASH SHA256=${PROTAL_ZSTD_SHA256}
        PREFIX ${PROTAL_STATIC_DEPS_DIR}/zstd
        SOURCE_SUBDIR build/cmake
        CMAKE_ARGS ${PROTAL_STATIC_DEPS_CMAKE_ARGS}
            -DZSTD_BUILD_STATIC=ON
            -DZSTD_BUILD_SHARED=OFF
            -DZSTD_BUILD_PROGRAMS=OFF
            -DZSTD_BUILD_TESTS=OFF
            -DZSTD_BUILD_CONTRIB=OFF
            -DZSTD_MULTITHREAD_SUPPORT=ON
        BUILD_BYPRODUCTS ${PROTAL_STATIC_DEPS_DIR}/lib/libzstd.a)

ExternalProject_Add(protal_static_deflate_build
        URL ${PROTAL_LIBDEFLATE_URL}
        URL_HASH SHA256=${PROTAL_LIBDEFLATE_SHA256}
        PREFIX ${PROTAL_STATIC_DEPS_DIR}/libdeflate
        CMAKE_ARGS ${PROTAL_STATIC_DEPS_CMAKE_ARGS}
            -DLIBDEFLATE_BUILD_STATIC_LIB=ON
            -DLIBDEFLATE_BUILD_SHARED_LIB=OFF
            -DLIBDEFLATE_BUILD_GZIP=OFF
            -DLIBDEFLATE_BUILD_TESTS=OFF
        BUILD_BYPRODUCTS ${PROTAL_STATIC_DEPS_DIR}/lib/libdeflate.a)

add_library(protal_static_zstd STATIC IMPORTED GLOBAL)
set_target_properties(protal_static_zstd PROPERTIES IMPORTED_LOCATION ${PROTAL_STATIC_DEPS_DIR}/lib/libzstd.a)
add_dependencies(protal_static_zstd protal_static_zstd_build)

add_library(protal_static_deflate STATIC IMPORTED GLOBAL)
set_target_properties(protal_static_deflate PROPERTIES IMPORTED_LOCATION ${PROTAL_STATIC_DEPS_DIR}/lib/libdeflate.a)
add_dependencies(protal_static_deflate protal_static_deflate_build)
