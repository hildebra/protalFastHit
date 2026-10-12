# zstd for the static binaries (protal_static, simulate_metagenomes_static), built in the build
# tree, with -DPROTAL_STATIC_FETCH_DEPS=ON (as `just static` sets): the system then needs no
# libzstd.a. Defines the imported static library `protal_static_zstd` and its headers' directory
# PROTAL_STATIC_DEPS_INCLUDE_DIR. (ISA-L is lib/isa-l, built for every binary: lib/isa-l.cmake.)
#
# The release tarball is downloaded when the build first needs it (the sha256 is checked):
#   https://github.com/facebook/zstd/releases/download/v1.5.7/zstd-1.5.7.tar.gz
# Without network access, download it elsewhere and pass its path as -DPROTAL_ZSTD_URL=... . To
# update: change the version and sha256 here.
#
# zstd is built with its own CMakeLists.txt, as a static library only, with protal's C compiler and
# no -march: it chooses its SSE/AVX2 code at run time, so the static binaries still run on any
# x86-64. It keeps its multithreading support (protal compresses with several workers).

include(ExternalProject)
# Extracted files get the time of the extraction, so a new tarball rebuilds them.
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(PROTAL_ZSTD_URL "https://github.com/facebook/zstd/releases/download/v1.5.7/zstd-1.5.7.tar.gz"
        CACHE STRING "zstd 1.5.7 release tarball (URL or local path) for the static binaries")
set(PROTAL_ZSTD_SHA256 eb33e51f49a15e023950cd7825ca74a4a2b43db8354825ac24fc1b7ee09e6fa3)

set(PROTAL_STATIC_DEPS_DIR ${CMAKE_BINARY_DIR}/static-deps)
set(PROTAL_STATIC_DEPS_INCLUDE_DIR ${PROTAL_STATIC_DEPS_DIR}/include)
file(MAKE_DIRECTORY ${PROTAL_STATIC_DEPS_INCLUDE_DIR})
# Release whatever protal's build type: the library is only linked, never debugged here.
# lib, not lib64 (GNUInstallDirs picks lib64 on Red Hat-like systems).
ExternalProject_Add(protal_static_zstd_build
        URL ${PROTAL_ZSTD_URL}
        URL_HASH SHA256=${PROTAL_ZSTD_SHA256}
        PREFIX ${PROTAL_STATIC_DEPS_DIR}/zstd
        SOURCE_SUBDIR build/cmake
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_INSTALL_PREFIX=${PROTAL_STATIC_DEPS_DIR}
            -DCMAKE_INSTALL_LIBDIR=lib
            -DCMAKE_INSTALL_INCLUDEDIR=include
            -DZSTD_BUILD_STATIC=ON
            -DZSTD_BUILD_SHARED=OFF
            -DZSTD_BUILD_PROGRAMS=OFF
            -DZSTD_BUILD_TESTS=OFF
            -DZSTD_BUILD_CONTRIB=OFF
            -DZSTD_MULTITHREAD_SUPPORT=ON
        BUILD_BYPRODUCTS ${PROTAL_STATIC_DEPS_DIR}/lib/libzstd.a)

add_library(protal_static_zstd STATIC IMPORTED GLOBAL)
set_target_properties(protal_static_zstd PROPERTIES IMPORTED_LOCATION ${PROTAL_STATIC_DEPS_DIR}/lib/libzstd.a)
add_dependencies(protal_static_zstd protal_static_zstd_build)
