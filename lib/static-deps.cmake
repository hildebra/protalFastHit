# zstd and ISA-L for the static binaries (protal_static, simulate_metagenomes_static), built in
# the build tree, with -DPROTAL_STATIC_FETCH_DEPS=ON (as `just static` sets): the system then
# needs no libzstd.a or libisal.a. Defines the imported static libraries `protal_static_zstd`
# and `protal_static_isal`, and their headers' directories PROTAL_STATIC_DEPS_INCLUDE_DIR (zstd)
# and PROTAL_STATIC_ISAL_INCLUDE_DIR. Without an ISA-L on the system, protal and the tests link
# this one too.
#
# The release tarballs are downloaded when the build first needs them (the sha256 is checked):
#   https://github.com/facebook/zstd/releases/download/v1.5.7/zstd-1.5.7.tar.gz
#   https://github.com/intel/isa-l/archive/refs/tags/v2.32.1.tar.gz
#   https://www.nasm.us/pub/nasm/releasebuilds/2.16.03/nasm-2.16.03.tar.xz (only without a nasm of
#   2.14.01 or later on the PATH: ISA-L's x86-64 code is nasm assembly)
# Without network access, download them elsewhere and pass their paths as -DPROTAL_ZSTD_URL=...,
# -DPROTAL_ISAL_URL=... and -DPROTAL_NASM_URL=... . To update: change the versions and sha256 here.
# (ISA-L publishes no release tarballs; GitHub's tag archive is the release's source.)
#
# zstd and ISA-L are built with their own CMakeLists.txt, as static libraries only, with protal's
# compilers and no -march: they choose their SSE/AVX2/AVX-512 code at run time, so the static
# binaries still run on any x86-64. zstd keeps its multithreading support (protal compresses with
# several workers). ISA-L is built without its tests (which need zlib), programs and zlib shim.

include(ExternalProject)
# Extracted files get the time of the extraction, so a new tarball rebuilds them.
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(PROTAL_ZSTD_URL "https://github.com/facebook/zstd/releases/download/v1.5.7/zstd-1.5.7.tar.gz"
        CACHE STRING "zstd 1.5.7 release tarball (URL or local path) for the static binaries")
set(PROTAL_ZSTD_SHA256 eb33e51f49a15e023950cd7825ca74a4a2b43db8354825ac24fc1b7ee09e6fa3)
set(PROTAL_ISAL_URL "https://github.com/intel/isa-l/archive/refs/tags/v2.32.1.tar.gz"
        CACHE STRING "ISA-L 2.32.1 source tarball (URL or local path) for the static binaries")
set(PROTAL_ISAL_SHA256 d9f7179ab0e14a3db9b610fac22793854a1435e8423ec9ce07f4cbedc5f92f5e)
set(PROTAL_NASM_URL "https://www.nasm.us/pub/nasm/releasebuilds/2.16.03/nasm-2.16.03.tar.xz"
        CACHE STRING "nasm 2.16.03 source tarball (URL or local path), built when the PATH has no nasm 2.14.01 or later")
set(PROTAL_NASM_SHA256 1412a1c760bbd05db026b6c0d1657affd6631cd0a63cddb6f73cc6d4aa616148)

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

# nasm for ISA-L: the PATH's if it is 2.14.01 or later (ISA-L's minimum), else built here (a C
# compiler and make are all it needs, ~30 s).
find_program(PROTAL_NASM_EXECUTABLE nasm)
set(PROTAL_NASM_VERSION "")
if (PROTAL_NASM_EXECUTABLE)
    execute_process(COMMAND ${PROTAL_NASM_EXECUTABLE} -v OUTPUT_VARIABLE PROTAL_NASM_V ERROR_QUIET)
    if (PROTAL_NASM_V MATCHES "version ([0-9]+\\.[0-9]+(\\.[0-9]+)?)")
        set(PROTAL_NASM_VERSION ${CMAKE_MATCH_1})
    endif()
endif()
if (PROTAL_NASM_VERSION AND NOT PROTAL_NASM_VERSION VERSION_LESS 2.14.01)
    message(STATUS "nasm for ISA-L: ${PROTAL_NASM_EXECUTABLE} (${PROTAL_NASM_VERSION})")
    set(PROTAL_ISAL_NASM ${PROTAL_NASM_EXECUTABLE})
    set(PROTAL_ISAL_DEPENDS "")
else()
    message(STATUS "nasm for ISA-L: none of 2.14.01 or later on the PATH, nasm 2.16.03 is built")
    # nasm's Makefile remakes `configure` with autotools when configure.ac looks newer, so its files keep the
    # tarball's times (CMake before 3.24 keeps them anyway); it is built in its source tree, which a new tarball
    # replaces, so nothing old is kept.
    set(PROTAL_NASM_TIMES "")
    if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
        set(PROTAL_NASM_TIMES DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    endif()
    ExternalProject_Add(protal_static_nasm_build
            URL ${PROTAL_NASM_URL}
            URL_HASH SHA256=${PROTAL_NASM_SHA256}
            ${PROTAL_NASM_TIMES}
            PREFIX ${PROTAL_STATIC_DEPS_DIR}/nasm
            BUILD_IN_SOURCE 1
            CONFIGURE_COMMAND ./configure --prefix=${PROTAL_STATIC_DEPS_DIR}/nasm-install CC=${CMAKE_C_COMPILER}
            BUILD_COMMAND make nasm
            INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${PROTAL_STATIC_DEPS_DIR}/nasm-install/bin
                COMMAND ${CMAKE_COMMAND} -E copy nasm ${PROTAL_STATIC_DEPS_DIR}/nasm-install/bin/nasm
            BUILD_BYPRODUCTS ${PROTAL_STATIC_DEPS_DIR}/nasm-install/bin/nasm)
    set(PROTAL_ISAL_NASM ${PROTAL_STATIC_DEPS_DIR}/nasm-install/bin/nasm)
    set(PROTAL_ISAL_DEPENDS protal_static_nasm_build)
endif()

# ISA-L in a prefix of its own: protal_lib may use it too (src/CMakeLists.txt), without zstd's headers.
set(PROTAL_STATIC_ISAL_DIR ${PROTAL_STATIC_DEPS_DIR}/isal-install)
set(PROTAL_STATIC_ISAL_INCLUDE_DIR ${PROTAL_STATIC_ISAL_DIR}/include)
file(MAKE_DIRECTORY ${PROTAL_STATIC_ISAL_INCLUDE_DIR})
ExternalProject_Add(protal_static_isal_build
        URL ${PROTAL_ISAL_URL}
        URL_HASH SHA256=${PROTAL_ISAL_SHA256}
        PREFIX ${PROTAL_STATIC_DEPS_DIR}/isa-l
        DEPENDS ${PROTAL_ISAL_DEPENDS}
        CMAKE_ARGS ${PROTAL_STATIC_DEPS_CMAKE_ARGS}
            -DCMAKE_INSTALL_PREFIX=${PROTAL_STATIC_ISAL_DIR}
            -DCMAKE_ASM_NASM_COMPILER=${PROTAL_ISAL_NASM}
            -DBUILD_SHARED_LIBS=OFF
            -DISAL_BUILD_TESTS=OFF
            -DISAL_BUILD_PERF_TESTS=OFF
            -DISAL_BUILD_ISAL_SHIM=OFF
            -DISAL_BUILD_IGZIP_CLI=OFF
        BUILD_BYPRODUCTS ${PROTAL_STATIC_ISAL_DIR}/lib/libisal.a)

add_library(protal_static_zstd STATIC IMPORTED GLOBAL)
set_target_properties(protal_static_zstd PROPERTIES IMPORTED_LOCATION ${PROTAL_STATIC_DEPS_DIR}/lib/libzstd.a)
add_dependencies(protal_static_zstd protal_static_zstd_build)

add_library(protal_static_isal STATIC IMPORTED GLOBAL)
set_target_properties(protal_static_isal PROPERTIES IMPORTED_LOCATION ${PROTAL_STATIC_ISAL_DIR}/lib/libisal.a)
add_dependencies(protal_static_isal protal_static_isal_build)
