# ISA-L (igzip), protal's gzip library, built for every binary (protal and the static binaries alike)
# as the imported static library `protal_isal` (libisal.a), with its headers in
# PROTAL_ISAL_INCLUDE_DIR (isa-l/*.h, as protal includes them; rapidgzip includes them without the
# isa-l/ prefix). No ISA-L of the system is used.
#
# lib/isa-l is ISA-L 2.32.1 from
#   https://github.com/intel/isa-l/archive/refs/tags/v2.32.1.tar.gz
#   (sha256 d9f7179ab0e14a3db9b610fac22793854a1435e8423ec9ce07f4cbedc5f92f5e)
# without its doc/, examples/, programs/, tests/ and tools/ directories, with the changes rapidgzip
# (lib/rapidgzip.cmake) needs from its fork of ISA-L (https://github.com/mxmlnkn/isa-l, as rapidgzip
# 0.16.0 vendors it), marked by comments starting "protal":
#  1. igzip/igzip_inflate.c, include/igzip_lib.h: stopping points (inflate_state's points_to_stop_at,
#     stopped_at, tmp_out_stopped_at and btype; the fork's commit 18660e0, "Add stopping points to
#     isal_inflate"). isal_inflate stops where asked (after a stream header, a block header, a block, a
#     stream), so that rapidgzip can hand the rest of a chunk to ISA-L at a deflate block. Upstream had
#     since added a second check for corrupt data after the copy into tmp_out; it now ends as the fork's
#     do, through the restore of the input. With no stopping point asked for (isal_inflate_init sets
#     none) isal_inflate inflates as upstream's.
#  2. igzip/igzip_inflate.c, include/igzip_lib.h, igzip/huff_codes.h: set_codes,
#     set_and_expand_lit_len_huffcode, make_inflate_huff_code_lit_len and make_inflate_huff_code_dist
#     public, their counts const, and struct huff_code in igzip_lib.h (the fork's 68027bb and part of its
#     3d91c4d): rapidgzip's own deflate decoder builds its Huffman tables with them.
#  3. igzip/igzip_inflate.c, include/igzip_lib.h: isal_inflate_set_dict's dictionary const (part of the
#     fork's 3d91c4d).
# Deflate is untouched: protal's gzip output (BGZF .sam.gz, simulated reads) is as with 2.32.1. The
# fork's commits did not apply to 2.32.1 as patches (upstream has reformatted igzip_inflate.c since the
# fork's base); to update, replace lib/isa-l with the new release, less those directories, and carry the
# three changes over.
#
# Upstream's CMakeLists.txt builds it in the build tree (a static library, position-independent for
# protal's PIE binary, without tests, shim or programs), with protal's C compiler and no -march: ISA-L
# chooses its SSE/AVX2/AVX-512 code at run time, so protal still runs on any x86-64. Its x86-64 code is
# nasm assembly: the PATH's nasm if it is 2.14.01 or later (ISA-L's minimum; Ubuntu and conda: nasm),
# else nasm 2.16.03, downloaded when the build first needs it (sha256 checked) and built here (a C
# compiler and make are all it needs, ~30 s). Without network access, pass its tarball:
#   -DPROTAL_NASM_URL=<path of nasm-2.16.03.tar.xz>
# (from https://www.nasm.us/pub/nasm/releasebuilds/2.16.03/nasm-2.16.03.tar.xz).

include(ExternalProject)
# Extracted files get the time of the extraction, so a new tarball rebuilds them.
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(PROTAL_NASM_URL "https://www.nasm.us/pub/nasm/releasebuilds/2.16.03/nasm-2.16.03.tar.xz"
        CACHE STRING "nasm 2.16.03 source tarball (URL or local path), built when the PATH has no nasm 2.14.01 or later")
set(PROTAL_NASM_SHA256 1412a1c760bbd05db026b6c0d1657affd6631cd0a63cddb6f73cc6d4aa616148)

set(PROTAL_ISAL_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/isa-l)
set(PROTAL_ISAL_DIR ${CMAKE_BINARY_DIR}/isa-l-install)
set(PROTAL_ISAL_INCLUDE_DIR ${PROTAL_ISAL_DIR}/include)
# The imported target's include directory must exist before ISA-L is built.
file(MAKE_DIRECTORY ${PROTAL_ISAL_INCLUDE_DIR}/isa-l)

# nasm: the PATH's if it is 2.14.01 or later, else built here.
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
    set(PROTAL_NASM_DIR ${CMAKE_BINARY_DIR}/nasm)
    ExternalProject_Add(protal_nasm_build
            URL ${PROTAL_NASM_URL}
            URL_HASH SHA256=${PROTAL_NASM_SHA256}
            ${PROTAL_NASM_TIMES}
            PREFIX ${PROTAL_NASM_DIR}
            BUILD_IN_SOURCE 1
            CONFIGURE_COMMAND ./configure --prefix=${PROTAL_NASM_DIR}/install CC=${CMAKE_C_COMPILER}
            BUILD_COMMAND make nasm
            INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${PROTAL_NASM_DIR}/install/bin
                COMMAND ${CMAKE_COMMAND} -E copy nasm ${PROTAL_NASM_DIR}/install/bin/nasm
            BUILD_BYPRODUCTS ${PROTAL_NASM_DIR}/install/bin/nasm)
    set(PROTAL_ISAL_NASM ${PROTAL_NASM_DIR}/install/bin/nasm)
    set(PROTAL_ISAL_DEPENDS protal_nasm_build)
endif()

# Release whatever protal's build type: the library is only linked, never debugged here. lib, not lib64
# (GNUInstallDirs picks lib64 on Red Hat-like systems). BUILD_ALWAYS: a change in lib/isa-l is rebuilt (the
# step checks its own files; installing an unchanged header leaves its time, so protal is not rebuilt).
ExternalProject_Add(protal_isal_build
        SOURCE_DIR ${PROTAL_ISAL_SOURCE_DIR}
        PREFIX ${CMAKE_BINARY_DIR}/isa-l
        DEPENDS ${PROTAL_ISAL_DEPENDS}
        BUILD_ALWAYS 1
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_INSTALL_PREFIX=${PROTAL_ISAL_DIR}
            -DCMAKE_INSTALL_LIBDIR=lib
            -DCMAKE_INSTALL_INCLUDEDIR=include
            -DCMAKE_ASM_NASM_COMPILER=${PROTAL_ISAL_NASM}
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            -DBUILD_SHARED_LIBS=OFF
            -DISAL_BUILD_TESTS=OFF
            -DISAL_BUILD_PERF_TESTS=OFF
            -DISAL_BUILD_ISAL_SHIM=OFF
            -DISAL_BUILD_IGZIP_CLI=OFF
        BUILD_BYPRODUCTS ${PROTAL_ISAL_DIR}/lib/libisal.a)

add_library(protal_isal STATIC IMPORTED GLOBAL)
set_target_properties(protal_isal PROPERTIES
        IMPORTED_LOCATION ${PROTAL_ISAL_DIR}/lib/libisal.a
        INTERFACE_INCLUDE_DIRECTORIES ${PROTAL_ISAL_INCLUDE_DIR})
add_dependencies(protal_isal protal_isal_build)
message(STATUS "ISA-L: lib/isa-l (2.32.1 with rapidgzip's changes), built into ${PROTAL_ISAL_DIR}")
