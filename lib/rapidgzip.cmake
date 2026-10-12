# rapidgzip, for protal the header-only (INTERFACE) target `protal_rapidgzip`: a gzip file that is
# not BGZF inflated on several threads (src/IO/ParallelGzip.cpp, the only file that includes it).
#
# lib/rapidgzip is rapidgzip 0.16.0 (Maximilian Knespel; MIT or Apache-2.0: LICENSE-MIT,
# LICENSE-APACHE), the src/core, src/filereader, src/huffman, src/indexed_bzip2 (its chunk fetcher
# includes the bzip2 decoder) and src/rapidgzip directories of its source distribution, unmodified:
#   https://files.pythonhosted.org/packages/95/9a/d94edac485ade88fbee6864d057eae8a5363bf734da5760f4e99f7a02d94/rapidgzip-0.16.0.tar.gz
#   (sha256 8b124f29bc12de4249ab81e83e5ad35e67742a1a8ff4acb61b74c0d9fda1c14e)
# It is built as its own CMake builds librapidgzip, with LIBRAPIDARCHIVE_WITH_ISAL (ISA-L inflates the
# rest of a chunk once its window is known, and computes the CRC-32s; lib/isa-l carries the changes
# rapidgzip needs from its ISA-L fork) and without rpmalloc (its allocator; FasterVector falls back to
# malloc). Its zlib paths (gzip/zlib.hpp, used without ISA-L) are compiled, not called;
# lib/rapidgzip-zlib/zlib.h gives them protal's zlib-ng. To update: replace the five directories with
# the new release's, and compare lib/isa-l with its ISA-L fork (src/external/isa-l of the source
# distribution).

find_package(Threads REQUIRED)
add_library(protal_rapidgzip INTERFACE)
# SYSTEM: its headers' warnings are not protal's.
target_include_directories(protal_rapidgzip SYSTEM INTERFACE
        ${CMAKE_CURRENT_LIST_DIR}/rapidgzip/src
        ${CMAKE_CURRENT_LIST_DIR}/rapidgzip-zlib
        ${PROTAL_ISAL_INCLUDE_DIR}/isa-l)
target_compile_definitions(protal_rapidgzip INTERFACE LIBRAPIDARCHIVE_WITH_ISAL)
target_link_libraries(protal_rapidgzip INTERFACE protal_isal zlib-ng Threads::Threads)
