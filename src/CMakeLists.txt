project(protal_lib)
# C++ standard, OpenMP and the per-build-type flags (optimisation, NDEBUG) come from the top level.

# CHECK DEBUG
if(DEFINED ENV{GCC} AND NOT "$ENV{GCC}" STREQUAL "")
    message(STATUS "Using GCC-prefixed binutils from Conda: $ENV{GCC}")
    message(STATUS "Using GCC Prefix: ${GCC_PREFIX}")
    message(STATUS "PREFIX PATH: ${PREFIX}")
    SET(CMAKE_AR  "$ENV{GCC}-ar")
    SET(CMAKE_RANLIB  "$ENV{GCC}-ranlib")
else()
    message(STATUS "Using system default ar/ranlib")
    set(CMAKE_AR "ar")
    set(CMAKE_RANLIB "ranlib")
endif()


SET(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> rcs <TARGET> <LINK_FLAGS> <OBJECTS>")
SET(CMAKE_C_ARCHIVE_FINISH   true)

# Prefer shared libraries (avoid forcing static)
# set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")

message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")



#set(CMAKE_PREFIX_PATH "/usr/users/QIB_fr017/fritsche/miniconda3/envs/gcconda/" ${CMAKE_PREFIX_PATH})

# This creates the target and must be specified before using target_include etc
set(PROTAL_LIB_SOURCES
        RunProtal.h
        Options.h
        Build.h
        Classify.h
        Utilities/Benchmark.h
        Utilities/Utilities.h
        Utilities/ProgressBar.h
        Utilities/xxhash64.h
        Utilities/BinaryClassifierEvaluator.h
        Utilities/LineSplitter.h
        Utilities/Compressor.h
        Utilities/RunStatus.h
        Utilities/BuildInfo.h
        Utilities/BuildInfo.cpp
        Utilities/ReferenceFingerprint.h
        SequenceUtils/SeqReader.h
        SequenceUtils/FastaBatches.h
        SequenceUtils/KmerIterator.h
        SequenceUtils/KmerUtils.h
        SequenceUtils/Minimizer.h
        Hash/Seedmap.h
        Hash/KmerPutter.h
        Hash/KmerLookup.h
        Profiling/Profiler.h
        Profiling/Strain.h
        Profiling/ProfilerDefinitions.h
        Profiling/SampleContext.h
        SequenceUtils/GenomeLoader.h
        SequenceUtils/PackedSequence.h
        IO/AlignmentOutputHandler.h
        Alignment/FastAlignment.h
        Alignment/InternalReadAlignment.h
        Alignment/ScoreAlignments.h
        Alignment/AlignmentUtils.h
        Core/CoreBenchmark.h
        IO/FastxReader.h
        IO/FastxReader.cpp
        IO/ThreadedGzStream.h
        IO/Bgzf.h
        IO/SamFile.h
        IO/SamHandler.h
        IO/BufferedOutput.h
        IO/BufferedOutput.cpp
        Core/SeedingStrategy.h
        Core/ChainingStrategy.h
        Core/ChainAnchorFinder.h
        Core/AlignmentStrategy.h
        Core/LongReads.h
        VarkitInterface/VarkitInterface.h
        Taxonomy/Taxonomy.h
        Taxonomy/Taxonomy.cpp
        Statistics.h
        Constants.h
        ReadType.h
        SNPs/SNP.h
        SNPs/SNPUtils.h
        SequenceUtils/SequenceRange.h
        SequenceUtils/SequenceRange.cpp
        SequenceUtils/SequenceRangeHandler.h
        SequenceUtils/SequenceRangeHandler.cpp
        SequenceUtils/Variant.h
        SequenceUtils/VariantHandler.h
        IO/TaxonStatisticsOutput.h
        Utilities/Matrix.h
        Alignment/WFA2Wrapper2.h
        RandomForest/MetagenomeSimulator.h
        RandomForest/MetagenomeSimulator.cpp
        RandomForest/MetagenomeTypes.h
        RandomForest/CommunityProfileDesigner.h
        RandomForest/CommunityProfileDesigner.cpp
        RandomForest/ArtIlluminaWrapper.h
        RandomForest/ArtIlluminaWrapper.cpp
)

add_library(protal_lib
        ${PROTAL_LIB_SOURCES})

add_library(protal_lib_static STATIC
        ${PROTAL_LIB_SOURCES})

target_include_directories(protal_lib
        PUBLIC
        ../lib
        IO
        Hash
        Core
        VarkitInterface
        Alignment
        Utilities
        Profiling
        SequenceUtils
        RandomForest
        Taxonomy
        SNPs)

target_include_directories(protal_lib_static
        PUBLIC
        ../lib
        IO
        Hash
        Core
        VarkitInterface
        Alignment
        Utilities
        Profiling
        SequenceUtils
        RandomForest
        Taxonomy
        SNPs)


set_target_properties(protal_lib PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(protal_lib_static PROPERTIES LINKER_LANGUAGE CXX)

target_link_directories(protal_lib
        PUBLIC
        ../IO
        ../lib
        ../Hash
        ../Core
        ../VarkitInterface
        ../Alignment
        ../Utilities
        ../Profiling
        # ../Profiler
        ../SequenceUtils
        ../Taxonomy
        ../SNPs
        ../cPMML
)

target_link_directories(protal_lib_static
        PUBLIC
        ../IO
        ../lib
        ../Hash
        ../Core
        ../VarkitInterface
        ../Alignment
        ../Utilities
        ../Profiling
        ../SequenceUtils
        ../Taxonomy
        ../SNPs
        ../cPMML
)

# target_include_directories(protal_lib
#         PUBLIC
#         ../../lib
#         ../IO
#         ../Hash
#         ../Core
#         ../VarkitInterface
#         ../Alignment
#         ../Utils
#         ../Profiling
#         ../SequenceUtils
#         ../Taxonomy
#         ../SNPs)


# zlib-ng (vendored in lib/zlib-ng, built static): gzip that is not BGZF, gzstream. No system zlib.
include(${PROJECT_SOURCE_DIR}/../lib/zlib-ng.cmake)

# zstd: compressed database files (index.prx.zst, reference.fna.zst). Ubuntu: libzstd-dev,
# conda: zstd.
find_path(ZSTD_INCLUDE_DIR zstd.h REQUIRED)
find_library(ZSTD_LIBRARY NAMES zstd REQUIRED)
get_filename_component(ZSTD_LIBRARY_DIR ${ZSTD_LIBRARY} DIRECTORY)
message(STATUS "zstd: ${ZSTD_LIBRARY} (${ZSTD_INCLUDE_DIR})")
target_include_directories(protal_lib PUBLIC ${ZSTD_INCLUDE_DIR})

# libdeflate: gzip as protal writes it (BGZF: .sam.gz, the simulator's reads) and reads BGZF input.
# Ubuntu: libdeflate-dev, conda: libdeflate.
find_path(LIBDEFLATE_INCLUDE_DIR libdeflate.h REQUIRED)
find_library(LIBDEFLATE_LIBRARY NAMES deflate REQUIRED)
get_filename_component(LIBDEFLATE_LIBRARY_DIR ${LIBDEFLATE_LIBRARY} DIRECTORY)
message(STATUS "libdeflate: ${LIBDEFLATE_LIBRARY} (${LIBDEFLATE_INCLUDE_DIR})")
target_include_directories(protal_lib PUBLIC ${LIBDEFLATE_INCLUDE_DIR})

# The static binaries link zstd and libdeflate statically: the system's libzstd.a and libdeflate.a
# (by name, next to the libraries found above), or with PROTAL_STATIC_FETCH_DEPS (`just static`)
# releases downloaded and built in the build tree (lib/static-deps.cmake).
option(PROTAL_STATIC_FETCH_DEPS "Download and build zstd and libdeflate for the static binaries" OFF)
if (PROTAL_STATIC_FETCH_DEPS)
    include(${PROJECT_SOURCE_DIR}/../lib/static-deps.cmake)
    message(STATUS "Static binaries: zstd and libdeflate built in ${PROTAL_STATIC_DEPS_DIR}")
    # BEFORE: their headers, not the system's, which a conda include directory may also hold.
    target_include_directories(protal_lib_static BEFORE PUBLIC ${PROTAL_STATIC_DEPS_INCLUDE_DIR})
    add_dependencies(protal_lib_static protal_static_zstd_build protal_static_deflate_build)
    set(PROTAL_STATIC_ZSTD protal_static_zstd)
    set(PROTAL_STATIC_DEFLATE protal_static_deflate)
else()
    target_include_directories(protal_lib_static PUBLIC ${ZSTD_INCLUDE_DIR} ${LIBDEFLATE_INCLUDE_DIR})
    target_link_directories(protal_lib_static PUBLIC ${ZSTD_LIBRARY_DIR} ${LIBDEFLATE_LIBRARY_DIR})
    set(PROTAL_STATIC_ZSTD zstd)
    set(PROTAL_STATIC_DEFLATE deflate)
endif()

target_include_directories(protal_lib PUBLIC ../lib/robin)
target_include_directories(protal_lib PUBLIC ../lib/tsl)
target_include_directories(protal_lib_static PUBLIC ../lib/robin)
target_include_directories(protal_lib_static PUBLIC ../lib/tsl)

# Add-Subdirectory with ${CMAKE_CURRENT_SOURCE_DIR}/


include(${PROJECT_SOURCE_DIR}/../lib/wfa2-lib.cmake)
add_subdirectory(${PROJECT_SOURCE_DIR}/../lib/gzstream gzstream)
# protal scores taxa within its own threads (one per sample); cPMML's default of 4 threads per
# score would oversubscribe the cores and ignore -t.
set(THREADS 1)
add_subdirectory(${PROJECT_SOURCE_DIR}/../lib/cPMML cPMML)

target_link_libraries(protal_lib PRIVATE isa_baseline)
target_link_libraries(protal_lib_static PRIVATE isa_baseline)

# Create an external library target
#add_library(cPMML STATIC IMPORTED)
#
## Set the library location (will be set after build)
#set_target_properties(cPMML PROPERTIES
#        IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/src/cPMML/libcPMML.a
#)
## Ensure the library is built before linking
#add_dependencies(cPMML protal_lib)

#find_package(cPMML REQUIRED)


#set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")  # Set to OFF for static build


#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#find_library(CPMML_LIBRARY
#        NAMES libcPMML_lib.a
#        PATHS
#        ${CMAKE_BINARY_DIR}
#        NO_DEFAULT_PATH
#)

#set(CPMML_LIBRARY "${CMAKE_BINARY_DIR}/src/cPMML/libcPMML.a")



message("ProjectSource: ${PROJECT_SOURCE_DIR}")
message("CPMML_LIBRARY: ${CPMML_LIBRARY}")
message("CMAKE PREFIX PATH: ${CMAKE_PREFIX_PATH}")
message("Find cPMML: ${PROJECT_SOURCE_DIR}/../")
message("Build dir: ${CMAKE_BINARY_DIR}/src/cPMML/")
message("CPMML: ${cPMML}")


# Add this before target_link_libraries
#find_library(RT_LIBRARY rt)
#message("RT LIB: ${RT_LIBRARY}")

# Necessary to link to stdc++fs (filesystem experimental support)
target_link_libraries(protal_lib
        LINK_PUBLIC
        #${RT_LIBRARY}
        cPMML
        stdc++fs
        wfa_lib
        gzstream_lib
        zlib-ng
        ${ZSTD_LIBRARY}
        ${LIBDEFLATE_LIBRARY}
)

target_link_libraries(protal_lib_static
        LINK_PUBLIC
        cPMML
        stdc++fs
        wfa_lib
        gzstream_lib_static
        zlib-ng
        ${PROTAL_STATIC_ZSTD}
        ${PROTAL_STATIC_DEFLATE}
)

message("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")


message("Ranlib: ${CMAKE_RANLIB}")
message("Ranlib: ${CMAKE_AR}")
message("GCC: $ENV{GCC}")
