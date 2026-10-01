cmake_minimum_required(VERSION 3.22)
project(protal VERSION 0.7.1)

# An optimised build unless another build type is asked for. -DCMAKE_BUILD_TYPE=Debug gives an
# unoptimised build with debug info and assert() checks enabled.
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release, Debug, RelWithDebInfo)" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CONDA_PREFIX: $ENV{CONDA_PREFIX}, PREFIX: $ENV{PREFIX}")

# god-awful hack because it seems to not run correct tests to determine this:
set(__CHAR_UNSIGNED___EXITCODE 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_FLAGS}")

# The generated version header lives in the build tree. Older builds wrote it into src/, where it
# would shadow the current one, so a stale copy there is removed.
configure_file(protal_config.h.in ${PROJECT_BINARY_DIR}/generated/protal_config.h)
include_directories(${PROJECT_BINARY_DIR}/generated)
if (EXISTS ${PROJECT_SOURCE_DIR}/src/protal_config.h)
	file(REMOVE ${PROJECT_SOURCE_DIR}/src/protal_config.h)
	message(STATUS "Removed stale generated src/protal_config.h")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(OPENMP_FLAG "")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set (OpenMP_C "${CMAKE_C_COMPILER}")
	set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
	set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
	set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
	set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
	set (OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
	set (OpenMP_libomp_LIBRARY "omp")
	set (OpenMP_libgomp_LIBRARY "gomp")
	set (OpenMP_libiomp5_LIBRARY "iomp5")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lomp")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	if (CMAKE_DISABLE_FIND_PACKAGE_OpenMP)
		message(STATUS "OpenMP disabled by configuration")
	else()
		find_package(OpenMP REQUIRED)
		set (OPENMP_FLAG "-fopenmp")
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} ${OPENMP_FLAG}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${OPENMP_FLAG}")
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
endif()

# Optimisation level and NDEBUG come from the build type (CMAKE_CXX_FLAGS_<CONFIG>), not from here.
set(PROTAL_COMMON_COMPILE_OPTIONS -pthread)
if (CMAKE_DISABLE_FIND_PACKAGE_OpenMP)
	# Avoid failing builds when OpenMP pragmas are present but support is disabled.
	add_compile_options(-Wno-error=unknown-pragmas)
else()
	list(APPEND PROTAL_COMMON_COMPILE_OPTIONS -Werror=unknown-pragmas)
endif()
if (OPENMP_FLAG)
	list(APPEND PROTAL_COMMON_COMPILE_OPTIONS ${OPENMP_FLAG})
endif()

add_library(isa_baseline INTERFACE)
target_compile_options(isa_baseline INTERFACE -march=x86-64 -mtune=generic)

add_library(isa_avx2 INTERFACE)
target_compile_options(isa_avx2 INTERFACE -march=x86-64-v3 -mtune=generic)

include_directories(src/ lib/)

add_executable(protal
        src/main.cpp)
target_compile_options(protal PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})

add_executable(protal_static
        src/main.cpp)
target_compile_options(protal_static PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})
target_link_options(protal_static PRIVATE -static -static-libgcc -static-libstdc++)
set_target_properties(protal_static PROPERTIES OUTPUT_NAME "protal_${PROJECT_VERSION}_static")

add_executable(protal_avx2
        src/main.cpp)
target_compile_options(protal_avx2 PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})

add_executable(simulate_metagenomes
        src/simulate_metagenomes_main.cpp)
target_compile_options(simulate_metagenomes PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})

add_executable(simulate_metagenomes_static
        src/simulate_metagenomes_main.cpp)
target_compile_options(simulate_metagenomes_static PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})
target_link_options(simulate_metagenomes_static PRIVATE -static -static-libgcc -static-libstdc++)


# Link to cmake file in src
target_link_libraries(protal PRIVATE protal_lib dl isa_baseline)
target_link_libraries(protal_static PRIVATE protal_lib_static dl isa_baseline)
target_link_libraries(protal_avx2 PRIVATE protal_lib dl isa_avx2)
target_link_libraries(simulate_metagenomes PRIVATE protal_lib dl isa_baseline)
target_link_libraries(simulate_metagenomes_static PRIVATE protal_lib_static dl isa_baseline)

add_subdirectory(src/)

# Unit tests (GoogleTest): cmake -DPROTAL_BUILD_TESTS=ON ..., then ctest --test-dir <build dir>
option(PROTAL_BUILD_TESTS "Build the GoogleTest unit tests" OFF)
if (PROTAL_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests/)
endif()
