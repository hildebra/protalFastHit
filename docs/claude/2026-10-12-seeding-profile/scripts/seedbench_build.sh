#!/usr/bin/env bash
# Builds seedbench (experiment only) in a protal source tree SRC with protal's own flags and libraries: the file is
# copied to SRC/misc/seedbench.cpp and a target appended to SRC/CMakeLists.txt (once), then BUILD is configured as
# protal's Release build and the target built.
#
#   bash seedbench_build.sh SEEDBENCH_CPP SRC BUILD [JOBS]
set -euo pipefail
cpp=$1 src=$2 build=$3 jobs=${4:-4}
mkdir -p "$src/misc"
cp "$cpp" "$src/misc/seedbench.cpp"
if ! grep -q "add_executable(seedbench" "$src/CMakeLists.txt"; then
  cat >> "$src/CMakeLists.txt" <<'EOF'

# seedbench: experiment only (docs/claude, the seeding profile), not part of protal.
add_executable(seedbench misc/seedbench.cpp)
target_compile_options(seedbench PRIVATE ${PROTAL_COMMON_COMPILE_OPTIONS})
target_link_libraries(seedbench PRIVATE protal_lib dl isa_baseline)
if (TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(seedbench PRIVATE OpenMP::OpenMP_CXX)
endif()
EOF
fi
if [ ! -f "$build/build.ninja" ] && [ ! -f "$build/Makefile" ]; then
  generator=(); command -v ninja > /dev/null && generator=(-G Ninja)
  cmake -S "$src" -B "$build" "${generator[@]}" -DCMAKE_BUILD_TYPE=Release ${EXTRA_CMAKE:-}
fi
cmake --build "$build" -j "$jobs" --target seedbench
ls -la "$build/seedbench"
