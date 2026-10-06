#!/bin/bash
# ISA-L 2.32.1 CMake build: static only, Release, no tests/perf/shim/CLI.
# nasm is NOT on PATH; it is handed over with CMAKE_ASM_NASM_COMPILER.
cd ~/isal
echo "nasm on PATH: $(command -v nasm || echo none)"

# 1. Without the compiler variable: does configure fail?
rm -rf isal-cmake-nonasm && mkdir isal-cmake-nonasm
( cd isal-cmake-nonasm && cmake ../isa-l-2.32.1 -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
    -DISAL_BUILD_TESTS=OFF -DISAL_BUILD_PERF_TESTS=OFF -DISAL_BUILD_ISAL_SHIM=OFF > cmake.log 2>&1; echo "configure without nasm: rc=$?"; grep -iE 'nasm|error' cmake.log | head -8 )

# 2. With the compiler variable
rm -rf isal-cmake-build && mkdir isal-cmake-build && cd isal-cmake-build
t0=$(date +%s.%N)
cmake ../isa-l-2.32.1 -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
    -DISAL_BUILD_TESTS=OFF -DISAL_BUILD_PERF_TESTS=OFF -DISAL_BUILD_ISAL_SHIM=OFF -DISAL_BUILD_IGZIP_CLI=OFF \
    -DCMAKE_ASM_NASM_COMPILER=$HOME/isal/prefix/bin/nasm \
    -DCMAKE_INSTALL_PREFIX=$HOME/isal/prefix > cmake.log 2>&1
echo "configure rc=$?"
t1=$(date +%s.%N)
cmake --build . -j6 > build.log 2>&1
echo "build rc=$?"
t2=$(date +%s.%N)
cmake --install . > install.log 2>&1
echo "install rc=$?"
t3=$(date +%s.%N)
awk -v a=$t0 -v b=$t1 -v c=$t2 -v d=$t3 'BEGIN{printf "configure %.1f s, build -j6 %.1f s, install %.1f s\n", b-a, c-b, d-c}'
grep -A12 'configuration summary' cmake.log
grep -i nasm CMakeCache.txt | grep -v '^//' | head
grep -c '\.asm' build.log
cat install.log
