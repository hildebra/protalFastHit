#!/bin/bash
# ASan/UBSan Debug build of protal at the strain-build commit, in ~/audit6/io/src (never in ~/strain-build).
set -eu
mkdir -p ~/audit6/io
rsync -a --exclude /build --exclude '/build*/' --exclude /data ~/strain-build/src/ ~/audit6/io/src/
cd ~/audit6/io/src
git log --oneline -1 || true
FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"
cmake -S . -B build-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPROTAL_BUILD_TESTS=ON \
  -DCMAKE_CXX_FLAGS="$FLAGS" -DCMAKE_C_FLAGS="$FLAGS" -DCMAKE_EXE_LINKER_FLAGS="$FLAGS" > ~/audit6/io/asan_configure.log 2>&1
nice -n 5 cmake --build build-asan -j 2 --target protal protal_tests > ~/audit6/io/asan_build.log 2>&1
echo BUILD_DONE
ls -la build-asan/bin build-asan/src 2>/dev/null | head
find build-asan -maxdepth 3 -name protal -type f -o -maxdepth 3 -name protal_tests -type f
