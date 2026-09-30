#!/bin/bash
# ASan/UBSan build (asserts on) of a COPY of the strain-build source, with PROTAL_LR_DEBUG
# diagnostics patched in (patch_debug.py; only active with that environment variable).
set -eu
S=$(dirname "$0")
W=~/audit6/longreads
mkdir -p $W/src
rsync -a --delete --exclude /build ~/strain-build/src/ $W/src/
python3 $S/patch_debug.py $W/src
if [ ! -f $W/build-asan/Makefile ]; then
cmake -S $W/src -B $W/build-asan -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O1 -g" -DCMAKE_C_FLAGS_RELWITHDEBINFO="-O1 -g" \
    -DMY_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
    -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined" \
    -DPROTAL_BUILD_TESTS=ON > $W/asan_configure.log 2>&1
fi
nice -n 10 cmake --build $W/build-asan -j 2 --target protal > $W/asan_build.log 2>&1
echo PROTAL_DONE
[ -n "${WITH_TESTS:-}" ] && nice -n 10 cmake --build $W/build-asan -j 2 --target protal_tests >> $W/asan_build.log 2>&1
echo BUILD_DONE
ls -la $W/build-asan/protal $W/build-asan/tests/protal_tests 2>&1
