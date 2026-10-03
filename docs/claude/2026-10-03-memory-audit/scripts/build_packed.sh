#!/usr/bin/env bash
# Builds the packed index change and its baseline in WSL from Linux-side copies: WORK/base is the
# committed HEAD, WORK/src is HEAD with this change's files laid over it (other sessions' uncommitted
# work in the checkout stays out). Then the unit tests of the change.
SRC=${SRC:-/mnt/c/Users/hildebra/Documents/locDev/protal}
WORK=${WORK:-$HOME/protal-pack}
FILES="src/Hash/Seedmap.h src/Hash/KmerLookup.h src/SequenceUtils/GenomeLoader.h src/RunProtal.h tests/test_PackedIndex.cpp"
mkdir -p "$WORK" && cd "$WORK" || exit 1
for tree in base src; do
  if [ ! -d "$tree/src" ]; then
    mkdir -p "$tree" && (cd "$SRC" && git archive --format=tar HEAD) | tar -x -C "$tree" || exit 1
  fi
done
for f in $FILES; do cp "$SRC/$f" "src/$f" || exit 1; done
grep -q test_PackedIndex src/tests/CMakeLists.txt || sed -i 's/^        test_Index.cpp$/        test_Index.cpp\n        test_PackedIndex.cpp/' src/tests/CMakeLists.txt
cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > cmake.log 2>&1 || { tail -20 cmake.log; exit 1; }
cmake --build build --target protal protal_tests simulate_metagenomes -j "$(nproc)" > build.log 2>&1 || { grep -E "error" -A4 build.log | head -80; exit 1; }
if [ ! -x base/build/protal ]; then
  cmake -S base -B base/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF > base-cmake.log 2>&1 || { tail -20 base-cmake.log; exit 1; }
  cmake --build base/build --target protal -j "$(nproc)" > base-build.log 2>&1 || { grep -E "error" -A3 base-build.log | head -40; exit 1; }
fi
echo "built: $(ls -la build/protal base/build/protal | awk '{print $5, $9}')"
ctest --test-dir build -R "PackedIndex|IndexLoad|IndexHeader|UniqueKmers|GeneWindows" --output-on-failure 2>&1 | tail -25
