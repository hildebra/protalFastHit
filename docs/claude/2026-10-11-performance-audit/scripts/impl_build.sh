#!/usr/bin/env bash
# The implementation's WSL build: ~/impl/src = git archive of HEAD (refreshed when HEAD changes) with this session's
# changed files copied over (rsync --checksum, so that ninja sees every change); ~/impl/build (Release) and, with
# REF=1, ~/impl/ref (HEAD unchanged, for output comparisons). Targets: protal, protal_tests (TARGETS overrides).
# Pinned to cores 0-3.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/impl
targets=${TARGETS:-protal protal_tests}
files=(src/Build.h src/Options.h src/Profiling/Profiler.h src/RunProtal.h src/SequenceUtils/AncestrySites.h src/SequenceUtils/ColumnWeights.h
       src/SequenceUtils/StrainAlleles.h src/Utilities/Database.h src/Utilities/BinaryTable.h src/Utilities/Zstd.h tests/test_ColumnWeights.cpp tests/test_StrainAlleles.cpp tests/test_Database.cpp
       tests/test_ReadEvidence.cpp tests/e2e/test_protal_e2e.py src/main.cpp ${EXTRA:-})
mkdir -p "$W"
launcher=()
command -v ccache > /dev/null && launcher=(-DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DCMAKE_C_COMPILER_LAUNCHER=ccache)
pin=(taskset -c 0-3 nice -n 5)
fresh() {  # dir: a tree of HEAD there, rebuilt from scratch if HEAD changed
  local dir=$1
  if [ "$(cat "$dir/commit" 2> /dev/null)" != "$(cat "$A/archives/head.commit")" ]; then
    rm -rf "$dir"; mkdir -p "$dir/src"
    tar -xf "$A/archives/head.tar" -C "$dir/src"
    cp "$A/archives/head.commit" "$dir/commit"
  fi
}
build() {  # dir name
  local dir=$1 name=$2 start
  start=$(date +%s)
  if "${pin[@]}" cmake -S "$dir/src" -B "$dir/build" -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON "${launcher[@]}" > "$dir/configure.log" 2>&1 &&
     "${pin[@]}" cmake --build "$dir/build" -j 4 --target $targets > "$dir/build.log" 2>&1; then
    echo "OK $name in $(( $(date +%s) - start )) s"
  else
    echo "FAIL $name"; grep -E 'error|Error' -A4 "$dir/build.log" | head -60
    return 1
  fi
}
fresh "$W/work"
for f in "${files[@]}"; do
  mkdir -p "$(dirname "$W/work/src/$f")"
  rsync -a --checksum --no-times "$REPO/$f" "$W/work/src/$f"
done
build "$W/work" work || exit 1
if [ "${REF:-0}" = 1 ]; then fresh "$W/ref"; build "$W/ref" ref; fi
