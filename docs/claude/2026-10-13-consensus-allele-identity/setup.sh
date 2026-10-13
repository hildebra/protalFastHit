#!/bin/bash
# The clones of the comparison in ~/idtest, each a clone of the checkout at 49cd583, built (protal, simulate_metagenomes;
# new: every target, for the tests):
#   base    as committed
#   idonly  + the allele identities (the code of f2b08ba)
#   new     + the allele identities and the consensus tie-break (minimax_tiebreak.patch)
# (In the run, new also carried the map-folder test fix of f2b08ba.) usage: setup.sh [clone ...]     (WSL; 4 cores)
set -eu
CLONES=${*:-base new idonly}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/idtest
CODE="src/Profiling/Profiler.h scripts/model_features.py scripts/test_model_pmml.py tests/test_StrainAlleles.cpp tests/test_Parsing.cpp"
mkdir -p $W/logs
for v in $CLONES; do
  D=$W/$v
  [ -d $D/.git ] || git clone --quiet $REPO $D
  cd $D
  git fetch --quiet origin
  git checkout --quiet -f 49cd583
  git checkout --quiet -- .
  case $v in
    idonly) git show f2b08ba -- $CODE | git apply ;;
    new) git show f2b08ba -- $CODE | git apply; git apply $HERE/minimax_tiebreak.patch ;;
  esac
  git diff --name-only | xargs -r touch  # newer than any object built before
  echo "$v: $(git log -1 --format=%h) + $(git diff --name-only | paste -sd' ')"
  if [ ! -f build/build.ninja ]; then
    TESTS=OFF; [ $v = new ] && TESTS=ON
    cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=$TESTS > $W/logs/configure_$v.log
  fi
  target="--target protal simulate_metagenomes"; [ $v = new ] && target=""
  taskset -c 0-3 nice -n 5 cmake --build build -j 4 $target > $W/logs/build_$v.log 2>&1 || { tail -40 $W/logs/build_$v.log; exit 1; }
  ./build/protal --version | tail -1
done
# Then the reduced database's inputs: ../2026-10-13-mock-local-test/make_inputs.sh <a clone> ~/idtest/inputs; then run_all.sh.
