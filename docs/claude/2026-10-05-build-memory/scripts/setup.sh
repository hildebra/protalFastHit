#!/usr/bin/env bash
# Builds protal from a commit in WSL from a Linux-side copy (git archive), and converts the 900-species
# tuning world once: setup.sh [COMMIT=HEAD] [TREE_NAME=base]
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal
WORK=$HOME/buildmem
commit=${1:-HEAD}; name=${2:-base}
mkdir -p $WORK && cd $WORK || exit 1
if [ ! -d $name/src ]; then mkdir -p $name && (cd $SRC && git archive --format=tar $commit) | tar -x -C $name || exit 1; fi
cmake -S $name -B $name/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF > $name-cmake.log 2>&1 || { tail -20 $name-cmake.log; exit 1; }
nice cmake --build $name/build --target protal -j 6 > $name-build.log 2>&1 || { grep -E "error" -A3 $name-build.log | head -40; exit 1; }
ls -la $name/build/protal
if [ ! -d db900.pristine ]; then
  python3 $name/scripts/mini_db/gtdb_to_protal_db.py --gtdb ~/tune/world --outdir db900.pristine -t 6 > db900.convert.log 2>&1 || { tail -20 db900.convert.log; exit 1; }
fi
du -sh db900.pristine
