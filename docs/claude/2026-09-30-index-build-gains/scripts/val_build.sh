#!/bin/bash
# 56217e1 with the build speed-ups only (mine.tar): build with tests, all tests, then index.prx and
# unique_kmers.tsv against 56217e1's binary on the tuning world at 25, 50, 100% of its species.
set -e
S=$(cd "$(dirname "$0")" && pwd)  # this folder; head56.tar and mine.tar (the changed files) next to it, see README
B=$HOME/protal-gains; OLD=$HOME/protal-uniq/build/protal
if [ "$1" != "--skip-build" ]; then
  rm -rf $B; mkdir -p $B/src
  tar -xf $S/head56.tar -C $B/src
  tar -xf $S/mine.tar -C $B/src
  cmake -S $B/src -B $B/build -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $B/configure.log 2>&1
  cmake --build $B/build -j 8 --target protal simulate_metagenomes protal_tests > $B/make.log 2>&1 || { grep -m 20 -B2 -A8 "error" $B/make.log; exit 1; }
  echo "built"
  ctest --test-dir $B/build 2>&1 | tail -n 3 | head -n 1
  cd $B/src
  python3 -m unittest scripts/mini_db/test_mini_db.py 2>&1 | tail -n 1
  PROTAL=$B/build/protal bash scripts/mini_db/build_mini_db.sh data/mini_db > $B/minidb.log 2>&1 && echo "mini db ok"
  PROTAL_TEST_DB=$B/src/data/mini_db/protal_db PROTAL=$B/build/protal SIMULATE=$B/build/simulate_metagenomes \
      python3 -m unittest tests/e2e/test_protal_e2e.py 2>&1 | tail -n 3
fi
W=/tmp/scale; mkdir -p $W
R=$B/src
if [ ! -f $W/in100/reference.fna ]; then
  cut -f2 ~/tune/release_p/*_taxonomy_r226.tsv | sed 's/.*;s__/s__/' | sort -u > $W/species.txt
  python3 - $W <<'EOF'
import random, sys
w = sys.argv[1]
species = open(f"{w}/species.txt").read().split("\n")[:-1]
random.Random(1).shuffle(species)
for keep in (50, 25):
    open(f"{w}/exclude{keep}.txt", "w").write("\n".join(species[len(species) * keep // 100:]) + "\n")
EOF
  for p in 100 50 25; do
    ex=""; [ $p != 100 ] && ex="--exclude_species $W/exclude$p.txt"
    python3 $R/scripts/mini_db/gtdb_to_protal_db.py --gtdb ~/tune/release_p --outdir $W/in$p -t 8 $ex 2> $W/convert$p.log
  done
  echo "inputs converted"
fi
build() {  # label binary size threads
  local D=$W/db_$1; rm -rf $D; cp -r $W/in$3 $D
  /usr/bin/time -f "%e s, %M KB" -o $W/time_$1 $2 --build --no_profile --no_compress -t $4 --db $D \
      --reference $D/reference.fna --full_reference $D/full_reference.fna > $W/log_$1 2>&1
  echo "$1: $(cat $W/time_$1), index $(md5sum < $D/index.prx | cut -c1-8), unique_kmers $(md5sum < $D/unique_kmers.tsv | cut -c1-8)"
  grep -E " took |Uniqueness check:|Distance-two" $W/log_$1 | grep -v "Preload\|Run protal" | sed 's/^/    /'
  rm -rf $D
}
for p in 25 50 100; do
  build old_p${p}_t8 $OLD $p 8
  build new_p${p}_t8 $B/build/protal $p 8
  build new_p${p}_t1 $B/build/protal $p 1
done
build old_p100_t1 $OLD 100 1
