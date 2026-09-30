#!/bin/bash
# Peak RSS of protal --build on the tuning world: build_memory.sh THREADS  (the build packs database.protal and
# removes the separate files, so the conversion is repeated for every run)
set -e
cd ~/protal-mem
t=$1
rm -rf db900b
python3 ~/protal/scripts/mini_db/gtdb_to_protal_db.py --gtdb ~/tune/world --outdir ~/protal-mem/db900b > db900b.convert.log 2>&1
/usr/bin/time -v -o build_t$t.time ./build/protal --build --no_profile -t $t --db db900b --reference db900b/reference.fna --full_reference db900b/full_reference.fna > build_t$t.log 2>&1
grep -E "Maximum resident|Elapsed" build_t$t.time
