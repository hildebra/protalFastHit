#!/bin/bash
set -eu
W=~/audit6/longreads
mkdir -p $W
cd $W
if [ ! -d mini_db ]; then
  mkdir mini_db
  cp ~/strain-build/mini_db/protal_db/* mini_db/
  ~/strain-build/bin/protal --unpack_db --db mini_db/database.protal -t 2 > unpack.log 2>&1 || { tail -20 unpack.log; exit 1; }
fi
ls -la mini_db
tail -5 unpack.log
grep -c . mini_db/*.tsv || true
for f in mini_db/*.xml; do echo $f; head -c 600 $f | grep -o 'protal:placeholder' || echo "not placeholder"; done
