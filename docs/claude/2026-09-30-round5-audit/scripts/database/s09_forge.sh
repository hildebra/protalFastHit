#!/usr/bin/env bash
# Adversarial forged directories/seek tables. Each should be rejected cleanly (no crash, no huge alloc).
set -u
P=~/strain-build/bin/protal
W=~/audit6/database
S=$(dirname "$0")
SRC=$W/bundle_only/database.protal
mkdir -p $W/forged
# open a forged bundle via --unpack_db (just needs Bundle::Open to succeed/fail)
try(){ # name
  local f=$W/forged/$1.protal
  python3 $S/forge.py $SRC $1 $f >/dev/null 2>$W/forged/$1.forgeerr || { echo "[$1] FORGE FAILED"; cat $W/forged/$1.forgeerr; return; }
  /usr/bin/timeout 60 /usr/bin/env bash -c "ulimit -v 4000000; '$P' --unpack_db --db '$f' --unpack_dir $W/forged/out_$1 -t 2" >$W/forged/$1.log 2>&1
  local rc=$?
  echo "[$1] exit=$rc  $(grep -iE 'error|invalid|not a|cannot|corrupt|traversal|file name|twice|tile|cover|version|missing' $W/forged/$1.log | head -2 | tr '\n' '|')"
  # did anything escape to a traversal/absolute path?
  [ -e $W/forged/escapee.map ] && echo "   !!! escapee.map CREATED outside unpack dir"
  [ -e /tmp/escapee.map ] && echo "   !!! /tmp/escapee.map CREATED (absolute traversal)"
}
for c in version traversal absolute dupe notile badcount hugedir hugecontent; do try $c; done

echo "=== a random regular file given as --db"
head -c 100000 /dev/urandom > $W/forged/random.bin
/usr/bin/timeout 30 "$P" --unpack_db --db $W/forged/random.bin --unpack_dir $W/forged/out_rand -t 2 >$W/forged/rand.log 2>&1
echo "[random.bin] exit=$? $(grep -iE 'error|neither|single-file|cannot' $W/forged/rand.log | head -2 | tr '\n' '|')"

echo "=== a plain reference.fna.zst (non-bundle seekable) given as --db"
zstd -q -f -3 $W/folder/reference.fna -o $W/forged/reference.fna.zst
/usr/bin/timeout 30 "$P" --unpack_db --db $W/forged/reference.fna.zst --unpack_dir $W/forged/out_ref -t 2 >$W/forged/ref.log 2>&1
echo "[reference.fna.zst] exit=$? $(grep -iE 'error|neither|single-file|cannot' $W/forged/ref.log | head -2 | tr '\n' '|')"
