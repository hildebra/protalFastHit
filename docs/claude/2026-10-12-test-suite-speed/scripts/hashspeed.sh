#!/bin/bash
export PATH=~/micromamba/envs/protal-db-build/bin:$PATH
W=/tmp/hashspeed; rm -rf $W; mkdir -p $W; cd $W
P=~/testspeed/build/protal
$P --unpack_db --db ~/testspeed/src/data/mini_db/protal_db/database.protal --unpack_dir $W/f -t 4 > /dev/null 2>&1
mkdir raw; for f in f/*; do ln -s $W/$f raw/$(basename $f); done
$P --decompress_db --db raw -t 4 > /dev/null 2>&1
for py in python3 /usr/bin/python3; do
taskset -c 0-3 $py - <<'PY'
import hashlib, zlib, time, sys
path = "raw/index.prx"
def run(name, fn):
    t = time.time(); v = fn(); print(f"{sys.version.split()[0]} {name}: {time.time() - t:.2f} s", flush=True)
def blocks():
    with open(path, "rb") as f:
        yield from iter(lambda: f.read(1 << 24), b"")
def h(algo):
    x = hashlib.new(algo)
    for b in blocks(): x.update(b)
    return x.hexdigest()
def crc():
    c = 0
    for b in blocks(): c = zlib.crc32(b, c)
    return c
def readonly():
    n = 0
    for b in blocks(): n += len(b)
    return n
run("read only", readonly)
run("sha256", lambda: h("sha256"))
run("sha1", lambda: h("sha1"))
run("blake2b", lambda: h("blake2b"))
run("md5", lambda: h("md5"))
run("crc32", crc)
PY
done
rm -rf $W
