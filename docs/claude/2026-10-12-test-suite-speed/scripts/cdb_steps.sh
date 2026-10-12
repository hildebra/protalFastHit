#!/bin/bash
# The steps of the e2e CompressedDatabaseTest, timed one by one, on the mini database.
export PATH=~/micromamba/envs/protal-db-build/bin:$PATH
P=~/testspeed/build/protal
DB=~/testspeed/src/data/mini_db/protal_db/database.protal
W=/tmp/cdb_steps; rm -rf $W; mkdir -p $W; cd $W
TS="taskset -c 0-3 nice -n 5"
t() { local s=$(date +%s.%N); "$@" > /dev/null 2>&1; echo "$(echo "$(date +%s.%N) - $s" | bc | cut -c1-5) s: ${*:0:12}"; }
t $TS $P --unpack_db --db $DB --unpack_dir $W/files -t 4
for k in raw seekable seekable1 single bundle bundle1; do mkdir $k; for f in files/*; do n=$(basename $f); case $n in index.prx*|reference.fna*) ;; *) ln -s $W/$f $k/$n;; esac; done; done
ln -s $W/files/index.prx.zst raw/index.prx.zst; ln -s $W/files/reference.fna raw/reference.fna 2>/dev/null || ln -s $W/files/reference.fna.zst raw/reference.fna.zst
ls -la raw | head -5
echo "== decompress_db (raw copy)"; t $TS $P --decompress_db --db raw -t 4
ls -la raw/index.prx
echo "== sha256 of the raw index (python hashlib)"; t $TS python3 -c "import hashlib,sys; h=hashlib.sha256(); f=open('raw/index.prx','rb'); [h.update(b) for b in iter(lambda: f.read(1<<24), b'')]; print(h.hexdigest())"
echo "== sha1 (python hashlib)"; t $TS python3 -c "import hashlib,sys; h=hashlib.sha1(); f=open('raw/index.prx','rb'); [h.update(b) for b in iter(lambda: f.read(1<<24), b'')]; print(h.hexdigest())"
echo "== cmp raw against itself (read twice)"; t $TS cmp raw/index.prx raw/index.prx
echo "== zstd -1 -T0 of the raw index (single)"; t $TS zstd -q -f -1 -T0 raw/index.prx -o single/index.prx.zst
for k in seekable seekable1 bundle bundle1; do ln -s $W/raw/index.prx $k/index.prx; ln -s $W/raw/reference.fna $k/reference.fna; done
echo "== compress_db --no_bundle level 3"; t $TS $P --compress_db --no_bundle --db seekable -t 4 --compress_level 3 --compress_frame_mb 1
echo "== compress_db --no_bundle level 1"; t $TS $P --compress_db --no_bundle --db seekable1 -t 4 --compress_level 1 --compress_frame_mb 1
echo "== compress_db bundle level 3"; t $TS $P --compress_db --db bundle -t 4 --compress_level 3 --compress_frame_mb 1
echo "== compress_db bundle level 1"; t $TS $P --compress_db --db bundle1 -t 4 --compress_level 1 --compress_frame_mb 1
ls -la seekable/index.prx.zst seekable1/index.prx.zst bundle/database.protal bundle1/database.protal
for k in seekable bundle; do
  mkdir rt_$k; for f in $k/*; do ln -s $(readlink -f $f) rt_$k/$(basename $f); done
  echo "== round trip decompress_db of $k"; t $TS $P --decompress_db --db rt_$k -t 4
done
rm -rf $W
