#!/usr/bin/env bash
# Determinism / byte-identity / atomicity of compress <-> decompress <-> bundle. At most 2 threads.
set -u
P=~/strain-build/bin/protal
W=~/audit6/database
S=$(dirname "$0")
Q(){ "$@" >/tmp/rt.log 2>&1; echo "  exit=$? ($*)"; grep -iE 'error|fail|warn|byte-identical|kept|differ' /tmp/rt.log | head -6 | sed 's/^/    /'; }

# reads for profiling exercises
python3 $S/mkreads.py $W/folder/reference.fna $W/reads sa >/dev/null

echo "=== 1. decompress folder to raw (index.prx + reference.fna)"
rm -rf $W/raw; mkdir -p $W/raw
for f in $W/folder/*; do ln -sf "$(realpath "$f")" $W/raw/$(basename "$f"); done
Q $P --decompress_db --db $W/raw -t 2
ls -la $W/raw | grep -E 'index.prx|reference.fna'

echo "=== 2. compress raw -> bundle, twice, same params; must be byte-identical (determinism)"
rm -rf $W/rtA $W/rtB; mkdir -p $W/rtA $W/rtB
for f in $W/raw/index.prx $W/raw/reference.fna $W/raw/reference.map $W/raw/internal_taxonomy.dmp $W/raw/unique_kmers.tsv $W/raw/model_pe.xml; do
  ln -sf "$(realpath "$f")" $W/rtA/$(basename "$f"); ln -sf "$(realpath "$f")" $W/rtB/$(basename "$f");
done
Q $P --compress_db --db $W/rtA -t 2 --compress_level 3 --compress_frame_mb 1
Q $P --compress_db --db $W/rtB -t 2 --compress_level 3 --compress_frame_mb 1
if cmp -s $W/rtA/database.protal $W/rtB/database.protal; then echo "  BYTE-IDENTICAL bundles across two runs"; else echo "  DIFFER: two --compress_db runs are NOT deterministic"; cmp $W/rtA/database.protal $W/rtB/database.protal; fi

echo "=== 3. decompress the bundle back to raw; compare index.prx & reference.fna to the originals"
rm -rf $W/rtC; mkdir -p $W/rtC
ln -sf "$(realpath $W/rtA/database.protal)" $W/rtC/database.protal
Q $P --decompress_db --db $W/rtC -t 2
for n in index.prx reference.fna reference.map internal_taxonomy.dmp unique_kmers.tsv model_pe.xml; do
  if cmp -s $W/rtC/$n $W/raw/$n; then echo "  $n: identical to raw"; else echo "  $n: DIFFERS from raw"; fi
done

echo "=== 4. leftover .partial after a simulated crash is ignored by Locate (folder still loads)"
touch $W/folder/index.prx.zst.partial $W/folder/database.protal.partial
Q $P --db $W/folder -1 $W/reads/sa_R1.fq -2 $W/reads/sa_R2.fq --prefix sa -o $W/out_partial -t 2 --no_qcmsa --no_strains
rm -f $W/folder/index.prx.zst.partial $W/folder/database.protal.partial
ls $W/out_partial 2>/dev/null | head
