#!/usr/bin/env bash
# Does a forged seek-table decompressed_size trigger a huge allocation, or is it caught first?
# Case 1: forge a small member (reference.map) written by protal (known content size) -> --decompress_db
#         loads it (FrameStreambuf/ParallelRead). Expect the getFrameContentSize guard to catch it.
set -u
P=~/strain-build/bin/protal
W=~/audit6/database
S=$(dirname "$0")
SRC=$W/bundle_only/database.protal

# forge reference.map frame (index 51 in the 0-based frame list = member reference.fna? re-check via dbinfo)
python3 - "$SRC" <<'PY'
import sys, struct, subprocess
sys.path.insert(0, "$(dirname "$0")")
import forge
data, frames = forge.parse(sys.argv[1])
ver, mem = forge.members(data, frames)
for name, first, nf in mem:
    print(name.decode(), "first_frame=", first)
PY

# reference.map is a single frame; from dbinfo it's frame index 52. Forge that frame's content size to ~4GB.
python3 - <<PY
import sys
sys.path.insert(0, "$S")
import forge
data, frames = forge.parse("$SRC")
forge.emit("$W/forged/mapbig.protal", data, frames, seek_overrides={52: 0xFFFFFFF0})
print("forged reference.map (frame 52) content-size -> 4GB")
PY

echo "=== load the forged bundle with a 2 GB address-space cap; a 4GB alloc would fail with bad_alloc"
/usr/bin/timeout 120 /usr/bin/env bash -c "ulimit -v 2000000; '$P' --decompress_db --db $W/forged/mapbig.protal -t 2" >$W/forged/mapbig.log 2>&1
echo "exit=$?"
grep -iE 'error|bad_alloc|holds|seek table says|invalid|corrupt|reference.map|out of memory|terminate' $W/forged/mapbig.log | head -6
# clean up the raw files a decompress may have written next to it
rm -f $W/forged/index.prx $W/forged/reference.fna $W/forged/reference.map $W/forged/internal_taxonomy.dmp $W/forged/unique_kmers.tsv $W/forged/model_pe.xml
