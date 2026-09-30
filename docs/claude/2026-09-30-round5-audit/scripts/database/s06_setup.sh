#!/usr/bin/env bash
# Build working databases from the mini DB: a folder of separate files, plus a fresh bundle copy.
set -eu
P=~/strain-build/bin/protal
W=~/audit6/database
rm -rf $W/folder $W/bundle_only
mkdir -p $W/folder $W/bundle_only
cp ~/strain-build/mini_db/protal_db/database.protal $W/bundle_only/database.protal
# Unpack into a folder (index.prx.zst + raw reference.fna + text + model)
$P --unpack_db --db $W/bundle_only/database.protal --unpack_dir $W/folder -t 2 2>&1 | sed 's/^/[unpack] /'
echo "== folder contents"; ls -la $W/folder
echo "exit-of-unpack captured above"
