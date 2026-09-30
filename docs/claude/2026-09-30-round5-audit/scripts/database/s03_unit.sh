#!/usr/bin/env bash
# Runs the existing unit tests of the database container, zstd and index codec (no large index load).
set -u
T=~/strain-build/src/build/tests/protal_tests
mkdir -p ~/audit6/database/work
cd ~/audit6/database/work
OMP_NUM_THREADS=2 "$T" --gtest_filter='Database.*:Zstd.*:ZstdSeekable.*:IndexCodec.*' 2>&1 | tail -40
echo "exit=$?"
echo "== mini bundle members (zstd -l on a copy)"
cp -f ~/strain-build/mini_db/protal_db/database.protal ~/audit6/database/work/mini.protal
zstd -lv ~/audit6/database/work/mini.protal 2>&1 | head -20
