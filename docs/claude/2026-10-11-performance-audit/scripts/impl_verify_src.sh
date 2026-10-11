#!/usr/bin/env bash
S=$HOME/impl/work/src; REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
for f in src/Utilities/Zstd.h src/Utilities/Database.h src/RunProtal.h tests/test_Database.cpp; do
  if cmp -s "$S/$f" "$REPO/$f"; then echo "current: $f"; else echo "STALE: $f"; fi
done
grep -c ReadFramesPrefix $S/src/Utilities/Zstd.h $S/src/Utilities/Database.h
ls -la --time-style=+%T $HOME/impl/work/build/protal $HOME/impl/asan 2>/dev/null | head -3
