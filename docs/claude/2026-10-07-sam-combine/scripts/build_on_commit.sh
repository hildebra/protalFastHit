#!/usr/bin/env bash
# build_on_commit.sh TREE FILE... - copies this session's changed FILEs (paths in the checkout) into TREE, an extracted
# `git archive` with a build folder (verify_commit.sh), and rebuilds it (docs/claude/2026-10-07-sam-combine). The
# checkout holds other sessions' uncommitted work, so only the named files are taken. 4 cores.
set -eu
TREE=$1; shift
SRC_WIN=${SRC_WIN:-/mnt/c/Users/hildebra/Documents/locDev/protal}
for f in "$@"; do
  mkdir -p "$TREE/$(dirname "$f")"
  cp "$SRC_WIN/$f" "$TREE/$f"
  touch "$TREE/$f"
done
cd "$TREE"
taskset -c 0-3 nice -n 5 ninja -C build -j4 > build.log 2>&1 || { echo "build failed"; grep -E "error|Error" build.log | head -30; exit 1; }
grep -E "warning" build.log | grep -v "lto-wrapper\|serial compilation" | head -10
echo "built $(ls -la build/protal)"
