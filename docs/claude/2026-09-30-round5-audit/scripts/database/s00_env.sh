#!/usr/bin/env bash
# Survey of the environment for the database-container audit.
set -u
echo "== binaries"; ls -la ~/strain-build/bin/ 2>&1 | head -20
echo "== tests"; find ~/strain-build -maxdepth 4 -name protal_tests -type f 2>/dev/null
echo "== mini db"; ls -la ~/strain-build/mini_db/protal_db 2>&1
echo "== tune"; ls -la ~/tune/H2/protal_db 2>&1 | head; ls ~/genes_study/tuneH2 2>&1 | head
echo "== git"; cd ~/strain-build/src 2>/dev/null && git log --oneline -1 2>&1; git status --short | head
echo "== tools"; which zstd python3 gdb valgrind cmake ninja 2>&1; python3 -c 'import zstandard' 2>&1
nproc; free -g; df -h ~ | tail -1
