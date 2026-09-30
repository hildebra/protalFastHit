#!/bin/bash
set -u
mkdir -p ~/audit6/io
ls -la ~/strain-build/ ~/strain-build/bin ~/strain-build/mini_db ~/strain-build/mini_db/protal_db 2>&1 | head -60
ls ~/strain-build/src/build/tests/ 2>&1 | head
cd ~/strain-build/src && git log --oneline -1 2>/dev/null
which art_illumina zstd gzip bgzip samtools python3 pigz
nproc; uptime
~/strain-build/bin/protal --help 2>&1 | head -120
