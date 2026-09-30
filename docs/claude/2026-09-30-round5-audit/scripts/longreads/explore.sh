#!/bin/bash
set -u
mkdir -p ~/audit6/longreads
ls -la ~/strain-build/ ~/strain-build/bin ~/strain-build/mini_db ~/strain-build/mini_db/protal_db 2>&1 | head -80
echo ---
ls ~/audit5/world/ ~/audit5/world/protal_db 2>&1 | head -40
ls ~/audit5/world/gtdb_r226/simulation/ 2>&1 | head
nproc; uptime
cd ~/strain-build/src && git log --oneline -1 2>/dev/null || true
ls ~/strain-build/src | head -50
