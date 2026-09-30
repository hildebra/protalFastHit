#!/bin/bash
set -u
ls -la ~/strain-build/mini_db/ ~/strain-build/mini_db/protal_db
ls ~/strain-build/mini_db/protal_db/* | head -40
cut -f2 ~/audit5/world/gtdb_r226/simulation/genomes.tsv | sort | uniq -c | head -30
f=$(sed -n 2p ~/audit5/world/gtdb_r226/simulation/genomes.tsv | cut -f3)
zcat "$f" | grep -c '>'
zcat "$f" | awk '/^>/{if(n)print n; n=0; next}{n+=length($0)}END{print n}' | head
ls ~/audit5/world/gtdb_r226/
ls ~/audit5/
