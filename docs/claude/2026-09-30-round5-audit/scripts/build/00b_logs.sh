#!/bin/bash
# Show the logs of the prebuilt tree (read-only) for timings and settings.
for f in ~/strain-build/build.log ~/strain-build/e2e.log ~/strain-build/mini_db.log; do echo "=== $f"; cat "$f"; done
echo "=== configure.log (head)"; head -40 ~/strain-build/configure.log
echo "=== mini_db/build.log tail"; tail -5 ~/strain-build/mini_db/build.log
ls -la ~/strain-build/mini_db/gtdb_r226 ~/strain-build/mini_db/gtdb_r226/simulation
head -5 ~/strain-build/mini_db/gtdb_r226/simulation/genomes.tsv
