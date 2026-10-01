#!/bin/bash
# Keep the previous binaries, then build the current checkout.
set -e
source "$(dirname "$0")/env.sh"
mkdir -p $PERF_DIR/bin-old
[ -e $PERF_DIR/bin-old/protal_avx2 ] || cp $PERF_DIR/build-rel/protal_avx2 $PERF_DIR/build-rel/simulate_metagenomes $PERF_DIR/bin-old/ 2>/dev/null || true
bash $here/build.sh
