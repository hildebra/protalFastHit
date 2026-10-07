#!/usr/bin/env bash
# idle_check.sh TREE CHANGED_TAR DELETED_LIST STAGE [E2E_FILTER...] - brings TREE (an extracted archive with a build folder)
# to a newer commit (the files of CHANGED_TAR, less DELETED_LIST) plus this session's staged files (STAGE), rebuilds
# protal only and runs the e2e tests matching the filters (docs/claude/2026-10-07-sam-combine). Everything at idle CPU
# priority (SCHED_IDLE) and one build job, for a machine whose cores others need.
set -u
TREE=$1; TAR=$2; DELETED=$3; STAGE=$4; shift 4
IDLE="chrt --idle 0 nice -n 19"
cd "$TREE" || exit 1
tar -xf "$TAR"
while IFS= read -r f; do [ -n "$f" ] && rm -f "$f"; done < "$DELETED"
(cd "$STAGE" && find . -type f) | while IFS= read -r f; do cp "$STAGE/$f" "$TREE/$f"; touch "$TREE/$f"; done
$IDLE ninja -C build -j1 protal > build_protal.log 2>&1 || { echo "build failed"; grep -E "error" build_protal.log | head -20; exit 1; }
echo "built $(ls -la build/protal | awk '{print $6, $7, $8}')"
filters=()
for k in "$@"; do filters+=(-k "$k"); done
PROTAL_TEST_DB=$HOME/samcombine/db/database.protal PROTAL=$TREE/build/protal SIMULATE=$TREE/build/simulate_metagenomes \
  $IDLE python3 -m unittest tests/e2e/test_protal_e2e.py -v "${filters[@]}" 2>&1 | grep -E " \.\.\. |^Ran|^OK|^FAILED|Error" | tail -30
echo "done"
