#!/usr/bin/env bash
# The builds of the report: base (48e8cd0), base-reorder, work (the change), work-reorder.
HERE=$(cd "$(dirname "$0")" && pwd)
for spec in "base-reorder 48e8cd0 reorder" "work tree" "work-reorder tree reorder"; do $HERE/build.sh $spec || exit 1; done
