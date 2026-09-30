#!/usr/bin/env bash
cd ~/audit6/database/fuzz || exit 1
ASAN_OPTIONS=abort_on_error=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ./fuzz
echo "fuzz-exit=$?"
