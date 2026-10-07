#!/usr/bin/env bash
# massif_pair.sh - heap profiles of two w900 samples with strain MSAs, by the binary of 12b059b and by the changed one
# (docs/claude/2026-10-07-sam-combine), one after the other.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
bash "$HERE/massif.sh" "$HOME/samcombine/bin_before/protal" before_strains 2
bash "$HERE/massif.sh" "$HOME/samcombine/bin_trim_protal" packed_strains 2
echo "pair done"
