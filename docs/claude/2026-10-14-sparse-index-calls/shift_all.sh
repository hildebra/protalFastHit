#!/usr/bin/env bash
# feature_shift.pl for the read lengths (150 against 100 and 75) and the seeding (ungated against the grid's sparse index).
S=$(cd "$(dirname "$0")" && pwd)
W=~/strainfn
echo "### sanity: g0 re-profiled against len150 (other binary, same reads)"; perl $S/feature_shift.pl $W/run_cf_g0/profiles $W/run_len150/profiles 0.05 | head -8
echo "### 150 -> 100"; perl $S/feature_shift.pl $W/run_len150/profiles $W/run_len100/profiles
echo "### 150 -> 75"; perl $S/feature_shift.pl $W/run_len150/profiles $W/run_len75/profiles
echo "### seeding: g0 -> g4"; perl $S/feature_shift.pl $W/run_cf_g0/profiles $W/run_cf_g4/profiles
