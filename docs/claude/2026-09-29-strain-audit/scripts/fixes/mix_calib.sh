#!/bin/bash
S=$(cd "$(dirname "$0")" && pwd)
for l in "$@"; do echo "######## $l"; python3 $S/mix_calib.py $l; done
