#!/bin/bash
python3 $(dirname "$0")/variants_cmp.py "$@" | column -t -s $'\t'
