#!/bin/bash
python3 $(dirname "$0")/tree_summary.py "$@" | column -t -s $'\t'
