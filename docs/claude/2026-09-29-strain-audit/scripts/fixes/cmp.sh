#!/bin/bash
python3 $(dirname "$0")/compare_steps.py "$@" | column -t -s $'\t'
