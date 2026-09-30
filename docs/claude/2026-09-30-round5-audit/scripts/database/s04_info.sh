#!/usr/bin/env bash
set -u
S=$(dirname "$0")
python3 $S/dbinfo.py ~/audit6/database/work/mini.protal
