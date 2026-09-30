#!/bin/bash
S=$(dirname "$0")
cd ~/audit6/longreads
python3 $S/check_world.py ~/strain-build/mini_db/gtdb_r226 ~/audit6/longreads/mini_db
head -3 mini_db/reference.map
