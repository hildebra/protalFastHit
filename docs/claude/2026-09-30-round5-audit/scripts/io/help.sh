#!/bin/bash
~/strain-build/bin/protal --help 2>&1 | sed -n '120,400p' | grep -E '^\s+(-[a-zA-Z], )?--'
