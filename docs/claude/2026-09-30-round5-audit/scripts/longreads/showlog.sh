#!/bin/bash
# Prints a protal log from "Start parallel" on, long lines cut.
f=$1
awk '/Start parallel/{on=1} on' $f | cut -c1-${2:-220} | head -${3:-60}
echo "---- lines longer than 300 chars: $(awk 'length($0)>300' $f | wc -l)"
