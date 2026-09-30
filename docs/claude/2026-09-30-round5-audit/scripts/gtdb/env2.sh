#!/bin/bash
exec > ~/audit6/gtdb/env2.out 2>&1
find / -xdev \( -path /proc -o -path /mnt \) -prune -o -type d -name sklearn -print 2>/dev/null | head
ls ~/.local/lib 2>/dev/null
ls -d ~/*/
python3 -c 'import scipy; print(scipy.__version__)'
