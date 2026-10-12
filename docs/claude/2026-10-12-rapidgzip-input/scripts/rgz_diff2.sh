#!/bin/bash
# Files present in both upstream 2.32.1 and the fork that differ, with line counts.
set -uo pipefail
cd ~/rgz
U=up/isa-l-2.32.1
cd fork
for f in $(find . -type f | sort); do
  if [ -f ../$U/$f ]; then
    if ! cmp -s $f ../$U/$f; then
      n=$(diff ../$U/$f $f | grep -c '^[<>]')
      echo "differs: $f ($n lines)"
    fi
  else
    echo "fork only: $f"
  fi
done
cd ..
diff -u $U/include/igzip_lib.h fork/include/igzip_lib.h > igzip_lib.diff
diff -u $U/igzip/igzip_inflate.c fork/igzip/igzip_inflate.c > igzip_inflate.diff
wc -l igzip_lib.diff igzip_inflate.diff
