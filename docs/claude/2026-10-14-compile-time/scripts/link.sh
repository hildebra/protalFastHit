#!/bin/bash
# Link time of protal from the one main.cpp.o and from the split prototype's seven objects.
W=~/ctime
B=$W/build
D=$W/split
cd $B
LINK=$(cat $D/link_cmd.txt)
ORIG=${LINK/CMakeFiles\/protal.dir\/src\/main.cpp.o/$D/seq_original.o}
ORIG=${ORIG/-o protal /-o $D/protal_orig }
SPLIT=${LINK/CMakeFiles\/protal.dir\/src\/main.cpp.o/$D/seq_main.o $D/seq_split_align.o $D/seq_split_profile.o $D/seq_split_strains.o $D/seq_split_buildmaint.o $D/seq_split_buildindex.o $D/seq_split_options.o}
SPLIT=${SPLIT/-o protal /-o $D/protal_split2 }
for i in 1 2; do
    /usr/bin/time -f "link one object: %e s" taskset -c 2 bash -c "$ORIG"
    /usr/bin/time -f "link seven objects: %e s" taskset -c 2 bash -c "$SPLIT"
done
ls -la $D/protal_orig $D/protal_split2 | awk '{print $5, $9}'
$D/protal_orig --version | head -1; $D/protal_split2 --version | head -1
size $D/protal_orig $D/protal_split2
