#!/bin/bash
# Link without LTO: cPMML's own sources (ninja -t commands cPMML) compiled again with -fno-lto into an archive.
W=~/ctime
B=$W/build
D=$W/split
L=$W/link
rm -rf $L
mkdir -p $L
cd $B
ninja -t commands cPMML | grep -E ' -c [^ ]*\.cc( |$)' > $L/cpmml_cmds.txt
echo "cPMML sources: $(wc -l < $L/cpmml_cmds.txt)"
i=0
while IFS= read -r cmd; do
    i=$((i+1))
    c=${cmd/ -flto/ -fno-lto}
    c=$(echo "$c" | sed "s# -o [^ ]*# -o $L/cpmml_$i.o#")
    taskset -c 2 bash -c "$c" 2>/dev/null || echo "compile $i failed"
done < $L/cpmml_cmds.txt
ar rcs $L/libcPMML_nolto.a $L/cpmml_*.o
ls -la $L/*.a | awk '{print $5, $9}'
LINK=$(cat $D/link_cmd.txt)
ORIG=${LINK/CMakeFiles\/protal.dir\/src\/main.cpp.o/$D/seq_original.o}
NOLTO=${ORIG/src\/cPMML\/libcPMML.a/$L/libcPMML_nolto.a}
NOLTO=${NOLTO/-o protal /-o $L/protal_nolto }
for i in 1 2 3; do
    /usr/bin/time -f "link, cPMML without LTO: %e s" taskset -c 0-3 bash -c "$NOLTO"
done
$L/protal_nolto --version | head -1
