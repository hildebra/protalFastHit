#!/bin/bash
# main.cpp.o with parts of protal::Run cut out (#ifndef in a patched copy of RunProtal.h): each part's exclusive
# share of the -O3 compile. The patched header goes into its own include directory searched first.
W=~/ctime
B=$W/build
R=$W/res
P=$R/patched
mkdir -p $P
H=$W/src/src/RunProtal.h
# Bottom to top, so that the line numbers stay valid.
sed -e '4298i #ifndef NO_STRAINS' -e '4298a #endif' \
    -e '4247i #ifndef NO_PROFILE' -e '4247a #else\n            Profiles profiles;\n#endif' \
    -e '4225i #ifndef NO_ALIGN' -e '4238a #endif' \
    -e '4221i #ifndef NO_PROFILE' -e '4221a #endif' \
    -e '4214i #ifndef NO_PROFILE' -e '4214a #endif' \
    -e '3838i #ifndef NO_BUILDMAINT' -e '3880a #endif' \
    -e '3735s/^    size_t/    static size_t/' \
    -e '2625s/^    TaxidList/    static TaxidList/' \
    -e '2398s/^    Profiles/    static Profiles/' \
    -e '2338s/^    Profiles/    static Profiles/' \
    $H > $P/RunProtal.h
diff $H $P/RunProtal.h | head -60
CMD=$(grep -o '"command": "[^"]*protal.dir/src/main.cpp[^"]*"' $B/compile_commands.json | head -1 | sed 's/^"command": "//; s/"$//')
BASE=${CMD% -o *}
SRC=$W/src/src/main.cpp
# -iquote $P: "RunProtal.h" from main.cpp resolves to the patched copy (the includer's own directory is searched
# before -iquote only for the includer's directory... so compile a copy of main.cpp placed next to the patched header).
cp $SRC $P/main.cpp
# Syntax check of every variant first (fast fail)
ALL="-DNO_ALIGN -DNO_PROFILE -DNO_STRAINS -DNO_BUILDMAINT"
taskset -c 0-3 nice -n 5 $BASE $ALL -fsyntax-only $P/main.cpp && echo "syntax ok (all cut)"
run() {  # name core flags...
    local name=$1 core=$2; shift 2
    /usr/bin/time -f "$name: %e s wall, %U s user, %M KB" taskset -c $core nice -n 5 $BASE "$@" -ftime-report -c $P/main.cpp -o $R/v_$name.o 2> $R/ft_$name.txt
    tail -1 $R/ft_$name.txt
}
uptime
run full 0 &
run no_align 1 -DNO_ALIGN &
run no_profile 2 -DNO_PROFILE &
run no_strains 3 -DNO_STRAINS &
wait
uptime
run no_buildmaint 0 -DNO_BUILDMAINT &
run core 1 $ALL &
run noclones 2 -DPROTAL_NO_CLONES &
run O2 3 -O2 &
wait
uptime
for f in $R/ft_*.txt; do echo "$(basename $f) $(grep -E '^ TOTAL' $f)"; done
