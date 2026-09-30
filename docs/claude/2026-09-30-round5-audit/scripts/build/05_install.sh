#!/bin/bash
# `just install` into a temporary prefix (a fake nproc caps the justfile's -j$(nproc) at 2), then
# the installed launcher: which binary it runs, PROTAL_NO_AVX2, a CPU without x86-64-v3 (a fake
# grep on $PATH, since the launcher reads /proc/cpuinfo with grep), a missing AVX2 binary, no
# binary at all, a symlinked launcher.
set -u
A=~/audit6/build
P=$A/prefix
mkdir -p $A/fakebin
printf '#!/bin/sh\necho 2\n' > $A/fakebin/nproc; chmod +x $A/fakebin/nproc
export PATH=$A/fakebin:$PATH
cd $A/src
rm -rf $P
rm -rf "$A/src/prefix="; just install prefix=$P > $A/install_documented_form.log 2>&1
echo "documented form 'just install prefix=\$P' rc=$?"; ls -d "$A/src/prefix="*/*/*/*/*/*/bin 2>&1; ls $P 2>&1; grep -E '^Installed|command not found' $A/install_documented_form.log
rm -rf "$A/src/prefix="; { time just install $P ; } > $A/install.log 2>&1; echo "just install (positional) rc=$?"
tail -4 $A/install.log
echo "== installed files"; ls -la $P/bin; find $P -path '*share*' | head -3; echo "(share/ entries above, if any)"
echo "== launcher: which binary (sh -x)"
sh -x $P/bin/protal --version 2>&1 | grep -E '^\+ exec|^protal'
echo "== PROTAL_NO_AVX2=1"; PROTAL_NO_AVX2=1 sh -x $P/bin/protal --version 2>&1 | grep -E '^\+ exec|^protal'
echo "== PROTAL_NO_AVX2=0 (non-empty: also baseline)"; PROTAL_NO_AVX2=0 sh -x $P/bin/protal --version 2>&1 | grep -E '^\+ exec'
echo "== CPU without AVX2 (fake grep prints a flags line without avx2)"
mkdir -p $A/nov3; printf '#!/bin/sh\necho "flags\t\t: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr sse sse2 ssse3 sse4_1 sse4_2 popcnt"\n' > $A/nov3/grep; chmod +x $A/nov3/grep
PATH=$A/nov3:$PATH sh -x $P/bin/protal --version 2>&1 | grep -E '^\+ exec|^protal'
echo "== CPU with avx2 but without abm (lzcnt): must fall back"
printf '#!/bin/sh\necho "flags\t\t: fpu sse sse2 avx avx2 bmi1 bmi2 fma f16c movbe popcnt"\n' > $A/nov3/grep
PATH=$A/nov3:$PATH sh -x $P/bin/protal --version 2>&1 | grep -E '^\+ exec'
echo "== no /proc/cpuinfo readable (grep fails) -> baseline"
printf '#!/bin/sh\nexit 2\n' > $A/nov3/grep
PATH=$A/nov3:$PATH sh -x $P/bin/protal --version 2>&1 | grep -E '^\+ exec'
echo "== AVX2 binary hidden"
mv $P/bin/protal_avx2 $P/bin/protal_avx2.hidden
sh -x $P/bin/protal --version 2>&1 | grep -E '^\+ exec|^protal'
echo "== no binary at all"
mv $P/bin/protal_baseline $P/bin/protal_baseline.hidden
(PATH=/usr/bin:/bin $P/bin/protal --version; echo "rc=$?") 2>&1
mv $P/bin/protal_avx2.hidden $P/bin/protal_avx2; mv $P/bin/protal_baseline.hidden $P/bin/protal_baseline
echo "== symlinked launcher"
mkdir -p $A/linkdir; ln -sf $P/bin/protal $A/linkdir/protal
sh -x $A/linkdir/protal --version 2>&1 | grep -E '^\+ exec'
echo "== exit status passes through (bad option)"
$P/bin/protal --no_such_option > /dev/null 2>&1; echo "rc=$?"; $P/bin/protal_baseline --no_such_option > /dev/null 2>&1; echo "baseline rc=$?"
echo "== prefix with a space"
rm -rf "$A/with space"
just install "$A/with space" > $A/install_space.log 2>&1; echo "rc=$?"; tail -3 $A/install_space.log; ls -d "$A/with" "$A/with space" $A/space 2>&1; ls "$A/with space/bin" 2>&1 | head -3
rm -rf "$A/with space" "$A/with" $A/space ~/audit6/build/src/space 2>/dev/null
echo "== just --list"; just --list 2>&1 | head -60
echo DONE
