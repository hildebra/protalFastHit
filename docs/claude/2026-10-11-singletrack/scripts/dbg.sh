#!/usr/bin/env bash
# dbg.sh: a scratch copy of the work tree's WFA2-lib with the Singletrack backtrace's state printed where it fails,
# linked into repro2 (the failing pair, ends-free with the read's end free, wf-adaptive on), next to the classic CIGAR.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad
R=$HOME/stk/repro; rm -rf $R/inc; mkdir -p $R/inc; cp -r $HOME/stk/work/lib/wfa2-lib $R/inc/wfa2-lib
B=$R/inc/wfa2-lib/wavefront/wavefront_backtrace.c
# print the state before the exit
perl -0pi -e 's/(    fprintf\(stderr,"\[WFA::Backtrace\] I\?\/D\?-Beginning backtrace error\\n"\);\n    fprintf\(stderr,">%\.\*s\\n",pattern_length)/    fprintf(stderr,"STATE v=%d h=%d score=%d k=%d offset=%d in_m=%d l=%d cigar=%s\\n",v,h,score,k,(int)offset,(int)in_mmatrix,l,cigar->operations+cigar->begin_offset+1);\n$1/' $B
grep -c STATE $B
cd $R/inc/wfa2-lib && rm -f *.o
for f in $(grep -oE '[a-z_/0-9]+\.c' $HOME/stk/work/lib/wfa2-lib.cmake) bindings/cpp/WFAligner.cpp; do
  case $f in *.cpp) g++ -O2 -std=c++17 -I. -c $f -o $(basename $f).o ;; *) gcc -O2 -I. -c $f -o $(basename $f).o ;; esac
done
ar rcs libwfadbg.a *.o
cat > $R/repro2.cpp <<'EOF'
#include <climits>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include "WFA2Wrapper2.h"
using namespace protal;
static std::string Read(const char* p) { std::ifstream in(p); std::string s; std::getline(in, s); return s; }
int main(int argc, char** argv) {
    std::string const ref = Read(argv[1]), read = Read(argv[2]);
    WFA2Wrapper2 classic(4, 6, 2, 0, WfaTraceback::Classic);
    classic.Alignment(read, ref, 0, static_cast<int>(read.size()), 0, 0, INT_MAX);
    std::printf("classic score %d\n%s\n", -classic.GetAlignmentScore(), classic.Cigar().c_str());
    std::fflush(stdout);
    WFA2Wrapper2 single(4, 6, 2, 0, WfaTraceback::Singletrack);
    single.Alignment(read, ref, 0, static_cast<int>(read.size()), 0, 0, INT_MAX);
    std::printf("singletrack score %d\n%s\n", -single.GetAlignmentScore(), single.Cigar().c_str());
}
EOF
H=$HOME/stk/work
g++ -O2 -std=c++20 -I$H/src -I$H/src/IO -I$H/src/Hash -I$H/src/Core -I$H/src/SequenceUtils -I$H/src/Alignment -I$H/src/Utilities \
  -I$H/src/Taxonomy -I$H/lib -I$H/lib/robin -I$H/lib/tsl -I$R/inc -I$R/inc/wfa2-lib -I$H/lib/gzstream -I$H/build/zlib-ng -I$H/lib/zlib-ng \
  -I$H/build/generated $R/repro2.cpp $R/inc/wfa2-lib/libwfadbg.a -o $R/repro2 2>&1 | grep error | head
$R/repro2 $R/pattern.txt $R/text.txt 2>&1 | grep -v "^[<>]" | cut -c1-2000
