#!/usr/bin/env bash
# variants.sh NAME FIXES ARCH: a copy of the work tree's WFA2-lib in ~/wfasimd/libs/NAME/wfa2-lib with FIXES applied
# (comma-separated), compiled with protal's flags (-fopenmp -O3 -DNDEBUG) for -march=ARCH, and wfareplay built against
# it as ~/wfasimd/libs/NAME/wfareplay. Fixes:
#   bounds   the vectorizable loops of wavefront_heuristic.c (wf-adaptive's distances, X-drop's scores) read lo/hi into
#            locals first: GCC could not compute their trip count (the stores may alias the wavefront struct)
#   term     the ends-free extend kernel calls wavefront_termination_endsfree (FORCE_NO_INLINE) only for a diagonal that
#            reached the end of the text or the pattern, the only cases in which it can stop the alignment
set -uo pipefail
NAME=$1; FIXES=${2:-}; ARCH=${3:-x86-64}
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfasimd
H=$HOME/wfasimd/src; V=$HOME/wfasimd/libs/$NAME; L=$V/wfa2-lib
rm -rf $V; mkdir -p $V; cp -r $H/lib/wfa2-lib $L
case ",$FIXES," in *,bounds,*)
  perl -0pi -e 's/  PRAGMA_LOOP_VECTORIZE\n  for \(k=wavefront->lo;k<=wavefront->hi;\+\+k\) \{/  const int k_lo = wavefront->lo, k_hi = wavefront->hi;\n  PRAGMA_LOOP_VECTORIZE\n  for (k=k_lo;k<=k_hi;++k) {/g' $L/wavefront/wavefront_heuristic.c
  echo "bounds: $(grep -c 'const int k_lo = wavefront->lo' $L/wavefront/wavefront_heuristic.c) loops" ;;
esac
case ",$FIXES," in *,term,*)
  perl -0pi -e 's/if \(wavefront_termination_endsfree\(wf_aligner,mwavefront,score,k,offset\)\)/if ((WAVEFRONT_H(k,offset) >= wf_aligner->sequences.text_length || WAVEFRONT_V(k,offset) >= wf_aligner->sequences.pattern_length) \&\& wavefront_termination_endsfree(wf_aligner,mwavefront,score,k,offset))/g; s/if \(endsfree \&\& wavefront_termination_endsfree\(wf_aligner,mwavefront,score,k,offset\)\)/if (endsfree \&\& (WAVEFRONT_H(k,offset) >= wf_aligner->sequences.text_length || WAVEFRONT_V(k,offset) >= wf_aligner->sequences.pattern_length) \&\& wavefront_termination_endsfree(wf_aligner,mwavefront,score,k,offset))/g' $L/wavefront/wavefront_extend_kernels.c
  echo "term: $(grep -c 'WAVEFRONT_H(k,offset) >= wf_aligner->sequences.text_length' $L/wavefront/wavefront_extend_kernels.c) call sites" ;;
esac
case ",$FIXES," in *,noavxext,*)
  # WFA2's own AVX2/AVX-512 extend kernels off (they end ends-free alignments on other diagonals; 2026-10-02 report),
  # so that -march=x86-64-v3 vectorizes the rest without changing the extend kernel
  perl -pi -e 's/^\s*#if __AVX2__ &&\s+__BYTE_ORDER == __LITTLE_ENDIAN/  #if 0 \/\* protal: portable extend kernels \*\//' $L/wavefront/wavefront_extend_kernels.c
  echo "noavxext: $(grep -c 'protal: portable extend kernels' $L/wavefront/wavefront_extend_kernels.c) dispatch sites" ;;
esac
cd $L; mkdir -p obj
for f in $(sed -n '/set(WFA2_SOURCES/,/)/p' $H/lib/wfa2-lib.cmake | grep -oE '[a-z_0-9]+/[a-z_0-9]+\.c'); do
  gcc -D_FILE_OFFSET_BITS=64 -I. -Iwavefront -Iutils -fopenmp -O3 -DNDEBUG -march=$ARCH -c $f -o obj/$(basename $f .c).o || echo "FAIL $f"
done
g++ -D_FILE_OFFSET_BITS=64 -I. -fopenmp -O3 -DNDEBUG -march=$ARCH -c bindings/cpp/WFAligner.cpp -o obj/WFAligner.o
ar rcs $V/libwfa.a obj/*.o
g++ -O3 -march=x86-64 -mtune=generic -std=c++20 -fopenmp -I$V -I$L -I$H/src -I$H/src/IO -I$H/src/Hash -I$H/src/Core -I$H/src/SequenceUtils \
  -I$H/src/Alignment -I$H/src/Utilities -I$H/src/Taxonomy -I$H/lib/robin -I$H/lib/tsl -I$H/lib/gzstream -I$H/build/zlib-ng -I$H/lib/zlib-ng \
  -I$H/build/generated $S/wfareplay.cpp $V/libwfa.a -o $V/wfareplay 2>&1 | grep -E ' error' | head -5
ls -la $V/wfareplay | awk '{print "built " $NF}'
