#!/bin/bash
set -euo pipefail

rm -rf cmake-build-release
mkdir -p cmake-build-release "${PREFIX}/bin"

cmake \
	-DCMAKE_INSTALL_PREFIX="$PREFIX" \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_MAKE_PROGRAM=make \
	-DCMAKE_C_COMPILER="$GCC" \
	-DCMAKE_CXX_COMPILER="$GXX" \
	-G "Unix Makefiles" \
	-S "$SRC_DIR" \
	-B cmake-build-release

cmake --build cmake-build-release --target protal protal_avx2 simulate_metagenomes -- -j "${CPU_COUNT:-6}"

# Same layout as `just install`: the launcher is `protal` and picks protal_avx2 or the baseline build.
install -m 755 cmake-build-release/protal               "${PREFIX}/bin/protal_baseline"
install -m 755 cmake-build-release/protal_avx2          "${PREFIX}/bin/protal_avx2"
install -m 755 cmake-build-release/simulate_metagenomes "${PREFIX}/bin/simulate_metagenomes"
install -m 755 protal_launcher                          "${PREFIX}/bin/protal"
install -m 755 scripts/qcmsa.py                         "${PREFIX}/bin/qcmsa"
install -m 755 scripts/protal_map_utils                 "${PREFIX}/bin/protal_map_utils"
