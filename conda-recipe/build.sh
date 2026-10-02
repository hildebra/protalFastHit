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

cmake --build cmake-build-release --target protal simulate_metagenomes -- -j "${CPU_COUNT:-6}"

# Same layout as `just install`. protal is one binary for every x86-64 CPU: it runs x86-64-v3 code
# (AVX2) in its hot functions where the CPU has it.
install -m 755 cmake-build-release/protal               "${PREFIX}/bin/protal"
install -m 755 cmake-build-release/simulate_metagenomes "${PREFIX}/bin/simulate_metagenomes"
install -m 755 scripts/qcmsa.py                         "${PREFIX}/bin/qcmsa"
install -m 755 scripts/protal_map_utils                 "${PREFIX}/bin/protal_map_utils"
install -m 755 scripts/protal_profile_utils             "${PREFIX}/bin/protal_profile_utils"

# The database build and training scripts, in the layout they import each other from (their Python
# requirements are in envs/protal-db-build.yaml, not in this package):
#   python $CONDA_PREFIX/share/protal/scripts/build_gtdb_database.py --gtdb ... --outdir ...
mkdir -p "${PREFIX}/share/protal/scripts/mini_db"
for f in scripts/*.py scripts/random_forest.xml; do
	install -m 644 "$f" "${PREFIX}/share/protal/scripts/"
done
for f in scripts/mini_db/*.py scripts/mini_db/*.sh scripts/mini_db/markers_r226.tsv; do
	install -m 644 "$f" "${PREFIX}/share/protal/scripts/mini_db/"
done
chmod 755 "${PREFIX}"/share/protal/scripts/*.py "${PREFIX}"/share/protal/scripts/mini_db/*.py
