# Installing protal

protal runs on Linux on x86-64 CPUs. There is no macOS or native Windows version; on Windows it
builds and runs under WSL2 (see [Windows (WSL2)](#windows-wsl2)).

At run time protal needs `python3` for the strain MSA post-filter [qcmsa](qcmsa.md), which the
bioconda package brings. It compresses its outputs itself (zstd, and gzip with libdeflate), so
no external compressor is needed.

## bioconda (recommended)

```bash
conda create -n protal protal -c bioconda -c conda-forge
conda activate protal
protal --version
```

Then download a database, as described on the
[website](https://protal.earlham.ac.uk/main.php?site=documentation#download-the-database).

## Static binaries

The [GitHub releases](https://github.com/4less/protal/releases) have statically linked binaries,
for clusters where conda is not an option or where compute nodes differ from the node that
installed the software. qcmsa is `scripts/qcmsa.py` of the source (or run protal with
`--no_qcmsa`).

To build them yourself (needs the static libraries of zstd and libdeflate):

```bash
just static        # -> build/protal_<version>_static, build/simulate_metagenomes_static
```

They are compiled for plain x86-64 (SSE2), so they run on any x86-64 CPU, and use AVX2 where the CPU
has it, as `protal` does (see [One binary for every CPU](#one-binary-for-every-cpu)).

## Building from source

Requirements: CMake 3.22 or later, a C++20 compiler with OpenMP (GCC 13 and 14 are tested),
zstd and libdeflate development files. On Ubuntu:

```bash
sudo apt-get install cmake ninja-build g++ libzstd-dev libdeflate-dev python3
```

All other libraries (zlib-ng, WFA2-lib, cPMML, gzstream, robin-map, ...) are in `lib/` and built
with protal; no system zlib is needed. Then:

```bash
git clone https://github.com/4less/protal.git
cd protal
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target protal simulate_metagenomes -j 8
./build/protal --version
```

A default build is `Release`; `-DCMAKE_BUILD_TYPE=Debug` builds without optimisation and with
`assert()` checks on.

| Target | Binary | |
|---|---|---|
| `protal` | `build/protal` | runs on any x86-64 CPU, with AVX2 where the CPU has it ([below](#one-binary-for-every-cpu)) |
| `protal_static` | `build/protal_<version>_static` | fully static build of `protal` |
| `simulate_metagenomes` | `build/simulate_metagenomes` | read simulator, see [simulation.md](simulation.md) |
| `simulate_metagenomes_static` | `build/simulate_metagenomes_static` | static simulator |
| `protal_tests` | `build/tests/protal_tests` | unit tests; needs `-DPROTAL_BUILD_TESTS=ON` and GoogleTest, see [testing.md](testing.md) |

The [justfile](../justfile) wraps these: `just baseline` (protal), `just simulate`,
`just build-all` (the two binaries that get installed), `just static`, and `just clear` to
delete the build trees.

### One binary for every CPU

`protal` is compiled for plain x86-64 (SSE2), so it runs on any x86-64 CPU. Where the CPU has AVX2
and the other x86-64-v3 instructions (Intel Core since Haswell, 2013; AMD since 2015), protal runs faster code
in the places where that was measured to pay. It chooses when it starts, and its output is the
same on every CPU.

- The syncmer scan and sequence packing have AVX2 kernels of protal's own.
- zlib-ng, libdeflate and zstd choose their own.
- The hot functions marked `PROTAL_CLONE_V3` (`src/Utilities/TargetClones.h`) are compiled twice,
  for x86-64 and for x86-64-v3, and the loader picks one (GCC's `target_clones`). This needs GCC
  12 or later on Linux with glibc. With other compilers, or with
  `-DCMAKE_CXX_FLAGS=-DPROTAL_NO_CLONES`, they are compiled once, for x86-64.

WFA2-lib is built for plain x86-64, because its AVX2 kernels were not faster in protal's use
([report](claude/2026-10-02-wfa-avx2/README.md)).

### Installing a source build

```bash
just install                    # into $HOME/.local/bin; `just install /opt/protal` into /opt/protal/bin
```

This rebuilds first and installs the same layout as the conda package:

| Installed as | From |
|---|---|
| `protal` | the build |
| `simulate_metagenomes` | the simulator |
| `qcmsa` | `scripts/qcmsa.py`, the strain MSA post-filter |
| `protal_map_utils` | `scripts/protal_map_utils`: `generate` a map from read folders, `merge` or `flatten` maps, `validate` one |
| `protal_profile_utils` | `scripts/protal_profile_utils`: merge profiles into one abundance table (see the [website](https://protal.earlham.ac.uk/main.php?site=documentation#species-profiles)) |

The conda package also installs the database build and training scripts under
`share/protal/scripts/` (the layout of `scripts/`); their Python requirements are not part of the
package, see [Tools to build a database](#tools-to-build-a-database).

Up to version 0.7.2, `protal` was a launcher that ran one of two builds, `protal_baseline` or
`protal_avx2`; `just install` removes those two from the prefix.

protal finds qcmsa next to its own binary, then on `$PATH`, then through the environment variable
`PROTAL_QCMSA_SCRIPT`, and finally as `scripts/qcmsa.py` of a source checkout; `--qcmsa_script`
names it directly.

### A conda package from the checkout

`conda-recipe/` builds a conda package from the local checkout (the version in `meta.yaml` must
match the one in `CMakeLists.txt`):

```bash
conda install conda-build
conda build conda-recipe -c conda-forge -c bioconda --output-folder conda-build
conda create -n protal_local -c "file://$PWD/conda-build" -c conda-forge -c bioconda protal
```

## Tools to build a database

Building and training a database ([building-a-database.md](building-a-database.md)) needs more than
profiling: compilers for protal from the checkout, ART for the simulations, Python with
numpy, pandas, scikit-learn and joblib for the training, and NCBI's `datasets` for strain genomes.
`envs/protal-db-build.yaml` is a conda environment with all of them:

```bash
conda env create -f envs/protal-db-build.yaml      # or: micromamba create -f envs/protal-db-build.yaml
conda activate protal-db-build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
cmake --build build --target protal simulate_metagenomes -j 16
```

protal is built from the checkout because bioconda's protal 0.6.0a predates the single-file
database, `--add_model` and the read types the workflow uses.

## Windows (WSL2)

protal builds and runs in WSL2 with Ubuntu, with the packages listed above. Build from a copy of
the source on the Linux file system (for example under `~/`), not from `/mnt/c/...`: Windows file
systems are case-insensitive, so cPMML's `#include "options.h"` finds protal's `src/Options.h`
and the build fails. Keep databases on the Linux file system too: reading and writing them
through `/mnt/c` is slow.

## Checking the installation

```bash
protal --version
protal --help          # --full_help adds the alignment and developer options
just example           # from a source checkout: builds a tiny database, profiles simulated reads and checks the result
```

`just example` needs no download; see [`examples/mini_db/`](../examples/mini_db/README.md).
