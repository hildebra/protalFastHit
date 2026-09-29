# protal documentation

How to install protal, run it and read its output is described on the website:
https://protal.earlham.ac.uk/main.php?site=documentation (installation, database download, a
worked four-sample example, the map file, and every output file). The pages here cover the rest.

| Page | For | What it covers |
|---|---|---|
| [installation.md](installation.md) | users, packagers | bioconda, static binaries, building from source with CMake, `just` or conda-build, the AVX2 launcher, what gets installed |
| [running.md](running.md) | users | where outputs go in each mode, exit codes, reruns and `--profile_only`, options the website does not list |
| [database-files.md](database-files.md) | users, database builders | the files of a database, the single-file `database.protal`, compression, converting and unpacking |
| [building-a-database.md](building-a-database.md) | database builders | a database from a GTDB release, reduced marker sets, the complete build-and-train workflow |
| [model-training.md](model-training.md) | database builders, developers | the presence model: what protal checks, training data, features, trainers, installing a new model |
| [qcmsa.md](qcmsa.md) | users | the strain MSA post-filter: what it removes, its parameters, re-filtering existing output |
| [simulation.md](simulation.md) | developers, benchmarking | `simulate_metagenomes`: mock communities, strain sharing, provenance and exact replays |
| [testing.md](testing.md) | developers | unit tests, end-to-end tests, the mini database and example, CI, the strain test harness |
| [claude/](claude/README.md) | developers | audits, reviews and benchmarks written with Claude Code |

The example in [`examples/mini_db/`](../examples/mini_db/README.md) builds a tiny database,
profiles simulated reads with it and checks the result, without any download.
