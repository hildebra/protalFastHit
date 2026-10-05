# protal documentation

How to install protal, run it and read its output is described on the website:
https://protal.earlham.ac.uk/main.php?site=documentation (installation, database download, a
worked four-sample example, the map file, and every output file). The pages here cover the rest.

| Page | For | What it covers |
|---|---|---|
| [installation.md](installation.md) | users, packagers | bioconda, static binaries, building from source with CMake, `just` or conda-build, one binary for every CPU, what gets installed |
| [running.md](running.md) | users | where outputs go in each mode, exit codes, reruns and `--profile_only`, options the website does not list |
| [version_changes.md](version_changes.md) | users, developers | what changed from 0.6.0a to 0.7.6, and every version's detection, abundance, strain, speed and memory benchmarks on the same data |
| [database-files.md](database-files.md) | users, database builders | the files of a database, the single-file `database.protal`, compression, converting and unpacking |
| [building-a-database.md](building-a-database.md) | database builders | download, build and train a database in one command, for one or several GTDB releases, full or reduced; the steps one by one |
| [model-training.md](model-training.md) | database builders, developers | the presence model: what protal checks, training data, features, trainers, installing a new model |
| [features.md](features.md) | users, developers | every feature of the presence model: since which version, what it measures, its importance per read type at GTDB r226, and the situations it matters for |
| [qcmsa.md](qcmsa.md) | users | the strain MSA post-filter: what it removes, its parameters, re-filtering existing output |
| [simulation.md](simulation.md) | developers, benchmarking | `simulate_metagenomes`: mock communities, strain sharing, provenance and exact replays |
| [testing.md](testing.md) | developers | unit tests, end-to-end tests, the mini database and example, CI, the strain test harness |
| [claude/](claude/README.md) | developers | audits, reviews and benchmarks written with Claude Code |

The example in [`examples/mini_db/`](../examples/mini_db/README.md) builds a tiny database,
profiles simulated reads with it and checks the result, without any download.
