# protal documentation

How to install protal, run it and read its output is described on the website:
https://protal.earlham.ac.uk/main.php?site=documentation (installation, database download, a
worked four-sample example, the map file, and every output file). The pages here cover the rest.

| Page | For | What it covers |
|---|---|---|
| [installation.md](installation.md) | users, packagers | bioconda, static binaries, building from source, one binary for every CPU, the tools to build a database |
| [running.md](running.md) | users | where outputs go, read files and read types, reruns and `--profile_only`, exit codes, options the website does not list, memory |
| [strains.md](strains.md) | users | strain MSAs: which samples, genes and reads enter them, strains in long-read samples, filtering with qcmsa, trees |
| [databases.md](databases.md) | database builders | the files of a database and `database.protal`; building and training one from GTDB in one command, for several releases, reduced, or step by step; the presence model and how to train, check and install it |
| [features.md](features.md) | users, developers | every feature of the presence model: since which version, what it measures, its importance per read type at GTDB r226, and when it matters |
| [development.md](development.md) | developers | unit and end-to-end tests, mini databases, `simulate_metagenomes`, CI, the strain test harness |
| [version_changes.md](version_changes.md) | users, developers | what changed from 0.6.0a to 0.7.6, and every version's detection, abundance, strain, speed and memory benchmarks |
| [claude/](claude/README.md) | developers | audits, reviews and benchmarks written with Claude Code |

The example in [`examples/mini_db/`](../examples/mini_db/README.md) builds a tiny database,
profiles simulated reads with it and checks the result, without any download.
