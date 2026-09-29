# Protal

Protal profiles bacterial and archaeal communities from paired-end shotgun metagenomes, and
resolves strains across samples. It aligns reads to the marker genes of GTDB (the prebuilt
database covers GTDB r226), decides with a random forest which species are present, and writes
species abundances and per-species multiple sequence alignments for strain phylogenies.

- Website and user documentation: https://protal.earlham.ac.uk
  ([documentation](https://protal.earlham.ac.uk/main.php?site=documentation),
  [downloads](https://protal.earlham.ac.uk/main.php?site=downloads),
  [FAQ](https://protal.earlham.ac.uk/main.php?site=faq))
- Further documentation, for what the website does not cover: [`docs/`](docs/README.md)

## Quick start

protal runs on Linux (x86-64). Install it from bioconda, download the database, and profile a
sample:

```bash
conda create -n protal protal -c bioconda -c conda-forge
conda activate protal

wget https://protal.earlham.ac.uk/data/dbs/protal-db-r226-0.5.1a.tar.gz
tar -xzf protal-db-r226-0.5.1a.tar.gz          # -> protal_0.5.1_r226/

protal --db protal_0.5.1_r226 -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
    --prefix sample --outdir results -t 8
```

`results/sample.profile` lists the species found and their relative abundances. Many samples
are run with a map file (`--map`, see `protal --map_help`); strain MSAs are written for species
found in at least two samples. The website's
[documentation](https://protal.earlham.ac.uk/main.php?site=documentation) walks through a
four-sample example and describes every output file.

## Documentation

| Page | What it covers |
|---|---|
| [Installation](docs/installation.md) | bioconda, static binaries, building from source, the `protal` launcher |
| [Running protal](docs/running.md) | details beyond the website: where outputs go, exit codes, reruns, less common options |
| [Database files](docs/database-files.md) | `database.protal`, compression, converting and unpacking a database |
| [Building a database](docs/building-a-database.md) | a custom database from a GTDB release, marker subsets, the full build-and-train workflow |
| [Training the model](docs/model-training.md) | the presence model: what protal expects, training data, features, trainers |
| [Strain MSA filtering (qcmsa)](docs/qcmsa.md) | the post-filter of the strain MSAs and how to re-filter without re-running protal |
| [Simulating metagenomes](docs/simulation.md) | `simulate_metagenomes`: mock communities, strain sharing, reproducible replays |
| [Testing and development](docs/testing.md) | unit and end-to-end tests, the mini database, CI, the strain test harness |
| [Audits and benchmarks](docs/claude/README.md) | reports written with Claude Code |

## Citation

Joachim Fritscher, Anthony Duncan, Falk Hildebrand. Protal: Ultra-fast metagenomic profiling and
strain-resolved analysis. bioRxiv. https://doi.org/10.64898/2026.08.03.742433

## License

GNU GPL version 2 (GPL-2.0-only), see [LICENSE](LICENSE). Developed in the [Hildebrand lab](http://falk.science) at the
Earlham Institute and the Quadram Institute.
