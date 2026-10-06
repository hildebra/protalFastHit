# Protal

Protal profiles bacterial and archaeal communities from shotgun metagenomes (paired-end,
single-end, PacBio HiFi and Nanopore reads), and resolves strains across samples. It aligns reads to the marker genes of GTDB (the prebuilt
database covers GTDB r226), decides with a model of decision trees which species are present, and writes
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
| [Installation](docs/installation.md) | bioconda, static binaries, building from source |
| [Running protal](docs/running.md) | details beyond the website: outputs, read types, reruns, exit codes, less common options |
| [Strains](docs/strains.md) | strain MSAs, filtering them with qcmsa, building trees |
| [Databases](docs/databases.md) | database files; building and training a database from GTDB; the presence model |
| [The model's features](docs/features.md) | every feature, its importance per read type, and when it matters |
| [Development](docs/development.md) | tests, mini databases, `simulate_metagenomes`, CI |
| [What changed](docs/version_changes.md) | 0.6.0a to 0.7.8: changes and benchmarks of every version |
| [Audits and benchmarks](docs/claude/README.md) | reports written with Claude Code |

## Citation

Joachim Fritscher, Anthony Duncan, Falk Hildebrand. Protal: Ultra-fast metagenomic profiling and
strain-resolved analysis. bioRxiv. https://doi.org/10.64898/2026.08.03.742433

## License

GNU GPL version 2 (GPL-2.0-only), see [LICENSE](LICENSE). Developed in the [Hildebrand lab](http://falk.science) at the
Earlham Institute and the Quadram Institute.
