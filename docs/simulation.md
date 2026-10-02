# Simulating metagenomes

`simulate_metagenomes` draws mock communities from a table of genomes and simulates paired-end
Illumina reads for them with [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art)
(`art_illumina`). It is used to test protal, to benchmark it, and to generate the training data of
the presence model ([model-training.md](model-training.md)). It needs `art_illumina`.

Build it with the other binaries ([installation.md](installation.md)):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target simulate_metagenomes     # or: just simulate
```

## Input

A tab-separated genome table with three columns and no header: genome name (accession), GTDB
taxonomy string, and the path to the genome FASTA (`.gz` works); optionally a fourth, the genome's
length (its letters outside header lines). Without the lengths the simulator reads every genome of
the table for its length at the start of each run (~20 ms per genome of GTDB's size), whichever it
simulates. `scripts/build_gtdb_database.py` writes one for a GTDB release, with the lengths
(`OUTDIR/genomes.tsv`), and the synthetic releases of the mini database have one in
`simulation/genomes.tsv`.

## A run

```bash
./build/simulate_metagenomes \
  --genome_table genomes.tsv \
  --output_dir sims/ \
  --samples 3 \
  --sample_prefix sim \
  --total_read_pairs 100000 \
  --species_per_sample 15 \
  --distribution power_law \
  --strains_per_species "0.4,0.2" \
  --seed 1 -t 8 \
  --protal_metafile sims/protal
```

Reads are simulated per genome, concatenated per sample into `reads/<sample>_R1.fq.gz` and
`reads/<sample>_R2.fq.gz` (BGZF, compressed in process as each genome's reads are appended: no uncompressed copy is written), and the composition is recorded:

| Output | |
|---|---|
| `reads/` | the reads of each sample |
| `manifest.tsv` | one row per (sample, genome): read pairs, relative abundance, vertical coverage, the FASTA it came from and the `art_seed` ART used for it |
| `manifests/<sample>.tsv` | the same rows, one file per sample |
| `abundance_matrix.tsv` | relative abundance of each species (rows) in each sample (columns) |
| `run_params.tsv` | the command line, the RNG seed (recorded even when `--seed` was not given) and the ART settings |
| `protal.meta`, `protal_goldstd/<sample>.profile_truth` | with `--protal_metafile DIR`: a protal map for the reads, with outputs under `DIR` and a `PROFILE_TRUTH` column pointing at the true species of each sample |
| `plots/<sample>.png` | with `--plot_png`: abundance bar plots (needs `Rscript`) |

With `protal --db DB --map sims/protal.meta`, protal prints true and false positives and false
negatives per sample and writes `<profile>.truth_annotated`.

## Options

`simulate_metagenomes --help` lists them all.

| Option | Default | |
|---|---|---|
| `-n, --samples` | 1 | number of samples |
| `--sample_prefix` | `sample` | sample names are `<prefix>_<n>` |
| `--total_read_pairs` | 100000 | read pairs per sample |
| `--species_per_sample` | 10 | a number or an inclusive range, e.g. `20-80` |
| `--distribution` | `poisson_lognormal` | abundance model: `power_law` (`--alpha`), `negative_binomial` (`--nb_r`, `--nb_p`), `poisson_lognormal` (`--pln_mu`, `--pln_sigma`) |
| `--strains_per_species` | none | probabilities of a 2nd, 3rd, ... strain of a species, e.g. `0.4,0.2,0.1` |
| `--include_species` | | species to put in every sample |
| `--genus`, `--taxon` | | fixed designs: `g__A:10,g__B:2` species from those genera, `d__Archaea:10` from any taxon; the rest of the sample is drawn from all species |
| `--pick_random_demand_if_fail` | off | if `--genus`/`--taxon` ask for more species than exist, cap and fill randomly instead of failing |
| `--strain_sharing_file` | | cross-sample strain sharing, below |
| `--seed` | random | RNG seed |
| `--read_length`, `--fragment_mean`, `--fragment_stdev` | 150, 350, 50 | ART read and fragment sizes |
| `--sequencer` | `HS25` | ART error profile |
| `--extra_art_args` | | passed to ART, e.g. `"--qprof1 q1 --qprof2 q2"` |
| `--art_path` | on `$PATH` | (`--pigz_path` is accepted and ignored: the reads are compressed in process) |
| `-t, --threads` | 1 | samples simulated at a time, each on one thread (ART and the compression); the samples are the same for any number (their designs and ART seeds are drawn first, in order) |
| `--test` | off | write the design, manifests and truth, but no reads |
| `--keep_tmp` | off | keep the reads of each genome |
| `-v, --version` | | the version and the commit it was built from (as `protal --version`) |

### Strain sharing across samples

`--strain_sharing_file` takes a tab-separated file (`#` lines are comments) with one species per
row, to test strain tracking:

| Column | |
|---|---|
| `SPECIES` | the species, as in the genome table |
| `SAMPLE_FRACTION` | fraction of all samples that contain the species (0-1) |
| `N_STRAINS` | distinct strains of it drawn across all samples |
| `MIN_OCCURRENCE` | each strain appears in at least this many samples |
| `MIN_VCOV` | optional: minimum vertical coverage per strain and sample (0 = none) |
| `CONSPECIFIC_STRAINS` | optional: probabilities of a 2nd, 3rd, ... of these strains in one sample, as in `--strains_per_species` |

## Reproducing a simulated dataset

A manifest fully describes a dataset's composition (genomes, read pairs, ART seeds), so
replaying one does not depend on reproducing the community-design RNG:

```bash
./build/simulate_metagenomes \
  --from_manifest sims/manifest.tsv \
  --output_dir sims_replay/
```

The replay reproduces the reads byte for byte, provided it uses the same `art_illumina`
build and the same ART settings. The manifest stores the ART seeds but not the settings,
so `--read_length`, `--fragment_mean`, `--fragment_stdev`, `--sequencer` and
`--extra_art_args` still come from the command line; the replay warns about each one that
differs from the original run's `run_params.tsv`, and also checks `--read_length` against
the depths the manifest implies. All sampling options (`--distribution`,
`--species_per_sample`, `--seed`, ...) are ignored. A per-sample manifest replays just that
sample.

Manifests written before the `fasta_path` and `art_seed` columns existed still replay:
pass `--genome_table` so the genomes can be resolved by name, and the composition and
per-genome depth are reproduced exactly while the reads themselves are fresh
realizations. The replay's own manifest carries seeds, so it is exactly reproducible
from then on.

## Without ART

For quick tests, `scripts/mini_db/simulate_reads.py` draws error-prone read pairs straight from
the genomes of a mock community (substitutions only, no ART needed) and writes the truth table
next to them; the [mini database example](../examples/mini_db/README.md) uses it.
