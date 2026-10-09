# What changed from protal 0.6.0a to 0.7.9

The shipped 0.6.0a (tag `0.6.0a`, commit `014f4a9`, August 2026) and the 0.7 series (0.7.0 on
2026-09-30 to 0.7.9 on 2026-10-07, branch `audit-fixes`) are compared here: what changed, and how
every version scores on the same simulated data in detection, abundance, strains, speed and memory.
The numbers come from the benchmark and audit reports in [`docs/claude/`](claude/README.md), which
hold the commands, scripts and raw tables; this page only collects them. Nothing here was run again
for this page.

| Version | Commit | Date | Commits since the one before | In short |
|---|---|---|---|---|
| 0.6.0a | `014f4a9` | 2026-08-10 | | paired-end reads only, one model, raw database files |
| 0.7.0 | `4b21427` | 2026-09-30 | 98 | the audits' fixes; single-end and long reads; the single-file database; twice as fast; a database build-and-train pipeline |
| 0.7.1 | `1c11a00` | 2026-10-01 | 30 | mates paired across neighbouring genes; the reads' other candidates as model features; depth margin 0.08; genes at 2 bits per base |
| 0.7.2 | `fd1e719` | 2026-10-02 | 21 | divergence and conservation features; one sample profiled on all threads; `gene_congeners.tsv` |
| 0.7.3 | `b146951` | 2026-10-02 | 17 | one binary for every CPU; gene neighbours from every genome; strain rows from long reads; HiFi training reads |
| 0.7.4 | `8b2e376` | 2026-10-03 | 23 | species calls judged against the sample's congeners (relatives, distance, read EM); long reads through their chain; prefetched lookups |
| 0.7.5 | `b892133` | 2026-10-03 | 17 | against false positives: the sample's depth and the reads' divergence as default features, suspect gene copies, species priors; the index packed in memory |
| 0.7.6 | `13aac36` | 2026-10-05 | 17 | GTDB-scale speed (paired-end run 131 → 58 s, HiFi 104 → 25 s on r226); outputs independent of record order; species priors opt-in; reduced marker sets; in-silico strains and denser depths in training |
| 0.7.7 | `f2bafd9` | 2026-10-05 | 33 | GTDB-scale speed and memory (paired-end run 58 → 42 s, HiFi 25 → 15.6 s, peak 38.1 → 34.2 GB on r226): AVX2 flex scan, unaligned reads counted in the SAM header (`--write_unmapped_reads`), the index packed key by key as it loads; read EM capped at 100 sweeps; `.fq.zst` input |
| 0.7.8 | `5de3224` | 2026-10-06 | 16 | GTDB-scale speed (paired-end run 42 → 33 s, HiFi 15.6 → 12.3 s on r226, the same outputs): AVX2 tie masks in the seed lookup, flat gene tables, the k-mer screen from packed genes and reads, only shared seeds sorted; gzip with ISA-L (paired-end runs no longer wait for their input); `--add_model` in place; gradient-boosted models by default |
| 0.7.9 | `ffbedb3` | 2026-10-07 | 37 | against false positives in complex communities: the reads' consistency, the genes' shape, the database's neighbourhood, the sample's complexity and the ancestry sites (which side a read takes where the species differs from its congeners) as default features, untested at r226 until the next build; strain MSAs packed (159 → 37 MB per dense sample), spilled to disk and merged over runs; builds keep the reads behind every model error (`--error-reads`, `--share-logs` with the ancestry report), calibrate the in-silico strains' dN/dS on real strains, vary the scenarios and choose the knob by bootstrap; `simulate_metagenomes` makes Illumina, Ultima and long reads itself and streams large samples; `--build` resumes, packs the index straight into `database.protal`; every test suite in CI |

## Changes since 0.6.0a

### Reads and input

- Single-end (`se`), PacBio (`pb`) and Nanopore (`ont`) reads besides paired-end, each profiled with
  a model of its own (`--read_type`, or `READ_TYPE` in a map file). Long reads are aligned per gene,
  reads over 65 kb in chunks; since 0.7.4 through every link of their chain (31-38% fewer alignment
  instructions, the same F1). A single read file's type is told from its first reads (0.7.3).
- Pipes, gzip, BGZF and FASTA input. Input is inflated outside the reader lock (libdeflate and a
  vendored zlib-ng; since 0.7.8 ISA-L, 1.6-1.8x as fast on single-member gzip). Unusable inputs fail
  their sample instead of passing silently.

### The database

- One file, `database.protal` (seekable zstd, the index in columns), checked against its reference;
  107 MB instead of 3.4 GB for the 765-species tuning world, built 3x faster. The old raw files still
  work and can be converted both ways ([databases.md](databases.md#the-files-of-a-database)).
- `--add_model` installs a read type's model into the database (several at once since 0.7.3).
- Index format 2 with correct unique flags, built in parallel and the same for any thread count.
  No k-mer with an ambiguous base enters the index.
- Per-database tables that the profiling uses: the genes' within-species conservation
  (`gene_conservation.tsv`, 0.7.1), how the genes differ between congeners (`gene_congeners.tsv`,
  0.7.2), gene neighbours from every genome with lines for a species' own gene order
  (`gene_neighbours.tsv`, 0.7.1, in the database file since 0.7.3), gene copies near-identical to
  another genus's (`suspect_copies.tsv`, 0.7.5) and species priors from GTDB's clusters
  (`species_priors.tsv`, 0.7.5).
- In memory the reference genes are held at two bits per base (0.7.1; 19.5 GB → 5.1 GB at GTDB
  r226 size) and the index's values packed at 42 bits (0.7.5; 35 → 27 GB at r226). The files are
  unchanged by either.
- `gene_table.bin` (0.7.6): `reference.map` and `unique_kmers.tsv` in binary, with the reference's
  fingerprint, as a member of `database.protal` beside the text tables. A run loads it on all threads
  without parsing (r226 size, six threads: 1.76 → 0.43 s). `--build` writes it; `--compress_db` adds
  it to an older single-file database; without it the text tables are read as before.
- Reduced databases (0.7.6): a subset of the marker genes (`--build_gene_subset`, the converter's
  `--genes`) now reaches every build phase. `scripts/rank_genes.py` picks the most distinctive genes
  with a share reserved for each domain, so that archaea keep genes
  ([databases.md](databases.md#building-a-database)).

### Profiles and strains

- Abundance from a species' own reads: `--depth_identity_margin` (default 0.08 since 0.7.1; 0.7.0's
  0.04 dropped 6-26% of the bases of strains 1-4.5% from the reference and biased abundances), a
  fragment's bases counted once where the mates overlap (0.7.3), foreign genes left out of depth
  and MSAs (`--keep_foreign_genes`, 0.7.3).
- Mates paired across neighbouring genes, and a mate rescued on the neighbour (0.7.1).
- Species calls: the presence model's features grew from the reads' counts and identities to the
  reads' other candidates (MAPQ, congener fits; 0.7.1), their divergence beyond what their base
  qualities expect and the conservation of the genes they hit (0.7.2), the sample's congeners
  (relatives and distance features, a read EM over the reads' alternatives; 0.7.4), the sample's
  depth and the evidence before the MAPQ filter (0.7.5). Species priors (`+priors`) are opt-in.
- Variants called per fragment with strand tests, deletions and unbiased depth; MSAs chosen by
  `--msa_knob`; strain rows from long reads with phasing (`--no_phasing` turns it off; 0.7.3).
  0.6.0a wrote 0-based partition files, which put every gene boundary one column off in IQ-TREE.
- Failures are reported and give a non-zero exit; outputs are written crash-safe.
- Outputs no longer depend on the order of the SAM's records (0.7.6), which multi-threaded alignment
  varies: sums of doubles are exact (fixed point), and ties break by taxid. Before, repeated runs of
  one build could differ in the last digit of a value or two
  ([order independence](claude/2026-10-05-order-independence/README.md)).
- The per-taxon files `misc/<taxon>.statistics.tsv` are written only with `--taxon_statistics`
  (0.7.6). They took 15 s of a 107 s r226 run on a network file system; the per-sample profile files
  hold the same numbers.
- After 0.7.9 (2026-10-08): the sample's composition. `species_priors.tsv` gives each species' genome
  size (GTDB's genomes, corrected by CheckM) and its marker genes' share of it; the aligner writes the
  reads it read into the SAM header; a profile ends with the share of the genomes its species do not
  explain (`?`), and its abundances are shares of all genomes (`--no_unknown_share`: of the called
  species, as before). `<profile>.composition` gives the share of the reads explained, the average
  genome size and its quantiles, and how many species the rest would be
  ([running.md](running.md#what-the-called-species-explain-the-unknown-share)). Older databases get the
  sizes with `gtdb_to_protal_db.py --priors_only` and `--add_tables`.

### Speed

- Short reads aligned from their anchors' exact matches, with flanks of up to three mismatches
  aligned without WFA2 (0.7.3 to 0.7.5); branch-free and AVX2 syncmer scans from 2-bit codes;
  k-mer lookups prefetched within the read and one read ahead (0.7.4).
- SAMs compressed in the alignment threads (`.sam.zst` by default), parsed as views and
  decompressed on threads of their own (0.7.2). One sample is profiled on all threads (0.7.2) and
  the profiling pass itself got cheaper (coverage as a difference array, the read EM on indices;
  0.7.4 and 0.7.5).
- One binary for every CPU: 13 hot functions are also compiled for x86-64-v3 (AVX2) and chosen at
  run time (0.7.3). Huge pages for the index where the kernel allows them.
- The database is loaded and the index built in parallel; since 0.7.6 the start-up's parts load
  side by side (the index beside the genome preload, the small tables each on a thread;
  `--sequential_load` restores the old order) and the gene tables are added per genome on all
  threads.
- At GTDB scale (0.7.6): a k-mer screen refuses a candidate before WFA2 when read and gene window
  share too few k-mers for any alignment within the score budget (exact; `--no_alignment_screen`
  turns it off). A long read's window ends at its chain's last link, so most long reads align through
  their chain. The failed candidates of the `ZF` tag are counted in vectors by taxid, not in a hash
  map per chunk; at r226 the paired-end SAM holds 45M unmapped records naming 140,000 taxa.
- Every run prints per sample the reads, the candidates tried, refused and aligned, and the records
  written. Every stage, the start-up and the teardown are timed (0.7.6).

### Training a database

- `build_gtdb_database.py` builds and trains a database from a GTDB release: held-out clades and
  species, other strains, an independent test set, one model per read type, resumable reruns
  ([databases.md](databases.md#building-a-database), [databases.md](databases.md#the-presence-model)).
  0.6.0a shipped one model trained on older databases.
- After 0.7.9 (2026-10-08): the hold-out keeps species complexes whole (congeners within 0.01 on
  the references, `species_clouds.tsv` from `protal --write_species_neighbours` before anything is
  built; the training table's `meta_novel_distance` and the report's errors by that distance), and
  the foreign scan tiles each database's full reference (every species alike, every copy listed)
  instead of the genomes the samples are drawn from, which had told the r226 v17 models the
  simulation's species ([features.md](features.md#the-foreign-features-leak)); after the r226 v18
  build (no leak, 0.001 of AUC) the build scans only with `--foreign-rates`; and the ancestry sites
  are the consensus of the congeners (nine in ten of three or more compared, the nearest congener
  alone with fewer) rather than the nearest congener's differences, the comparison is chained
  across indels instead of ending at the first, and the indels the congeners share against the
  species are sites too, two more features (80 in the default set)
  ([features.md](features.md#which-side-the-reads-take-ancestry-079)). The build's `--outdir` holds
  `protal_db/`, `model_logs/`, `logs/`, `work/` and `console.log`, each file once: ~118 files at
  r226 instead of ~11,900 (the in-silico strains on the samples' disk, one error-read table per read
  type; [outputs](databases.md#outputs), [report](claude/2026-10-08-build-outputs.md)). Strain
  alleles: `--build` keeps up to 4 alleles of each species' other genomes per gene
  (`strain_alleles.tsv`, edits of the representative's copy, each nearer it than the nearest
  congener's), a short read's candidates are scored with their species' best allele
  (`--no_allele_scores`: not), and two `alleles` features join the default set (82); the GTDB
  build takes the alleles from half of the genomes by a hash of the accession and simulates strains
  from the other half only, so that no simulated strain is its own allele
  ([features.md](features.md#the-strain-alleles-alleles-2026-10-08)).
- After 0.7.9 (2026-10-09): the species' polymorphic sites (where one of its alleles differs). In the alignment, a
  read's difference that any of the species' alleles has counts as a match, and one at a polymorphic site with another
  base counts as half a difference. These shifts decide only the reads protal is unsure about (a candidate of another
  species within 3 mismatches). Three `polymorphic` features join the default set (85): of the polymorphic sites a
  taxon's reads cover, the shares with a known allele's base and with another, and the species' base at the ancestry
  sites fixed within the species against all of them
  ([features.md](features.md#the-polymorphic-sites-polymorphic-2026-10-09)). The in-silico strains put as many of their
  marker substitutions at their nearest congener's sites, with its base, as the real strains have there
  (`species_clouds.tsv` is now made before them; [report](claude/2026-10-09-r226-v19/README.md), section 7).
- 0.7.3 to 0.7.5 made the training data more realistic (PacBio HiFi reads with qualities by length
  instead of pbsim3's quality-0 reads; a mixed design of lognormal sigma 1.3 and 2.0 communities;
  30% of the species held out; deeper design points) and the build faster (one simulation queue,
  deep samples in chunks, the training database on `--scratch`, both collections profiled in one
  protal run).
- `--version` names the commit a binary was built from, which the build script checks.
- 0.7.6 changed the training defaults and design:
  - The species priors are opt-in (`+priors`); their gain at r226 was the cluster size alone, a rule
    the simulation cannot test.
  - Every species with one genome gets an in-silico strain to be simulated from: a codon-aware
    mutated copy of its representative, as far from it as the real strains are from theirs
    (`scripts/insilico_strains.py`, `--insilico-strains`).
  - Training depths of 2,000, 50,000 and 200,000 read pairs are added; r226 v10's test errors sat
    between the old points.
  - Databases of several GTDB releases, full and reduced, in two phases
    (`download_gtdb_releases.py`, `build_gtdb_releases.py`).
  - Simulated worlds can evolve each gene at the real r226 gene's speed
    (`simulate_gtdb_release.py --gene_rates r226`).
  - `trace_relatives.py` follows reads by their contig, so the build's gene-conservation trace works
    on GTDB's genomes.
  ([r226 v10 evaluation](claude/2026-10-04-r226-v10-evaluation/README.md),
  [mini database genes](claude/2026-10-05-mini-database-genes/README.md),
  [trace](claude/2026-10-05-trace-relatives.md)).

## Benchmarks

All versions were run on one simulated world: 900 species (12% archaea, 135 in no database),
marker genes in operon-like clusters, genes of different conservation
([v0.7.1 benchmark](claude/2026-10-01-v071-benchmark/README.md)). Each 0.7 version used the
databases and models its own pipeline trained on that world; 0.6.0a used its own databases and the
model it ships. "full" means every sampled species is in the database, "missing" that the training
database lacks 290 species (349 for 0.7.5, whose holdout design changed, so its missing rows are
not like-for-like). 0.7.4 never had a database of its own; 0.7.5 holds its changes. Machine: WSL2,
Intel Core Ultra 7 258V, 6 threads, 23 GB. The rows for 0.7.0 to 0.7.2 are the
[v0.7.3 benchmark](claude/2026-10-02-v073-benchmark/README.md)'s, those for 0.6.0a, 0.7.3 and
0.7.5 the [v0.7.5 benchmark](claude/2026-10-03-v075-benchmark/README.md)'s; the reruns of older
versions reproduced their calls to the fourth digit.

0.7.6 was not rerun on this world. Its short-read changes leave the alignments and profiles
identical to 0.7.5's (checked on 500k and 5M pairs), apart from the last digit of the few values
that depended on the SAM's record order. Long reads differ where a window's end moved (on the test
reads 51 → 56 of 60 now align through their chain; one HiFi read of the r226 sample). Its training
changes show in the GTDB models below.

### Detection and abundance, short reads

Means over 26 paired-end samples (1,000 to 5M pairs) and 24 single-end samples; F1 of the species
reported at knob 0.5, false positives per sample, Bray-Curtis dissimilarity of the abundances to
the truth (0 is identical), recall over all present archaea:

| reads | database | version | F1 | precision | recall | FP per sample | Bray-Curtis | archaea recall |
|---|---|---|---|---|---|---|---|---|
| pe | full | 0.6.0a | 0.634 | 0.993 | 0.556 | 1.00 | 0.345 | 0.242 |
| pe | full | 0.7.0 | 0.903 | 0.976 | 0.869 | 2.96 | 0.131 | 0.736 |
| pe | full | 0.7.1 | 0.908 | 0.980 | 0.873 | 2.12 | 0.087 | 0.747 |
| pe | full | 0.7.2 | 0.911 | 0.982 | 0.875 | 1.88 | 0.085 | 0.749 |
| pe | full | 0.7.3 | 0.913 | 0.982 | 0.876 | 1.35 | 0.084 | 0.742 |
| pe | full | **0.7.5** | **0.915** | **0.988** | 0.876 | **1.27** | 0.085 | 0.745 |
| pe | missing | 0.6.0a | 0.629 | 0.982 | 0.560 | 1.50 | 0.357 | 0.225 |
| pe | missing | 0.7.0 | 0.895 | 0.963 | 0.865 | 2.00 | 0.131 | 0.697 |
| pe | missing | 0.7.1 | 0.900 | 0.967 | 0.869 | 1.54 | 0.093 | 0.704 |
| pe | missing | 0.7.2 | **0.906** | 0.977 | 0.869 | 1.08 | 0.088 | 0.700 |
| pe | missing | 0.7.3 | 0.902 | 0.963 | 0.872 | 1.46 | 0.089 | 0.708 |
| pe | missing | 0.7.5 (349 species out) | 0.902 | 0.966 | 0.872 | 1.77 | 0.098 | 0.706 |
| se | full | 0.7.0 | 0.874 | 0.975 | 0.830 | 2.88 | 0.145 | 0.682 |
| se | full | 0.7.1 | 0.888 | 0.981 | 0.844 | 1.92 | 0.107 | 0.696 |
| se | full | 0.7.2 | 0.888 | 0.982 | 0.843 | 1.75 | 0.107 | 0.700 |
| se | full | 0.7.3 | **0.893** | 0.984 | 0.847 | 1.08 | 0.106 | 0.708 |
| se | full | **0.7.5** | **0.893** | **0.987** | 0.847 | 1.08 | 0.106 | 0.702 |
| se | missing | 0.7.0 | 0.866 | 0.962 | 0.826 | 2.04 | 0.146 | 0.648 |
| se | missing | 0.7.1 | 0.881 | 0.973 | 0.836 | 1.42 | 0.113 | 0.648 |
| se | missing | 0.7.2 | 0.884 | 0.977 | 0.839 | 1.33 | 0.109 | 0.648 |
| se | missing | 0.7.3 | **0.886** | 0.978 | 0.838 | 1.21 | 0.111 | 0.652 |
| se | missing | 0.7.5 (349 species out) | 0.881 | 0.962 | 0.846 | 1.71 | 0.113 | 0.647 |

0.6.0a profiles paired-end reads only. The gain from 0.6.0a to 0.7.0 is the model's: 0.7.0 with
0.6.0a's model detects as 0.6.0a does, so it is what a database trained for its own reference gives
([0.7 against 0.6.0a](claude/2026-09-30-v07-vs-v06/README.md)). The step from 0.7.0 to 0.7.1 in
Bray-Curtis (0.131 → 0.087) is the depth margin 0.08. Per sample, 0.7.1 over 0.7.0 is +0.0055 F1
(significant), 0.7.2 over 0.7.1 +0.0022 (significant), 0.7.3 over 0.7.2 +0.0020 and 0.7.5 over
0.7.3 +0.0019 (both within noise); 0.7.3 over 0.6.0a is +0.279 (+0.185, +0.370).

By depth, F1 / Bray-Curtis with the full database (from the
[v0.7.1 benchmark](claude/2026-10-01-v071-benchmark/README.md); the later versions change these
rows by at most a few thousandths):

| read pairs | 0.6.0a | 0.7.0 | 0.7.1 |
|---|---|---|---|
| 1,000 | 0.242 / 0.588 | 0.749 / 0.192 | 0.755 / 0.155 |
| 10,000 | 0.602 / 0.444 | 0.949 / 0.094 | 0.954 / 0.063 |
| 500,000 | 0.964 / 0.040 | 0.984 / 0.076 | 0.984 / 0.028 |
| 5,000,000 | 0.986 / 0.043 | 0.980 / 0.090 | 0.992 / 0.030 |

Archaea in samples of 10,000 pairs: 0.6.0a found 3%, the 0.7 versions 82-83%; at 500,000 pairs 63%
against 98-100%. Species at 0.01-0.1% abundance in 10,000 pairs: 8% against 84-85%.

### Long reads

F1 (false positives per sample) on 4 PacBio HiFi and 4 Nanopore samples of 3 and 90 Mb each, made
by 0.7.3's collector (HiFi reads with qualities by length, which is what real HiFi reads carry;
0.7.0 to 0.7.2 were trained on pbsim3 reads with every base quality 0):

| reads | database | 0.7.0 | 0.7.1 | 0.7.2 | 0.7.3 | 0.7.5 |
|---|---|---|---|---|---|---|
| pb (HiFi) | full | 0.818 (3.00) | 0.824 (1.38) | 0.817 (2.75) | 0.826 (0.62) | 0.826 (**0.12**) |
| pb (HiFi) | missing | 0.792 (2.50) | 0.798 (1.75) | 0.793 (2.62) | 0.800 (0.88) | **0.811** (**0.12**) |
| ont | full | 0.865 (1.25) | 0.870 (0.88) | 0.872 (0.38) | 0.872 (0.62) | 0.873 (0.38) |
| ont | missing | 0.843 (2.38) | 0.854 (1.75) | 0.845 (2.50) | 0.858 (1.12) | 0.860 (0.75) |

0.7.5 over 0.7.3 on PacBio with species missing is +0.0113 per sample (significant), with
0.75 fewer false positives and closer abundance. The long-read F1 is lower than the short-read one
because half the samples are 3 Mb, where the depth decides.

### Strains

Twelve species of 8 strains each, one strain per sample at 2-40x, each version against its own
database; SNPs of the strain MSAs after qcmsa against the true genotypes, trees by IQ-TREE against
the true topology:

| reads | version | MSA columns (raw / qcmsa) | SNP precision | SNP recall | SNP F1 | miscalled SNPs per 100 kb | true topology |
|---|---|---|---|---|---|---|---|
| pe | 0.6.0a | 114,125 / 108,841 | 1.000 | 0.768 | 0.869 | 0.0 | 4 of 12 |
| pe | 0.7.0 | 114,640 / 114,416 | 0.998 | 0.937 | 0.967 | 2.0 | 10 of 12 |
| pe | 0.7.1 to 0.7.5 | 114,695 / 114,695 | 0.998 | 0.939 | 0.967 | 1.9 | 11 of 12 |
| pb | 0.7.0 to 0.7.2, 0.7.3 and 0.7.5 `--no_phasing` | 115,220 / 114,112 | 0.951-0.953 | 0.867-0.868 | 0.907-0.908 | 45.8-47.9 | 6-7 of 12 |
| pb | 0.7.5 (phasing on) | 115,171 / 114,064 | 0.948 | 0.861 | 0.903 | 50.9 | 8 of 12 |
| ont | 0.7.0 to 0.7.5 | 114,681 / 114,110 | 0.994 | 0.889 | 0.939 | 5.7 | 8 of 12 |

0.6.0a left many true SNPs as N at low depth. On the strain audit's world (8 species, 1-50x), its SNP
recall at covered true SNPs was 0.57 at 1x and 0.65 at 2x against 0.993 and 0.996 for 0.7.0, and its
qcmsa removed more cells the deeper the sample (0.78 of the positions called at 50x against 0.96).
The PacBio samples hold one strain each, so 0.7.3's and 0.7.5's phasing, which split 8 of 96
sample rows, is wrong there and is scored without those rows.

### Speed

Mean wall / CPU seconds per run at `-t 6`, both databases, no qcmsa, from the v0.7.5 benchmark
(0.7.0 to 0.7.2 from the v0.7.3 benchmark, whose deep runs shared the machine with a build, so
their 5M-pair and long-read rows are inflated; 0.6.0a, 0.7.3 and 0.7.5 ran in turn on an otherwise
idle machine). Repeated runs vary by 10-20%.

| reads | depth | 0.6.0a | 0.7.0 | 0.7.1 | 0.7.2 | 0.7.3 | 0.7.5 |
|---|---|---|---|---|---|---|---|
| pe | 1,000 pairs | 3.1 / 3 | 1.1 / 4 | 1.1 / 4 | 1.1 / 4 | 0.9 / 3 | 1.0 / 3 |
| pe | 10,000 pairs | 3.9 / 4 | 1.1 / 4 | 1.1 / 4 | 1.0 / 4 | 1.3 / 4 | 1.5 / 4 |
| pe | 500,000 pairs | 13.2 / 36 | 7.6 / 32 | 6.5 / 25 | 5.1 / 25 | 6.1 / 28 | 5.7 / 25 |
| pe | 5M pairs | 59.8 / 276 | 76.0 / 339 | 81.0 / 311 | 58.4 / 307 | 30.5 / 167 | 29.7 / 158 |
| se | 500,000 reads | - | 4.0 / 16 | 3.7 / 14 | 2.9 / 12 | 3.0 / 13 | 3.3 / 13 |
| pb | 90 Mb | - | 10.1 / 54 | 11.1 / 60 | 10.6 / 59 | 6.0 / 32 | 4.7 / 25 |
| ont | 90 Mb | - | 10.4 / 54 | 10.4 / 54 | 10.0 / 55 | 9.1 / 49 | 7.2 / 39 |

Like-for-like comparisons on a quiet machine, from the reports that made them:

| comparison | data | 0.6.0a or older | newer | source |
|---|---|---|---|---|
| 0.6.0a → 0.7.0 | 5M pairs 2x150, 150 species, tuning world | 100 s (310 CPU s) | 49.1 s (187 CPU s) | [0.7 against 0.6.0a](claude/2026-09-30-v07-vs-v06/README.md) |
| 0.6.0a → 0.7.0 | 500k pairs 2x150 | 10.2 s | 6.0 s | same |
| 0.6.0a → 0.7.0 | strain run A, 42 samples with qcmsa | 116 s | 36.6 s | same |
| 0.6.0a → 0.7.0 | index load, 5M-pair run | 9.4 s | 1.7 s | same |
| 0.6.0a → 0.7.5 | 5M pairs, operon world | 59.8 s | 29.7 s | v0.7.5 benchmark |
| 0.6.0a → 0.7.5 | 500k pairs, operon world | 13.2 s | 5.7 s | v0.7.5 benchmark |
| 0.7.3 → 0.7.5 | 90 Mb Nanopore / PacBio HiFi | 9.1 / 6.0 s | 7.2 / 4.7 s | v0.7.5 benchmark |

At GTDB scale: 0.7.5 against 0.7.6 on the r226 database (143,614 species, 27 GB file) with real reads
(a gut paired-end sample, a PacBio HiFi barcode; the read counts were not recorded), AMD EPYC 9634 at
32 threads, warm runs, medians of three
([performance at GTDB scale](claude/2026-10-04-performance-gtdb-scale/README.md)). The 0.7.6 column
is the fourth cluster run (`8c37dd6`, wall time; the stage times from the third, `7693f28`), before
`gene_table.bin` and the concurrent start-up:

| | paired-end 0.7.5 | paired-end 0.7.6 | HiFi 0.7.5 | HiFi 0.7.6 |
|---|---|---|---|---|
| wall | 131.5 s | **58.5 s** (−55%) | 103.7 s | **24.6 s** (−76%) |
| aligning | 61.9 s | 38.3 s (the k-mer screen) | 83.7 s | 11.8 s (screen and chain end) |
| profiling | 36.0 s | 6.4 s (failed candidates by taxid) | 3.4 s | 0.5 s |
| per-taxon statistics files | 14.9 s | 0 (`--taxon_statistics`) | 1.6 s | 0 |
| gene tables | (untimed) | 3.4 s (6.0 s before the parallel load) | (untimed) | 3.4 s |
| untimed | 22-26 s | under 1 s | 9 s | under 1 s |

The remaining fixed cost per run is about 11 s (13.8 s in the third run): the gene tables, the
genome preload (3.2 s), the index (3.2 s warm, 131 s cold from NFS) and the small tables (1.4 s).
`gene_table.bin` and the concurrent start-up should take it lower; that run has not been made. Seeding (20.6 s per thread) is what is
left of the paired-end alignment.

Where the time went, by performance round (callgrind instruction counts do not depend on the
machine's load; wall times are one thread unless said):

- [Round 1](claude/2026-09-29-performance-profiling/README.md) (before 0.7.0): gzip was inflated
  inside the reader lock (1M pairs at 8 threads: 21 s gzipped, 8.5 s plain), the syncmer test was
  recomputed per window (40% of the loop on a realistic mix), the 3 GB key map was read through
  4 KB pages (huge pages: seeding 5.9 → 3.7 s), the model's parser threw 229k exceptions at start.
  Fixed in 0.7.0: threaded gzip reader, branch-free and AVX2 syncmers, anchored alignment, zstd SAM
  output, zlib-ng, parallel index passes, huge pages, cPMML patches.
- [Round 2](claude/2026-09-30-performance-round2/README.md) (just before 0.7.1, against a build of
  the day before 0.7.0): aligning 1M realistic pairs at one thread 18.0 → 9.5 s, at six 5.6 → 3.2 s;
  then three cheap fixes (no per-read clock
  reads, a table reverse complement, no per-anchor read copy) took another 19% of the wall time, the
  windowed gene decode and syncmers from 2-bit codes a few percent more.
- [Round 3](claude/2026-10-02-performance-round3/README.md) (0.7.3 → 0.7.4): flanks of one
  mismatch without WFA2 and a faster index decoder (−7.7% of a 100k-pair run's instructions);
  long reads through their chain (−31% to −38% of the long-read alignment's instructions, 15-21%
  less CPU on 90 Mb samples). [Prefetching](claude/2026-10-02-prefetch/README.md) the k-mer
  lookups: about 4-5% of alignment. The [one binary](claude/2026-10-02-one-binary/README.md):
  −3.5% instructions on short reads against plain x86-64, and the launcher is gone.
- [Review](claude/2026-10-03-performance-review/README.md) and
  [pass 2](claude/2026-10-03-performance-pass2/README.md) (0.7.4 → 0.7.5): the profiling stage had
  grown to 30-40% of a short-read run with the congener features (5M pairs at six threads 34.0 s,
  of which profiling 12.6 s, the read EM 7.5 s on one thread); the EM on indices and cheaper
  congener sketches took it back (profiling −65% at one thread), then coverage as a difference
  array and 2-3 mismatch flanks without WFA2 (profiling −16.7% of instructions, aligning −2.4%).
  At 0.7.5 a 500k-pair run at one thread is 13.0 s: aligning 9.4 s, profiling 2.0 s, loading 1.9 s.

The [multithreading audit](claude/2026-10-01-multithreading-audit/README.md) found the locks
cheap at 0.7.1 (the reader lock 1.3-1.5% of thread time at six threads, the SAM writer nothing);
what did not scale was one sample profiled on one thread, fixed in 0.7.2.

### Memory

Peak RSS of a run on the operon world (765-species database with gene neighbours, 6 threads):

| version | peak RSS | source |
|---|---|---|
| 0.6.0a | 3.21-3.39 GB | v0.7.1, v0.7.2 and v0.7.5 benchmarks |
| 0.7.0 | 3.48-3.56 GB | same |
| 0.7.1 | 3.44-3.61 GB | same |
| 0.7.2 | 3.44-3.49 GB | same |
| 0.7.3 | 3.45-3.48 GB | v0.7.5 benchmark |
| 0.7.5 | 3.42-3.45 GB | v0.7.5 benchmark |

At this size the index's fixed 3.2 GB key map is most of the memory, so the versions hardly differ
(0.7.0 takes about 0.15 GB more than 0.6.0a, the later versions a little less than 0.7.0). The changes
that matter are at GTDB scale, where the index values and the genes dominate. Breakdown of a full
GTDB r226 run (16.6M genes, 2.90 G index entries), estimated from measured per-value and per-gene
costs except where said ([memory profiling](claude/2026-09-30-memory-profiling/README.md),
[the 2-bit gene store](claude/2026-09-30-gene-store/README.md),
[memory audit](claude/2026-10-03-memory-audit/README.md)):

| part | 0.6.0a and 0.7.0 | 0.7.1 to 0.7.4 (2-bit genes, `Gene` 112 → 32 bytes) | 0.7.5 (index values packed at 42 bits) |
|---|---|---|---|
| index values | 35 GB (measured on the r226 v7 index: 35.0) | 35.0 | **27.3** (computed from the measured slot count) |
| key map | 3.2 | 3.2 | 3.2 |
| reference genes | 17 (1 byte per base) | **4.3** | 4.3 |
| gene tables | 1.9 | **~0.8** | ~0.8 |
| **run, about** | **~58-59 GB** | **~43 GB** | **~35 GB** |

The 59 GB is the figure the website's download page still gives for 0.6. Measured on r226 runs of
0.7.5 (the four cluster runs above, both read types): **38.0 GB** peak RSS, of which the index
27.3 + 3.2 GB as loaded and the genes, tables and run buffers 7.5 GB. 0.7.6 has not been measured
there; its changes add no table in memory. On the operon world the packed index saved 70 MB of
213 MB of values with byte-identical outputs. A database build at r226 holds the index and the reference
(peak 88 GiB with both the training and the final index built at once in the first r226 run; about
50 GB with `--one-build-at-a-time`). Reads, threads and depth add little (≤ 20 MB per thread;
profiling ~0.1 GB per million alignments, kept per sample only in `--map` runs).

### Databases: size and build time

The tuning world (765 species, `reference.fna` 95 MB), built from the same converted files at
6 threads ([0.7 against 0.6.0a](claude/2026-09-30-v07-vs-v06/README.md)):

| build | wall | CPU | peak RSS | on disk |
|---|---|---|---|---|
| 0.6.0a, raw files | 40.3 s (31.0 s at 1 thread) | 78.5 s | 3.47 GB | 3.82 GB |
| 0.7, raw files (`--no_compress`) | 13.1 s | 27.9 s | 3.65 GB | 3.80 GB |
| 0.7, `database.protal` at zstd level 3 | 4.9 s | 21.1 s | 4.35 GB | 129 MB + full reference |
| 0.7, `database.protal` at level 19 (the default) | 101.2 s | 180.5 s | 4.88 GB | 107 MB + full reference |

0.6.0a was slower on 6 threads than on one (its uniqueness check got slower with threads). The 0.7
build is bound by writing the index when it writes raw files, and by zstd at level 19.

GTDB r226 (143,614 species) with `build_gtdb_database.py`, on the HPC
([r226 build evaluation](claude/2026-10-02-r226-build-evaluation/README.md),
[build profiling](claude/2026-10-03-build-profiling-r226/README.md),
[r226 v9](claude/2026-10-03-r226-v9-evaluation/README.md)):

| run | code | threads | whole workflow | notes |
|---|---|---|---|---|
| v1 | `ff57266` (in 0.7.1) | 16 | 8:54 h | conversion 2:46, training database 0:40, final database (level 19) 1:20 in the background, training-data simulations 1:52 (paired-end) + 1:23 (long reads), protal on 480 samples 0:14, test set 1:01, training 6 min, `--add_model` x4 22 min; peak 88 GiB; `database.protal` 20.45 GB (its index 18.05 GB) |
| v5 | `27423c6` (0.7.4) | 64 | 1:55 h after the holdout | the collector waited ~40 min for the long-read simulations, which ran last and alone |
| v8 | `a3e397d` | 52 | 2:17 h from the end of the conversion | |
| v9 | `26b065c` (0.7.5's defaults) | 52 | 2:02 h from the end of the conversion | with 30% of the species held out, a 30M-pair point and the mixed design; index build 33.5 min, of which the suspect-copy scan 11.7 |

Between v1 and v5 the conversion no longer pipes 86 GB through one zstd process, the gene
neighbours take 9.7x less CPU, the simulations run on the node's scratch, and the collector profiles
all samples of a collection in one protal run. The 0.7.5 simulation queue (long reads largest first
and in chunks, beside the paired-end points) halved a six-core test collection's wall time
(6:05 → 3:11) with the same reads. 0.6.0a had no build-and-train pipeline, and its index build was
never run at r226 size; from the tuning-world builds it would take about three times 0.7.0's.

### The models at GTDB scale

The r226 runs trained a model per read type and scored it on an independent test set of the same
run. The test sets differ between runs (species held out 20% → 30%, abundances sigma 1.3 → 2.0,
deeper points), so the F1 columns are not strictly comparable; the trend is what they show. F1 of
the independent test set at knob 0.5 (false positives per sample where reported):

| run | code | pe | se | pb | ont | source |
|---|---|---|---|---|---|---|
| v1 | `ff57266` (0.7.1 era) | 0.948 (4.0) | 0.947 (3.9) | 0.886 (4.2) | 0.904 (5.9) | [r226 build evaluation](claude/2026-10-02-r226-build-evaluation/README.md) |
| v5 | `27423c6` (0.7.4) | 0.9615 | 0.9570 | 0.9714 | 0.9704 | [v5/v6](claude/2026-10-03-r226-v5-v6-training/README.md) |
| v7 | `a3e397d`, divergence features, knob curve | 0.9654 (2.50) | 0.9606 (2.91) | 0.9764 (0.93) | 0.9737 (1.00) | [v7/v8](claude/2026-10-03-r226-v7-v8-evaluation/README.md) |
| v8 | `a3e397d`, + the sample's depth | 0.9671 (1.64) | 0.9685 (1.46) | 0.9763 (0.64) | 0.9730 (0.78) | same |
| v9 | `26b065c` (0.7.5 defaults, `+priors`), harder test set | **0.9699** | 0.9668 | 0.9727 | 0.9741 | [v9](claude/2026-10-03-r226-v9-evaluation/README.md) |
| v10 | `281a4ba` (after 0.7.5: the default set without priors, as 0.7.6 trains), same samples as v9 | 0.9630 (2.15) | 0.9596 (2.26) | 0.9675 (0.90) | 0.9689 (1.10) | [v10](claude/2026-10-04-r226-v10-evaluation/README.md); [features.md](features.md) |

From v1 to v9 the paired-end false positives per test sample fell from 4.0 to about 1.6 at the same
or higher sensitivity; v1 was precision-limited, with false positives growing with depth (0.2 at
1,000 pairs to 8.6 at 500,000). The v9 ablation attributes +0.007 to +0.009 of paired-end and
single-end F1 to the species priors (mostly their cluster-size singleton rule, which still wants a
check on real data), which is why they are opt-in since `281a4ba` (0.7.6). v10 is the same samples
without them: 0.007 lower on the paired-end test set, all of it absent one-genome species called and
strains of multi-genome species missed. v11, the first build with 0.7.6's training design (in-silico
strains of one-genome species, training depths of 2,000, 50,000 and 200,000 pairs), will show
whether the priors' gain was the simulation's. The PacBio and Nanopore rows between v1 and v5 also
reflect the change from quality-0 to HiFi training reads and the knob curves.

## What these benchmarks do not cover

Accuracy on real samples (real samples were only timed, at r226); a 0.6.0a run at GTDB scale (its
memory and index build there are estimates); 0.7.6 on the benchmark world (by its checks the same
calls as 0.7.5 for short reads) and its memory at r226; strain mixtures and long-read strains of
HiFi reads; 0.7.4 on its own.
Every number above is from one laptop (6 threads) or one HPC node, on simulated worlds whose strains
are 0.4-4% from their references at the markers.

## Sources

The version benchmarks: [0.7 against 0.6.0a](claude/2026-09-30-v07-vs-v06/README.md),
[0.6.0a, 0.7.0 and 0.7.1](claude/2026-10-01-v071-benchmark/README.md),
[0.7.2](claude/2026-10-02-v072-benchmark/README.md), [0.7.3](claude/2026-10-02-v073-benchmark/README.md),
[0.7.5](claude/2026-10-03-v075-benchmark/README.md); for 0.7.6 the
[performance at GTDB scale](claude/2026-10-04-performance-gtdb-scale/README.md) and the
[r226 v10 evaluation](claude/2026-10-04-r226-v10-evaluation/README.md). Everything else is linked where it is used; the
full list is in [`docs/claude/README.md`](claude/README.md). The condensed changelog is the header
comment of `src/main.cpp`.
