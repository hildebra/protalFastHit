# Running protal: details beyond the website

The [website](https://protal.earlham.ac.uk/main.php?site=documentation#usage) explains how to run
protal on one sample or on many with a map file, and what each output file holds. This page adds
what it leaves out. `protal --help` lists the common options, `protal --full_help` all of them, and
`protal --map_help` the map file format.

## Where the outputs go

| Output | `-1 -2 --prefix P -o DIR` | `--map` |
|---|---|---|
| alignments | `DIR/P.sam.zst` | `#SAM_OUTPUT_DIR/<SAM>`, default `#OUTPUT_DIR/alignments/` (without a `SAM` column `<PREFIX>.sam.zst`) |
| profile and its logs | `DIR/P.profile`, `DIR/P.profile.log`, ... | `#PROFILE_OUTPUT_DIR/<PROFILE>`, default `#OUTPUT_DIR/profiles/` |
| strain MSAs and tables | `DIR/strains/` | `#STRAIN_OUTPUT_DIR`, default `#OUTPUT_DIR/strains/` |
| coverage, SNP counts, statistics | `DIR/misc/` | `#MISC_OUTPUT_DIR`, default `#OUTPUT_DIR/misc/` |

- `-o` replaces a map's `#OUTPUT_DIR`, and `--profile_dir` moves the profiles of either mode.
- Without `--prefix`, the prefix is the longest common prefix of the two read file names, which
  needs both files in the same folder; for single-end reads it is the file name without its
  FASTQ/FASTA and compression extensions (`S1.fq.gz` gives `S1`).
- In a map, `#SAMPLEID` (the first column) is the sample name in MSA rows, logs and statistics.
  Every row needs a value in every column the header declares, columns are separated by tabs,
  and no two samples may share a SAM or profile file; protal checks all of this before it starts.
- The SAM is an intermediate file, and is zstd-compressed by default: names protal picks end in
  `.sam.zst` (`--sam_format gz` gives `.sam.gz`, `--sam_format sam` plain `.sam`). In a map, the
  `SAM` column's name chooses the format: `.sam.zst` zstd, `.sam.gz` gzip, any other name plain SAM.
  The alignment threads compress as they write, so no tool is needed. zstd is the faster choice:
  for 1M pairs of a marker-rich sample at 8 threads, aligning and writing `.sam.zst` took 22–38%
  less time than `.sam.gz`, for a 12% smaller file, and profiling reads it faster;
  `zstdcat S1.sam.zst` or `zstd -dc` decompresses it (samtools does not read zstd). A `.sam.gz` is
  BGZF (gzip blocks of 64 KB, as `bgzip` writes), which `zcat`, `gzip -d` and samtools read.
  `protal_map_utils generate` (and `merge --use-sampleid`) and `simulate_metagenomes
  --protal_metafile` write `.sam.zst` names; `protal_map_utils --gzip` writes `.sam.gz`, `--nogzip`
  `.sam`. A rerun also takes a plain `P.sam` from an earlier run as the SAM of `P.sam.zst`.
- The SAM header (`@SQ`) lists the genes that the alignments name, not every gene of the database
  (the full r226 database has millions); `--full_sam_header` lists every gene, as protal did before.
- `<sam>.err` lists the reads whose alignment does not fit the database (a gene it lacks, a
  position past a gene's end, bases that differ from the gene). They are left out of the
  profile, and protal warns with their number.
- `misc/` also receives `P_seedsizes_histogram.tsv`, `P_anchorsizes_histogram.tsv` and
  `P_runtime.tsv`, diagnostics of the seeding and alignment stages. `P_runtime.tsv` has one row
  per stage (reading, k-mers, seeding and its steps, alignment, output): the seconds spent in it
  summed over threads, the number of threads, and the seconds per thread that `--verbose` prints.

Strain MSAs are written for the species that pass the model in at least two samples, each with a
row for every sample in which the species passes. A run of one sample therefore writes no MSAs.
`--msa_species s__Genus_species,...` writes MSAs for the named species only. A species passes
with a model probability of at least `--knob` in the profiles and of at least `--msa_knob`
(default: `--knob`) for the MSAs, so by default both hold the same samples. An MSA takes the genes
with reads in its samples. For a species with relatives in the database, whose genes have long
unique k-mers (a 15-mer core shared with a relative) in 90% of cases or more, the genes without
any are left out: a relative may share them unchanged, and its reads would then show as a second
strain.

`misc/unreported_species.tsv` lists the species that a sample's profile leaves out, their
probability below `--knob`, although their own reads are strong evidence that they are present:
1x or more depth from reads within `--depth_identity_margin` of their best ones, reads on 90% of
their genes, best reads 98% identical to the reference or more, and at most half of their
aligned bases from reads of lower identity. protal warns when it lists any. In the model's test
sets on the toy database this held for 9% of the present species below the knob and for none of
14,390 absent ones. Columns: sample, species, taxid, score, own_depth, hit_gene_fraction,
top_identity, low_identity_share, and passes_msa_knob (whether the sample enters the species'
MSA).

## Single-end reads

`-1 reads.fq` without `-2` aligns single-end reads (100-300 bp, e.g. Illumina single-end runs).
In a map, a single-end sample has `-` in the `SECOND` column, and a map without a `SECOND` column
has only single-end samples; one run may mix single-end and paired-end samples, and their strain
MSAs are joint.

- Each read is aligned on its own, with no mate to recover anchors from or to pair alignments
  with. Its candidates are written as unpaired records (FLAG 0 or 16, and 256 for the ones after
  the first), each alignment once; the first carries the read's MAPQ, which compares its best and
  second-best candidate as for a mate aligned alone.
- Single-end samples are profiled with a model trained on single-end reads: the database's
  `model_se.xml`, or `--model_se`. The model of paired-end reads does not fit them (its features,
  such as hits and MAPQ, are counted per mate of a pair). `--model` alone replaces both models.
  protal checks before aligning that each kind of sample in the run has its model.
- `--profile_only` tells a single-end SAM from its records (no `0x1` flag) and profiles it with
  the single-end model.

## Reruns and profiling existing alignments

protal skips the alignment of a sample whose SAM file already exists and profiles that SAM;
`--force` aligns again. SAMs are written under a temporary `.partial` name and renamed when
complete, so an interrupted run never leaves a truncated SAM that a rerun would reuse. A truncated
`.sam.gz` or `.sam.zst`, or a SAM aligned against another database (its `@SQ` genes missing or of
another length), stops with an error.

`--profile_only a.sam,b.sam.gz,c.sam.zst` profiles existing SAM files without loading the index. The
prefixes come from `--prefix` (one per file) or from the SAM names, and the outputs go to `-o`,
or next to each SAM without it. This is the quick way to try another `--knob`, model or
`--depth_identity_margin`.

## Exit status

| Exit code | Meaning |
|---|---|
| 0 | every sample and output was written |
| 1 | the run finished, but at least one sample or output failed; the errors are listed at the end (`protal finished with N error(s)`) |
| other | protal stopped before aligning: invalid options, missing or inconsistent input files, a database or model that fails its checks. The message says what is wrong |

Inputs, the database files and the model are all checked before the index is loaded, and all
missing files are reported at once, so a bad run fails in seconds rather than after hours of
alignment. Workflow managers can rely on a non-zero status.

## Options the website does not list

| Option | Default | |
|---|---|---|
| `-t, --threads` | 1 | threads for alignment (which also compresses the SAM), database loading and profiling. Set it: the default is one thread. While aligning, each read file is also decompressed by a thread of its own (two for paired reads); BGZF files (`bgzip`, `simulate_metagenomes`) decompress about 2x faster than other gzip files (libdeflate against zlib-ng) |
| `--knob` | 0.5 | detection threshold, 0 to 1 (checked). Choose it on data like yours; see [model-training.md](model-training.md) |
| `--depth_identity_margin` | 0.04 | a read counts towards a species' abundance, and towards its strain MSA rows, only if its identity is at most this far below that of the species' best reads (98th percentile). Reads of relatives the database lacks still count for detection, not for depth or strains. 1 lets every read count |
| `--model` | `model.xml` of the database (`model_se.xml` for single-end samples) | a PMML file, or the name of another model in the database folder (`<name>.xml`); for all samples unless `--model_se` is given. protal checks the model before aligning, see [model-training.md](model-training.md) |
| `--model_se` | `--model`, else `model_se.xml` of the database | the model of single-end samples, given as `--model` |
| `--sam_format` | `zst` | the format of the SAM files protal names: `zst` (`.sam.zst`), `gz` (`.sam.gz`) or `sam`; a map's `SAM` names keep their own ending |
| `--no_strains` | off | no MSAs or SNP tables. Variants are still called, since the model uses them, so profiles are the same with and without it |
| `--msa_knob` | `--knob` | model probability a sample's species needs for the sample to enter the species' strain MSA. By default the MSA holds the samples whose profile reports the species; lower it to add samples the profile leaves out |
| `--msa_min_hcov` | 1000 | minimum non-N, non-gap bases for a sample's sequence to stay in an MSA; passed to qcmsa as `--reapply-hcov` |
| `--msa_min_depth` | 1 | reads a position needs to be written in an MSA, else `-`. Where its reads all show one allele, that many suffice; a second allele (an IUPAC code) needs `--snp_min_cov` reads, and a position whose reads disagree otherwise is `N`. 2 is the behaviour before the 2026-09-29 strain audit |
| `--snp_max_alleles` | 3 | alleles encoded as an IUPAC ambiguity code in the MSA: 1 = only the top allele, 2 = two-allele mixtures (R, Y, ...), 3 = also three-allele mixtures (B, H, ...) |
| `--qcmsa_script` | | the qcmsa executable; see [installation.md](installation.md#installing-a-source-build) for how protal finds it otherwise |
| `--preload_genomes_off` | off | read reference genes on demand instead of loading `reference.fna`: less memory, slower. Needs the database as separate files, see [database-files.md](database-files.md) |
| `--verbose` | off | more progress output |

The SNP filters (`--snp_min_cov`, `--snp_min_phred_sum`, `--snp_min_mean_qual`, `--snp_min_af`,
`--snp_no_strand`) are described on the website. Base qualities are read with the standard
Phred+33 offset. `--snp_min_af` defaults to 0.15 (0.2 for ONT reads): an allele needs
that share of the reads with a base at a position to be called or to enter an IUPAC code. The
profiles do not depend on it.

### Alignment options

Shown by `--full_help`. The defaults are what the model and the profiler's MAPQ cutoffs were
calibrated with; change them for experiments, not for production profiles.

| Option | Default | |
|---|---|---|
| `-c, --align_top` | 3 | candidate anchors aligned per read, best first. MAPQ depends on it: with `-c 1` far more reads look unique |
| `-m, --max_out` | 1 | alignments written per read. Profiles are the same with more, SAMs larger |
| `-u, --max_key_ubiquity` | 256 | skip seeds whose best-matching flex-key occurs more often than this |
| `-s, --max_seed_size` | 128 | seeding stops at this many seeds, if `-w` lookups have succeeded |
| `-w, --min_successful_lookups` | 4 | successful core k-mer lookups needed before `-s` stops seeding |
| `-a, --max_score_ani` | 0.9 | give up an alignment once it diverges below about this identity |
| `-x, --x_drop` | 1000 | X-drop of the alignment (WFA2), added to its adaptive pruning; 0 turns it off. The default prunes nothing in practice (outputs identical to `-x 0`); `-x 50` loses a few alignments and changes MAPQs |

### Developer options

Also shown by `--full_help`: `--build` and its options ([building-a-database.md](building-a-database.md)),
the database conversions `--compress_db`, `--unpack_db`, `--decompress_db`
([database-files.md](database-files.md)), `--profile_truth` for the training dump
([model-training.md](model-training.md)), `--benchmark_alignment` (checks alignments against the
`taxid_geneid` encoded in simulated read names), `--mapq_debug_output`, `--full_sam_header` (every
gene in the SAM header, see above), and `--whole_read_alignment`.

Short reads are aligned from their anchor's exact matches: WFA aligns the read left and right of
them (and between them), each part anchored at a match, instead of the whole read into the gene
window. Alignments come out as good as before by the aligner's scoring (in tests never worse; 0.5%
of SAM records differ, mostly in where a gap sits among equally good places, and about 0.15% of
pairs get a MAPQ a few points apart) and faster, since anchors on relatives' genes are given up on
sooner. `--whole_read_alignment` aligns every read as a whole, as protal did
before, e.g. to reproduce earlier results; long reads are always aligned as a whole.

## Environment variables

| Variable | |
|---|---|
| `PROTAL_DB_PATH` | the database, when `--db` is not given |
| `PROTAL_QCMSA_SCRIPT` | the qcmsa executable, when it is not next to protal or on `$PATH` |
| `PROTAL_NO_AVX2` | set to anything: the `protal` launcher runs the baseline build even on AVX2 CPUs |

## Memory

protal keeps the index and the reference genes in memory: about 59 GB for the full r226 database
and 12 GB for the reduced one ([downloads](https://protal.earlham.ac.uk/main.php?site=downloads)).
It prints the machine's total memory at start. The index is read at random, one lookup per k-mer,
so protal asks Linux for transparent huge pages for it; the usual setting (`madvise` in
`/sys/kernel/mm/transparent_hugepage/enabled`) grants them, and seeding is about a third faster.
With THP set to `never` protal uses normal pages. The profiling stage streams each SAM and, once a
sample's outputs are written, keeps only what the strain MSAs need: the variants and read ranges
of the species that pass the model (nothing with `--no_strains`).
