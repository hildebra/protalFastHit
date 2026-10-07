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
| coverage, SNP counts, timers | `DIR/misc/` | `#MISC_OUTPUT_DIR`, default `#OUTPUT_DIR/misc/` |

- `-o` replaces a map's `#OUTPUT_DIR`, and `--profile_dir` moves the profiles of either mode.
  Without `-o`, a `-1 -2` run writes into the current folder (runs without `-o` aborted up to 0.7.8).
- Without `--prefix`, the prefix is the longest common prefix of the two read file names, which
  needs both files in the same folder; for single-end reads it is the file name without its
  FASTQ/FASTA and compression extensions (`S1.fq.gz` gives `S1`).
- In a map, `#SAMPLEID` (the first column) is the sample name in MSA rows, logs, statistics and
  file names (`misc/<sample>_runtime.tsv`), so it must be safe as a file name on Linux and macOS: no
  `/` or `:`, no control characters, not starting with `-`, not `.` or `..`, UTF-8, at most 200
  bytes (and no whitespace with strain MSAs). Every row needs a value in every column the header
  declares, columns are separated by tabs, and no two samples may share a SAM or profile file;
  protal checks all of this before it starts, and lists every unsafe sample ID with its line.
- The SAM is an intermediate file, and is zstd-compressed by default: names protal picks end in
  `.sam.zst` (`--sam_format gz` gives `.sam.gz`, `--sam_format sam` plain `.sam`). In a map, the
  `SAM` column's name chooses the format: `.sam.zst` zstd, `.sam.gz` gzip, any other name plain SAM.
  The alignment threads compress as they write, so no tool is needed. zstd is the faster choice:
  for 1M pairs of a marker-rich sample at 8 threads, aligning and writing `.sam.zst` took 22–38%
  less time than `.sam.gz`, for a 12% smaller file, and profiling reads it faster;
  `zstdcat S1.sam.zst` or `zstd -dc` decompresses it (samtools does not read zstd). Its records
  go straight into the file, behind room left for the header (16 KB to 1 MB, from the database's
  gene count and the size of the read files), which a skippable zstd frame pads; the room the header
  does not need is never written, a hole where the file system has them, so `ls -l` can show up to
  that much more than `du`. A `.sam.gz` or `.sam` collects its records in a temporary file and copies
  them behind the header at the end, as a `.sam.zst` does when its header needs more room. A `.sam.gz` is
  BGZF (gzip blocks of 64 KB, as `bgzip` writes), which `zcat`, `gzip -d` and samtools read.
  `protal_map_utils generate` (and `merge --use-sampleid`) and `simulate_metagenomes
  --protal_metafile` write `.sam.zst` names; `protal_map_utils --gzip` writes `.sam.gz`, `--nogzip`
  `.sam`. A rerun also takes a plain `P.sam` from an earlier run as the SAM of `P.sam.zst`.
  `protal_map_utils generate` names a sample after its read files without the mate number, an `R`
  before it and Illumina's chunk number (`x_R1.fq` is `x`, `x_S1_L001_R1_001.fastq.gz` is `x_S1_L001`).
  `protal_map_utils merge` keeps the SAMs the merged maps' runs wrote, by their absolute paths, so a
  run of the merged map profiles them instead of aligning again
  ([strains.md](strains.md#strain-msas-over-several-runs)).
- The SAM header (`@SQ`) lists the genes that the alignments name, not every gene of the database
  (the full r226 database has millions); `--full_sam_header` lists every gene, as protal did before.
- `<sam>.err` lists the reads whose alignment does not fit the database (a gene it lacks, a
  position past a gene's end, bases that differ from the gene). They are left out of the
  profile, and protal warns with their number.
- `misc/` also receives `P_seedsizes_histogram.tsv`, `P_anchorsizes_histogram.tsv` and
  `P_runtime.tsv`, diagnostics of the seeding and alignment stages. `P_runtime.tsv` has one row
  per stage (reading, k-mers, seeding and its steps, alignment and its k-mer screen, output): the seconds spent in it
  summed over threads, the number of threads, and the seconds per thread that `--verbose` prints.
  The stages that run for every read are timed on every 61st call only (reading the clock costs
  10-15% of the alignment time otherwise), so their seconds are estimates: the mean timed interval
  times the number of calls. Since 0.7.6 it also has the profiling steps' wall-clock times (reading the
  SAM, record evidence and sample context, read EM, congener distances, SNPs, scoring, writing), which
  protal prints in one line per sample too; a rerun that profiles the sample again replaces them.
- Every run prints per sample how its reads went: the reads, those with an anchor, the candidate
  alignments tried, those the k-mer screen refused, those aligned, and the records written. A second
  line ("seeding:") counts the k-mer lookups, the index blocks they scanned and their sizes, the
  seeds, the seeds that share their taxon and gene with another seed of their read (the only ones
  anchors are made of), the lookups dropped as too ubiquitous (more tied entries than
  `--max_key_ubiquity`) and the anchors. At the end it times loading, freeing memory and every
  stage, so that a run's wall time adds up.
- `misc/<taxon>.statistics.tsv` (a taxon's coverage, reads, ANI and MAPQ in each sample) is written
  only with `--taxon_statistics` since 0.7.6: on a GTDB-sized database a sample has reads on
  thousands of taxa, and that many small files took 15 s of a run on a network file system. The
  profile files hold the same numbers per sample.

Strain MSAs (which samples, genes and reads enter them, long-read strains, qcmsa, trees) are
described in [strains.md](strains.md).

`misc/unreported_species.tsv` lists the species that a sample's profile leaves out, their
probability below the sample's knob, although their own reads are strong evidence that they are present:
1x or more depth from reads within `--depth_identity_margin` of their best ones, reads on 90% of
their genes, best reads 98% identical to the reference or more, and at most half of their
aligned bases from reads of lower identity. protal warns when it lists any. In the model's test
sets on the toy database this held for 9% of the present species below the knob and for none of
14,390 absent ones. Columns: sample, species, taxid, score, own_depth, hit_gene_fraction,
top_identity, low_identity_share, and passes_msa_knob (whether the sample enters the species'
MSA).

Species called in several samples of a run on a read or two, always beside the same more abundant
relative, point to a reference artefact rather than to the species: a gene copy the database's
species carries from another lineage, which every present organism with that gene puts a perfect
read on. `scripts/recurrent_calls.py OUTPUT_DIR` reads a run's `*.profile.log` files and lists such
species with their companion and the rank the two share (`recurrent_calls.tsv`; `--max-hits`,
`--min-samples`, `--ratio` set what counts as thin, recurring and more abundant). `--build` flags
such copies in advance where another genus's copy is in the database (`suspect_copies.tsv`,
`--keep_suspect_copies`); the script finds what it could not ([report](claude/2026-10-03-false-positive-anatomy/README.md)).

The samples of one study share species, and a species present in most of them is likelier present in
the one at hand than the model's training prior says. `scripts/prevalence_calls.py OUTPUT_DIR` reads a
run's `*.profile.log` files and moves each taxon's probability in a sample by the odds ratio of its
prevalence in the other samples to the base rate (a Bayesian update; the sample itself never counts),
then calls at `--knob` (`prevalence_calls.tsv`; `--profiles DIR` writes adjusted `.profile` files). By
default (`--direction up`) only the boosts of prevalent species are applied; `--direction both` adds the
cuts of species seen nowhere else, the full update, which is right only when the samples do share
species: on 42 simulated samples drawn independently it removed a quarter of the calls, nearly all
true, while the default changed nothing. Use `--check` to see how far the probabilities moved and how
many calls flipped, and validate on a study with a truth first ([report](claude/2026-10-03-false-positive-fixes/README.md)).

## Read files

Read files are FASTQ or FASTA, plain, gzip, BGZF or zstd, told apart by their first bytes (not
their names); bases are read as uppercase. Other compressions (bzip2, xz) are read through a pipe,
which protal reads as a file, as it does process substitution: `-1 <(xz -dc a_R1.fq.xz) -2 <(xz -dc
a_R2.fq.xz) --prefix a` (give `--prefix`: the pipe's name, `/dev/fd/63`, names no sample). A gzip
file may hold several members, one after the other (`cat a.fq.gz b.fq.gz`, also after BGZF ones),
and zero bytes of padding; a zstd file several frames (`cat a.fq.zst b.fq.zst`, `pzstd`'s and the
seekable format's skippable frames, any `zstd --long` window). Each file is decompressed in a
thread of its own, so the alignment threads only copy reads that are already decompressed.

A read file that is missing, a directory or not readable stops protal before the index is loaded.
A sample fails (exit 1, no SAM file) when its reads are not FASTQ or FASTA, when a record is
incomplete, has no `+` line or has more or fewer qualities than bases, when a gzip or zstd file is
truncated or corrupt (whatever follows a gzip member must be another member: a damaged member
header is an error, not the end of the file; a zstd file must not end inside a frame and holds
zstd frames only), when paired files hold different numbers of reads,
or when a single-end read is longer than 1,000 bp (below). A sample without reads, and pairs
whose mates' names differ (other than in a last `1` and `2`, as in `r/1` and `r/2`), get a
warning: the files of a pair may not be in the same order.

## The kind of reads

`--read_type` (or a map's `READ_TYPE` column) names the kind of a sample's reads: `pe`, `se`, `pb`
or `ont`. Without it, a sample with a second read file has paired-end reads, and for a single read
file protal looks at its first 200 reads (at most 5 Mb):

- At most a quarter of them longer than 1,000 bp: single-end reads. A longer read later in the
  file stops the sample (exit 1), as any long read does in a sample given as `se`, whose first
  100 reads protal checks before aligning.
- Else long reads: ONT reads if most of them are named as MinKNOW and dorado name reads (a UUID;
  MinKNOW adds `runid=`), PacBio reads if named by movie and ZMW (`m64011_190830_220126/123/ccs`).
  Else their median read quality decides, each read's quality being that of its bases' mean error
  probability: PacBio from Q25 on (HiFi reads are Q30 and better), ONT below (about Q12-25). Long
  reads without qualities named neither way are taken for PacBio reads: FASTA, or Q0 at every
  base, which means unknown (pbsim3 writes it for its PacBio reads).

protal says what it found, e.g. `Note: Sample s: ONT reads (the first 200 reads up to 48.1 kb
long, median 6.3 kb, median read quality Q18.4 (PacBio from Q25, ONT below) in s.fq.gz), aligned
and profiled as such`, and the sample is aligned and profiled with that kind's settings and model.
Reads in a pipe are not looked at (they would be lost to the alignment): give `--read_type` for
long reads read through a pipe.

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
`--force` aligns again. Its reads are profiled as the kind the SAM's header names (protal writes
it there), not as the read files suggest; a kind given with `--read_type` or `READ_TYPE` wins.
protal warns when the two differ. SAMs are written under a temporary `.partial` name and renamed
when complete, so an interrupted run never leaves a truncated SAM that a rerun would reuse. A truncated
`.sam.gz` or `.sam.zst`, or a SAM aligned against another database (its `@SQ` genes missing or of
another length), stops with an error.

`--profile_only a.sam,b.sam.gz,c.sam.zst` profiles existing SAM files without loading the index, and
builds strain MSAs over all of them, also when they come from several runs
([strains.md](strains.md#strain-msas-over-several-runs)). An item with a wildcard, such as
`'runs/*/alignments/*.sam.zst'`, stands for the SAM files it matches (quoted or not). The prefixes
come from `--prefix` (one per file) or from the SAM names (from their folders where names repeat),
and the outputs go to `-o`, or without it the profiles next to each SAM and `strains/` and `misc/`
into the current folder. This is also the quick way to try another `--knob`, model or
`--depth_identity_margin`. `--profile_only` stops before it starts where an earlier run's results
would be overwritten (a profile it would write, or `species.tsv` in its strain folder): give another
`-o`, or `--force` to write them again. Arguments that follow no option stop protal, except SAM files after
`--profile_only`: an unquoted `-1 *_1.fq` once aligned only its first file.

Since 0.7.6 a profile does not depend on the order of the SAM's records, which multi-threaded
alignment varies from run to run: repeated runs on the same reads and database give the same
profile, to the last digit.

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
| `-t, --threads` | 1 | threads for alignment (which also compresses the SAM), database loading and profiling. Set it: the default is one thread. While aligning, each read file is also decompressed by a thread of its own (two for paired reads): gzip, BGZF or not, with ISA-L at about 1 GB/s of FASTQ, zstd with libzstd; at many threads a paired-end run can still wait for its input. Samples are profiled in parallel, the largest SAM first, each on threads in proportion to its SAM's share of all the samples' bytes (at least one; a single sample on all), with the same results as on one thread. A sample profiled on several threads also has its SAM read by a thread of its own; from 8 threads on, the zstd frames of a `.sam.zst` that protal wrote are decompressed ahead by `-t`/6 more threads (2 to 8) |
| `--knob` | 0.5 | detection threshold, 0 to 1. A model with a knob curve over the sample's depth uses the curve unless `--knob` is given; models trained with the depth as a feature (the default) have none ([databases.md](databases.md#how-protal-calls-species)). Choose it on data like yours |
| `--fdr` | off | report each sample's species while their expected share of false calls stays at or below this, for a model with calibrated calls; at r226 it called slightly below the knob curve, so it is off ([databases.md](databases.md#calls-at-a-target-share-of-false-calls)) |
| `--singleton_congener` | 0 | veto a single-fragment species beside a congener of at least this many fragments when its read looks like the congener's; 0: no rule ([databases.md](databases.md#the-singleton-rule)) |
| `--keep_suspect_copies` | off | count reads on the database's suspect gene copies (near-identical to another genus's copy: contamination or a transferred gene) as evidence; by default they are left out, and the log says how many per sample ([databases.md](databases.md#2-build-the-index)) |
| `--depth_identity_margin` | 0.08 | a read counts towards a species' abundance only if its identity is at most this far below that of the species' best reads (98th percentile). Reads of relatives the database lacks still count for detection, not for depth. The margin is the same on every gene unless `--gene_conservation` scales it. The default counts the reads of strains up to about 5% from the reference; 0.04, the earlier default, dropped up to 40% of the reads of strains 3–5% away and undercounted those strains ([report](claude/2026-09-30-depth-margin-stress/README.md), [scaled per gene](claude/2026-10-01-gene-scaled-margin/README.md)). 1 lets every read count. The depth counts each base a fragment covers on a gene once: where a pair's mates overlap there, the second mate adds only what the first did not cover, as the strain MSA counts them ([report](claude/2026-10-02-fragment-depth/README.md)) |
| `--gene_conservation` | `none` | scale `--depth_identity_margin` per gene by its conservation factor: `db` for the database's, or a file of `geneid<TAB>factor`; a gene's margin is then 0.03 plus the rest times its factor. `none` did best over three simulated worlds ([report](claude/2026-10-01-gene-scaled-margin/README.md)). The model's conservation features use the factors either way |
| `--model` | `model.xml` of the database (`model_se.xml` for single-end samples) | a PMML file, or the name of another model in the database folder (`<name>.xml`); for all samples unless `--model_se` is given. protal checks the model before aligning, see [databases.md](databases.md#the-presence-model) |
| `--model_se` | `--model`, else `model_se.xml` of the database | the model of single-end samples, given as `--model` |
| `--sam_format` | `zst` | the format of the SAM files protal names: `zst` (`.sam.zst`), `gz` (`.sam.gz`) or `sam`; a map's `SAM` names keep their own ending |
| `--no_strains` | off | no MSAs or SNP tables (profiles are the same). The other strain options are in [strains.md](strains.md#options) |
| `--keep_foreign_genes` | off | keep foreign genes (below) in their taxon's depth and strain MSAs; by default they are left out of both. `.profile.genes.log` lists them either way |
| `--preload_genomes_off` | off | read reference genes on demand instead of loading `reference.fna`: less memory, slower. Needs the database as separate files, see [databases.md](databases.md#the-files-of-a-database) |
| `--verbose` | off | more progress output |

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
| `-x, --x_drop` | 1000 | X-drop of the alignment of short reads (WFA2), added to its adaptive pruning; 0 turns it off. The default changes no short-read alignment in tests (outputs identical to `-x 0`); `-x 50` loses a few alignments and changes MAPQs. Long reads (`pb`, `ont`) are aligned without X-drop: over their gene-long windows even 1000 lost the own species' alignment of genes an ONT read ends in |
| `--no_mate_guidance` | off | paired-end reads: do not let a mate that is sure of its alignment (MAPQ 20 or more) guide the other when they did not align together (below) |
| `--no_alignment_screen` | off | align every candidate with WFA2. By default a candidate whose read and gene window share too few k-mers for any alignment within the score budget is refused before WFA2 (exact: it changes no alignment; at r226 it took the paired-end alignment from 62 to 39 s). For measuring |
| `--long_read_budget E` | 0 | long reads: once a candidate of a gene has aligned, its other candidates get the best's budget plus E edits. Saves 10-25% of PacBio aligning but changes MAPQs and the alternatives of such genes, so the models would need retraining; off until benchmarked |
| `--no_gene_neighbours` | off | do not use the database's gene neighbours (below): no mates looked for past their gene's end, no pairs over two neighbouring genes, no genes looked for next to a long read's genes, and the profile's `adjacent_*` features 0. A database without `gene_neighbours.tsv` works as with it |

A read is one organism, so its parts are given one taxon. When the mates of a pair aligned to two
genes, or a long read to several, protal compares the taxa of their candidates on the parts both
have a candidate on, and gives every part its alignment of the winning taxon (a long read's genes
from at least half of them); a part without an alignment of that taxon (a gene the database lacks
for that species), and a long read's gene that another taxon's gene fits clearly better (MAPQ 4 or
more between them: a chimeric read, or a homolog of a gene that lies elsewhere on the read), keeps
its own best alignment but is written with MAPQ 0, so that it counts for no other taxon (long reads:
tagged `ZR:i:2`; `ZR:i:1` marks a gene whose alignment or MAPQ the read's taxon changed). When only one mate of a pair aligned, or the other has no candidate of the first one's
taxon, a mate that is sure of its alignment guides the other: the other mate's own anchor of that
taxon is aligned (only a read's best anchors are), or the other mate is looked for on the first
one's gene where the fragment can reach (up to 1,000 bases on its side), placed by its 12-mers,
and aligned in part if it runs past the gene's end. The log reports how many mates each way found.

Marker genes often lie next to each other (the ribosomal protein operons, rpoB and rpoC, ...), so a
fragment or a long read can span two of them. A database built with gene neighbours
(`gene_neighbours.tsv`, [databases.md](databases.md#gene-neighbours)) knows,
for each clade, how often each gene's end faces which other gene's end in its species' genomes and
how far apart they are. For a taxon, a pairing's share in its family leans on the clades above it, the
more the fewer species the family has
([databases.md](databases.md#gene-neighbours)): ends that face each other in
20% or more are expected, in 5% or less unlikely. With it:
- mates on two genes of one taxon that each run towards the end facing the other gene (expected
  neighbours), as one
  fragment of at most 1,000 bases would, are paired like mates on one gene: written as a proper
  pair (flag 2) on two references (RNEXT the other gene, TLEN 0);
- a guiding mate whose fragment reaches past its gene's end also looks for the other mate on the
  gene the clade has there, over the stretch the fragment reaches after the gap;
- a long read that goes on past the end of a gene of its taxon is searched, where the clade puts
  that gene's neighbour, for the neighbour if no part of the read has it: a candidate of it beyond
  the longest is aligned, and a gene too divergent to be seeded is placed by its 12-mers (one that
  was aligned already and did not pass is not tried again);
- the profile gets `adjacent_expected_share`, `adjacent_unlikely_share` and `adjacent_support`: of
  the genes next to each other on a taxon's reads, the shares that are expected and unlikely
  neighbours in the taxon's clade, and the mean frequency of their pairings there
  ([features.md](features.md#gene-neighbours-adjacency-071-and-073));
- a gene whose reads' neighbouring genes are mostly unlikely neighbours in its taxon's clade (4 or
  more such links judged, more than half of them unlikely) is *foreign*: its reads come from another
  genome, one that carries a copy of it among other genes (a transferred gene, a relative's homolog in
  another gene order, a contaminant of the reference). It is left out of the taxon's depth (it then
  counts as a gene the taxon lacks) and of its strain MSA rows, unless `--keep_foreign_genes`;
  `.profile.genes.log` lists every gene's `LinksJudged`, `LinksUnlikely` and `Foreign`, and the log
  says how many a sample has. A gene with one end in context and one not is not foreign, and with
  the species' own lines in the table ([databases.md](databases.md#gene-neighbours))
  the genes next to each other on its own genome's reads are expected, so its own genes are not taken
  for foreign ones;
- a long-read sample's strains are phased with reads cut between unlikely neighbours
  ([strains.md](strains.md#strains-in-long-read-samples)).

The log says what the database has (`Gene neighbours: ... rules of ... clades`) and how many
fragments and genes were found so. Gene order is much the same across bacteria, so this helps to
align and pair reads; it does not tell a species from its congeners.

Each read's best record carries `ZA:Z:<taxid>:<edits more>,...`: its other candidates in other
taxa (at most 5 edits more than the best, or `*`), from which the profiler counts the reads that a
congener fits as well ([features.md](features.md#the-reads-other-candidates-071)). A read's first record carries
`ZF:Z:<taxid>,...` when the read seeded on taxa strongly enough to be aligned against them (the
align-top anchors) but did not align to them; the reads that seeded on taxa and aligned nowhere are
counted per taxon in one header line (`@CO protal failed candidates of unaligned reads: <taxid>:<reads>,...`),
and nothing else is written for unaligned reads. `--write_unmapped_reads` writes a minimal unmapped
record (flag 4, no sequence) with the tag for each of them instead, as protal did up to 0.7.6 (on
GTDB r226 these were 95% of a paired-end sample's records, which the profiler read only to count
them); so does `--full_sam_header`, whose header is written before the reads. A map's
`UNMAPPED_READS` column chooses per sample: `write` (the records) or `count` (the header line; `-`
leaves it to `--write_unmapped_reads`); the database build uses it to keep the non-hits of its
samples ([databases.md](databases.md#the-reads-behind-the-errors)). The profiler counts
per taxon the reads that failed on it (`failed_candidate_rate`): a relative the database lacks seeds on its nearest species and fails there, a present species' reads align
([report](claude/2026-10-03-false-positive-fixes/README.md)). Profiling a SAM of an older protal, without these, gives the feature 0.

### Developer options

Also shown by `--full_help`: `--build` and its options ([databases.md](databases.md#building-a-database)),
the database conversions `--compress_db`, `--unpack_db`, `--decompress_db`
([databases.md](databases.md#the-files-of-a-database)), `--profile_truth` for the training dump
([databases.md](databases.md#training-data)), `--benchmark_alignment` (checks alignments against the
`taxid_geneid` encoded in simulated read names), `--mapq_debug_output`, `--full_sam_header` (every
gene in the SAM header, see above), `--write_unmapped_reads` (an unmapped record for each read that
aligned nowhere, see above), `--whole_read_alignment` and `--profile_after_alignment` (both below),
`--sequential_load` (load the database's parts one after another, as before 0.7.6, instead of the
index beside the genome preload and the small tables each on a thread; the run is the same either way).

Short reads are aligned from their anchor's exact matches: WFA aligns the read left and right of
them (and between them), each part anchored at a match, instead of the whole read into the gene
window. Alignments come out as good as before by the aligner's scoring (in tests never worse; 0.5%
of SAM records differ, mostly in where a gap sits among equally good places, and about 0.15% of
pairs get a MAPQ a few points apart) and faster, since anchors on relatives' genes are given up on
sooner.

Long reads (PacBio, ONT) are aligned from their anchor's exact matches too. Their indels shift the
diagonal along a gene, so the alignment goes through every link of the chain, aligning the bases
between two links end to end, and since few of the index's k-mers survive their errors (a long
read's seeds cover a third of its gene or less), it first looks for more exact matches (12 bases,
unique near the diagonal the neighbouring link leads to) before, between and after the links. On
the long-read benchmark that left the profiles as they were (the same F1, Bray-Curtis between the
two 0.0005 or less on 90 Mb samples), took 15-21% less CPU, and gave alignments of the same
penalty, placing indels differently where several places are as good (in 38-46% of the records of
simulated ONT and older PacBio reads, 8% of HiFi reads').
`--whole_read_alignment` aligns every read as a whole, as protal did before, e.g. to reproduce
earlier results.

`--profile_ahead` (several samples): a sample is profiled as soon as its SAM file is complete,
while the next sample's reads are aligned, by a worker on a quarter of the threads; the profiling
stage after the alignment takes the rest (the last sample, and any the worker is behind with) on
all threads. A sample's profile depends on its SAM and the database only, so the profiles are the
same either way. It can only pay where the profiling stage leaves cores idle; on a 6-core machine
with the alignment on every core it gained nothing, so it is off by default.

## Environment variables

| Variable | |
|---|---|
| `PROTAL_DB_PATH` | the database, when `--db` is not given |
| `PROTAL_QCMSA_SCRIPT` | the qcmsa executable, when it is not next to protal or on `$PATH` |

## Memory

protal keeps the index and the reference genes in memory: 36.5 GB peak for a run on the full r226
database, measured on real paired-end and HiFi samples (2026-10-05, `2d809cf`, 32 threads), and about
34 GB expected since the index load frees its threads' buffers as it ends and the index's values take
41 bits (26.6 GB;
its key map, 3.2 GB; the genes at two bits per base, the tables and the run's buffers, the rest;
[databases.md](databases.md#the-database-in-memory)) and correspondingly less for the
reduced one ([downloads](https://protal.earlham.ac.uk/main.php?site=downloads); the figures there,
59 and 12 GB, predate the 2-bit genes and the packed index).
It prints the machine's total memory at start and, after loading the index, the memory it takes. The index is read at random, one lookup per k-mer,
so protal asks Linux for transparent huge pages for it; the usual setting (`madvise` in
`/sys/kernel/mm/transparent_hugepage/enabled`) grants them, and seeding is about a third faster.
With THP set to `never` protal uses normal pages. The profiling stage streams each SAM and, once a
sample's outputs are written, keeps only what its strain MSA rows need, packed: the alleles and
coverage of the genes of the species that enter MSAs (nothing with `--no_strains`; every taxon's
numbers with `--taxon_statistics`), 37 MB for a dense sample of 1M pairs from 60 species, 159 MB up
to 0.7.8; with `--strain_spill DIR` it goes to a file there and the strain stage reads it back per
species, so the memory of a run over many samples does not grow with them
([strains.md](strains.md#strain-msas-over-several-runs)).
