# Running protal: details beyond the website

The [website](https://protal.earlham.ac.uk/main.php?site=documentation#usage) explains how to run
protal on one sample or on many with a map file, and what each output file holds. This page adds
what it leaves out. `protal --help` lists the common options, `protal --full_help` all of them, and
`protal --map_help` the map file format.

## Where the outputs go

| Output | `-1 -2 --prefix P -o DIR` | `--map` |
|---|---|---|
| alignments | `DIR/P.sam.zst` | with a `SAM` column `#SAM_OUTPUT_DIR/<SAM>`, default `#OUTPUT_DIR/alignments/`; without one `#OUTPUT_DIR/<PREFIX>.sam.zst` |
| profile and its logs | `DIR/P.profile`, `DIR/P.profile.log`, ... | with a `PROFILE` column `#PROFILE_OUTPUT_DIR/<PROFILE>`, default `#OUTPUT_DIR/profiles/`; without one `#OUTPUT_DIR/<PREFIX>.profile` |
| strain MSAs and tables | `DIR/strains/` | `#STRAIN_OUTPUT_DIR`, default `#OUTPUT_DIR/strains/` |
| coverage, SNP counts, timers | `DIR/misc/` | `#MISC_OUTPUT_DIR`, default `#OUTPUT_DIR/misc/` |

- `-o` replaces a map's `#OUTPUT_DIR`, and `--profile_dir` moves the profiles of either mode.
- protal makes an output folder when it writes into it. A map without `SAM` and `PROFILE` columns leaves no
  `alignments/` or `profiles/`, and a run with `--no_strains` no `strains/`. Up to 2026-10-13 a map run made all four
  folders, and left the unused ones empty.
  Without `-o`, a `-1 -2` run writes into the current folder (runs without `-o` aborted up to 0.7.8).
  A relative output folder (`-o` or `#OUTPUT_DIR`) is relative to the folder protal runs in, and a
  relative `#SAM_OUTPUT_DIR`, `#PROFILE_OUTPUT_DIR`, ... lies in it. Up to 2026-10-08 a map run put a
  relative output folder twice in front of the folders the map did not name (`alignments`,
  `profiles`, `strains`, `misc`): `-o out` wrote `out/out/alignments` (an absolute folder, and folders
  the map names, were right).
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
  A sample's lanes in one folder (`x_S1_L001_...`, `x_S1_L002_...`) become one sample, `x_S1`, with
  its files listed in lane order ([Read files](#several-files-of-one-sample)); a sample of one lane
  keeps its name, and `--split-lanes` makes a sample of each lane, as `generate` did up to 2026-10-12.
  `protal_map_utils merge` keeps the SAMs the merged maps' runs wrote, by their absolute paths, so a
  run of the merged map profiles them instead of aligning again
  ([strains.md](strains.md#strain-msas-over-several-runs)).
- The SAM header (`@SQ`) lists the genes that the alignments name, not every gene of the database
  (the full r226 database has millions); `--full_sam_header` lists every gene, as protal did before.
- `misc/<sample>.err` lists the reads whose alignment does not fit the database (a gene it lacks, a
  position past a gene's end, bases that differ from the gene), as SAM records. They are left out of
  the profile, and protal warns with their number; a sample without any gets no file. Up to 0.7.8 it
  was `<sam>.err`, next to the SAM, which `--profile_only` may read from another run's folder.
- `misc/` also receives `P_seedsizes_histogram.tsv`, `P_anchorsizes_histogram.tsv` and
  `P_runtime.tsv`, diagnostics of the seeding and alignment stages. `P_runtime.tsv` has one row
  per stage (reading, k-mers, seeding and its steps, alignment and its k-mer screen, for long reads the indel bound (which ONT reads run), output): the seconds spent in it
  summed over threads, the number of threads, and the seconds per thread that `--verbose` prints.
  The stages that run for every read are timed on every 61st call only (reading the clock costs
  10-15% of the alignment time otherwise), so their seconds are estimates: the mean timed interval
  times the number of calls. Since 0.7.6 it also has the profiling steps' wall-clock times (reading the
  SAM, record evidence and sample context, read EM, congener distances, SNPs, scoring, writing), which
  protal prints in one line per sample too; a rerun that profiles the sample again replaces them. Since
  2026-10-12 a last row, "Profiling: after the profile", times what follows (the unreported species,
  freeing the sample's read data and evidence, spilling its strain evidence, returning freed memory),
  each part in the log line `Sample <name> after its profile: ...`. Freed memory goes back to the system
  (`malloc_trim`) after each sample's profile but the run's last, and not while another sample's is
  doing it: at r226 it took 0.4-0.5 s of a paired-end sample.
- `misc/cpu.tsv` (since 2026-10-09) has a row per stage of the run, one after the other: the start-up
  (loading the database's tables), `index and genomes` (the rest of the index load and the genome
  preload; since 2026-10-12, before which the first sample's row held them), each sample's alignment
  (`aligning <sample>`, up to its SAM written), the
  profiling stage and the rest of the run, with its wall-clock seconds, the CPU seconds the whole
  process used meanwhile (all threads, user and system time) and `-t`. CPU over wall is how many cores
  the stage kept busy. With `--profile_ahead` the profiling worker's time falls into the samples' rows.
- Every run prints per sample how its reads went: the reads, those with an anchor, the candidate
  alignments tried, those the k-mer screen refused, those the indel bound refused (ONT reads, in place of the k-mer
  screen, since 2026-10-09; 0 for other reads), those aligned, and the records written. A second
  line ("seeding:") counts the k-mer lookups, the index blocks they scanned and their sizes, the
  seeds, the seeds that share their taxon and gene with another seed of their read (the only ones
  anchors are made of), the lookups dropped as too ubiquitous (more tied entries than
  `--max_key_ubiquity`) and the anchors. At the end it times loading, freeing memory and every
  stage, so that a run's wall time adds up. Since 2026-10-11 a run that aligns loads the tables only profiling
  reads (gene conservation, suspect copies, species priors and neighbours, congener gaps, foreign rates, ancestry
  sites, and the strain alleles and column weights unless short reads' alignments are scored with them) beside the
  alignment: "Loading the taxonomy, models and tables" times what the alignment waits for, and the line "Loading the
  profiling's tables beside the alignment took ... s; the profiling waited ... s for them" (from the start-up to
  the last of them loaded, and the profiling's wait) comes with their messages when the first sample is
  profiled (at r226 the column weights alone had held every run back by 12 s;
  [report](claude/2026-10-11-performance-audit/README.md)). A table nothing in the run reads stays unread: the
  NCBI taxids without `--ncbi_taxids`, the suspect copies with `--keep_suspect_copies`, the column weights and
  ancestry sites with `--no_site_weights`, the gene neighbours with `--no_gene_neighbours`, the profiling's tables
  with `--no_profile`, and there also the strain alleles and column weights unless short reads are scored with
  them (not with `--no_allele_scores`, nor for long reads).
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
With `--profiles`, a sample whose own `.profile` ends with the unknown share keeps it, and the adjusted calls share the rest.

## What the called species explain: the unknown share

Since 2026-10-08 protal tells, per sample, how much of the sample the species it calls explain. A
species' depth on its marker genes, times its genome size (`species_priors.tsv`,
[databases.md](databases.md#species-priors)), is the bases its whole genome put into the sample.
Their sum over the called species, against the bases read (a pair's overlapping bases counted once,
as in the depth), is the share they explain. The rest came from organisms the profile does not
name: species the database lacks or the model did not call, plasmids, viruses, eukaryotes, host. At
the called species' average genome size, the rest is a number of genomes, and it ends the profile
as a last line, the unknown share:

```
GCF_000005845.2	d__Bacteria;...;s__Escherichia coli	0.412
...
?	?	0.183
```

The species' abundances are then shares of all genomes sequenced (of cells), adding up to 1 with the
`?` line; `--no_unknown_share` gives shares of the called species' summed depth as before. A
profile has no `?` line when the database has no genome sizes (converted before 2026-10-08; the run
says so when it loads it), or when the SAM does not say how many reads were read (aligned by protal
before 2026-10-08, or with `--full_sam_header`). A sample without called species gets `? ? 1`. Where
the called species explain more than was read (their depth or their sizes overestimated, or reads
of the marker genes alone, as in the tests), the share is 0.
`protal_profile_utils merge` puts the `?` row last, with NA for a profile without one, and
`scripts/prevalence_calls.py --profiles` keeps a sample's `?`.

The numbers behind it are in `<profile>.composition` (one header line and one line per sample;
`protal_profile_utils composition --input '*.profile.composition'` makes one table of many), written
with `--no_unknown_share` too, with NA where a number is not known:

| Column | |
|---|---|
| `ScannedFragments`, `ScannedReads`, `ScannedBases` | what the aligner read: pairs (or single reads), reads, bases |
| `FragmentBaseShare` | the aligned reads' fragment bases per read base: below 1 by the bases a pair's mates share |
| `CalledSpecies`, `SpeciesWithGenomeSize` | the species called, and those the database gives a size (one without counts at the average size) |
| `AttributedShare`, `AttributedFragments` | the share of the bases read (and as many fragments) that the called species' depth times genome size explains; above 1 where depths or sizes are overestimated |
| `AverageGenomeSize`, `GenomeSizeSD` | the called species' genome sizes weighted by their depth, which counts cells: the genome of the average cell (the average genome size of MicrobeCensus, here from the called species), and the spread around it |
| `GenomeSizeP10` to `GenomeSizeP90` | the depth-weighted 10%, 25%, 50%, 75% and 90% quantiles of the called species' genome sizes |
| `GenomeEquivalents`, `UnknownGenomeEquivalents` | the called species' depths summed (genomes sequenced at 1x), and the bases left over divided by the average genome size |
| `UnknownShare` | the profile's `?`: the unknown genome equivalents over all |
| `MedianDepth`, `LowestDepth` | of the called species |
| `MissingSpeciesAtMedianDepth`, `MissingSpeciesAtLowestDepth` | the unknown genome equivalents divided by the median and by the lowest depth of the called species: how many species the rest would be if each were as abundant as a typical called species, or as the least abundant one. Species below what protal detects would be more |

`.profile.log` gives each taxon's `GenomeSize` and `GenomeFragments` (the fragments its depth implies
over its whole genome), and the run prints a line per sample, e.g. `Sample S1: 81.7% of 1000000
fragments explained by the 45 species called; average genome size 3.42 Mb (10-90%: 2.10-5.31 Mb);
unknown 18.3% of genomes, about 12 species at the median depth (210 at the lowest)`. The aligner
writes what it read into the SAM header (`@CO protal scanned reads: fragments=F reads=R bases=B`),
so `--profile_only` and reruns compute the same.

The estimate leans on two assumptions: that a called species' genome in the sample has its GTDB
species' mean size (a strain's own size can differ by 10-20%), and that the organisms the profile
does not name have the called species' average genome size. Host or eukaryotic reads count as
unknown genomes of that size too, so in a host-rich sample remove host reads first or read the
`AttributedShare` rather than the `?`.

## Read files

Read files are FASTQ or FASTA, plain, gzip, BGZF or zstd, told apart by their first bytes (not
their names); bases are read as uppercase. Other compressions (bzip2, xz) are read through a pipe,
which protal reads as a file, as it does process substitution: `-1 <(xz -dc a_R1.fq.xz) -2 <(xz -dc
a_R2.fq.xz) --prefix a` (give `--prefix`: the pipe's name, `/dev/fd/63`, names no sample). A gzip
file may hold several members, one after the other (`cat a.fq.gz b.fq.gz`, also after BGZF ones),
and zero bytes of padding; a zstd file several frames (`cat a.fq.zst b.fq.zst`, `pzstd`'s and the
seekable format's skippable frames, any `zstd --long` window). Each file is decompressed in a
thread of its own, so the alignment threads only copy reads that are already decompressed. A BGZF
file (as `bgzip` and, by default, bcl2fastq write) is inflated on 1 to 4 threads per file (one per 4
alignment threads and file: 4 per file for a paired-end run on 32 threads), about 1 GB/s of FASTQ
each. A gzip file of many members that is not BGZF (members one after the other, as tools that
compress a file in chunks write it, its second member within 16 MiB of its start) is read the same
way: protal cuts it where each member's header begins and inflates whole members on those threads. A
header found by chance in a member's data, or a damaged member, sends the rest of the file to one
thread, so such a file reads, and fails, as on one thread. A file read in pieces (below) goes on by
members where small members follow one of its members (a large gzip file and a file of many members
concatenated). Another gzip file (one member, as `gzip` and most pipelines write it) cannot be cut at block
boundaries, as its blocks refer back to the 32 KB before them: protal inflates it with ISA-L in pieces
on 1 to 8 threads per file (one per 4 alignment threads and file: 2 per file for a paired-end run on
16 threads, 4 on 32). Pieces of the file are decoded ahead with the 32 KB before them unknown (marked,
and filled in once the piece before is done) while the file's own thread decodes the others in order;
on FASTQ about 1.7, 2.0 and 2.4 GB/s on 1, 2 and 4 threads (on 4 cores), as a piece decoded ahead costs
about 2.6 times one decoded in order. A zstd file's frames are independent, as BGZF blocks are, and are
decompressed on as many threads as a BGZF file's: protal's simulator (a frame per batch of reads) and
`pzstd` write many frames, `zstd -T` one (read on one thread all the same); a frame whose header does not
say its size (`cat reads.fq | zstd`) and the rest of the file after it are read as a stream. A pipe is
read on one thread. A regular file's next 64 MB are asked to be read ahead as it is read, so that a
file read from a network file system for the first time arrives before it is needed. After each
sample's alignment the log says how each of its files was read and how long the alignment waited for
them, e.g. `Sample s1 input: s1_R1.fq.gz gzip of 22606 members, inflated by members on 4 threads; ...;
the 32 alignment threads spent 1.20 s each getting reads, 5% of the alignment's 24.20 s; 32 alignment
and 10 reading threads on 84 CPUs`. A quarter or more of the alignment spent getting reads ("the input
held the alignment back") means the input was the limit; `--verbose` also says how each file is read
as it is opened. For runs on many threads, BGZF is still the
fastest and the cheapest: `bgzip -@ 8 -c reads.fq > reads.fq.gz` (or `zcat reads.fq.gz | bgzip -@ 8 >
reads.bgz.fq.gz`) turns a gzip file into BGZF. `PROTAL_INFLATE_THREADS=N` in the environment sets the
threads per gzip or zstd file.

A read file that is missing, a directory or not readable stops protal before the index is loaded.
A sample fails (exit 1, no SAM file) when its reads are not FASTQ or FASTA, when a record is
incomplete, has no `+` line or has more or fewer qualities than bases, when a gzip or zstd file is
truncated or corrupt (whatever follows a gzip member must be another member: a damaged member
header is an error, not the end of the file; a zstd file must not end inside a frame and holds
zstd frames only), when paired files hold different numbers of reads,
or when a single-end read is longer than 1,000 bp (below). A sample without reads, and pairs
whose mates' names differ (other than in a last `1` and `2`, as in `r/1` and `r/2`), get a
warning: the files of a pair may not be in the same order.

### Several files of one sample

A sample sequenced in several lanes or runs is one sample, not one per lane: its depth, the
features that depend on the sample's size and the read EM need all its reads. Give its files in
one map cell, comma-separated, in the same order in `FIRST` and `SECOND` (they pair up by position):

```
#SAMPLEID	PREFIX	FIRST	SECOND
x	x	x_L001_R1.fq.gz,x_L002_R1.fq.gz	x_L001_R2.fq.gz,x_L002_R2.fq.gz
```

or on the command line with `--one_sample`, which makes all files of `-1` (and `-2`) one sample:
`-1 x_L001_R1.fq.gz,x_L002_R1.fq.gz -2 x_L001_R2.fq.gz,x_L002_R2.fq.gz --one_sample`. Without
`--prefix`, such a sample is named after the beginning its files' names share, without a lane
number (`x_S1_L001_R1_001` and `x_S1_L002_R1_001` give `x_S1`); the files must then be in one folder.
The files are read one after the other as one stream, each in its own format (a gzip lane beside a
zstd one is fine), each relative to the map's `#INPUT_DIR`; a file whose last line lacks its
newline gets one. The profile is the one the reads would give in one pair of files. A lane whose
two files hold different numbers of reads fails the sample ("the two files of one of its lanes hold
different numbers of reads"), even where the sample's totals agree, instead of pairing the next
lane's reads with the wrong mates. A file that is truncated or corrupt fails the sample as a single
file does, and the message names it. A file name cannot contain a comma. `protal_map_utils
generate` writes such cells for Illumina's lane files, and `validate`, `flatten` and `merge`
handle them.

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
| `-t, --threads` | 1 | threads for alignment (which also compresses the SAM), database loading and profiling. Set it: the default is one thread. While aligning, each read file is also decompressed by a thread of its own (two for paired reads): gzip with ISA-L at about 1 GB/s of FASTQ, zstd with libzstd, and a BGZF file on 1 to 4 threads (`-t`/4 per file; see [Read files](#read-files)); at many threads a paired-end run of one-member gzip files can still wait for its input. Samples are profiled in parallel, the largest SAM first, each on threads in proportion to its SAM's share of all the samples' bytes (at least one; a single sample on all), with the same results as on one thread. A sample profiled on several threads also has its SAM read by a thread of its own; from 8 threads on, the zstd frames of a `.sam.zst` that protal wrote are decompressed ahead by `-t`/6 more threads (2 to 8) |
| `--knob` | 0.5 | detection threshold, 0 to 1. A model with a knob curve over the sample's depth uses the curve unless `--knob` is given; models trained with the depth as a feature (the default) have none ([databases.md](databases.md#how-protal-calls-species)). Choose it on data like yours |
| `--fdr` | off | report each sample's species while their expected share of false calls stays at or below this, for a model with calibrated calls; at r226 it called slightly below the knob curve, so it is off ([databases.md](databases.md#calls-at-a-target-share-of-false-calls)) |
| `--fragment_count` | `auto` | how the model's `fragments` feature counts: `mapq`, the fragments the MAPQ and length filters keep; `soft`, those plus each dropped fragment whose best record is on the taxon by the share the read EM leaves the taxon ([features.md](features.md)); `auto`, as the database records it in `build_settings.tsv` (`mapq` when it does not say). The models were trained with their database's count: a run that counts otherwise gets a warning. The profile's fragments, depth and the sample's depth keep the filters' count either way |
| `--singleton_congener` | 0 | veto a single-fragment species beside a congener of at least this many fragments when its read looks like the congener's; 0: no rule ([databases.md](databases.md#the-singleton-rule)) |
| `--no_unknown_share` | off | abundances as shares of the called species' summed depth (adding up to 1), without the profile's last line `?`, the share of the genomes the called species do not explain ([above](#what-the-called-species-explain-the-unknown-share)) |
| `--ncbi_taxids` | off | two more columns on every profile line: the species' NCBI species taxid and NCBI's name of it, from the database's `ncbi_taxids.tsv` (the taxid most of the species' GTDB genomes have; [databases.md](databases.md#ncbi-taxids)), e.g. `GCF_003697165.2  d__Bacteria;...;s__Escherichia coli  0.412  562  Escherichia coli`. Most GTDB species have only an NCBI placeholder ("uncultured Lachnospiraceae bacterium"), which several species can share: the lineage column still tells them apart. `NA NA` for a species the table lacks, `? ?` on the `?` line. A database without the table stops the run before it aligns; one converted before 2026-10-11 gets it with `gtdb_to_protal_db.py --priors_only` and `--add_tables` |
| `--keep_suspect_copies` | off | count reads on the database's suspect gene copies (near-identical to another genus's copy: contamination or a transferred gene) as evidence; by default they are left out, and the log says how many per sample ([databases.md](databases.md#2-build-the-index)) |
| `--depth_identity_margin` | 0.08 | a read counts towards a species' abundance only if its identity is at most this far below that of the species' best reads (98th percentile). Reads of relatives the database lacks still count for detection, not for depth. The margin is the same on every gene unless `--gene_conservation` scales it. The default counts the reads of strains up to about 5% from the reference; 0.04, the earlier default, dropped up to 40% of the reads of strains 3–5% away and undercounted those strains ([report](claude/2026-09-30-depth-margin-stress/README.md), [scaled per gene](claude/2026-10-01-gene-scaled-margin/README.md)). 1 lets every read count. The depth counts each base a fragment covers on a gene once: where a pair's mates overlap there, the second mate adds only what the first did not cover, as the strain MSA counts them ([report](claude/2026-10-02-fragment-depth/README.md)) |
| `--gene_conservation` | `none` | scale `--depth_identity_margin` per gene by its conservation factor: `db` for the database's, or a file of `geneid<TAB>factor`; a gene's margin is then 0.03 plus the rest times its factor. `none` did best over three simulated worlds ([report](claude/2026-10-01-gene-scaled-margin/README.md)). The model's conservation features use the factors either way |
| `--model` | `model.xml` of the database (`model_se.xml` for single-end samples) | a PMML file, or the name of another model in the database folder (`<name>.xml`); for all samples unless `--model_se` is given. protal checks the model before aligning, see [databases.md](databases.md#the-presence-model) |
| `--model_se` | `--model`, else `model_se.xml` of the database | the model of single-end samples, given as `--model` |
| `--model_kind` | `auto` | which kind of the database's models profiles the samples: `gbm` (gradient-boosted trees) or `forest` (a random forest), for every read type or per read type (`pe:forest,pb:gbm`); `auto`: the database's default models. A GTDB build since 2026-10-13 holds both kinds at every gate (`model_pe.xml` and `model_pe_forest.xml`, `model_pe_c1.xml` and `model_pe_c1_forest.xml`, ...); each sample takes the model of the kind at the gate nearest its own. With `--profile_only` a run's SAMs are profiled again with the other kind. A read type without a model of the kind stops the run; `--model` and the like take precedence ([databases.md](databases.md#kinds-of-models)) |
| `--sam_format` | `zst` | the format of the SAM files protal names: `zst` (`.sam.zst`), `gz` (`.sam.gz`) or `sam`; a map's `SAM` names keep their own ending |
| `--scratch_dir` | none | a folder for the run's intermediate files, best on the machine's own disk (e.g. `$TMPDIR` on a cluster node) when the outputs go to a network file system. protal makes a folder of its own in it (`protal_<host>_<pid>_<time>`) and removes it at the end of the run, also when the run fails (a killed run leaves it). In it: each sample's strain evidence, read back per species ([strains.md](strains.md#strain-msas-over-several-runs); this replaces `--strain_spill`, which still works and means the same); with `--no_strains`, the SAM files protal names (with `-1`/`-2`, or a map without a `SAM` column), removed once their sample is profiled, as nothing after the run reads them (writing and closing a SAM on a network file system took 1-7 s of a 32-thread r226 run); the records of a `--sam_format gz` or `sam` SAM before its header is written; qcMSA's working files. Without `--no_strains` the SAMs stay in the output folder, for strain MSAs over several runs, and a SAM an earlier run left there is profiled as before |
| `--one_sample` | off | the files of `-1` (and `-2`) are one sample's, its lanes or runs, read one after the other ([Several files of one sample](#several-files-of-one-sample)); in a map, a `FIRST` and `SECOND` cell lists a sample's files, comma-separated |
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
| `--wfa_traceback` | `auto` | how WFA2 recovers each alignment's operations: `classic` (WFA2's high-memory mode, which keeps every wavefront) or `singletrack` (keeps only the match wavefronts and recomputes a gap's steps while tracing back; López-Villellas et al. 2025). Both give the same alignment scores, but among equally good alignments they may choose differently, so SAM records, and the features a model was trained on, can differ. `auto` aligns as the database was built: `--build` aligns with `singletrack` unless told otherwise (`classic` before protal 0.8.2) and records its traceback in `build_settings.tsv` ([databases.md](databases.md#2-build-the-index)); a database without the file was built with `classic`. A run that aligns otherwise than its database was built gets a warning ([report](claude/2026-10-11-singletrack/README.md)) |
| `--no_mate_guidance` | off | paired-end reads: do not let a mate that is sure of its alignment (MAPQ 20 or more) guide the other when they did not align together (below) |
| `--no_alignment_screen` | off | align every candidate with WFA2. By default a candidate whose read and gene window share too few k-mers for any alignment within the score budget is refused before WFA2 (exact: it changes no alignment; at r226 it took the paired-end alignment from 62 to 39 s); also turns off the indel bound below. For measuring |
| `--no_indel_bound` | off | ONT reads: screen the candidates by their k-mers, as other reads (and ONT reads before 2026-10-09). By default ONT reads are screened by their indel distance instead: a candidate whose read and gene window are too far apart by it for any alignment within the score budget is refused before WFA2. At ONT's identity floor of 0.85 the k-mer screen refuses next to nothing, while most failing candidates are unrelated sequence; PacBio reads keep the k-mer screen, as the candidates it passes mostly align and the bound cost more than it saved there. Exact: neither changes an alignment. The `P_runtime.tsv` row "Indel bound" times it. For measuring |
| `--flex_indels TYPES` | `none` | score the flex cells of every seed lookup (the 8 bases either side of the k-mer's 15-base core, against the read's) allowing one insertion or deletion in each flank, for reads of these types (`ont`, `pb,ont`, `all`; a map's optional `FLEX_INDELS` column, `1`, `0` or `-`, sets one sample's). Without it an indel in a flank shifts the bases beyond it, which then match only by chance: one next to the core costs ~6 of the flank's 8 bases, and a relative's cell may take the seed. With it a single indel costs one base, as a mismatch does; only an exact match still scores 16, so the unique seeds are the same. Each flank scores the most of its bases equal in place and, from its first differing base on, the read's bases one place in or out. Meant for ONT reads, whose errors are mostly indels, and for gated runs (`--sparse_first`, `--sparse_seeding`), where each lookup counts more; the gate's `--sparse_first_score` counts the bases so too. The scan takes ~3-4 times the instructions per cell. It changes seeds, records and the models' features, and the database's models are trained without it ([report](claude/2026-10-13-flex-indels/README.md)) |
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
`ZF:Z:<taxid>:<gene>,...` when the read seeded on taxa strongly enough to be aligned against them (the
align-top anchors) but did not align to them, each with the gene of its longest anchor (since
2026-10-06; `<taxid>` alone before); the reads that seeded on taxa and aligned nowhere are
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
Since 2026-10-06 a read's best record (each mate's, each long-read gene's) also carries
`ZN:i:<taxa>`: the taxa whose anchors were at least 0.8 as long as the read's longest, how many
species its seeds could not tell apart, also those it was never aligned against; the profiler takes
`seed_crowding`, the genes of `ZF` and the long reads' `ZR` as features against the false positives
of complex communities ([features.md](features.md#against-false-positives-in-complex-communities-consistency-shape-neighbourhood-2026-10-06)).
Since 2026-10-07 a short read's first record also carries `ZC:Z:<taxid>,...`: up to 8 taxa of that crowd it was never
aligned against (beyond `--align_top`), the strongest anchor first; the profiler counts them per taxon
(`untried_candidate_rate`). And a short read whose best alignment is divergent (identity below 0.99) tries up to
`--adaptive_candidates` (7) more anchors of its crowd, of taxa of that alignment's genus: a strain whose genes lie
between its own species' reference and a congener's ranks its own species below the congeners by its seeds, and
`--align_top` alone never aligned it there ([report](claude/2026-10-07-error-read-signatures/README.md), section 5).
`--adaptive_candidates 0` aligns as before; the log's `adaptive candidates:` line counts the extra alignments.
Since 2026-10-08, with a database that has strain alleles (`strain_alleles.tsv`, [databases.md](databases.md#2-build-the-index)),
each candidate alignment of a short read gets a site shift against its species' known alleles:
- the read's differences from the reference that an allele has count as matches;
- a difference at a site where the species' alleles vary, but none has the read's base, counts as half a difference;
- with column weights in the database (`column_weights.tsv`,
  [features.md](features.md#the-columns-conservation-weights-2026-10-10)), a remaining substitution at a column the
  family's genera change freely (below 3 nats of among-genus conservation) counts half a difference; a copy without
  alleles gets this part of the shift too. On simulated worlds whose every column is variable it cost 0.006 of F1
  ([report](claude/2026-10-09-site-weighted-evidence-plan/README.md), section 9); a GTDB database, whose families
  have conserved columns, has to show whether it pays. Between 2026-10-10 (when it became an option,
  `--weighted_site_shift`) and 2026-10-11 the option had no effect: the alignment threads' copies of the run's settings
  dropped it ([report](claude/2026-10-11-r226-v23/README.md), section 6). Fixed, it was neutral to +0.003 of F1 on a
  small real world (section 8 there) at ~1% of the alignment's instructions. It was on by default from 2026-10-11 to
  2026-10-12 and is off again: at r226 (v24b) it lowered the short reads' species depths by 1-2% without a measurable
  F1 gain, and on the small real world, with the same samples, it lowered 19-22% of the called species' depths
  ([report](claude/2026-10-12-r226-v24b/README.md), section 4). `--weighted_site_shift` turns it on. A database built
  in those two days (r226 v24b; its `build_metadata.tsv` says `weighted_site_shift on`) trained its models with it:
  profile it with `--weighted_site_shift` for the calls of its training. Until the
  evening of 2026-10-11 a run with `--no_profile` left the column weights unread, so its short reads were aligned
  without this part of the shift and its SAM differed from a profiled run's; now it reads them too, and the two
  write the same records ([report](claude/2026-10-11-performance-audit/README.md), follow-up).

Since 2026-10-09 the shifts decide only the reads protal is unsure about: a candidate of another species within 3
mismatches of the best by the references alone. Then a read of a known strain scores on its species as on the strain's
gene. A read one species fits clearly better keeps the references' scores. The records keep their alignment to the
reference (CIGAR, identity), while `AS` and MAPQ count the shifts. The log's `strain alleles:` line counts:
- the candidates on copies with alleles, and those with a shift;
- the unsure reads that took the shifts, and those whose best species they changed.

A database's models are trained with these scores; `--no_allele_scores` turns them off, to measure them
([features.md](features.md#the-strain-alleles-alleles-2026-10-08)).

### Developer options

Also shown by `--full_help`: `--build` and its options ([databases.md](databases.md#building-a-database)),
the database conversions `--compress_db`, `--unpack_db`, `--decompress_db`
([databases.md](databases.md#the-files-of-a-database)), `--write_species_neighbours FILE` (every species' nearest
congeners by their references, from a folder of separate files, without building: the build script's
`species_clouds.tsv`), `--profile_truth` for the training dump
([databases.md](databases.md#training-data)), `--benchmark_alignment` (checks alignments against the
`taxid_geneid` encoded in simulated read names), `--mapq_debug_output`, `--full_sam_header` (every
gene in the SAM header, see above), `--write_unmapped_reads` (an unmapped record for each read that
aligned nowhere, see above), `--whole_read_alignment` and `--profile_after_alignment` (both below),
`--sequential_load` (load the database's parts one after another, as before 0.7.6, instead of the
index beside the genome preload and the small tables each on a thread, the profiling's tables beside the
alignment; the run is the same either way).

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

`--stage1_probe` measures, in an ordinary run, what the first stage of a two-stage protal would keep:
for every lookup, its best flex score and whether its k-mer's core falls in a sparse sample of 1/c of
the cores (c = 1 to 32), tallied per read (a pair as one) for the rules "at least m lookups scoring T
or more" over all reads and those with records, per taxon, and per mate the best score of its first
lookups. It writes `misc/<sample>_stage1_probe.tsv`, `_stage1_probe_taxa.tsv` and
`_first_lookups_probe.tsv` and a log line; the SAM and profile are the same as without it
([report](claude/2026-10-12-stage1-probe/README.md)). Since 2026-10-13 it also tallies the rules for
a stride pick, every c-th core of each mate in read order (spread evenly over the read, where the hash
sample leaves gaps by chance), into `_stage1_probe_stride.tsv` and the taxa table's `s<c>_t<T>_m<m>`
columns, and gives each rule the share of the lookups a gate with it would make (`lookups_left`)
([report](claude/2026-10-13-stride-pick/README.md)).

`--sparse_first C` acts on such a rule: a read (pair; a long read over all its chunks) is looked up
by the 1/C of its k-mer cores in the sample first, and given up, before the rest of its seeding and
any alignment, unless at least `--sparse_first_hits` (1) of those lookups have a flex cell sharing
`--sparse_first_score` (14) or more of the 16 flank bases with the read. A read kept is seeded and
aligned as without it. At GTDB r226 the probe found C = 2 keeps 9.4% of a gut sample's pairs with
99.6% of those that have records, and C = 8 18% of a HiFi sample's reads with 99.2%. Reads given up
count no failed candidates, which the models' `failed_candidate_rate` uses, so the models must be
trained on runs gated as the run is ([report](claude/2026-10-12-sparse-first/README.md)). A database
built with `--sparse_first` records each read type's gate in its `build_settings.tsv`, and runs on it
gate as it says (`auto`, the default; a database without one: off). `--sparse_first` takes one gate for
all read types (`0` for off, or 1, 2, 4, 8, 16, 32) or read types and theirs (`pe:1,se:1,pb:8,ont:1`;
those not named keep the database's), and a map's optional `SPARSE_FIRST` column one sample's (`-`:
its read type's); `--sparse_first_score` and `--sparse_first_hits` likewise default to the database's
(else 14 and 1). The run says the gates it uses beside the database's, and warns when a sample is
gated otherwise. `--stage1_probe` and an explicit gate exclude each other; with `--stage1_probe` the
database's gates are not followed.

`--sparse_first_pick` (since 2026-10-13) chooses which 1/C of a read's cores the gate looks up:
`hash`, the cores whose hash falls in 1/C of its range (the same cores in every read, so a read's
sample has gaps by chance; the only pick before), or `stride`, every C-th core of each mate (each
chunk of a long read) in read order, spread evenly over the read. On local test data the stride kept
far more of the pairs with records at the same C (C = 8: 99.4% against 92.0%) for a few more lookups,
and at a recall of 99% needed C = 8 where the hash needed C = 2
([report](claude/2026-10-13-stride-pick/README.md)); at C = 1 both look up every core. At GTDB r226
the stride lost 2.5-17 times fewer present species than the hash at equal or better F1, and it is the
default since 2026-10-13 ([report](claude/2026-10-13-r226-v27-c-gates.md)). Like the gates, `--build`
records the picks (`build_settings.tsv`; `stride` where not given), runs follow the database's
(`auto`, the default; `stride` where it records none, `hash` for a database that records its gates but
no pick, built on 2026-10-12/13), `--sparse_first_pick` takes one pick for all read types or read types
and theirs (`pe:stride,pb:hash`), a map's optional `SPARSE_FIRST_PICK` column one sample's (`hash`,
`stride`, or `-`: its read type's), and a sample picked otherwise than its database is warned.
`--sparse_seeding` stands for an index of the hash pick's cores: a sample it would seed with the
stride pick stops the run (give `--sparse_first_pick hash` with it).

A database may hold models trained at other gates than its own (since 2026-10-13: `model_pe_c1.xml`,
`model_PB_c32.xml`, ..., which `build_gtdb_database.py --model-gates` trains and
`--add_model FILE --read_type TYPE --model_gate C` stores). A sample is then profiled with the model
trained at the gate nearest its own (on the scale of log2 C, off next to C 1; between two, the lower),
and is warned only if that gate is not its own. The run names each model it loads with its gate
("Model of paired-end reads gated at 1: model_pe_c1.xml ..."), and the gates line says which other gates
the database has models of. With `--profile_only` the gate (`--sparse_first`, else the database's)
only chooses the model. With `--model_kind forest` (or `gbm`) the gates are those the database has
models of that kind at, and each sample takes that kind's model of the gate nearest its own
(`model_pe_c1_forest.xml`, ...; [databases.md](databases.md#kinds-of-models)).

`--sparse_seeding` (with a gate) seeds the reads the gate keeps from its 1/C of the cores alone, as an
index holding only those cores would, and logs after the index load what such an index would hold
of the loaded one's entries, values and keys for each C. It measures what a smaller index for gated
runs would cost in records and calls; the whole index is still loaded, so it saves no memory itself.
Short reads at C = 1 are seeded as without it ([report](claude/2026-10-13-sparse-seeding/README.md)).

**The grid pick and the sparse index.** The stride and hash picks choose a read's cores on the read, so any key may be one a gate looks up and
the whole index must be in memory. The grid pick (`--sparse_first_pick grid`, since 2026-10-13, on a
database built with grid levels: `Index features: ... grid levels`) chooses them on the reference
instead: `--build` keeps one core per window of C consecutive cores of each gene at each gate C, the one
a read of the species most likely carries intact (by the species' strain alleles, else by how
conserved its columns are within the family's genera) and the least shared with other genes
([databases.md](databases.md#2-build-the-index), step 5). A read cannot tell which of its cores
those are, so it looks up every core and only the grid's entries count: about every C-th core of the
stretch it covers, as with the stride. A read kept is seeded from those entries alone, so it has about
1/C of its seeds; a short read's gene with a single seed anchors on it (two seeds of a gene on one
diagonal otherwise, as in every other run), since one grid hit fixes the read's place on its gene. On the
mini database that turned the paired-end runs' 3.2% fewer records at C 4 into 0.05%, for 1-5% more
candidate alignments than an ungated run; `--no_single_seed_anchors` anchors on two seeds only, for
comparison. Because only
those entries are ever used, a run whose every sample is gated with the grid pick (C 2 or more)
loads only them (`--sparse_index auto`, the default: the entries of grid level log2 C of its smallest
gate, about 1/C of the values; the key map keeps its size), and the index line says so ("a sparse
index: the entries of grid level 2 or more, N of the file's M, on K of its J cores with entries"). On the whole index (`--sparse_index
off`, or a run whose samples are not all grid-picked) the grid only gates: the same reads are given up,
and a read kept is seeded from every core, a gene anchored on two seeds, as in other runs (the run's
"Grid pick:" line says which). `--sparse_index emulate` loads the whole index and seeds the
grid-picked samples as the sparse index does (the same records as `auto`; for comparison, and for the
grid rows of `build_gtdb_database.py`, whose collections' runs load the whole index); `on` stops where
it cannot load a sparse index. A database without
grid levels stops a grid-picked run (exit 2). The models must be trained on runs picked so, as for
the other picks (`build_gtdb_database.py --test-grid` scores the test set with it beside the database's
pick; [report](claude/2026-10-13-gene-grid-index.md)).

`--compact_keys` (since 2026-10-14) keeps the index's key map for the cores with values alone: a
presence bit for each of the 2^30 cores with ranks (153 MB) and 3 bytes a used core, instead of
16-bit offsets for all of them (3.2 GB, of which r226 uses 9.2%). The load builds it as it decodes, so
the direct key map is never held: on the reduced r226 database a run peaks at 0.8 GB instead of 3.4
GB, and at r226 the key map would take ~0.45 GB (less on a grid's sparse index). The records are the
same; a lookup reads a second line of the key map. Since 2026-10-14 that line is fetched with the read's
keys where those come a read ahead, the map is on huge pages and its ranks are counted with the CPU's
popcount: on the reduced r226 database, one thread, its seeding takes 56 ns a k-mer for paired-end reads
and 38 for HiFi, against 55 and 47 with the direct key map (before: 86 and 73, against 78 and 58), and
a grid's sparse index loads faster with it than without (0.65 against 0.96 s for the strain world's at
C 4; [report](claude/2026-10-14-sparse-load-key-map/README.md),
[first report](claude/2026-10-13-compact-key-map.md)). Off by default until measured at r226.

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
its key map, 3.2 GB, ~0.45 GB with `--compact_keys` (above); the genes at two bits per base, the tables and the run's
buffers, the rest;
[databases.md](databases.md#the-database-in-memory)) and correspondingly less for the
reduced one ([downloads](https://protal.earlham.ac.uk/main.php?site=downloads); the figures there,
59 and 12 GB, predate the 2-bit genes and the packed index).
It prints the machine's total memory at start and, after loading the index, the memory it takes. The index is read at random, one lookup per k-mer,
so protal asks Linux for transparent huge pages for it; the usual setting (`madvise` in
`/sys/kernel/mm/transparent_hugepage/enabled`) grants them, and seeding is about a third faster.
With THP set to `never` protal uses normal pages. The profiling stage streams each SAM and, once a
sample's outputs are written, drops its record evidence (since 2026-10-11: a per-taxon entry for every
taxon with a record or a failed candidate, ~140 MB per sample at r226, which a run of many samples kept to
its end) and keeps only what its strain MSA rows need, packed: the alleles and
coverage of the genes of the species that enter MSAs (nothing with `--no_strains`; every taxon's
numbers with `--taxon_statistics`), 37 MB for a dense sample of 1M pairs from 60 species, 159 MB up
to 0.7.8; with `--scratch_dir DIR` it goes to a file in the run's folder there and the strain stage reads it back per
species, so the memory of a run over many samples does not grow with them
([strains.md](strains.md#strain-msas-over-several-runs)).
