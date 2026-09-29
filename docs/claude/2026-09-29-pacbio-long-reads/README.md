# PacBio long reads: native alignment, chunking, and splitting reads compared

- **Date**: 2026-09-29.
- **Code**: branch `pacbio-long-reads` (worktree `.claude/worktrees/pacbio-long-reads`), based on
  `audit-fixes` at `995c4f1` plus the single-end support described in
  [2026-09-29-single-end-and-long-reads.md](../2026-09-29-single-end-and-long-reads.md) (in the
  main checkout's `docs/claude/`); nothing committed yet. `docs/` is not in git on `995c4f1`, so
  this branch carries only this report; the user documentation for PacBio reads is at its end, to
  go into `docs/` with the branch, and the report's line in `docs/claude/README.md` too.
- **Machine**: WSL Ubuntu 24.04 (gcc 13), 8 threads.
- **Data**: mini databases of `scripts/mini_db/build_mini_db.sh` (3 species, 3 genomes each; the
  first genome of a species is the reference, the others strains), with the default divergences
  (species 3.5% from their genus, strains 0.5%) and with close relatives
  (`--species_divergence 0.005 --strain_divergence 0.002`, species about 1% apart). Reads simulated
  from the genomes of `examples/mini_db/community.tsv` (Mockella alpha 0.5 from a strain, Mockella
  beta 0.3 from its reference, Fakibacter gamma 0.2 from a strain).

## Question

Integrate PacBio (HiFi) reads so that they pass through protal without problems; split reads
longer than the ~65 kb that seeding's 16-bit read positions allow; and evaluate whether it is
cleaner to split long reads by default, or to align the long read as a whole.

## Answer

Align long reads whole ("native"), and split only reads longer than the 16-bit limit, into
overlapping chunks. Splitting every read into short pieces is less code but worse wherever it
matters: with species about 1% apart, pieces lose a quarter to a third of the marker bases to
ambiguity and misestimate abundances 5-6 times worse (table 2), and they cannot use the read's
other genes to settle an ambiguous one, which native alignment can (87% of its remaining ambiguous
bases, section "Next"). Chunking every read by default would change nothing for HiFi reads (15-25
kb) and adds boundaries that must be handled anyway, so chunks are used only above the limit.

## Design as implemented

A PacBio sample (`--read_type pacbio`, or `pacbio` in a map's new `READ_TYPE` column) is aligned
by `LongReadAligner` (`src/Core/LongReads.h`):

1. **Seeding** of the whole read with all its k-mers (the short-read stop after `-s` seeds would
   leave most of a long read unseeded). A read longer than 65,000 bp is seeded in chunks of at most
   65,000 bp that overlap by the database's longest gene plus margins (`ChunkRead`); a gene hit
   belongs to the chunk whose core (split at the middle of each overlap, and reaching beyond the
   read's ends) holds the gene's centre, so each gene lies wholly in one chunk and is counted once.
2. **Segments**: the anchors (one per gene and taxon, from the existing anchor finder) are grouped
   by where their gene lies on the read: genes that cover half of the shorter one's place are one
   segment, a gene and its homologs in other taxa. Adjacent genes of an operon are separate segments.
3. **Windowed alignment**: per segment, the `--align_top` longest anchors are aligned with the
   existing WFA alignment (`SimpleAlignmentHandler::AlignAnchor`), but only over the read's window
   on the gene (the gene's place plus 100 bp + 5% of its length on each side), not the whole read.
4. **Ranking and MAPQ** per segment: bitscore as for short reads, each alignment once, MAPQ of the
   best against the second best (`MAPQv2`).
5. **SAM** (`ProtalLongReadOutputHandler`): the best segment's best hit is the primary record, the
   other segments' are supplementary (0x800), each with its segment's MAPQ; alternatives (with
   `-m` > 1) are secondary (0x100). All records are hard-clipped to their aligned bases (SEQ and QUAL
   are the aligned part, `nH` the rest of the read), so a 15 kb read with ten genes writes about
   ten genes' worth of sequence, not ten reads' worth. The header names the reads
   (`@CO<tab>protal reads: PacBio`), which `--profile_only` uses to pick the model.
6. **Profiling**: `ReadSamGroups` starts a read's group anew at each supplementary record, so each
   segment is profiled as a read. Nothing else in the profiler changed.
7. **Model**: `model_pacbio.xml` of the database, or `--model_pacbio` (or `--model`); packed by
   `--build` and `--compress_db` like `model_se.xml`. Every kind of read in a run needs its model,
   checked before aligning.

Also: reads given as short single-end reads but longer than 1,000 bp (the first 100 are checked)
stop the run before aligning, with a message pointing to `--read_type pacbio`; the per-sample read
type (`ReadType` in `src/ReadType.h`: paired-end, single-end, PacBio) replaces the single-end flag in
`Options`, `RunProtal` and the model loading; anchors' and alignments' unique-seed counts are 16 bit
(were 8, which a whole-gene anchor overflows).

## Tests

- Unit (`just test`): 118 of 118 pass. New in `tests/test_LongReads.cpp`: the chunk layout, and
  that every stretch up to the overlap, including one cut by a read end, lies in the chunk owning
  it; segment grouping; ranking and MAPQ; primary, supplementary and secondary records with hard
  clips and SEQ in reference orientation; read back as one group per segment; a skipped
  inconsistent hit. In `test_Parsing.cpp`: the `READ_TYPE` column, PacBio from the SAM header, the
  model per read type.
- End-to-end (`tests/e2e/test_protal_e2e.py`): 73 of 73 pass, 8 of them new in `PacBioTest`:
  synthetic long reads of 3-8 reference genes between random DNA, 0.1% errors, half the reads
  starting and half ending inside a gene, one read of 150 kb (chunked); every gene found once, at
  95% or more of its length; records hard-clipped and consistent with the read; one primary per
  read; profiles; `--profile_only` choosing the PacBio model from the header; the missing-model and
  long-reads-as-short-reads checks; a map with `READ_TYPE`.

## Evaluation: native vs split

`run.sh` (scripts in this folder) simulates 20 Mb of HiFi-like reads (`simulate_hifi.py`: lengths
N(15,000, 3,000) within 5-30 kb, 0.1% errors: 40% substitutions, 30% 1 bp insertions, 30% 1 bp
deletions), profiles them natively and, cut into 150, 250 and 1,000 bp pieces (`split_reads.py`),
as single-end reads. The database's `model.xml` stands in for `model_pacbio.xml` and
`model_se.xml`, so detection calls are not meaningful, but alignments and depths are.
`evaluate.py` counts the records the profiler takes (not secondary, MAPQ >= 4, more than 50 aligned
bases): their reference bases on a gene of the read's own species or of another, against the reads'
bases in marker genes of their genome (`marker_positions.tsv`). Abundance error is half the L1
distance of the VCov shares from the cell abundances.

    B=$HOME/protal-lr-build bash run.sh                      # default mini DB, in $B/mini_db
    bash close_db.sh && MINI=$B/mini_db_close bash run.sh    # close relatives

**Table 1: default divergences** (species 3.5% from their genus), 1,347 reads:

| | align (s) | SAM records | recovered | misassigned | MAPQ < 4 (bp) | abundance error |
|---|---|---|---|---|---|---|
| native | 1.3 | 8,986 | **0.984** | 1.28% | 0 | 0.008 |
| split 150 bp | 0.9 | 66,858 | 0.966 | 1.01% | 11,495 | 0.009 |
| split 250 bp | 0.6 | 43,699 | 0.974 | 1.10% | 2,274 | 0.008 |
| split 1000 bp | 0.8 | 17,434 | 0.981 | 1.27% | 725 | 0.007 |

Species this far apart are unambiguous even in 150 bp pieces, and all approaches do well. The
misassigned bases are about the same in absolute terms (94-122 kb). They come from marker genes a
strain has but its species' reference lacks (`--marker_loss`), which align to a relative whatever
the read length; short pieces leave more of them unaligned.

**Table 2: close relatives** (species 0.5% from their genus, about 1% apart), 1,337 reads:

| | recovered | misassigned | MAPQ < 4 (bp) | abundance error | alpha / beta / gamma (truth 0.5 / 0.3 / 0.2) |
|---|---|---|---|---|---|
| native | **0.784** | **2.9%** | 1.77 M | **0.028** | 0.472 / 0.327 / 0.201 |
| split 150 bp | 0.580 | 4.5% | 3.51 M | 0.151 | 0.349 / 0.325 / 0.326 |
| split 250 bp | 0.600 | 3.9% | 3.41 M | 0.158 | 0.342 / 0.341 / 0.318 |
| split 1000 bp | 0.670 | 3.5% | 2.83 M | 0.157 | 0.343 / 0.381 / 0.277 |

With close relatives, a piece of a gene often fits both species equally (MAPQ < 4, dropped), the
more so the shorter it is; what is dropped is mostly Mockella alpha and beta, whose shares then
collapse towards each other. Whole-gene alignments separate them far better. SAM files: native
20 MB, split 25-36 MB for the same reads. Alignment times are comparable (0.5-1.3 s here).

**A bug the evaluation found.** The first native run recovered 0.960 in table 1's setting, less than
1 kb pieces. `diagnose.py` showed 393 of 430 missed genes to be genes the read starts or ends in:
their centre lies beyond the read, outside every chunk's core, so they were dropped. The first and
last chunk's cores now reach beyond the read's ends (fixed before the tables above were made, and
covered by unit and e2e tests).

## Next

- **Read-level consensus.** The native run of table 2 still drops 1.77 M bp as ambiguous. On the
  same reads (`rescue.py`), 87% of these bases lie on reads whose confident genes (MAPQ >= 4) all
  name one species, and that species was the read's true one for all of them. Giving such
  ambiguous segments to that species would lift recovery from 0.78 towards 0.95. It changes what
  MAPQ means for those records (the evidence is the read's other genes), so it is a decision for
  the model's design, not made here.
- **`model_pacbio.xml`**: to be trained on long-read simulations (pbsim3 or similar); until then
  PacBio samples need `--model_pacbio`. The features change meaning (a "hit" is a gene segment, and
  one read gives many correlated ones).
- **FASTA input** writes no base qualities, and the profiler skips records without them, in every
  mode (short reads too). HiFi reads often come as FASTA or unaligned BAM: a constant quality, or
  BAM input, would be needed.
- **ONT**: the chunking handles the length; chaining (within 6 bp of the first seed's diagonal),
  scoring and the SNP filters still assume few indels (see the single-end report's evaluation).

## Follow-up: read-type names as in the database build

The read types and model files above were renamed to those of branch `zstd-compression`, whose
database build writes them: `--read_type pe|se|pb` (instead of `short|pacbio`; without it, `pe`
with a second read file and `se` without), `model_pe.xml` (or `model.xml` of older databases),
`model_se.xml` and `model_PB.xml` (instead of `model_pacbio.xml`), `--model_pb` (instead of
`--model_pacbio`), and the SAM header line `@CO<tab>protal read type: pb`. Behaviour is unchanged.

## Follow-up: genes of a read vote; FASTA qualities

**Read-level consensus** (`TaxonOfRead`, `SettleByRead` in `src/Core/LongReads.h`), as proposed
under "Next". After a read's segments are aligned, each confident segment (MAPQ >= 4) votes for the
taxon of its best hit; a taxon with at least two thirds of the votes is the read's taxon. An
ambiguous segment (MAPQ < 4) takes the hit of that taxon as its best if its best cannot be told
apart from it (the best's MAPQ against it is below 4); if the taxon's anchor at that place was not
among the `--align_top` aligned, it is aligned first. The settled segment gets the lowest MAPQ of
the segments that voted for the taxon, and its record is tagged `ZR:i:1`. protal reports how many
ambiguous gene hits were settled.

Votes count equally. The first version weighted them by MAPQ, which settled 7.6% of the settled
bases wrongly (`settled.py`): a marker gene that the read's species lacks in the database (the mini
GTDB's `--marker_loss`) aligns to a relative with no competitor and MAPQ ~140, and outvoted five of
the species' own genes of MAPQ 5-15. With equal votes, none of 1.75 Mb settled was wrong.

The native rows of table 2 (close relatives) with the consensus; table 1 is unchanged (no
ambiguous genes there):

| native | recovered | misassigned | MAPQ < 4 (bp) | abundance error | alpha / beta / gamma |
|---|---|---|---|---|---|
| without consensus | 0.784 | 2.9% | 1.77 M | 0.028 | 0.472 / 0.327 / 0.201 |
| votes weighted by MAPQ | 0.947 | 3.8% | 0.10 M | 0.020 | 0.491 / 0.320 / 0.189 |
| **equal votes** | **0.969** | **2.4%** | 0.02 M | **0.015** | 0.507 / 0.307 / 0.185 |

2,249 of 2,280 ambiguous gene hits were settled. The misassigned bases (222 kb) are those without
consensus: genes the species lacks in the database, which no read-level rule can move.

**FASTA input** now gets a constant quality per read type: Q30 (`?`) for paired-end, single-end and
PacBio reads, filled in by the readers (`FillMissingQuality`), so that the profiler takes their
records; before, it skipped every record without qualities, in every mode.

Tests: unit 122 of 122 (consensus votes and settling, the `ZR` tag, FASTA qualities of single-end
and paired readers); end-to-end 75 of 75 (FASTA input of single-end and PacBio reads, the log line).

## Documentation to add with the branch

`docs/running.md`, a section after "Single-end reads":

> ## PacBio long reads
>
> `--read_type pb -1 reads.fq.gz` aligns PacBio long reads (HiFi); in a map, `pb` in a
> `READ_TYPE` column (`pe`, `se`, or `-` for pe or se by the `SECOND` column), with `-` as
> `SECOND`. Each read is seeded whole, each marker gene it spans is aligned over the read's window
> on that gene, and each gene hit gets its own MAPQ against the homologs of other taxa at the same
> place. A read's hits are written as one primary and supplementary (0x800) records, hard-clipped
> to the aligned bases; the profiler counts each as a read. Reads longer than 65,000 bp are seeded
> in overlapping chunks (protal says how many). A gene hit that fits several taxa equally is
> settled by the read's other genes: if two thirds of the read's confident genes name one taxon,
> the gene takes that taxon's hit (tagged `ZR:i:1`). PacBio samples are profiled with the
> database's `model_PB.xml`, or `--model_pb`. Reads of more than 1,000 bp given as short reads stop
> the run with a hint to `--read_type pb`. Reads without qualities (FASTA) are taken as Q30.

In `docs/running.md`'s options table: `--read_type` (default: `pe` with `-2`, else `se`) and
`--model_pb` (default `--model`, else `model_PB.xml`). `docs/database-files.md`: `model_pe.xml`,
`model_se.xml` and `model_PB.xml`, `model.xml` for older databases. `docs/model-training.md`: the
PacBio model is trained on long-read samples, as `model_se.xml` on single-end ones. The website is
out of date for single-end and PacBio reads. `docs/claude/README.md`: this report's line in the
table.
