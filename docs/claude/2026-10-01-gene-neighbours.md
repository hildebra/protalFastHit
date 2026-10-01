# Gene neighbours in the database: how hard it is to record them, and what a run could do with them

Implemented the same day, from the genomes already downloaded, with uses 1, 2, 3 and 5:
[the implementation and its experiment](2026-10-01-gene-neighbours-run/README.md).

- **Date**: 2026-10-01.
- **Code**: branch `audit-fixes` at `c575aa5`, with the working tree's uncommitted mate guidance
  (`src/Core/MateGuidance.h`) of another session. Nothing was implemented for this report.
- **Data**: nothing measured. The file listing of GTDB r226 (`data.gtdb.ecogenomic.org`,
  `releases/release226/226.0/genomic_files_reps/`) and the synthetic release's
  `simulation/marker_positions.tsv` (`~/audit2/gtdb_r226` in WSL). The gene layouts below are from
  the literature on *E. coli* and other bacteria, recalled, not measured; experiments E1 and E2 at
  the end measure them.

## The question

In simulations a read pair can have its mates on two genes of its species that are neighbours in the
genome, and a PacBio or Nanopore read spans several. Should `--build` record, for each gene, which
marker gene is its neighbour and how far away it is (when within ~3 kb)? Stored in the database
first and not used; how could a run use it, and does it fit in memory?

## What a read across two genes does today

- **Pairing**: `JoinAlignmentPairs` (`src/Classify.h:369`) pairs a candidate of mate 1 with one of
  mate 2 only on the same taxon and gene, in opposite orientations. Mates on two genes become two
  single-mate candidates (`ScorePairedAlignment`: own score − 2, against (sum) / 2.2 for a pair).
- **Output**: `WriteSplitMates` (`src/IO/AlignmentOutputHandler.h:619`) writes such mates at their
  best alignments of the pair's consensus taxon (`SettleMatesByPair`), paired but not properly
  paired, TLEN 0.
- **Mate guidance** (uncommitted): a sure mate guides the other to its taxon from the other's own
  anchor of the taxon (this catches a mate on a neighbouring gene when it has an anchor), else by a
  k-mer diagonal search on the guiding mate's gene only, clipped at the gene's end
  (`kRescueMaxFragment` 1000).
- **Long reads** (`src/Core/LongReads.h`): each gene a read hits is a segment, aligned on its own;
  the read's consensus taxon settles the segments. The order of the segments on the read and the
  distances between them are not used.
- **Profiler**: `linked_share` counts the reads with two or more records on the taxon (both mates,
  or two genes of a long read), whether or not the genes are neighbours. A mate-concordance feature
  gained +0.0008 F1 for paired-end reads
  ([alignment features](2026-10-01-alignment-features/README.md)).

## How often marker genes are neighbours

Much of bac120 lies in a few conserved clusters. In *E. coli*, roughly:

- the S10-spc-alpha ribosomal protein operons: L3, L4, L2, L22, S3, L16, L24, S8, L6, S5, L15, SecY,
  S11, S4, RpoA, L17, 16 markers in about 12 kb, with 0 to ~900 bp between consecutive markers
  (often a non-marker ribosomal protein between them, sometimes overlapping ends);
- SecE, NusG, L11, L1, L10, RpoB, RpoC (7 markers, ~11 kb);
- RimP, NusA, IF-2, RbfA, TruB, PNPase; RpsB, Tsf, PyrH, Frr; L13-S9; PheS-PheT; RuvA-RuvB;
  AtpG-AtpD; Rnc-Era; RimM-TrmD; IF3-L20; L21-Obg; DnaA-DnaN-(RecF)-GyrB; MraY-MurD-MurC; S6-L9;
  DnaX-RecR; GrpE-DnaK in Firmicutes.

That is about half of the 120 markers with another marker within 3 kb. Archaea have the same
ribosomal super-operon. Two consequences:

1. **It matters most on the shortest genes.** 26 bac120 markers and 18 ar53 markers have HMM models
   of at most 150 aa (`scripts/mini_db/markers_r226.tsv`; the genes are somewhat longer, and a
   domain model such as PNPase's much shorter than its gene), most of them ribosomal proteins of
   ~300-450 bp. That is shorter than or as long as a typical 300-500 bp fragment, so most fragments
   that touch such a gene have a mate mostly outside it: on the neighbour or on DNA between the genes.
2. **Gene order is conserved across bacteria**, so a neighbour table helps align and pair reads; it
   does not tell congeners apart (a missing congener's reads pair across the same genes). It is not a
   lever against the false positives of [the false-positive report](2026-10-01-false-positives-v2-test/README.md).

The synthetic releases (`simulate_gtdb_release.py`) shuffle the markers of each species and space them
evenly, ~1.2 kb apart at the default 150 kb genome: pairs there never span two genes (long reads do),
and the order is not conserved between species. Testing any use needs an option for operon-like
clusters (conserved order, 0-300 bp gaps).

## Is it easy?

### Storing it: yes

The same path as `gene_conservation.tsv`: a file name constant, an entry in the bundle's sources
(`src/Build.h:256`), a loader and table in `GenomeLoader`, `--unpack_db` writing it out, checks at
load, unit tests, the mini database. A database without the file loads as before and the feature is
off. About a day, and with nothing using it, no output changes.

### Getting the coordinates: this is the work

The build's inputs carry no coordinates. The real r226 marker files hold the accession and nothing
else in each header (checked by the user on `TIGR00580.fna`: `>RS_GCF_001027105.1`, ...), and the
records are in genome order, not genomic order, so nothing about position can be read from them.
A second file with coordinates is needed; each marker is found in it by its exact sequence:

| Source | Covers | Cost | Notes |
|---|---|---|---|
| GTDB's Prodigal gene calls of the representatives, `gtdb_proteins_aa_reps_r226.tar.gz` | every species | 87.5 GB download, one streamed pass (hours on a node) | Prodigal's headers normally hold `# start # end # strand`; check GTDB kept them. A marker is one of the genome's Prodigal genes: match its protein (the marker tarball's `faa/`) by exact sequence. Gives contig and coordinates, not the DNA between genes. A contig's length is not in the file; the last gene on it bounds it |
| Whole genomes the build already downloads to simulate from (`download_gtdb.py`: ~8,000 species' representatives and ~12,000 strains) | ~6% of ~140k species | none extra: an exact search of each genome's marker sequences (GTDB's are that genome's own gene calls), both strands | Gives contig lengths, the DNA between genes and the strains' variation of the distances |
| `gtdb_genomes_reps_r226.tar.gz` (`--rep_genomes gtdb`) | every species | 127.3 GB | As above for all species |
| Clade-level table from the genomes above | every species, approximately | kilobytes | Per gene pair and phylum or class: share of genomes where the two are neighbours, median distance, spread. A species without its own coordinates takes its nearest clade's |

Mapping coordinates to protal's genes: the converter keeps a genome's first copy of a marker and
numbers the genes, so the coordinates must be those of that copy (an exact match finds it) and go
through the same accession → taxid and marker → gene id maps. `--exclude_species` and `--from_db`
must carry the file; `--build_gene_subset` and `scripts/subset_genes.py` change which genes are
hittable, so neighbours should be derived among the indexed genes at load, not frozen at build.

## What to store: positions, not only the nearest neighbour

Store each reference gene's place in its representative genome, and derive neighbours when the
database loads:

- `gene_positions.tsv`: `taxid  geneid  contig  start  end  strand`, plus each contig's length
  (or whether its ends are known), contig numbered within the species.
- The producers (Prodigal headers, genome search, the simulator's `marker_positions.tsv`, which
  already has `accession marker contig start end strand`) write one interchange format; the converter
  turns it into `gene_positions.tsv`.

Why positions rather than "nearest neighbour within 3 kb": paired-end reads need neighbours within
the longest fragment (~1 kb), long reads any two genes on one contig within a read's length (10-50
kb); the cutoff and the gene subset can then change without a rebuild; and the distance between any
two genes, not only adjacent ones, is a subtraction. Each gene end must also be one of three states
(neighbour, no marker within the cutoff, contig ends within the cutoff): many representatives are
fragmented MAGs, and a contig end is no evidence that a gene has no neighbour.

"Average distance" only arises with more than one genome per species (strains) or per clade: for the
representative, whose genes protal holds, one distance per neighbour is the right value. The strains'
spread (from the downloaded genomes) can widen the tolerance later.

## In memory

Small, and off the hot path:

| Layout | Bytes | At GTDB size (~16.6M genes) |
|---|---|---|
| Positions: start (u32), contig (u16), strand and end states (u16), in a vector beside each `Genome`'s genes | 8 per gene | ~130-200 MB (dense up to the highest gene id) |
| Derived neighbours: per gene end, neighbour gene (u16), distance (i16, overlaps negative), strand relation and state | 2 × 4 per gene | another ~130-200 MB, or computed from a per-species order instead |
| Clade-level table only | a few hundred bytes per clade | under 1 MB |

A run's memory at GTDB size is tens of GB (59 GB before the 2-bit gene store, extrapolated in
[memory profiling](2026-09-30-memory-profiling/README.md); about 40-45 GB with it), so this is under 1%. Not inside `Gene`: it is 32 B with no spare bytes; its counts and length leave
spare bits (lengths and counts fit 20 bits), but packing neighbours there ties two unrelated things
together for ~130 MB. Loading: a TSV of ~16.6M lines parses in seconds on `-t` threads, as
`unique_kmers.tsv` does; ~100 MB compressed in `database.protal`.

Lookups come only on rare events: a fragment whose mates are on two genes of one taxon, a guided
mate near a gene's end, the segments of a long read. The code already holds the `Genome`
(`genomes.GetGenome(taxid)`), so a lookup is an index into its vector: cache behaviour and table
layout do not matter.

## Uses in a run, by value against risk

1. **Mate rescue across the gene's end** (alignment, low risk). In mate guidance step 2, when the
   reach of the fragment (`kRescueMaxFragment` from the guiding mate) runs past the gene's end, go
   on into the neighbour at that end: past the gap, in the orientation the strand relation gives,
   with the window widened by the distance's uncertainty. Only fragments guidance already handles
   pay for it, and it only adds alignments where today there are none. Measured by a new
   `MateGuidanceCounts` field (rescued on the neighbour) and the bases the strains' SNPs gain on the
   short ribosomal genes.
2. **Pairs across neighbouring genes** (alignment, the largest change). In `JoinAlignmentPairs`,
   also pair a mate on gene A with a mate on gene B of the same taxon when B is A's neighbour at
   that end, the orientations fit the strand relation, and the implied fragment
   (rest of A + gap + start on B) is within the longest fragment. Such a fragment becomes one
   candidate with the pair's score, so the taxon that explains both mates ranks above one that
   explains one; it is written properly paired, with a fragment length in genome coordinates, and
   `WriteSplitMates` is left for discordant splits. The check costs a lookup only on same-taxon,
   different-gene combinations. Risk: it changes candidate ranking, MAPQ and the records the
   profiler sees (`low_mapq_share`, `linked_share`, MAPQ features), so all presence models need
   retraining and the independent test sets need to confirm it.
3. **Long reads: expected genes and colinearity** (alignment, moderate risk). (a) From a segment
   on gene A and the taxon's positions, predict where on the read genes B, C, ... lie; a predicted
   gene without a segment (no seed: too divergent) gets a targeted alignment there, as mate rescue
   does. (b) In `SettleByRead`, score the segments' order and spacing on the read against the taxon's
   positions (tolerance for Nanopore indels); a segment that does not fit (chimera, wrong taxon) gets a
   lower MAPQ. Distances between strains differ by tens to hundreds of bases between genes, so this
   checks structure, it does not separate species.
4. **Strain phasing across genes** (profiler, larger work). Pairs and long reads that link two
   neighbouring genes link their SNPs, so `Strain.h` could phase haplotypes across genes rather
   than within one. Most useful for long reads.
5. **A presence-model feature** (profiler, low value). The share of a taxon's cross-gene links that fit
   its gene positions. Cheap, but congeners fit too, and the mate-concordance feature it would
   refine gained almost nothing.

## Risks to settle before implementing

- **Models**: uses 2 and 3 change what the profiler sees; retrain and test all read types. Use 1
  adds records only, and use 1 alone may already recover most of what matters.
- **Contig ends**: on fragmented representatives, an unknown end must not count as "no neighbour".
- **Strain variation of the gaps**: the representative's distances, with a tolerance; measured
  from strain genomes (E2).
- **Which neighbour**: derived among hittable genes at load, so reduced marker sets work.
- **Synthetic worlds**: without operon-like clusters, the mini database tests the plumbing only.

## Suggested order

- **E1, does it matter** (no code): in the paired-end samples simulated from real genomes for
  training, count per gene the fragments written by `WriteSplitMates` (records whose RNEXT is
  another gene of the taxon) and the mate guidance counts, and the gene pairs they join. This says
  which genes lose how much and gives the empirical adjacency of the read data.
- **E2, how conserved** (a Python script): place the markers in the downloaded genomes; per gene
  pair, how often they are neighbours across genera, families and phyla, and how much the distances
  vary within species. If adjacency and distance are conserved within a family, the clade-level
  table suffices and the 87.5 GB download is not needed.
- Then the interchange format, `gene_positions.tsv` and the in-memory table (no behaviour change),
  use 1, and use 2 or 3 if E1 shows enough fragments to gain.
