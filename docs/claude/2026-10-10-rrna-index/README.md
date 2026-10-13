# rRNA genes next to the marker genes: does GTDB have them, and could protal place reads of lineages it lacks?

2026-10-10, `audit-fixes` at `0778be5` (no protal code was changed or run; the scripts use only
`scripts/download_gtdb.py`'s genome categories and `scripts/mini_db/markers_r226.tsv`). Data: GTDB r226 (226.0),
`bac120_metadata_r226.tsv.gz`, `ar53_metadata_r226.tsv.gz` and `genomic_files_all/ssu_all_r226.fna.gz`, downloaded
from data.gtdb.ecogenomic.org and checked against the release's MD5SUM.txt, kept in WSL `~/GTDB/r226/`. Machine: the
WSL laptop, 1 core for the coverage script (18 s), 2 for the held-out test (9 min).

The questions: (1) does GTDB hold enough to give every species a 16S and a 23S rRNA gene; (2) the larger one: could
protal store an rRNA index next to the species marker genes, to detect lineages the marker genes do not cover at
species level (the reads behind the profile's `?` line)?

## Short answer

1. **Per species, no.** 59% of the 143,614 r226 species have a 16S hit on some genome, 40% a full-length one; the
   23S numbers are the same. GTDB distributes the 16S sequences (`ssu_all`, every hit of every genome) but **no 23S
   sequence at all**: the metadata gives each genome's 23S count, length and contig, so the 23S has to be cut out of
   the genomes (representatives: the 127 GB `gtdb_genomes_reps`; other genomes from NCBI). The gap is the MAGs: 76% of
   species have a MAG as representative and 19% of those have a full 16S (isolates 91%).
2. **Per lineage, mostly yes**, and that is what detecting a missing lineage needs: were a species missing from the
   database, 76% of species still have a congener with a full 16S, 96% a family member (section 3).
3. **Reads of a missing species or genus would align and land in the right genus or family; a missing family's
   reads would not** (held-out test on 600 species, section 4): 92-98% of their 150 bp 16S reads stay above protal's
   0.9 identity floor, and a missing genus's reads point to its family in a median 54% against 4% elsewhere; a
   missing family's reads still align (83%) but go to a wrong order about as often as to the right one.
4. **It does not fit protal's marker-gene index as "two more genes"** (section 6): rRNA is multi-copy, identical
   across congeners, has no unique k-mers, and its conserved seed cores exceed the index's 2048-value cap. It needs
   its own small index (< 1 GB), an LCA classifier (protal has none yet) and copy-number-corrected abundances, about
   2-3 weeks of work, best preceded by an offline check on real samples with a large `?` share (section 7).

## 1. What GTDB r226 holds

| file | content |
|---|---|
| `genomic_files_all/ssu_all_r226.fna.gz` (70 MB) | 1,053,664 16S hits (nhmmer, Rfam RF00177/RF01959, >= 200 bp, E <= 1e-6) of 478,089 genomes, all hits per genome, with GTDB taxonomy; agrees exactly with the metadata's `ssu_count` and `ssu_length` |
| `genomic_files_reps/{bac120,ar53}_ssu_reps_r226.fna.gz` (24 MB) | the longest 16S of each representative |
| metadata `ssu_*` | count, longest length, its contig and contig length, best SILVA and Greengenes hit |
| metadata `lsu_23s_*`, `lsu_5s_*` | count, longest length, its contig, best SILVA 23S hit; **no sequence file** |
| `auxillary_files/hq_mimag_genomes_r226.tsv` | genomes with 5S, 16S (>= 1200/900 bp) and 23S (>= 1900 bp) |

GTDB warns that the 16S sequences are not quality-checked and that a contaminating 16S can sit in a MAG.

## 2. Per species (`coverage_output.txt`, sections 2-5)

Full length = MIMAG high-quality lengths: 16S >= 1200 bp (bacteria) / 900 bp (archaea), 23S >= 1900 bp.

| species | 16S any | 16S full | 23S any | 23S full | both full | both full on one contig |
|---|---|---|---|---|---|---|
| representative | 53.8% | 36.6% | 53.8% | 35.3% | 31.7% | 23.4% |
| any genome of the species | 58.9% | 39.8% | 58.8% | 39.2% | 34.9% | 27.3% |

- Other genomes of the species add little (+3 points): 4,565 species gain a full 16S, 5,609 a full 23S.
- By the representative's category, full 16S (representative / any genome): isolate 91.0% / 93.7%, MAG 19.2% /
  22.5%, SAG 79.6% / 82.3%. rRNA operons are repeats; short-read assemblies collapse them into short contigs that
  binners leave out. In isolates, 19% of the 16S sit on a contig < 2 kb (the collapsed operon), in MAGs 1.6%.
- Species with 10+ genomes have a full 16S in 80%, singletons (94,012 species) in 32%.
- Phyla range from 15% (Gemmatimonadota, 99% MAGs) to 59% (Halobacteriota) for the representative's full 16S;
  Bacillota, the largest gut phylum, 24.5% / 29.0%.
- Where a hit exists, the median length is full (16S 1,490 bp, 23S 2,797 bp), but a quarter are fragments (16S p25
  843 bp).
- Contamination signals are small: the SILVA best hit's domain differs from GTDB's for 0.6% of representatives'
  16S and 0.5-0.8% of their 23S; 0.3-0.4% hit an organelle sequence. Genus-level misbinning would not show here.

## 3. Per lineage: what a read of a missing species would meet (section 2b)

| rank | clades | with a full 16S | full 23S | both | species whose clade has a full 16S in another species | 23S | both |
|---|---|---|---|---|---|---|---|
| genus | 29,405 | 51.9% | 51.9% | 45.4% | 75.8% | 75.6% | 72.4% |
| family | 5,932 | 66.3% | 66.5% | 60.5% | 96.0% | 95.9% | 95.2% |
| order | 2,164 | 72.9% | 74.0% | 68.6% | 98.8% | 98.8% | 98.6% |
| class | 638 | 76.6% | 77.6% | 73.5% | 99.7% | 99.7% | 99.6% |

Half the genera have no full 16S, but they are small: three quarters of the species sit in a genus where another
species has one, and nearly every species in a family that has one. Lineages GTDB has no genome for at all are
outside this count (section 7).

## 4. Held-out test: would 16S reads of a missing lineage be placed? (`signal_output.txt`)

Pool: one full-length 16S for each of the 57,085 species that have one (52,511 the representative's own). For 600
random pool species, the species', genus', family's and order's 16S were held out in turn; the held-out 16S was
compared with the rest, full length and as error-free 150 bp reads every 50 bp (~27 per 16S). Identity is edit-based
(Myers' bit-parallel algorithm); a read's candidates are the sequences sharing the most 12-mers with it, in its parent
clade, the two clades above, and overall (`rrna_signal.py` docstring). No aligner is installed in WSL; this is the
script's own.

**Nearest relative, full length** (median, p10-p90): a missing species is 98.9% (93.9-99.8) identical to a congener;
a missing genus 94.8% (88.4-98.0) to its family; a missing family 90.0% (84.3-94.2) to its order; a missing order
89.7% to its class. The classic thresholds hold on GTDB's taxonomy as medians but with wide spread: 54% of held-out
genera are >= 94.5% to some genus of their family, i.e. look like members of a known genus by 16S.

**Reads** (identity >= 0.90 is protal's alignment floor, `--max_score_ani`):

| held out | reads >= 0.97 | >= 0.90 (alignable) | right to the parent (A + B) | of which B: all best hits in one known sibling | right only higher (C) | misplaced (D) |
|---|---|---|---|---|---|---|
| species | 87% | 98% | 57% (genus) | 18% | 37% | 6% |
| genus | 59% | 92% | 55% (family) | 25% | 32% | 13% |
| family | 42% | 83% | 27% (order) | 19% | 40% | 33% |
| order | 38% | 79% | 33% (class) | 23% | 36% | 32% |

Placement columns are shares of the alignable reads, LCA of the hits within 1 edit of the best ("loose"; with the
best hits only, misplaced reads rise to 10/21/47/44%).

**Per held-out clade** (what a profile would aggregate): a missing genus has a median 54% of its reads placed in the
right family (p25-p75 25-78%) against 4% misplaced, 89% of missing genera get >= 3 such reads per 16S, and 18% have
more misplaced than right reads. A missing family: median 17% right, 21% misplaced, and 55% of families have more
misplaced than right reads.

So:

- **Unknown species of a known genus and unknown genera of a known family** are recoverable: their reads align (92-98%)
  and a clade's reads point to the right genus or family far more often than not.
- **Unknown families and orders** align (79-83%) but are placed no better than chance between the right order and a
  wrong one: report them at class or phylum level, or not at all. GTDB's genome-based taxonomy is not monophyletic in
  the 16S at these depths, and contaminating 16S in MAG representatives (the pool was not filtered) add misplacements.
- **B is the subtle case**: a quarter of a missing genus's reads have all their best hits in one known sibling genus
  (3% of its alignable reads identical to it; 59% of all its reads are >= 0.97 to their best hit, above where a
  per-read genus cap at 94.5% would act). In protal this is resolvable at the profile level: the
  marker genes say which species of that sibling are present; rRNA reads on a sibling genus with no called species
  are an unknown member of the family (or an uncalled species of the sibling).
- Not modelled: sequencing errors (they lower all identities by the error rate), intragenomic copy variants, the 23S
  (longer and more variable: more reads and somewhat more resolution), species without a 16S in the pool (60%).

## 5. Read yield and abundance (section 6)

- A median complete bacterial genome (9,888 species with a complete isolate genome) has 4 16S and 4 23S copies
  (p10-p90 1-8): 18.2 kb of rRNA, 0.41% of the genome (p10-p90 0.17-1.03%); archaea 1 copy, 0.26%. The bac120
  marker genes are ~123 kb, 3.0% of a median 4.15 Mb genome: rRNA gives about **7 times fewer reads** than the
  markers (30 times fewer for a single-copy genome). At 5 M fragments, a lineage at 0.1% of the cells (~5,000
  fragments) gives ~20 rRNA fragments at 4 copies, ~5 at one copy: detectable down to roughly 0.05-0.1% per lineage;
  the markers reach far lower for the species they know.
- Abundance needs the copy number. Only 12% of families (9% of genera) have a complete isolate genome to count it
  (they hold 65% of the species); MAG counts are useless (median 1, operons collapsed). Within a family the 10th and
  90th percentile species differ 1.5-fold (median family; 2.8-fold at the 90th percentile). A family-median copy
  number gives an abundance good to a factor 1.5-3; rrnDB-style phylogenetic prediction would be the refinement.
- The `?` line already measures the reads the profile leaves unexplained (`Composition.h`); rRNA reads would split
  that share by lineage, and their total is a check on it.

## 6. What it would take in protal

From the index and profiler as they are (`src/`, read for this report):

- **Seeds**: 31-base windows around a 15-base closed-syncmer core; a core with >= 2048 values is dropped at build
  (`Hash/Seedmap.h` `max_key_multiplicity`), and a lookup gives no seed when > 256 cells tie (`--max_key_ubiquity`).
  57,085 species' 16S share their conserved cores far more often than that, so in the marker index the conserved
  half of the gene would carry no seeds; that is acceptable (those reads only place at phylum level) only if the
  variable regions still seed, which a dedicated index with deduplicated references makes likely.
- **One sequence per (taxon, gene)**, a duplicate is an error (`GenomeLoader.h`); identical copies in several
  species are separate records, and a gene without unique k-mers is "not hittable" and skipped by the profiler. The
  16S is identical between many congeners: in the marker index it would be ignored or would add ambiguous anchors
  (`--align_top` 3 plus all ties, unbounded).
- **No LCA anywhere**: `IntTaxonomy::LCA` exists but is not called; reads are assigned to species only.

So, as its own component:

1. **Build** (GTDB script): 16S from `ssu_all` (one or a few per species, full length, filtered by agreement with
   its genus), 23S cut out of the representatives (and of one genome per species whose representative lacks it)
   with nhmmer/barrnap on the HPC (~60k genomes); identical sequences merged, each labelled with the LCA of its
   species; a family-level copy-number table from complete genomes. `rrna.fna` and its seed index in
   `database.protal`: <= 0.25 Gbp of sequence before merging (57k x 1.5 kb 16S + 2.9 kb 23S), ~50 M seed entries.
   Not a second `Seedmap` as it is: its key map is a dense 4^15 table, ~3 GB whatever it holds; a hashed or sorted
   key list for this small set keeps it well under 1 GB.
2. **Reads**: rRNA reads are ~0.4% of reads; seed them against the rRNA index (cheap, small), align the candidates
   with the same WFA2 path, keep all hits within one or two edits of the best, and assign the LCA, capped by identity
   (the classic 16S thresholds: species 98.7%, genus 94.5%, family 86.5%, order 82.0%; Yarza et al. 2014).
3. **Profile**: a read whose best hits include a called species' rRNA is that species'; the rest, aggregated per
   clade, become "unknown species of genus G" / "unknown genus of family F" lines that split the `?` share, with
   abundances from rRNA bases / (copy number x rRNA length), and a minimum of reads from several windows.
4. **Decoys**: chloroplast and mitochondrial 16S, eukaryotic 18S/28S and the host's rRNA (human 45S unit, mt 12S/16S):
   conserved regions of these align to bacterial rRNA, and plant, food or host-rich samples have many such reads.

## 7. Alternatives and open points

- **Lineages without any genome** (clades known only from 16S surveys) are not in GTDB at all; SILVA SSU/LSU Ref NR99
  or Greengenes2 (16S placed on a GTDB-based genome tree) would add them, with taxonomy that maps to GTDB only down to
  their nearest genome-bearing clade.
- **Protein-level marker search** (SingleM's approach: translated windows of conserved single-copy markers): 7 times
  more sequence per genome than rRNA, single copy (no copy-number correction), and amino acids stay alignable across
  genera and families where nucleotides (and protal's 0.9 identity floor) do not; but it needs a translated seed and
  aligner, a larger change than an rRNA index.
- **Relaxing protal's identity floor on the marker genes** is unlikely to reach new genera: genomes of different
  genera are typically below ~80% ANI (not measured here for the bac120 genes, whose conserved ribosomal proteins sit
  somewhat higher), and at 80% identity a 15-base exact seed core survives at 0.8^15 = 3.5% of positions. Measuring
  the marker genes' nucleotide identity between genera (from `bac120_marker_genes_reps`, 6.7 GB) would settle it.

**Suggested next step** (before any C++): run `rrna_signal.py`'s classifier, extended to raw reads, on two or three
real samples with a large `?` share (soil, an undersampled gut): screen the reads for rRNA, place them by the LCA
against the 57,085-species 16S pool, and see whether the "unknown" lineages they name add up to the `?` share. That
answers whether the rRNA layer explains the unknown fraction before the index is built.

## Files

- `rrna_coverage.py` -> `coverage_output.txt`: `python3 rrna_coverage.py ~/GTDB/r226` (18 s, 2.6 GB).
- `rrna_signal.py` -> `signal_output.txt`: `python3 rrna_signal.py ~/GTDB/r226 --sample 600 --workers 2` (seed 1,
  9 min on 2 cores, 2.6 GB).
- Both run with `~/micromamba/envs/protal-db-build/bin/python3` (pandas, numpy) in WSL; inputs in `~/GTDB/r226`.
