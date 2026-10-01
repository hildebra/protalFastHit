# Gene neighbours from the downloaded genomes: what was built, and what it does on an operon-like world

- **Date**: 2026-10-01.
- **Code**: branch `audit-fixes` at `4fd375e` with this work and the alignment-feature work of another session
  (mate guidance, read consensus, `ZA` tag; [report](../2026-10-01-alignment-features/README.md)), which the
  uses build on, committed together. Binaries built in WSL (`~/protal-nbr/build`) from the working tree.
- **Data**: a synthetic release ([`run.sh`](run.sh)): 80 species (`gtdb_like_lineages.py --species 80
  --archaea 0.1 --seed 3`), 2 genomes each, `simulate_gtdb_release.py --operons`, 400 kb of other DNA per
  genome, 5 contigs, strains 0.5-3% and species 1.5-6% apart, gene rates by category. The database leaves 8
  species of genera with congeners out. Samples from 30 genomes (22 species of the database, half of them
  from their other genome, a strain; the 8 left out), lognormal abundances: 300,000 read pairs
  (`simulate_reads.py`, 150 bp, insert 350 ± 50), PacBio HiFi (pbsim3 `ERRHMM-SEQUEL`, 15 kb, 99.9%) and
  Nanopore (`QSHMM-ONT-HQ`, 8 kb, 97%), each at depth 2 on average. Presence models: those of the V3 tuning
  build, trained without gene neighbours, for both arms.
- **Machine**: WSL Ubuntu 24.04, Core Ultra 7 258V; runs with 4 threads, callgrind with 1.
- **Commands**: [`run.sh`](run.sh) (world, database, samples, 3 alternated runs with and without
  `--no_gene_neighbours`, tables by [`design.py`](design.py)), [`callgrind.sh`](callgrind.sh) (instructions
  on 40,000 read pairs, 400 PacBio and 600 Nanopore reads).

The feasibility study is [the earlier report](../2026-10-01-gene-neighbours.md); this one implements its
uses 1, 2, 3 and 5 from the genomes the build already downloads.

## What was built

| Part | |
|---|---|
| `scripts/mini_db/gene_neighbours.py` | places each species' database genes in its representative genome (exact sequence, either strand), takes each gene end's next marker within `--max_gap` (3000) or none (an end within that of its contig's end says nothing), and counts per clade (family, order, class, phylum, domain): `gene_neighbours.tsv`, ten numbers per line; prints a summary of how alike each rank's clades are |
| `build_gtdb_database.py` | runs it after the converter on the genome table it simulates from (`--no-gene-neighbours` to skip; part of the conversion's stage key); `build_metadata.tsv` records it |
| `gtdb_to_protal_db.py` | a training copy (`--from_db`) keeps the file |
| `protal --build` | checks it (format, every clade in the taxonomy, every gene in the database) and packs it into `database.protal`; `--unpack_db` writes it out |
| `src/SequenceUtils/GeneNeighbours.h` | the table: rules sorted by (clade, gene, end); each species' lineage bound at load; `Partners` and `Assess` take the nearest clade with data on a gene end: a partner seen there is expected, one never seen with the end informative in 3 or more species unlikely, else the next clade up decides |
| `src/Core/AcrossGenes.h` | where a read past a gene's end lies on the neighbour; whether two mates are one fragment across facing ends |
| uses | 1: mate guidance looks past the guiding mate's gene end on the neighbour; 2: `JoinAlignmentPairs` pairs mates on two neighbouring genes (a proper pair on two references, RNEXT the other gene, TLEN 0); 3: long reads look for the gene the clade puts next to a read's gene where no segment has it; 5: `adjacent_expected_share`, `adjacent_unlikely_share` in the dump |
| `simulate_gtdb_release.py --operons` | markers in clusters of 1-6 genes 0-150 bp apart, the same in every species, each family breaking clusters up with `--operon_breaks` (0.25); the shuffled layout before it never had a pair span two genes |
| tests | `tests/test_GeneNeighbours.cpp` (11), mini-database `GeneNeighboursTest` (3: placements equal the simulator's, family counts recomputed independently), e2e `GeneNeighboursTest` (3: build checks and packs, a bad table stops the build, pairs across genes with and without the option) |

Left out of use 3: a MAPQ penalty for segments out of order on a read. The read's consensus taxon already
settles the segments, and gene order is shared by congeners, so it would change MAPQs without telling
species apart. The adjacency features are in the dump but not among the features the models train on (below).

## Results

Calls and records with and without gene neighbours (run 1 of 3; the gene neighbour counts of the other two are the same):

| | pe with | pe without | pb with | pb without | ont with | ont without |
|---|---|---|---|---|---|---|
| species called: true / false | 22 / 0 | 22 / 0 | 22 / 2 | 22 / 2 | 22 / 0 | 22 / 0 |
| fragments paired across two genes | 5,423 | 0 | | | | |
| guided mates found on the next gene | 328 | | | | | |
| fragments written as a proper pair | 57,917 | 52,990 | | | | |
| genes found where the neighbours put them | | | 9 | | 39 | |
| gene hits without one of the read's taxon (MAPQ 0) | | | 45 | 49 | 61 | 80 |
| primary records | 131,359 | 131,158 | 4,473 | 4,469 | 5,739 | 5,725 |
| aligned bases | 18,370,369 | 18,351,335 | 4,318,247 | 4,317,611 | 5,304,552 | 5,303,326 |

- **Pairs**: one fragment in ten that is written as a proper pair now spans two genes (5,423); before,
  their mates were written as split mates of the pair's consensus taxon, not properly paired. Mate guidance
  found 328 more mates on the next gene; 0.1% more bases aligned.
- **Long reads**: few genes are missing from a read's segments that the neighbours can find (9 on 2,196
  PacBio reads, 39 on 4,129 Nanopore reads); gene hits without one of the read's taxon fell by 8% and 24%.
- **Calls**: unchanged on this world (the models were trained without gene neighbours).

The adjacency features (profile-only runs of the same SAMs; means of present / absent taxa):

| | `adjacent_expected_share` | AUC present vs absent | `adjacent_unlikely_share` | AUC |
|---|---|---|---|---|
| pe | 0.989 / 0.854 | 0.39 | 0.011 / 0.016 | 0.72 |
| pb | 0.957 / 0.888 | 0.59 | 0.043 / 0.112 | 0.41 |
| ont | 0.972 / 0.904 | 0.54 | 0.028 / 0.096 | 0.46 |

As the feasibility study expected, they hardly separate present from absent taxa: the absent taxa's reads
come from congeners left out of the database, whose genes lie in the same order. An AUC below 0.5 means the
absent taxa score higher (many have one or two linked reads, all expected). Neither is among the normalised
features; a database whose clades differ more in gene order than this world's could change that, which the
GTDB build can test.

## Cost

Instructions (callgrind, one thread; the alignment stage is the OpenMP region's outlined function, the
profiler `ProfileSam`):

| | alignment with | without | | profiling with | without | |
|---|---|---|---|---|---|---|
| pe (40,000 pairs) | 4,340,077,104 | 4,315,383,630 | +0.57% | 621,528,101 | 611,122,395 | +1.7% |
| pb (400 reads) | 10,518,317,383 | 10,482,987,800 | +0.34% | 209,866,157 | 206,812,856 | +1.5% |
| ont (600 reads) | 10,348,344,311 | 10,323,420,243 | +0.24% | 232,034,059 | 228,736,746 | +1.4% |

Wall-clock runs of the samples (~2-3 s, most of it loading the 3 GB k-mer key map) differed by their
noise (±30%). Peak memory was the same (3.2 GB). The table of this world has 6,261 lines (160 kB, 20 kB with
zstd, one line per clade, gene end and partner, about 250 per clade); a GTDB build from ~8,000 species'
genomes spans some 3,000 clades, so roughly 1 M lines, ~35 MB in memory: an estimate, not measured.

## Found on the way

- Long reads' neighbour search first re-aligned genes that the read had seeded and that had failed `-a`:
  148 alignments for 4 hits on 400 PacBio reads, +5.6% alignment instructions. It now aligns a seeded gene's
  candidate only if none was aligned, and places only unseeded genes by their 12-mers: +0.3%. A prefilter by
  shared k-mers instead gained little and lost hits.
- The adjacency of a long read's genes was first judged on the taxon's records after the profiler's filters
  (a gene in between with MAPQ 0 left its neighbours "unlikely"), then on read positions the SAM reader had
  lost: it drops hard clips from the CIGAR. The profiler now judges every best record of a read, of any taxon,
  before the filters, and the reader keeps the hard clips (`SamEntry::m_hard_clip_start`, `_end`); a Python
  replica of the rule gives the same shares.
- Consecutive genes of a long read further apart than the table's `max_gap` are not judged.

## Next

- The GTDB build: `gene_neighbours.log` gives the share of genes placed (genomes from NCBI must be the
  accessions GTDB called its genes on) and how alike families are; models trained on a database with gene
  neighbours see the pairs and rescues the runs make.
- Test the adjacency features on that build's training data before adding them to the normalised set.
- The website does not describe `gene_neighbours.tsv`, `--no_gene_neighbours` or pairs on two references.
