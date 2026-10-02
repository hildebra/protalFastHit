# Gene neighbour frequencies from every genome, in the one database file

- **Date**: 2026-10-02.
- **Code**: branch `audit-fixes` at `fd1e719` (0.7.2) with this work uncommitted (the files listed below).
  Binaries built in WSL (`~/bprof/gnb/build`, Release, from an rsync copy of the working tree).
- **Data**: an operon world (`~/opw`, [`make_operon_world.sh`](scripts/make_operon_world.sh), the tuning
  world's design of [2026-09-29](../2026-09-29-model-training-tuning/README.md) with `--operons`): 900 species of
  GTDB-like lineages (8% archaea), 3 genomes each (150 kb of other DNA, 3 contigs, strains 0.2-2% and species
  1.5-6% apart), markers in operon-like clusters of 1-6 genes that each family breaks up with probability 0.25 and
  whose order on the chromosome each species shuffles; the release holds 765 of the species (2,295 genomes), 135
  are in no database.
  Builds by `build_gtdb_database.py` with gene neighbours from all 2,295 genomes and 238 species held out of the
  training database ([`opw_build.sh`](scripts/opw_build.sh)).
- **Machine**: WSL Ubuntu 24.04, Core Ultra 7 258V (6 cores in WSL), 23 GB.
- **Commands**: the scripts in [`scripts/`](scripts/): [`opw_build.sh`](scripts/opw_build.sh) (build and train),
  [`check_placements.py`](scripts/check_placements.py) (placements against the simulator's),
  [`rule_check.py`](scripts/rule_check.py) (how rules for sparse clades judge the genomes' own gene order),
  [`deciding_rank.py`](scripts/deciding_rank.py), [`feature_auc.py`](scripts/feature_auc.py),
  [`compare_features.sh`](scripts/compare_features.sh) (the models retrained with and without the adjacency
  features), [`alignment_compare.sh`](scripts/alignment_compare.sh) (test samples with and without
  `--no_gene_neighbours`).

The requests, in order: a per-family table of how often each gene lies next to which, for all bacterial and all
archaeal genes; can it guide read pairs and long reads; adjacency as evidence for detecting species, perhaps in the
random forest; the table in `database.protal`, with no other file needed; adjacency from k-mer traces of all
downloaded genomes, not only the representatives, the per-genome positions kept in the database and only the
per-clade frequencies loaded; complete genomes are circular; species representatives in one sequence are circular;
a clade with too few species falls back to higher ranks.

## What changed

| Part | |
|---|---|
| `scripts/mini_db/gene_neighbours.py` | places each species' database genes in **every** genome of the species at hand (the representative by exact sequence; other strains by their k-mer trace: the gene's 24-mers against the genome's 24-mers at every 16th position, the best band of diagonals, 10% of the k-mers and 3 hits at least); a sequence is **circular** if its header says it is a whole replicon (`complete genome`, `complete sequence`, `topology=circular`) or it is a species representative in one sequence (its last gene faces its first across the origin, and every end on it counts); a species counts once per gene end, with the partner most of its genomes show (the representative's on a tie); writes the per-clade **frequencies** (`gene_neighbours.tsv`: the species with the pairing, of those in which that end is informative) and the **positions** (`gene_positions.tsv`: accession, taxid, contig, its length, circular, gene, start, end, strand, exact or trace, k-mer share); `--from_positions` derives the frequencies from the positions alone (byte-identical) |
| `gtdb_to_protal_db.py --from_db` | the training database's copy derives its frequencies anew from the positions, without the genomes of the species it leaves out, and keeps their positions out |
| `protal --build` | checks `gene_positions.tsv` (eleven fields, ranges, every gene one of its species in the database; exit 8 at the first problem) and packs it into `database.protal` beside the frequencies; prints `Gene positions: N genes (T placed by their k-mer trace) in G genomes of S species, C read as circular` |
| a run | loads only `gene_neighbours.tsv`; the positions stay in the file, unread. A pairing's share is smoothed over the species' clades, from the top down (below); expected from 0.2, unlikely up to 0.05, rare between; mate pairing across genes, mate rescue on the neighbour and the long-read neighbour search follow expected partners only, the most common first |
| profile | `adjacent_support`: the mean smoothed share of the pairings on a taxon's reads across genes, (sum + 0.5) / (judged + 1); 0.5 without |
| trainer | `random_forest_cmdline.py --features normalized+adjacency`: the normalised features and the three adjacency features (`model_features.ADJACENCY_FEATURES`); the default stays `normalized` |
| `build_gtdb_database.py` | step line `the genes placed in N genomes (X exactly, Y by their k-mer trace) and their neighbours counted in ...`; `build_metadata.tsv` row `gene_positions` |
| tests | unit (`test_GeneNeighbours.cpp`: smoothed shares, sparse lineages, partners of all clades, `adjacent_support`); mini-db `GeneNeighboursTest` (strains placed within 15 bases of the simulator's positions, family lines recomputed independently from every genome with the species majority, `--from_positions` identical, the training copy's table without a species) and `CircularGeneNeighboursTest` (a single-sequence release: representatives circular, the last gene facing the first, strains not; headers; `neighbour_ends` on circular and linear contigs); e2e `GeneNeighboursTest` (positions packed and unpacked unchanged, a bad positions line stops the build, a run from `database.protal` alone in an empty folder profiles as the full folder does) |
| docs | `building-a-database.md` (the gene neighbours section rewritten), `database-files.md` (`gene_positions.tsv`), `running.md`, `model-training.md`, `testing.md` |

## Placing genes by their k-mer trace

On the operon world's 2,295 genomes (`gene_neighbours.py -t 5`: 18-23 s wall in the builds; 13.6 s wall and 49 CPU-s with `-t 4` alone),
against the simulator's positions ([`check_placements.py`](scripts/check_placements.py), start within 30 bases):

| genomes | placed | at its place | near (start within 30 bases, an indel shifts its end) | elsewhere | missed |
|---|---|---|---|---|---|
| representatives (765) | 85,901 exactly | 85,901 | 0 | 0 | 0 |
| other strains (1,530) | 335 exactly, 168,006 by trace | 115,621 | 52,720 | 0 | 3,682 |

The 3,682 missed are genes the representative lost (the simulator's 2% marker loss): the database has no gene of
that species to place. The first version placed 258 strain genes cut at a contig's start at the end of the
previous contig (the trace's diagonal lies before the contig); the contig is now that of the gene's middle.

With all genomes, an end the representative's contigs cut off is informative in another strain's assembly:
170,718 gene ends informative in the species, 90.2% of them with a marker within 3 kb (median gap 125 bases). At
family level the most common partner of a gene end is that of 70.6% of the family's species, 64.3% at order and
60.6% at domain level.

## One file

`database.protal` holds both tables; a run needs nothing else. The e2e test copies the bundle alone into an empty
folder, runs the paired-end sample from it and gets the same profile (same rows, values equal to 1e-9) and the
same number of fragments paired across genes as from the build folder, with nothing written beside it.
`--unpack_db` writes both tables out unchanged. Sizes on this world: `gene_neighbours.tsv` 258,590 lines (8.6 MB,
0.9 MB with zstd -19), `gene_positions.tsv` 254,242 lines (20.8 MB, 2.5 MB with zstd -19) in a 110 MB
`database.protal`; about 110 lines per genome, so some 3 M lines (~30-50 MB compressed) for 30,000 GTDB genomes.
The helpers `trace_relatives.py` and `db_gene_counts.py` still read `genome2tiid.tsv` and `gene_congeners.tsv`
beside the bundle; they are build and analysis scripts, not protal.

## Sparse clades

The rule first asked for, a clade with too few species left out and the next clade with 5 informative species
deciding, moves the decision for 17-19% of the species' gene ends from the family to a higher rank on this world
([`deciding_rank.py`](scripts/deciding_rank.py): family 80.9% for bacteria, 83.2% for archaea; order 13.4% and
12.0%). There it dilutes the family's own gene order: an arrangement that only a small family has is rare or
unlikely by its order's share, though every species of the family shows it. [`rule_check.py`](scripts/rule_check.py) judges the
genomes' own gene order with each rule, on the training database's table: true pairs (genes next to each other
within 3 kb), skip pairs (a gene and the one after its neighbour, as a read whose middle gene did not align) and
far pairs (genes 10 kb or more apart). Rules: *nearest*, the nearest clade with data (with fewer than 3
informative species, a pairing seen is expected and one not seen goes up; the rule before this work); *fallback*,
the nearest clade with 5 or more decides; *smoothed M*, each clade's share pulled towards its parent's by M
pseudo-species from the top down, (species + M x parent share) / (informative + M).

In sample (600 genomes of the training database's species), species whose family has the end informative in fewer
than 5 species:

| rule | true pairs expected | true pairs unlikely | skip pairs unlikely | far pairs unlikely |
|---|---|---|---|---|
| nearest | 98.4% | 1.5% | 91.9% | 98.6% |
| fallback | 54.2% | 25.0% | 91.6% | 98.6% |
| smoothed 5 | 85.9% | 1.5% | 93.5% | 98.9% |

All species (the strict fallback differs only where families are small):

| rule | true expected | true unlikely | skip unlikely | far unlikely | held-out species: true expected | skip unlikely |
|---|---|---|---|---|---|---|
| nearest | 67.7% | 10.9% | 86.6% | 97.9% | 58.7% | 87.9% |
| fallback | 60.8% | 14.6% | 86.6% | 97.9% | 57.8% | 88.0% |
| smoothed 2 | 65.4% | 11.7% | 87.6% | 98.1% | 60.0% | 88.4% |
| smoothed 3 | 65.2% | 12.8% | 88.0% | 98.2% | 60.5% | 88.6% |
| smoothed 5 | 63.8% | 13.4% | 88.9% | 98.3% | 60.8% | 89.1% |
| smoothed 10 | 61.1% | 14.2% | 88.8% | 98.3% | 61.3% | 89.1% |

(Held out: the 238 species the training database lacks, judged on their own lineage, whose clades do not hold
them. A third of true pairs are unlikely in every rule: the junctions between clusters, whose order the simulator
shuffles per species, so no clade predicts them.)

What detection needs: a database species' own reads judged on its lineage against the reads of a relative it lacks
landing on it (a held-out species' true pairs judged on the lineage of a training-database species of its genus,
else family; 235 pairs of species):

| rule | own true pairs expected | relative's expected | own unlikely | relative's unlikely |
|---|---|---|---|---|
| nearest | 72.8% | 63.9% | 5.0% | 30.5% |
| fallback | 65.0% | 62.2% | 6.4% | 30.5% |
| smoothed 2 | 70.6% | 63.8% | 5.9% | 30.8% |
| smoothed 3 | 70.5% | 63.8% | 7.1% | 31.2% |
| smoothed 5 | 69.7% | 63.8% | 7.9% | 31.5% |

The strict fallback nearly erases the difference between a species' own pairs and a relative's (2.8 points of
expected instead of 8.9). Smoothing falls back to the higher ranks in proportion to how sparse a clade is, keeps
most of that difference, and rejects skip and far pairs better than either. protal uses smoothing with 3
pseudo-species (`kPriorSpecies`): a family with the end informative in 3 species counts as much as the ranks
above it. If even the top clade with data has fewer than 5 (`kMinInformative`; only in tiny databases), a pairing
seen is expected and another unknown.

## Adjacency as evidence for the model

The larger build (`opw_build.sh ~/opw/b_gn2 --samples 8 --test-samples 6 --congeners 3`: 72 pe, 72 se, 24 pb,
24 ont training samples, 54/54/18/18 test samples, three species of a genus in every sample; 18.6 min, 4.7 GB peak)
with the smoothed rule. How well each adjacency feature alone separates present from absent taxa on its test set
([`feature_auc.py`](scripts/feature_auc.py); AUC over all taxa, and over those with reads across genes):

| reads | taxa present / absent (with links) | `adjacent_expected_share` | `adjacent_unlikely_share` | `adjacent_support` |
|---|---|---|---|---|
| pe | 1,068 / 3,250 (789 / 864) | 0.72 (0.43) | 0.56 (0.55) | 0.78 (0.73) |
| pb | 196 / 66 (190 / 34) | 0.67 (0.40) | 0.57 (0.39) | 0.58 (0.35) |
| ont | 232 / 353 (206 / 246) | 0.57 (0.46) | 0.48 (0.39) | 0.56 (0.46) |

Single-end reads have no reads across genes (all three constant). Most of the separation over all taxa is that
present taxa have reads across genes at all, which depth already says. Among taxa with such reads, only
`adjacent_support` (the frequency, not the binary verdict) separates paired-end taxa (0.73); for long reads the
absent taxa score higher. Present taxa's long reads have 10-12% of their consecutive genes judged unlikely: the
junctions between clusters, whose order the simulator shuffles per species, so that a species' own junctions are
rare in its family (in [`rule_check.py`](scripts/rule_check.py), 11-13% of a species' own true pairs).

Each read type's model retrained on the build's training table with and without the three features
([`compare_features.sh`](scripts/compare_features.sh), trainer seeds 1-3, the build's settings, `--evaluation
basic`), on its independent test set; mean ± standard deviation over the seeds:

| reads | features | F1, species held out | test F1 | test AP | test log loss | FP per test sample |
|---|---|---|---|---|---|---|
| pe | normalized | 0.9973 ± 0.0006 | 0.9597 ± 0.0023 | 0.9893 ± 0.0006 | 0.0749 ± 0.0013 | 0.00 |
| pe | + adjacency | 0.9974 ± 0.0005 | 0.9622 ± 0.0055 | 0.9898 ± 0.0018 | 0.0732 ± 0.0076 | 0.00 |
| se | normalized | 0.9945 ± 0.0008 | 0.9792 ± 0.0009 | 0.9908 ± 0.0008 | 0.0681 ± 0.0036 | 0.03 ± 0.02 |
| se | + adjacency | 0.9941 ± 0.0013 | 0.9795 ± 0.0003 | 0.9905 ± 0.0003 | 0.0685 ± 0.0014 | 0.02 ± 0.01 |
| pb | normalized | 0.9805 ± 0.0000 | 0.9915 ± 0.0015 | 0.9997 ± 0.0000 | 0.0377 ± 0.0018 | 0.06 |
| pb | + adjacency | 0.9790 ± 0.0013 | 0.9889 ± 0.0015 | 0.9998 ± 0.0001 | 0.0400 ± 0.0035 | 0.06 |
| ont | normalized | 0.9920 ± 0.0028 | 0.9833 ± 0.0033 | 0.9923 ± 0.0009 | 0.0924 ± 0.0099 | 0.09 ± 0.08 |
| ont | + adjacency | 0.9907 ± 0.0032 | 0.9848 ± 0.0038 | 0.9934 ± 0.0003 | 0.0795 ± 0.0021 | 0.07 ± 0.03 |

Single-end reads, which carry no adjacency at all, show the noise of a changed feature set (±0.0003-0.0004).
Paired-end +0.0025 test F1 and Nanopore +0.0015 (log loss 0.092 → 0.080) are within one or two standard
deviations; PacBio loses 0.0026. The first, smaller build with the strict fallback rule (4 samples per design point,
one seed) gave the same picture: pe +0.0012, se +0.0013 (the control), pb +0.0086, ont +0.0041 test F1. So the
adjacency features stay out of the normalised set the models train on; `--features normalized+adjacency` tries them
on a build whose clades differ in gene order as real ones do, which this simulator's (random junctions per species,
breaks per family) may not.

## Pairs and long reads

The 6 paired-end samples of 100,000 read pairs (their single-end halves too), 6 PacBio and 6 Nanopore samples of
30 Mb of the larger build's test set, profiled again against its training database with its trained models,
with gene neighbours and with `--no_gene_neighbours` ([`alignment_compare.sh`](scripts/alignment_compare.sh)):

| | with | without |
|---|---|---|
| fragments paired across two neighbouring genes (pe) | 12,545 (2.1% of the pairs) | 0 |
| guided mates found on the next gene (pe) | 877 | 0 |
| long-read genes found where the neighbours put them (pb, ont) | 446 | 0 |
| calls pe+se / pb / ont: TP, FP, FN | 266, 0, 2 / 75, 0, 0 / 77, 1, 0 | the same |
| wall time, peak memory (all 24 samples, 5 threads) | 64.7 s, 3.55 GB | 64.7 s, 3.52 GB |

So the table guides mate pairing (a pair across two genes is written as one proper pair, RNEXT the other gene),
mate rescue on the neighbouring gene and the long-read search for genes no segment of the read has, as in the
[earlier run](../2026-10-01-gene-neighbours-run/README.md); on these deep samples the calls do not change. The
binary of `5c7d64c` (the rule before this work) cannot read this build's index (format 0.7.2), so the counts were
not compared with the earlier rule on the same reads.

## The website

Out of date: it does not describe `gene_neighbours.tsv`, `gene_positions.tsv`, `--no_gene_neighbours`, pairs on
two references or the `adjacent_*` features (as before this work), and its list of database files lacks both
tables.

## Tests

The final working tree with the binaries of `~/bprof/gnb/build`: unit tests 273 passed, 1 skipped (as before);
the mini-database suite 41 passed; the e2e suite 122 passed on a mini database built afresh by
`build_mini_db.sh`. On the repository's old `data/mini_db/protal_db`, which has no `model_pe.xml`, 32 e2e tests
fail for that reason alone.

## Next

- A GTDB build: how many genes the k-mer traces place in real strains (other assemblies, other gene calls), how
  many representatives are one sequence, how sparse families are, and whether `--features normalized+adjacency`
  helps where gene order is real.
- A species' own gene order separated its reads from a relative's best (the *nearest* rule, in sample); the
  positions in the database could give a species-level table or feature, which the clade frequencies smooth away.
- The pairing and rescue counts of the earlier rule against the smoothed one on the same reads (needs the old
  rule in a 0.7.2 build).

## Follow-up (2026-10-02): the adjacency features in the default feature set

At the user's request the trainer trains on them by default: `model_features.DEFAULT_FEATURE_SET` is
`normalized+adjacency` (`random_forest_cmdline.py --features`, and `build_gtdb_database.py --features`, which passes
it on and records it in `build_metadata.tsv` as `classifier_features`). `--features normalized` leaves them out, to
test them on real data: a GTDB build with each, or one `--evaluation full` training, whose feature-set study now
compares `normalized`, `normalized+adjacency` and `all` with species held out. A training table without the
three columns (an older protal's dumps) stops the default with their names; `--features normalized` trains on it.
The results above, within noise on the operon world, are unchanged.
