# Strain rows from long reads, foreign genes, and a species' own gene order

- **Date**: 2026-10-02.
- **Code**: branch `audit-fixes`, this work committed on `f343113` (`d11381f`'s gene neighbours, with the HiFi reads,
  the one-binary work and the WFA2 flank change of other sessions on top; measured on `aeb7bf4` plus this work):
  `src/Profiling/Haplotypes.h` (new),
  `src/Profiling/Profiler.h`, `src/Profiling/Strain.h`, `src/RunProtal.h`, `src/Options.h`,
  `src/SequenceUtils/VariantHandler.h`, `scripts/mini_db/gene_neighbours.py`; tests `tests/test_Haplotypes.cpp` (new),
  `tests/test_GeneNeighbours.cpp`, `tests/e2e/test_protal_e2e.py`, `scripts/mini_db/test_mini_db.py`; docs
  `running.md`, `building-a-database.md`, `database-files.md`, `model-training.md`.
- **Request**: implement uses 2 (strain haplotypes across genes from long reads) and 6 (per-gene link evidence) of
  [2026-10-02-gene-neighbour-uses](../2026-10-02-gene-neighbour-uses.md); then, whether the same reads go into the
  MSA as into abundance and detection ([below](#which-reads-go-where)).
- **Data**: the strain benchmark of [2026-10-02-v072-benchmark](../2026-10-02-v072-benchmark/README.md) (WSL
  `~/bench071/strains`: 12 species of the v0.7.1 benchmark world, 8 strains each of known SNPs and tree) and its 0.7.2
  database; the operon world's build `~/opw/b_gn2` of
  [2026-10-02-gene-neighbour-frequencies](../2026-10-02-gene-neighbour-frequencies/README.md) (765 species, 135
  unknown ones, 238 held out of its training database).
- **Machine**: WSL Ubuntu 24.04, Core Ultra 7 258V (6 cores), shared with other sessions (load 5-12), so wall times
  are not compared; CPU times are.
- **Scripts**: [`scripts/`](scripts/) (below); data in WSL `~/bench071/phasing`, `~/fgw`, `~/fg`; build in
  `~/protal-hap`.

## Summary

- **Species lines in `gene_neighbours.tsv`.** The table was per clade (family up to domain), so where a species'
  gene order differs from its family's, its own neighbours read as unlikely, and per-gene evidence (use 6) would take
  its own genes for foreign ones. `gene_neighbours.py` now also writes lines of the species itself for every partner
  any of its genomes shows that its clades alone would not expect (smoothed share below 0.2). protal reads them with
  no change: a species' lineage already starts with the species. `--no_species_lines` leaves them out.
- **Foreign genes (use 6).** A gene whose reads' neighbouring genes are mostly unlikely neighbours in its taxon's
  clade (4 or more links judged, more than half unlikely) is foreign: its reads come from another genome. It is left
  out of the taxon's depth and strain MSA rows (`--keep_foreign_genes` keeps it) and listed in `.profile.genes.log`.
  In samples where a genome absent from the database carries a near copy of a database species' gene among other
  neighbours, the gene was flagged in 3 of 8 cases with paired-end reads, 6 of 8 with PacBio and Nanopore reads
  (where the copy's context differed from the donor's family at both ends); left out, its foreign alleles leave the
  species' MSA row. The depth hardly changes either way (the median over genes ignores one inflated gene). False
  flags on the operon world's 18 test samples: 1 of 4,707 judged genes of present species (281 with the clades'
  lines alone); calls unchanged.
- **Strain rows (use 2).** A PacBio or ONT sample whose reads show two or more strains of a species gets a row per
  strain in the strain MSA (`<sample>_hap1`, ...), each called from its strain's reads. On 96 species-in-mixture
  cases (two or three strains at 85:15 to 50:50, 2-40x), PacBio HiFi reads recover 52% of the strains as a row of
  their own (9% without phasing: the majority strain of a lopsided mixture), with SNP recall 0.54 against 0.26 and
  precision 0.973 against 0.948; Nanopore 27% (9%), recall 0.40 (0.27), precision 0.928 (0.913). No sample of one
  strain got strain rows. 7-14% more CPU. `--no_phasing` writes one row per sample, as before.

## A species' own gene order: species lines

[`gene_neighbours.py`](../../../scripts/mini_db/gene_neighbours.py) counts, per clade, how many of its species have a
gene end facing a partner; a species counts once, with the partner most of its genomes show. protal smooths a pairing's
share from the top clade down to the family and calls it expected from 0.2, unlikely up to 0.05. A species whose
order differs from its family's, at an end of a large family, then has its own pairing unlikely: on the synthetic
operon world (whose simulator rearranges every species' clusters at random) 11-13% of a species' own true pairs were
([2026-10-02-gene-neighbour-frequencies](../2026-10-02-gene-neighbour-frequencies/README.md)). For the existing uses
that only weakened a feature; for per-gene evidence it would flag the species' own genes.

`count_clades` now adds, after the clades' lines, lines of the species (clade: its taxid; species 1, informative 1)
for every partner of an end that any of its genomes shows, if the partner's smoothed share in its clades (computed
as protal computes it, `smoothed_share`) is below 0.2 in clades with enough species to judge (the top clade
informative in 5 or more). protal's lineage of a species starts with the species, so the species becomes the nearest
clade: its partner's share is at least (1 + 3 x 0)/4 = 0.25, expected, while its family's pairings keep 3/4 of their
share. The first version wrote only the species' majority partner; the scenario below showed that a representative
whose neighbour differs from its two other strains' (partner seen in 1 of 3 genomes) still had its own genome's
pairing unlikely. Every partner seen is written since.

| table (operon world, 2,295 genomes of 765 species) | lines | of them species lines |
|---|---|---|
| clades only (as built at `d11381f`) | 258,592 | 0 |
| species lines, majority partner | 311,294 | 52,702 |
| species lines, every partner seen (the default) | 319,535 | 60,943 |

That is many (the simulator reorders each species); in real genomes, whose order changes less within families, there
should be fewer. Each line costs ~36 bytes in a run. A training database's copy derives them anew without the
held-out species (`--from_positions`, `gtdb_to_protal_db.py --from_db`), as for the clades. Tests: three unit tests
of `count_clades` (`SpeciesLinesTest`: the odd species of six gets the three lines; none without them or in a
family of four; a strain's partners too) and one of `Table::Assess` (`ASpeciesOwnLineMakesItsOwnPairingExpected`).

## Foreign genes

### What

`MicrobialProfile::FinishLink` already judged, for the adjacency features, each two genes next to each other on a
read (a pair's mates on two genes; a long read's consecutive genes at most 3 kb apart on the read) in the clades of
each gene's taxon. It now also counts per gene (`GeneLinks`: judged, unlikely), from every read's best records before
the filters, as the features are. A gene is *foreign* (`Taxon::ForeignGene`) with 4 or more judged links of which more
than half are unlikely (smoothed share 0.05 or less): a gene with one end in context and one not (half) is not. A
foreign gene is left out of its taxon's depth (from the median over genes and from the expected length, so it counts
as a gene the taxon lacks; `low_identity_share` keeps it) and of the sample's strain MSA row, and not phased;
`--keep_foreign_genes` keeps it. `.profile.genes.log` has three new columns, `LinksJudged`, `LinksUnlikely` and
`Foreign`, and the log a line per sample with foreign genes. Paired-end reads judge only the 2% of pairs across two
genes, so a gene needs some depth to be judged at all; long reads judge nearly every gene end that has a marker
within 3 kb.

### False flags: the operon world's test samples

[`foreign_compare.sh`](scripts/foreign_compare.sh) and [`foreign_summary.py`](scripts/foreign_summary.py): the build's
6 paired-end (100,000 pairs), 6 PacBio and 6 Nanopore (30 Mb) test samples profiled against its training database
(which lacks the 238 held-out species, so absent relatives' reads are there) with its trained models: the database as
built (clade lines only), its table derived again with species lines (`--from_positions` on its own
`gene_positions.tsv`: 233,399 lines), with foreign genes kept and dropped. "Present taxa": in the sample; "absent
taxa": taxa with reads that are not (a held-out species' reads on its congener, mostly).

| table | reads | TP | FP | FN | genes judged (4+ links) of present taxa | foreign | genes judged of absent taxa | foreign |
|---|---|---|---|---|---|---|---|---|
| clade lines only | pe | 77 | 0 | 0 | 933 | 6 | 268 | 7 |
| | pb | 75 | 0 | 0 | 1,761 | 137 | 150 | 5 |
| | ont | 77 | 1 | 0 | 2,010 | 138 | 549 | 25 |
| species lines, majority partner (first version) | pe | 77 | 0 | 0 | 934 | 6 | 268 | 7 |
| | pb | 75 | 0 | 0 | 1,761 | 13 | 150 | 5 |
| | ont | 77 | 1 | 0 | 2,011 | 12 | 549 | 25 |
| **species lines, every partner** | pe | 77 | 0 | 0 | 935 | **0** | 268 | 7 |
| | pb | 75 | 0 | 0 | 1,761 | **1** | 150 | 5 |
| | ont | 77 | 1 | 0 | 2,011 | **0** | 549 | 25 |

With the clades' lines alone, 7-8% of the long-read genes of present species are flagged: the simulator's
rearranged junctions, which their families do not share. The species' majority partner cut that tenfold; the rest
were representatives whose neighbour differs from their other strains' (the diagnosis below). With every partner
the species' genomes show, 1 gene of 4,707 judged of present taxa is flagged. Calls do not change in any arm (the
dropped genes leave the depth: present taxa's depth changes by at most 0.4%).

### Where there are foreign genes: a transferred gene

[`foreign_world.py`](scripts/foreign_world.py), [`foreign_world.sh`](scripts/foreign_world.sh),
[`foreign_world_summary.py`](scripts/foreign_world_summary.py), [`foreign_world_others.py`](scripts/foreign_world_others.py):
on the operon world, 8 triples of a recipient O (an unknown species' representative, absent from the database), a
donor T (a database species of another phylum) and a marker A of both, chosen so that each end of A in O faces
another marker within 3 kb and those pairings are unlikely in T's clades, while T's own neighbours of A are expected.
O' is O with its gene A replaced by T's, 1% substituted, in O's orientation: a recent transfer. Sample hgt<i>: T at
10x, O' at 30x, 20 database species at 1-5x; ctl<i>: the same without O'. Paired-end (ART HSXt 2x150), PacBio
(`hifi_reads.py`) and Nanopore (pbsim3) reads, profiled against the full database with species lines (its table
derived again: 319,535 lines), gene A kept and dropped, with strain MSAs.

| reads | T's gene A foreign | gene A's depth in hgt / ctl (median of the 8) | T's depth, hgt against ctl: kept / dropped (median) | gene A's IUPAC positions in T's MSA row, kept |
|---|---|---|---|---|
| pe | 3 of 8 (9-30 links, all unlikely); the 5 others: 0-14 links, none unlikely | 37 / 9 | -0.004 / -0.005 | 2-18 (all 8); dropped: none in the 3 flagged |
| pb | 6 of 8 (51-102 links, 61-84% unlikely); the 2 others 18-21 links, none unlikely | 31 / 12 | -0.018 / -0.018 | 0-18; dropped: none in the 6 flagged |
| ont | 6 of 8 (36-56 links, 53-75% unlikely); one at 17 of 36 (47%), one 9 links, none | 27 / 9 | -0.016 / -0.019 | 2-16; dropped: none in the 6 flagged |

Where gene A is not flagged, its judged links are all expected: the recipient's reads that cross from A into its
neighbours were not judged (which reads they are was not traced). The recipient's reads triple the gene's depth, but
the median over genes ignores one inflated gene: T's depth is a median 2% from the control's either way (single
triples up to 14% and 34%, the hgt and control samples being separate draws of reads, kept or dropped alike), so
leaving the gene out matters for the MSA (its foreign alleles: 2-18 IUPAC positions in T's row), not for abundance. Other flags in these 16 samples: on species of the sample, 0 (pe), 0 (pb) and 23 (ont, all in hgt
samples, so where the recipient's reads are); on taxa that are not in the sample (holding the recipient's or other
genomes' reads), 19, 194 and 226, none of them called.

The first run of this scenario, with the majority-partner species lines, flagged 23 (pe), 118 (pb) and 146 (ont)
genes of species of the samples, controls as often as hgt samples. [`foreign_diagnose.py`](scripts/foreign_diagnose.py)
on one PacBio control: each flagged gene's reads joined it, at one end, to a gene its representative (the sample's
genome) has there, while the species' two other genomes have none within 3 kb, so the species' line said "none"; hence
lines for every partner. The same scenario showed the phasing hazard above: T's one gene carried O's reads, a second
"strain" on one gene, which split T into two rows until phasing needed 3 genes.

### Tests

Unit: `MicrobialProfile.AGeneWhoseNeighboursAreUnlikelyIsForeign` (4 unlikely links make a gene foreign, half
unlikely does not, 3 are too few; the depth with the gene left out); the species lines' tests above.

## Strain rows from long reads

### How

[`Haplotypes.h`](../../../src/Profiling/Haplotypes.h). A long-read sample's taxa keep, for each record they take,
what the read shows of the gene in the part the variant caller trusts (`ReadAlleles`: interval, SNP bases and
qualities, positions without a base), with the read, its place on the read, strand and divergence
(`haplotypes::ReadRecord`; `Taxon::SetKeepPhaseRecords`, only for pb and ont samples when strains are written). At
the strain MSA (`RunProtal.h` `PhaseSample`, `GetMSAForTaxon`), per sample:

1. **Sites**: the sample's multi-allelic positions as the MSA writes them (an IUPAC code: the call and another base
   pass the SNP filters; `MultiAllelicSites`), from the taxon's own reads (`--msa_identity_margin`), on the genes the
   MSA takes, not on foreign genes it drops.
2. **Reads**: each read's records in read order, cut where two consecutive genes within 3 kb are unlikely neighbours
   in the taxon's clade (a chimera; `UnlikelyNeighbours`), each part a unit of its alleles at the sites; records of
   more than the MSA's divergence left out.
3. **Blocks**: the sites that units link, directly or through other sites (union-find).
4. **Haplotypes per block** (`detail::Clustering`): seed and extend. A haplotype starts from the unit with the most
   sites (the most like the block's majority alleles), and takes, pass after pass, the units that share sites with
   it, agree at all but 20% of them, and fit it better than any other haplotype; when none joins, a unit that
   conflicts with every haplotype it shares sites with (2 sites or more, 1 in a one-site block) seeds the next (at most
   4). Units that conflict a little then join their best haplotype; a haplotype with fewer than 3 units or 5% of the
   block's gives them to the others, as does one that differs from another at fewer than 2 sites.
5. **Rows**: as many as the blocks most often have haplotypes (weighted by their sites); their shares pooled by rank
   over the blocks with that many. The anchor block (that many haplotypes, the most sites) gives the rows its
   haplotypes by rank; any other block's haplotypes join the rows by the multinomial likelihood of their units'
   counts (each row a haplotype, several rows one haplotype if the block has fewer), if the best join is 20 times as
   likely as the next (log odds 3), else the block is not phased. No rows unless the phased blocks span 3 genes or
   more (below).
6. **Rows' sequences**: each row is called from the records of its units in the phased blocks
   (`Phasing::row_records`, `HaplotypeItem`), as a sample's row from the sample's reads (`MSAItem`: coverage per
   strand, SNPs by base quality, `PostProcessSNPBin` with the sample's read type's filters), `-` where its reads do
   not reach. The rows get their own `.meta.tsv` lines and `.snp_stats.tsv` rows, so qcmsa treats them as rows.
   `<species>.haplotypes.tsv` lists every block of every long-read sample.

### Three rounds

The design changed twice on the benchmark below:

- **Round 1, split and refine**: one cluster split at the site with the strongest second allele, then every read
  moved to the cluster whose majorities it matched best (k-means), recursively. In a 50:50 sample the clusters came
  out as 100, 99 and 53 reads of a species of two strains, all of whose reads were the two strains' (checked by the
  reads' source genomes): on a block longer than a read (up to 37 genes), the refinement settled into a cluster of
  one strain's left half with the other's right half. Replaced by seed and extend, as read-based haplotype assembly
  grows haplotypes along overlapping reads (`ABlockLongerThanItsReadsIsPhasedAlongThem`).
- **Round 2, the sample's row with the strains' bases at the sites**: a strain row was the sample's row, with the
  haplotype's allele at each phased site. In the e2e test (strain 2 at 30%, 1.5% from strain 1) the minor row showed
  strain 1's base at 6% of the positions where the strains differ: positions where strain 2's allele failed the
  sample's filters (low coverage) were no sites, so both rows inherited strain 1's call. Replaced by rows called from
  their own reads.
- **Round 3, the rows from their own reads**, plus the 3-gene minimum, which the foreign-gene scenario called for:
  there a donor species' one gene carried another genome's reads, which looked like a second strain on one gene and
  split the species into two rows.

### Benchmark

[`make_mixtures.py`](scripts/make_mixtures.py): from the strain benchmark's 12 species of 8 strains, 16 samples:
`pure1..8` (strain j of every species) and 8 mixtures (strains 1+2, 3+4, 5+6, 7+8 at 70:30; 1+5, 2+6 at 50:50; 3+7+8
at 50:30:20; 4+5 at 85:15), each species at its benchmark coverage (2-40x, the mean of its 8 strains'), with the 15
background genomes of the benchmark's sample of that number; PacBio HiFi reads by `hifi_reads.py` (15 +- 3 kb, the
collector's setup) and Nanopore by pbsim3 (QSHMM-ONT-HQ, 8 +- 6 kb, 97%). [`run_phasing.sh`](scripts/run_phasing.sh):
the 16 samples per read type in one run against the 0.7.2 database, with and without `--no_phasing`, with qcmsa.
[`phasing_score.py`](scripts/phasing_score.py): per mixture sample and species, the positions where its strains
differ; each row matched to the strain whose bases it shows most there; a strain *recovered* if a row shows its base
at half of those positions or more with at most 10% of its calls another strain's; SNPs over all gene positions
against the reference: TP and FN every strain's SNPs by its best row (a strain without a row misses all), FP a row's
base that is neither the reference's nor its strain's; *placed*: the row nearest (p-distance) to the pure sample of
its strain.

| reads | phasing | strains recovered (of 204) | resolved (mean) | wrong calls (mean) | nearest pure sample is the strain's | rows qcmsa keeps | SNP precision | SNP recall | pure samples given strain rows |
|---|---|---|---|---|---|---|---|---|---|
| pb | off | 0.09 | 0.126 | 0.064 | 0.94 | 0.41 | 0.948 | 0.259 | 0 |
| pb | round 1 | 0.34 | 0.492 | 0.113 | 0.97 | 0.69 | 0.886 | 0.545 | 0 |
| pb | round 2 | 0.45 | 0.550 | 0.086 | 0.97 | 0.68 | 0.913 | 0.598 | 0 |
| pb | **on** | **0.52** | 0.552 | **0.050** | 0.96 | **0.86** | **0.973** | 0.541 | 0 |
| ont | off | 0.09 | 0.166 | 0.078 | 0.93 | 0.41 | 0.913 | 0.268 | 0 |
| ont | round 1 | 0.23 | 0.434 | 0.136 | 0.93 | 0.63 | 0.853 | 0.479 | 1 |
| ont | round 2 | 0.25 | 0.410 | 0.114 | 0.92 | 0.59 | 0.862 | 0.466 | 2 |
| ont | **on** | **0.27** | 0.386 | 0.097 | 0.95 | **0.75** | **0.928** | 0.398 | 0 |

(204 strains: 12 species x (6 two-strain samples + one of three) x their strains. "Placed" and "kept" count the
strains that have a row.) Strains recovered by design and the species' coverage in the sample, the final version:

| reads | design | strain's share | < 8x | 8-20x | > 20x |
|---|---|---|---|---|---|
| pb | 2 strains | 0.85 | 0.75 (4) | 0.33 (3) | 1.00 (5) |
| pb | 2 strains | 0.7 | 0.38 (16) | 0.92 (12) | 1.00 (20) |
| pb | 2 strains | 0.3 | 0.06 (16) | 0.92 (12) | 1.00 (20) |
| pb | 2 strains | 0.15 | 0.00 (4) | 0.00 (3) | 1.00 (5) |
| pb | 2 strains | 0.5 | 0.00 (16) | 0.33 (12) | 0.40 (20) |
| pb | 3 strains | 0.5 / 0.3 / 0.2 | 0.50 / 0 / 0 (4) | 1.00 / 0.33 / 0 (3) | 0.40 / 0.60 / 0.20 (5) |
| ont | 2 strains | 0.85 | 0.25 (4) | 0.00 (3) | 1.00 (5) |
| ont | 2 strains | 0.7 | 0.00 (16) | 0.25 (12) | 0.95 (20) |
| ont | 2 strains | 0.3 | 0.00 (16) | 0.00 (12) | 0.85 (20) |
| ont | 2 strains | 0.15 | 0.00 (4) | 0.00 (3) | 0.20 (5) |
| ont | 2 strains | 0.5 | 0.00 (16) | 0.00 (12) | 0.20 (20) |
| ont | 3 strains | 0.5 / 0.3 / 0.2 | 0 / 0 / 0 (4) | 0.33 / 0 / 0 (3) | 0.60 / 0.20 / 0.20 (5) |

Without phasing only the majority strain of 85:15 is recovered (its row is its strain's, the minority's alleles
failing the 15% frequency). With PacBio reads from 8x on, both strains of 70:30 are (92-100%); 50:50 mixtures stay
mostly unrecovered, as no read links the blocks and their shares cannot tell which strain of one block is which of
another (only the anchor block is phased); a strain at 20% or less needs depth. Nanopore needs more than 20x. The
rows' SNP recall is below round 2's because a row now has only its own reads (a `-` where they do not reach) instead of
the sample's calls, which were often the other strain's; precision is higher than without phasing.

Cost (CPU, all 16 samples per read type): PacBio 696 s against 650 s (+7%), Nanopore 1,041 s against 913 s (+14%);
peak memory 3.76 GB against 3.76 GB (pb), 4.19 against 3.99 GB (ont).

### Tests

Unit (`test_Haplotypes.cpp`, 9): rows across two linked genes; blocks no read links joined by 70:30 shares and not by
50:50 (only the anchor); a block longer than its reads phased along them; a second allele at one site (an error) no
strain; two strains of two genes no strains with the run's 3-gene minimum; three strains; own reads only and reads
cut between unlikely neighbours; `ReadAlleles` from `AddAlignment`; a row called from its own reads
(`HaplotypeItem`). End-to-end (`PhasingTest`, 2): a two-strain PacBio sample gets `pm_hap1` and `pm_hap2`, each
showing its own strain's base at 98% or more of the differing positions, with `.haplotypes.tsv` and `.meta.tsv`
lines; with `--no_phasing` one row with IUPAC codes and no `.haplotypes.tsv`.

## Which reads go where

Not the same reads: the three uses take nested sets, from the same records. Every one starts from a read's best
alignment (`BestOfGroup`: the primary record of each read, each mate or gene of a long read with its own), in
`Profiler::PrepareMAPQ` (`Profiler.h`); from there (code at this commit):

| step | presence model (detection) | abundance (depth) | strain MSA row |
|---|---|---|---|
| records before any filter | the "evidence" features only: `low_mapq_share`, `congener_fit_share`, `other_genus_fit_share`, `excess_*`, `conserved_fast_record_ratio`, `adjacent_*` (`NoteRecord`, `NoteLinkedRecord`), also records the filters below drop | no | no |
| MAPQ >= 4 and more than 50 aligned bases (`m_min_mapq`, `m_min_alignment_length`) | yes, all other features | yes | yes |
| genes without unique k-mers (not hittable: `CheckSam` kSkip) | left out (the evidence features above keep them) | left out | left out |
| reads by identity, against the taxon's best reads (`TopIdentity`, 98th percentile) | every read | within `--depth_identity_margin` 0.08 (`OwnIdentityThreshold`; per gene with `--gene_conservation db`) | within `--msa_identity_margin` 0.04 (`IdentityThreshold`; never per gene), in steps of 0.5% (`DivergenceBin`) |
| foreign genes (new) | kept, except in `depth` | left out | left out (and not phased) |
| genes with no long unique k-mer, for a species with relatives (90% rule, `SelectGenesForTaxon`) | kept | kept | left out |
| within a record | counts and identity of the whole record; the allele features (A, AF, RAF) from the variants as the strain container records them (as for the MSA, but of every read) | all its reference bases, both mates of a fragment in full (`MappedLength`) | the part the variant caller trusts (`AddAlignment`: ends with a cluster of mismatches trimmed), the second mate without the part the first covered, N and deleted positions without a base |
| samples | each sample alone | each sample alone | samples whose taxon passes `--msa_knob` (default: the sample's knob, so those that report it) |
| strain rows (new, long reads) | - | - | the reads with an allele at a multi-allelic site, in a phased block, cut at unlikely neighbours |

Whether that makes sense, one by one:

- **Detection takes the most reads, on purpose.** Its features have to see the reads that are not the taxon's own to
  tell a present species from one that holds a relative's reads: `low_identity_share`, the congener fits, the
  divergence beyond the base qualities and the conservation ratios are all about those reads, and the evidence
  features are counted before the MAPQ filter because a relative's reads fit several taxa (MAPQ near 0). Sensible.
- **Abundance takes the taxon's own reads (0.08), the MSA stricter ones (0.04).** Both are relative to the taxon's
  best reads in the sample, not to the reference, so a strain 5% from the reference keeps its reads in both. The
  margins were chosen separately: 0.08 counted strains' bases best ([2026-09-30](../2026-09-30-depth-margin-stress/README.md));
  for the MSA a relative's reads within 0.08 put 10 times the false allele calls into mixed samples' raw MSAs. A row
  with fewer reads is a smaller loss than false alleles, so the stricter margin makes sense. One consequence for the
  new strain rows: in a sample of a close and a far strain (say 0.5% and 5% from the reference), the best reads are
  the close strain's, and the far one's fall outside 0.04: abundance counts them, the MSA and phasing do not, so the
  far strain gets no row. Not tested here (the benchmark's strains are 0.4-1.2% deep); worth a world with such pairs.
- **`--gene_conservation db` breaks the nesting on the most conserved genes.** The depth margin becomes 0.03 + 0.05 x
  the gene's factor, below the MSA's 0.04 for factors under 0.2, so there the MSA keeps reads the depth does not. Off
  by default, and the depth rule is meant to be strict there; harmless, but `min(msa margin, gene's depth margin)`
  would keep the MSA's reads a subset of the depth's. Not changed.
- **Foreign genes leave the depth and the MSA but not the detection features** (only `depth` changes). The model was
  trained with all genes, and its adjacency features already carry the taxon-wide signal of reads from another
  context. Sensible as it is; dropping them from detection would need retraining and a test.
- **The MSA's gene filter (no long unique k-mer) is the MSA's only.** Such a gene a relative may share unchanged, so
  its reads, the relative's included, would show as a second strain; for the depth, the median over genes hardly
  feels one inflated gene. Sensible.
- **Overlapping mates count twice in the depth, once in the MSA.** The depth sums the reference bases of every record,
  so a fragment shorter than two reads adds its overlap twice: the depth (also the model's `depth` feature) of a
  library of short fragments is higher than of one of long fragments at the same coverage. Within a sample every
  taxon gets the same factor, so relative abundances do not change; across libraries the absolute depths, and the
  depth feature the model was trained on (350 bp fragments), do. Counting a fragment's covered reference once, as the
  strain container does, would remove that; a small change, but it shifts a model feature, so it wants a retraining
  and was left alone here. The untrusted ends (a few bases of some reads) count for depth and not for the MSA: negligible.
- **Strain rows take a subset of the MSA's reads**: a read without an allele at any site cannot be given a strain, so
  a row is `-` where only such reads reach. That is the price of calling each row from its own reads, which raised
  the rows' SNP precision above the unphased rows' ([above](#benchmark)).

So the rules differ where the uses differ; the two I would change if they matter on real data are the overlapping
mates in the depth (a cleaner fragment count, with retraining) and the far-strain case for phasing (a margin relative
to each haplotype's own best reads, which the clustering could provide).

Follow-up (2026-10-02): the depth now counts a fragment's bases once, as the MSA does
([2026-10-02-fragment-depth](../2026-10-02-fragment-depth/README.md)); on the operon world's 54 paired-end test
samples no call changed, the present taxa's depth fell by a median 0-1.2% by read setup.

## Tests run

On the final code (WSL build `~/protal-hap` from the working tree), with a mini database built afresh (its
`gene_neighbours.tsv` with species lines):

- unit: 285 tests, 284 passed, 1 skipped (`GeneSequence.TheBasesOutsideAWindowAreNotLeftToChance`, skipped before
  too). An earlier full run, under load from other sessions, failed `SamFile.AFileCutAtABlockOrFrameBoundaryIsIncomplete`
  once; it passed alone and in two full runs since (SAM I/O, which this work does not touch);
- end-to-end: 123 tests, all passed (the 2 new `PhasingTest` included);
- mini database: 42 tests, passed, 2 skipped (the GTDB build tests, which need ART and a training Python).

Result tables: [`results/`](results/) (`phasing_final.md`, `phasing_round1.md`, `phasing_round2.md`,
`phasing_final_rows.tsv`, `mixtures_truth.tsv`, `phasing_cost.txt`, `foreign_test_samples.md`,
`foreign_test_samples_majority_lines.md`, `foreign_scenario.md`, `foreign_scenario_other_flags.md`,
`foreign_scenario_triples.tsv`).

Commands, in WSL (environment `~/micromamba/envs/protal-db-build`):

```bash
python3 scripts/make_mixtures.py ~/protal-hap/src/scripts ~/bench071   # the mixtures' reads
bash scripts/run_phasing.sh && python3 scripts/phasing_score.py ~/bench071
bash scripts/foreign_world.sh ~/fgw 6                                     # the scenario (and foreign_world_others.py ~/fgw)
bash scripts/foreign_compare.sh ~/fg 6                                    # false flags on the test samples
python3 scripts/foreign_diagnose.py ~/fgw ctl1 pb                         # one sample's flagged genes, read by read
```

## The website

Out of date: it does not describe `--no_phasing`, `--keep_foreign_genes`, the strain rows `<sample>_hap<k>` and
`<species>.haplotypes.tsv`, the new `.profile.genes.log` columns (`LinksJudged`, `LinksUnlikely`, `Foreign`), or
`gene_neighbours.py --no_species_lines` (it does not describe the gene neighbours at all yet).

## Open

- 50:50 mixtures: without a read across blocks, phasing needs a link of another kind (a long read's coverage
  pattern, or the same strain in another sample of other shares).
- A strain far from the reference beside a close one: the MSA's identity margin is relative to the taxon's best
  reads (the close strain's), so the far strain's reads may fall outside it and get no row ([below](#which-reads-go-where)).
- The species lines on real genomes (GTDB r226): how many, and the false-flag rate there.
