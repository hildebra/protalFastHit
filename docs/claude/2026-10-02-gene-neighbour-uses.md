# What uses the gene neighbours, and what else could

2026-10-02, branch `audit-fixes` at `d11381f` (0.7.2 plus gene neighbours from every genome). This is from reading
the code; nothing was run. The numbers quoted come from
[2026-10-02-gene-neighbour-frequencies](2026-10-02-gene-neighbour-frequencies/README.md) and
[2026-10-01-gene-neighbours](2026-10-01-gene-neighbours.md).

## What uses them now

A run loads only `gene_neighbours.tsv`, the frequencies per clade (`RunProtal.h` `LoadGeneNeighbours`). It does not
load `gene_positions.tsv`. Every use follows *expected* partners only (smoothed share 0.2 or more), at most the two
most common, and only genes the species has in the database (`across_genes::Neighbours`).

| where | what | code |
|---|---|---|
| read mapping, pe | mates on two genes of one taxon whose facing ends are expected neighbours, with a fragment short enough, become one proper pair (RNEXT the other gene, TLEN 0) with the pair's score | `classify::JoinAlignmentPairs`, `across_genes::PairAcrossNeighbours` |
| read mapping, pe | a guided mate (sure mate at MAPQ 20 or more) is also looked for past its gene's end, on the neighbour there | `MateGuidance.h` step 2, `across_genes::NeighbourStretches` |
| read mapping, long reads | past each end of a read's gene, the clade's neighbours are looked for where they should lie on the read; an unseeded gene is anchored on its best diagonal | `LongReadAligner::AddNeighbourGenes` |
| profiling and the model | `adjacent_expected_share`, `adjacent_unlikely_share`, `adjacent_support`, counted from each read's best records before the MAPQ and length filters | `RecordEvidenceCollector::FinishLink`, `TaxonFeatures` |
| training | the three features are in the trainer's default set `normalized+adjacency` (`model_features.DEFAULT_FEATURE_SET`, `build_gtdb_database.py --features`); the training copy of the database derives its frequencies without the held-out species (`gtdb_to_protal_db.py --from_db`) | `scripts/model_features.py`, `scripts/mini_db/gtdb_to_protal_db.py` |

So "read mapping and RF training" is right, with two additions:

- **The model depends on the table beyond its three features.** Pairs across genes change candidate ranking and
  MAPQ, so `mean_mapq` and `low_mapq_share` also differ with and without the table. A model trained on a database
  with gene neighbours expects them in every one of those features.
- **Nothing checks that the model and the run agree.** With `--no_gene_neighbours`, or a database without the table,
  the three features fall to 0, 0 and 0.5 ("no reads across genes"). That is what a taxon with no links looks like,
  so a model that learned from them sees present taxa with paired-end or long reads as shallower than they are. The
  run says "not used" but gives no warning about the model (`ModelContractProblemInXml` checks only that the inputs
  exist). The effect is probably small (test F1 changed within noise with and without the features), but the
  mismatch goes unnoticed.

The gain so far: 2.1% of pairs paired across genes, 877 mates rescued on the neighbour, 446 long-read genes found,
calls unchanged on deep test samples. The features changed test F1 within noise (pe +0.0025, ont +0.0015, pb
−0.0026), on a simulator whose random junctions make 11-13% of a species' own true pairs look unlikely.

Not used by: the long-read read consensus (`SettleByRead`), strains and MSAs (`Strain.h`), depth and abundance, and
the database build beyond checking and packing the tables.

## Further uses

Ranked by value against effort.

### 1. Record in the model whether it was trained with gene neighbours (small, correctness)

Have the trainer write a PMML header `Extension` (as `protal_depth_knobs` does), for example
`protal_gene_neighbours`. protal would then warn when such a model runs with `--no_gene_neighbours` or on a database
without the table. `build_metadata.tsv` records the feature set already, but the run does not read it.

### 2. Strain haplotypes across genes from long reads (large, high value for HiFi)

`Strain.h` calls one consensus per gene and sample, so two strains of a species blend gene by gene. Nothing says
which alleles of gene A go with which of gene B. A long read that covers several neighbouring markers is a haplotype
across them. The bac120 set has 27 ribosomal proteins. About 15 of its markers lie in the S10-spc-alpha cluster
(L3 to L17, with `secY` and `rpoA`), about 14 kb in *E. coli*, so one HiFi read of 15-20 kb can carry an eighth of a
species' markers. With the table, a read's genes can be checked as one stretch of one genome (expected pairs, gaps within the
clade's range), and reads can be clustered by their alleles across genes into per-strain MSAs. Paired-end reads
across genes add links between two genes only, and 2% of pairs carry one: too few to phase alone. This fits the
PacBio HiFi work that is going on (`docs/claude/2026-10-02-pacbio-hifi-reads/`).

### 3. Gene order in the long-read consensus: chimeric reads (moderate)

`SettleByRead` gives a read's genes the read's consensus taxon. A gene keeps its own hit only if another taxon fits
it clearly better (MAPQ 4 or more). It does not look at gene order. Two consecutive segments whose pairing is unlikely
in the consensus taxon's clade, or whose gap on the read is far outside the clade's range, point to a chimera
(library chimeras are around 1% or more of Nanopore reads) or to a gene from somewhere else on the read. The read
could then be split into two consensus groups at that junction, rather than pulling the second part onto the first
part's taxon. This needs a tolerance: own junctions are unlikely 11-13% of the time in the simulator, and real
rearrangements happen. Test it on a world with real gene order (use 4).

### 4. Real gene order in the synthetic worlds (cheap, makes the tests fair)

`simulate_gtdb_release.py --operons` lays markers out in clusters with random junctions per species and breaks per
family. The adjacency features were judged on that layout. A GTDB build's `gene_positions.tsv` holds the real
marker order of thousands of genomes. The simulator could take its layout per clade from such a file, with each
species keeping its clade's order plus rare rearrangements. Uses 2, 3 and 6, and the default feature set, could then
be tested without real metagenomes.

### 5. Checking the representatives at build time (moderate, uses what the build already computes)

`gene_neighbours.py` places every representative gene in every other genome of its species, and records its contig
and neighbours. Many GTDB representatives are MAGs. A representative's gene that no other genome of the species
carries, or that sits on a short contig in a context its clade never shows, is likely contamination of the bin
(or a gene the species varies in). The build could list such genes in `build_metadata.tsv` or a QC table, and
optionally leave them out. A contaminant marker gives the species reads that belong to another taxon. That inflates
its gene presence and puts a foreign gene in its MSA.

### 6. Per-gene link evidence in profiling (uncertain value, cheap to test)

The adjacency counts are summed per taxon. Kept per gene, a gene whose links are mostly unlikely in the taxon's clade
gets its reads from another genomic context (a relative's homolog, a transferred gene, a contaminant). It could be
left out of the depth average or of the MSA, in the way the long-unique filter leaves genes out. Paired-end links are
few per gene at modest depth, so this is mainly for long reads. Whether it is worth doing can be measured on existing
SAM files (`--profile_only`) before any change.

### 7. A fragment length for pairs across genes (small)

Such pairs are written with TLEN 0. Each gene's reach to its facing end plus the clade's median gap gives an
estimate. Downstream tools would see a fragment length. Pairs on one gene and across genes together could also
estimate the sample's fragment lengths, in place of the fixed `kRescueMaxFragment` (1,000) of mate guidance.

### 8. Replication rate from gene positions (speculative, new output)

Coverage falls from the replication origin to the terminus in growing cells. Peak-to-trough methods (iRep, CoPTR)
estimate growth rate from that slope. The complete genomes in `gene_positions.tsv` give each marker's distance from
the origin (found by GC skew at build time), and that distance does not change under the inversions that are common
around the origin-terminus axis. A per-species distance per gene, stored at build time, would let the profiler fit
depth against it. Caveats: about 120 points per species, clustered near the origin (ribosomal operons sit there in
many fast growers), MAG representatives without a known origin, and depth noise at low coverage. This is research
work, not a fix. A side benefit: depth corrected for position would lower `depth_cv` for growing taxa.

## Not worth pursuing

- **Telling congeners apart by gene order.** Order is conserved within families and is clade-level in the table, so
  a relative's reads fit as well as the species' own (2026-10-01 report, use 5).
- **Single-end overhangs onto the neighbour.** The part of a read past a gene's end is mostly intergenic DNA, which
  the database does not hold.

## Follow-up (2026-10-02): uses 2 and 6 implemented

Strain haplotypes across genes from long reads (2) and per-gene link evidence (6, "foreign genes") are implemented,
together with species lines in `gene_neighbours.tsv` that use 6 needs (without them a species' own unusual gene order
reads as foreign): [2026-10-02-phasing-and-foreign-genes](2026-10-02-phasing-and-foreign-genes/README.md).
