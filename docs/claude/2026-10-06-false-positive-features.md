# Features against false positives in complex communities: what protal collects, and what it could

**Basis.** Code at `1750475` (audit-fixes, with the uncommitted performance changes of the working tree, which touch
no feature): `src/Profiling/Profiler.h` (`TaxonFeatures`, `RecordEvidenceCollector`), `src/IO/SamHandler.h`,
`src/IO/AlignmentOutputHandler.h`, `src/Alignment/AlignmentUtils.h`, `src/Core/{AlignmentStrategy,ChainAnchorFinder,
LongReads}.h`, `src/Hash/KmerLookup.h`, `src/SequenceUtils/GeneNeighbours.h`. Evidence from earlier reports:
[features.md](../features.md) (importances of r226 v10/v9), [the v13 soil errors](2026-10-06-r226-v13-soil/README.md),
[false-positive anatomy](2026-10-03-false-positive-anatomy/README.md), [the fixes](2026-10-03-false-positive-fixes/README.md),
[v9](2026-10-03-r226-v9-evaluation/README.md). No new runs: this is a reading of the code against those results, to
decide which features to build and test next.

## The false positives to aim at

From the v13 soil report: 95-99% of soil's false positives are database species of a genus that has a species the
database lacks in the sample. 59-76% are beside a novel species that is near-identical to them on the marker genes
(`excess_scaled_median` <= 0.015), the rest beside a more divergent one. The final model separates them no better on
its own training rows (in-sample F1 only +0.003-0.013 above species held out), so what is missing is information,
not fit. The marker genes are ~3% of a genome and among its slowest genes; a sister species more than 5% away
genome-wide can be under 1% away on them, and then nothing in a marker read tells it from a strain.

That splits the candidates three ways:

- **the reads' divergence and its shape**: separates the divergent sources (24-41% of soil's false positives), blind
  to the near-identical ones;
- **the reads' ambiguity among the database's taxa, and the database's crowding around the taxon**: partly sees the
  near-identical ones, because their taxa have close relatives in the database too;
- **sequence outside the marker CDS**: the only kind of evidence that can see the species boundary itself.

## Collected features that identify false positives

| rank | feature(s) | group | what separates | evidence | limit |
|---|---|---|---|---|---|
| 1 | `failed_candidate_rate` (ZF) | unfiltered | reads that seeded on the taxon and failed there: a relative lacking from the database | median 0.94 for absent taxa, 0.13 for present ones; 2nd feature of short-read models | near-identical sources rarely fail; long reads rarely fail |
| 2 | `lu_per_kb`, `lsu_per_kb`, `l(s)u_gene_rate{,2,3}` | normalized | the taxon's own unique k-mers on its reads | a quarter to a third of every model's importance | a near-identical sister carries them too; see "single-value cores" below |
| 3 | `excess_median`, `excess_high_share`, `excess_scaled_median`, `identity`, `top_identity` | normalized, divergence | read divergence beyond the base qualities, on the genome's scale | above 0.02 scaled excess almost no taxon is called | flat between 0.01 and 0.03, where real strains and sisters overlap |
| 4 | `sample_log_fragments` | depth | single perfect reads at depth | halved r226's false positives (pe 218 → 122) | identifies a scenario sample when depths are fixed (fixed by varied depths, 2026-10-06) |
| 5 | `relative_skew`, `relative_distance`, `relative_spill`, `relative_close_share` | distance | spill from an abundant close congener | +0.004 pe F1 at r226 | needs the congener in the sample; the soil sources are not in the database |
| 6 | `congener_fit_share`, `low_mapq_share`, `mean_mapq`, `em_own_share`, `em_kept_own_share` | normalized; EM ones opt-in (`relatives`) | reads that fit a congener as well | the only per-sample features that still separate soil's top false positives from true positives at the same divergence (0.13 vs 0.02, 0.23 vs 0.05, MAPQ 58 vs 84; AUC 0.80-0.88) | `em_*` are not in the default set; ZA is capped (below) |
| 7 | `su/lu/lsu_rate_ref`, `uniqueness` | ref (default since 2026-10-06) | a reference with close relatives in the database | `lu_rate_ref` 0.64 vs 0.77 at the same score; +0.002-0.003 soil with boosting | a property of the database, the same in every sample |
| 8 | `mate_lost_share` | divergence | a mate that fits nowhere on the reference | 0.009 pe importance | pe only |
| 9 | `gene_presence_ratio`, `hit_gene_fraction`, `depth_cv` | normalized | reads piled on few genes | 0.01-0.02 importance | near-identical sources cover the genes evenly (top soil FP: 77% of genes) |
| 10 | `third_position_share` | divergence | errors (all codon positions) from divergence (third) | top-3 for ONT | tells error from divergence, not strain from sister |

Of little or no use against false positives: the `adjacency` group (AUC 0.39-0.72, ±0.001 F1), the conservation
pattern ratios (`conserved_*`: the MAPQ filter erases the pattern), `RAF*`, `multiallelic_sites_per_kb`,
`other_genus_fit_share` (rarely non-zero), `linked_share`, `gene_dispersion`, `excess_conserved_fast_ratio`, and
`cluster_ani_radius` (95 for nearly every species). The priors' cluster size gains, but for a reason the simulation
manufactures ([v9](2026-10-03-r226-v9-evaluation/README.md)).

The `relatives` group's `em_own_share` is worth one more look: its ambiguity is what separates soil's top false
positives, but it left the defaults on the design's test set with a forest (v10 follow-up: within noise). It has not
been tried with boosting on soil.

## Query coverage and 5'/3' alignment misses

Neither is a feature, and as protal aligns, neither would carry much:

- **Short reads are aligned end to end on the gene.** A read's ends are free only where it overhangs the gene's end
  (`AlignmentStrategy.h:531-532`: the overhang plus 9 bases; `PostProcessAlignment` then soft-clips them). A divergent
  5' or 3' end inside the gene is aligned with its mismatches and lands in `identity` and the excess features; a read
  that does not fit within the score budget is dropped and becomes a failed candidate (ZF). The soft clips that
  exist measure where the read sits on the gene, not how well it fits. The profiler ignores S and H in every feature
  (`DifferencesAndAligned`, `ReadExcess`), which is right for that reason.
- **Long reads are cut to their marker genes.** Each gene's record is hard-clipped to the gene's window
  (`LongReads.h:294-300`); the read between and around the markers is aligned to nothing. A long read's "query
  coverage" is the marker share of the read: the genome's geometry.
- **What the overhang could tell, if it were compared to something.** The bases of a read past a gene's end are the
  most divergence-sensitive bases protal sees, intergenic sequence that evolves several times faster than the
  markers, and they are clipped unseen because the reference ends at the CDS. That is new feature A below.

## What the aligner computes and throws away

(From a survey of the seeding and alignment code; only what bears on false positives.)

- **ZR is written but never read.** For long reads `ZR:i:1` marks a record whose hit the read's consensus taxon set
  (or raised the MAPQ of), `ZR:i:2` a record that is clearly another taxon's than the read's consensus
  (`LongReads.h:729-730`). `SamFromTokens` parses ZU, ZT, ZA and ZF only. The `ReadConsensus` itself (taxid, parts,
  MAPQ) is never written.
- **ZA is truncated in crowded genera.** It lists at most 4 other taxa, at most 5 edits worse, on the best record
  only (`AlignmentOutputHandler.h:189-215`); and only candidates that were aligned: anchors beyond `align_top` (3 plus
  ties, `AlignmentStrategy.h:706`) are never attempted and appear neither in ZA nor in ZF. In a soil genus with ten
  close database species, `congener_fit_share` and the EM see a few of them.
- **Seed-level ambiguity is discarded.** The number of taxa and genes a read's seeds reach, their anchors'
  `total_length` (exact-match bases), seeds per anchor, the read's k-mers with no index hit: all only in run-wide
  totals (`ChainAnchorFinder.h:152-178`, `Classify.h:149-162`).
- **ZF does not say why a candidate failed**: the k-mer screen, WFA2's budget, the proxy-ANI floor and the output
  filter look the same (`AlignmentStrategy.h:563-580, 634-637, 717`).
- **MAPQ mixes paralogs with congeners.** Single-end MAPQ takes the second distinct candidate, which may be the
  same taxon on another gene (`AlignmentUtils.h:443-488`), so `low_mapq_share` counts paralog ties as ambiguity;
  ZA lists other taxa only, so `congener_fit_share` does not.
- **Single-value cores never count as unique.** A k-mer whose core has a single value in the index (no flex cells)
  is emitted with `unique` and `two` false (`KmerLookup.h:240-246`): the read's flex part cannot be checked there, so
  the most specific k-mers of the database never reach ZU/ZT, `lu_per_kb` or the gene rates. The alignment checks
  those bases afterwards; whether `GetUniqueKmerCounts` (the `*_rate_ref` side) counts them is worth checking, since
  observed and expected rates should count the same k-mers.

## New features, ranked

By what they could do for soil's false positives against what they cost. "Class" is the false positives they aim at:
near-identical (59-76%) or divergent sources (24-41%).

| | feature | class | what it measures | data | cost |
|---|---|---|---|---|---|
| **A** | **marker flanks**: `flank_excess`, `flank_unaligned_share`, `mate_in_flank_share` | near-identical | genome sequence beside each marker (±300 bases for short reads, a few kb for long reads), used to extend alignments, not seeded: the overhang's and the mate's divergence where a sister species differs by its genome's distance, not its markers' | `gene_positions.tsv` already places every gene in every genome at hand at build; flanks of all representatives need their genomes | the largest: build input (all representative genomes), 2.6-3.6 GB of 2-bit sequence at ±300 bases for 143k species of 120-168 genes (no index growth), ends-free extension into the flank, new record fields. The simulation is honest about it: reads come from whole genomes |
| **B** | **long-read consistency from ZR**: `consensus_set_share` (records whose hit the read's consensus set, `ZR:i:1`), `inconsistent_read_share` (the taxon's reads with a segment clearly another taxon's, `ZR:i:2`) | near-identical (pb, ont) | a novel species that is near-identical to S on some markers and to S2 on others splits a long read between them; a present taxon wins its genes on their own | written already; parse ZR | small: a tag parser and two counters |
| **C** | **long-read gene spacing**: `neighbour_gap_deviation` | near-identical (pb, ont) | the distance on the read between two consecutive marker genes against that species' own gap in its reference; intergenic indels differ between species more than between strains | gene positions per species are in `gene_positions.tsv` (not loaded by a run); the read's records give the read's gap | small to medium; only species with a genome at the build have positions |
| **D** | **seed-level crowding**: `crowded_read_share` | both | the share of the taxon's reads whose seeds reached ≥ N taxa with an anchor near the best one's exact-match length: ambiguity ZA cannot show past its 4 alternatives and `align_top` | computed in `ChainAnchorFinder`, discarded; a new integer tag | small |
| **E** | **ambiguity beyond the references' distance**: `congener_gap_deficit` | both | for each record with a congener in ZA, its edits more against those expected from the two references' distance on that gene (`CongenerDistances::Compared` has them gene by gene); separates crowded genera, where a true taxon's reads are ambiguous too, from reads that are not the taxon's | ZA + the pair distances already loaded for `relative_distance` | small; censored at ZA's 5 edits |
| **F** | **paired-read split**: `mate_split_share` | both (pe) | fragments whose mates' best records are on two species of the genus: one novel species split across two references | both mates' records are in the SAM | small |
| **G** | **gene-to-gene divergence heterogeneity**: `excess_scaled_gene_dispersion` | divergent; maybe near-identical | spread of the per-gene scaled divergence (robust CV, or the share of hit genes above twice the median): a mosaic of recombined and diverged genes against a strain's divergence in proportion to the genes' factors | per-record scaled excess is collected; needs the gene id beside it | small; sign unknown, the trees do not need one |
| **H** | **within-gene breadth ratio**: `breadth_ratio` | divergent, shallow | covered bases over those expected at the gene's depth (as inStrain's breadth / expected breadth): a relative aligns where the gene is conserved | per-position coverage is kept with the strain data (`CalculateCoverageVector2`), or from the records' intervals | small |
| **I** | **per-gene failures**: `failed_gene_share` | divergent | the share of the taxon's hit genes on which more reads failed than aligned; a sister's fast genes fail, a strain's do not | ZF names taxa only: add the gene, or count failures per gene family | medium (tag format) |
| **J** | **fixed and polymorphic sites**: `fixed_difference_rate` (scaled by the genes' factors), `polymorphic_site_rate` | divergent; mixtures | the consensus' divergence from the reference, free of sequencing error (long reads above all), and sites where two alleles are both common (a present taxon plus a relative's spill) | the variant handler has the allele counts | small to medium |
| **K** | **gene complementarity with a congener**: `congener_gene_overlap` | divergent | whether the taxon and its most ZA-linked congener hit the same gene families less often than their fragments predict: one novel species split across two references | hit gene families per taxon (≤ 120 bits) | small; needs both with several genes |
| **L** | **database neighbourhood**: congeners within 0.01 / 0.02 / 0.05 marker distance, the nearest one's distance | near-identical | crowding around the reference, beyond the k-mer uniqueness; a per-species constant computed at build | the build's sketches | small; bounded like any per-species prior (~+0.003 in the v13 propensity experiment) |

What these can and cannot do: B-F and L are the cheap ways at the near-identical class, and they can only see the
part of it where the novel species is not equally near all its database relatives. A is the one that can see the
boundary itself in every read that crosses a gene's end; it is the light version of v13's "genome-wide sketch"
recommendation (#4 there), staying within protal's marker-centred design, and would also help the divergent-strain
misses (a strain 2-3% from its reference on the genome matches its flanks; a sister does not). G-K mostly sharpen the
divergent class, which the excess features already handle well.

Already tried and not worth repeating for soil (v13 report): the sample's and genus's context (the sample's taxa,
its low-identity share and median identity, the taxon against its genus's best), the GTDB priors, a soil-only model,
cuts by distance or fragments, per-sample cuts.

## How to test them

1. **Keep two soil points' SAMs** of the next build (pe and pb, one hold-in and one hold-out sample each, with the
   truth of `sim/`): the v13 SAMs were on a node's local SSD.
2. **B, D-H, J and K offline**: compute them from the SAMs, the database and the per-record tags in a script, join
   them to the training table's rows, and refit with `soil_experiments.py` (boosting, depth rounded) as the v13
   follow-up did. B, E, F, G and K need only the SAM and the database; D needs the new tag (a protal change and a
   re-alignment of those samples).
3. **A first as a measurement, not a feature**: for the top false positives of `top_false_positives.tsv` and as many
   true strains at the same divergence, take the reads crossing a gene's end and align their overhang to the
   reference genome's flank (the build has the genomes of the simulated species). If flank divergence separates
   them where marker divergence does not, build it.
4. **The single-value cores**: count how many of a sample's anchors contain one, and whether `lu_genome` counts
   them; if many, flag the k-mer unique after alignment when its bases match, and retrain (all unique-k-mer features
   shift).

## Follow-up: the single-value cores

Asked whether k-mers whose core occurs once in the index should count as unique (ZU): yes, and they did not. A core
with fewer than two values (`Seedmap::m_flex_threshold`) has no flex cells, so the lookup compares only the read's
15-base core with the entry, not the 8 flex bases on either side; `GetFromLookup` therefore emitted such seeds with
`unique` and `two` false, though `--build`'s uniqueness check had read the entry's whole 31-mer back from its gene and
cleared the flag of any shared with another taxon (`Build.h`, "single entries read back from their genes"). On the
reference side `unique_kmers.tsv` counts them as `short_unique` (`su_genome`, `su_rate_ref`), apart from the long
ones that `lu_rate_ref` and the gene rates use, so read and reference sides agreed, but the reads' strongest
evidence of a species never counted. How many there are:

| database | short unique (single-value cores) | long unique | long unique, distance two |
|---|---|---|---|
| GTDB r226 v13 (`local/v13/index_and_package.log`) | 5,117,631 (0.24%) | 2,121,759,296 | 1,748,639,322 |
| r226 v13 training database | 5,654,769 (0.28%) | 1,984,166,534 | 1,639,161,168 |
| mini database (3 species, `build_mini_db.sh`) | 56,859 (85%) | 9,919 | 5,059 |

At GTDB scale nearly every 15-base core has several values, so the change moves the r226 features by a fraction of a
percent; a small database's unique k-mers are mostly short, and its `lu` features change entirely.

Now (branch `fp-features`): the lookup keeps the entry's flag for a single-value core and marks the seed `single`
(`LookupResult::single`, carried through the seed sort key); the anchor compares the seed's whole k-mer with the
gene, on the anchor's strand, and counts it in ZU and ZT only if it matches (`ChainAnchorFinder::CountUniques`,
`WholeKmerMatches`; a single-value core has no other value, so it is unique at distance two too). The reference
side counts them with the long ones (`Gene::HasWholeUniques`, `lu_genome` and `lsu_genome` plus `su`, the gene
rates' denominators), so the observed and expected rates count the same k-mers. Models must be retrained.

## Follow-up: B, D, E and F-L implemented

On branch `fp-features` (worktree `../protal-fpfeat`, from `cdce3c0`), at the user's request: the features go into the
database build and the profiling, and the default feature set. A (marker flanks) and C (long-read gene spacing) are
not part of it. What each became:

| | feature(s) | group | where it comes from |
|---|---|---|---|
| B | `read_consensus_share`, `read_inconsistent_share` | consistency | the long reads' `ZR` tag, written before but never read; now a field of `SamEntry` (`m_settled`), parsed and counted per record (`RecordEvidence::settled_by_read`, `settled_inconsistent`). A segment both settled and inconsistent is written `ZR:i:2` (both tags before) |
| D | `seed_crowding` | consistency | new tag `ZN:i:<taxa>` on each mate's primary record and each long-read segment's best record: the taxa with an anchor at least 0.8 as long (exact-match bases) as the read's longest (`SimpleAlignmentHandler::CrowdedTaxa`, `LongReadSegment::crowding`); the feature is the mean log2 over the taxon's records, summed as integers |
| E | `unexpected_congener_fit_share` | consistency | ZA's congener alternatives against the two references' distance from `species_neighbours.tsv` (L): P(Poisson(d x aligned) <= edits more) < 0.01 (`UnexpectedFit`); computed per record in `NoteRecord`, so no sketching at run time. A congener not in the taxon's list is taken at the list's farthest distance (or 0.15), a lower bound |
| F | `split_fragment_share` | consistency | the fragment's (pair's, long read's) best records: `FinishLink` now runs without gene neighbours too and counts, per taxon, its fragments and those with another species of its genus |
| G | `gene_divergence_dispersion` | shape | per-gene record sums (`RecordEvidence::gene_records`: records, differences, aligned, expected errors), Pearson chi-square around one factor-scaled divergence per degree of freedom (`DivergenceDispersion`) |
| H | `breadth_ratio` | shape | the strain data's coverage vectors of the kept records: covered over L(1 - e^-c) per gene (`Taxon::SiteRates`) |
| I | `failed_gene_share` | shape | `ZF` entries now `taxid:gene` (the gene of the taxon's longest attempted anchor; old `taxid` entries still parse, gene 0); records' ZF genes appended as 64-bit keys per chunk and counted once (`FoldFailedCandidates`), for taxa with records only: no hash map per chunk. The header counts of reads that aligned nowhere stay per taxon, so these failures are those of reads that aligned elsewhere |
| J | `fixed_difference_rate`, `polymorphic_site_rate` | shape | the strain data's alleles and coverage at sites of 4 reads or more (`Taxon::SiteRates`) |
| K | `congener_gene_overlap` | consistency | in `ApplySampleContext`: the congener the taxon's ambiguity classes name most often, the two taxa's hit genes among the genes both references have, (both + 0.5) / (expected + 0.5) (`GeneOverlap`) |
| L | `db_congeners_01`, `_02`, `_05`, `db_nearest_congener` | neighbourhood | new database table `species_neighbours.tsv` (`SpeciesNeighbours.h`), written by `--build` after the suspect-copy scan (`WriteSpeciesNeighbours`: every two species of a genus, `ReferenceSketch` + `SketchedTaxonDistance`, the run's own distance, in batches of 10,000 species), packed into `database.protal`, loaded by a run (`LoadSpeciesNeighbours`); -1 without it |

Order independence: the per-gene sums are integers and are walked in gene order; the partner of K breaks ties by
taxid; the site rates walk genes in order. The trainer's groups are `consistency`, `shape` and `neighbourhood`
(`model_features.py`), all three in `DEFAULT_FEATURE_SET`; the named sets add them one at a time for ablations.

**Tests.** C++ unit tests: 385 of 385 pass (13 new in `tests/test_ComplexCommunityFeatures.cpp`: the seed flag
through the sort key, ZF genes and their parsing, ZR/ZN round trip, crowding, the Poisson rule, the neighbours
table, the build's distance against the run's, dispersion, failed genes, gene overlap, and three profiles of SAMs
that check every new feature's value, on 1 and 3 threads). End to end on a mini database built with this code
(`species_neighbours.tsv`: 2 of its 3 species are congeners): 132 of 133 passed; the failure was
`MsaSampleSelectionTest`, whose three perfect read pairs the mini database's 0.6.0a forest now scores 0.55 (1 pair
0.34, 2 pairs 0.48), because that forest's unique-k-mer inputs grew with the single-value cores (85% of the mini
database's unique k-mers); the test now uses one pair, and passes. The mini-database Python tests, with
`GtdbBuildTest` (a synthetic GTDB release built, trained on the new default set, protal's scores checked for parity
with the trainer's): 70 of 70 pass, none skipped. The trainer's tests (`scripts/test_model_pmml.py`, with the feature
groups' new expectations): 37 of 37 pass. On simulated reads of the mini database every
paired-end record carries `ZN`, 89 `ZF` entries a gene, 13,040 of 13,080 records `ZU` > 0; PacBio reads: `ZN` on
all 418 records, `ZR` on 252.

**Not measured.** Nothing here says what the features are worth: the next r226 build's tables will (its training
logs the default set; `--features normalized+adjacency+distance+depth+divergence+unfiltered+ref` trains without
them for the comparison, and `--evaluation full` scores the named sets). What to watch: the build's
`Species neighbours` line (pairs compared; at r226 the largest genera set the time), the profile's memory (per-gene
record maps per taxon per chunk), and whether `breadth_ratio` and the site rates, which read every hit gene's
coverage vector when a taxon is scored, cost noticeable profiling time at r226.
