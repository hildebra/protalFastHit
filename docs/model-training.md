# Training protal's presence model

protal decides whether a species is present with a random forest (PMML, `model_pe.xml` in the
database for paired-end reads; `model.xml` in databases of earlier versions). It scores every species with reads, and reports those whose probability of `TRUE` is at
least `--knob` (0 to 1, default 0.5), or, for a model with knobs by sample depth
([below](#knobs-by-sample-depth)), at least the model's knob for the sample's depth unless `--knob`
is given, or, with `--fdr` and a model with calibrated calls ([below](#calls-at-a-target-share-of-false-calls)),
the sample's highest-scoring species while their expected share of false calls stays at the target.
With `--singleton_congener N` a species of a single fragment beside a congener of N fragments or more
is never reported when its read looks like the congener's ([the singleton rule](#the-singleton-rule));
off by default.

The model's features come from the reads and the database, so a model belongs to the kind of
database it was trained on. [building-a-database.md](building-a-database.md#build-and-train-in-one-command)
runs the whole pipeline below in one command; this page describes its parts.

Features are counted per read, and a pair of reads counts twice, so a model also belongs to the
kind of reads it was trained on. Single-end samples are profiled with `model_se.xml` of the
database (or `--model_se`; [running.md](running.md#single-end-reads)), a model trained the same way
on single-end samples, and PacBio and Nanopore samples with `model_PB.xml` and `model_ONT.xml`
(`collect_training_data.py --read_types`).

## What protal expects of a model

At start, before any alignment, protal checks that the model
- takes only inputs protal computes: the features of `TaxonFeatures` in
  `src/Profiling/Profiler.h`, the columns of the training dump below;
- predicts a field with the value `TRUE`, whose probability is the species' score;
- scores a species.

A model that fails one of these stops the run with exit code 2, naming the problem.

## Training data

With a truth file (`--profile_truth`, one per sample, comma-separated, or a `PROFILE_TRUTH` column
in the map), protal writes `<profile>.truth_annotated`. This file has one row per species with
reads: whether the species was present (`truth`), the current model's call and probability, and
every feature. The features are written exactly as the model is given them.

A truth file names one species present per line: by GTDB lineage in any tab-separated field
(`d__...;s__Genus species`) or by internal taxid in the first field. `simulate_metagenomes
--protal_metafile` writes such files and a map that uses them ([simulation.md](simulation.md)).

`collect_training_data.py` produces such data at scale. It simulates metagenomes with
`simulate_metagenomes` over a grid of read setups (length, ART profile, fragment size) and depths,
profiles them against a database, and joins all dumps into `training_data.tsv`:

    python3 scripts/collect_training_data.py --db DB --genome_table genomes.tsv -o training \
        --archaea 2 --species_per_sample 10-40 -t 8

With `--read_types`, it collects other reads of the same communities too: `se`, the paired-end
samples' first reads alone, profiled as single-end reads (no new simulation); `pb` and `ont`,
long reads of each paired-end sample's community (from its manifest; point i of the long reads
replays the communities of paired-end point i). The collector draws a long-read sample's reads as
they arise until they hold `--long_read_bases`: each read's genome by abundance times genome length,
its length from a gamma distribution of the setup's mean and SD (100 bp to 1 Mb, as pbsim3's), its
start uniform over the genome's contigs of 100 bases or more, cut where its contig ends, either
strand. One read is made of each: PacBio HiFi reads by `scripts/hifi_reads.py`, Nanopore reads by
pbsim3 (`--strategy templ`) with the setup's quality model, in one run per sample. pbsim3 does not
simulate HiFi reads (it simulates the subreads of several passes, for PacBio's `ccs`), and its reads of
one pass with an error model have quality 0 throughout, which made protal's `excess_median` about
−1 for every PacBio taxon until 2026-10-02. `hifi_reads.py` gives each read a quality by its length,
Q50 up to 5 kb, Q30 at 25 kb, linearly between, Q20 at 50 kb and less beyond (and normal noise of the
setup's SD, 3), errors mostly as indels in homopolymers, and base qualities that say how likely each
base is wrong (calibrated: a read's differences are what its qualities expect)
([report](claude/2026-10-02-pacbio-hifi-reads/README.md)). (pbsim3's own sampling, once per genome,
cut every contig's last read to that contig's share of the bases, so that shallow samples were mostly 100 bp reads, one per contig,
[report](claude/2026-10-01-build-profiling/README.md).) All samples are profiled in one protal run, each with its read
type's model and settings (the map's `READ_TYPE` column), and each read type gets its table:
`training_data.tsv` (pe), `training_data_se.tsv`, `training_data_pb.tsv`, `training_data_ont.tsv`.

| Option | Default | |
|---|---|---|
| `--db`, `--genome_table`, `-o` | required | database, simulator genome table (accession, GTDB taxonomy, FASTA path, optionally the genome's length, which spares the simulator reading every genome for it), output directory |
| `--protal`, `--simulator` | on `$PATH` | the binaries |
| `--samples` | 4 | samples per design point |
| `--read_pairs` | `1000,5000,20000,100000,500000` | depths, one design point each |
| `--read_setups` | `100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50` | length : ART profile : fragment mean : fragment SD, one design point each; the profile `file=R1.txt+R2.txt` (or `file=P.txt`) uses quality profiles of `art_profiler_illumina` instead of a built-in one |
| `--species_per_sample` | `5-30` | N or MIN-MAX |
| `--strains_per_species` | one strain | probabilities of a second, third, ... strain of a species, e.g. `0.3,0.1` |
| `--abundance` | the simulator's (Poisson-lognormal, sigma 1.3) | `lognormal:SIGMA`, `powerlaw:ALPHA` or `negbin:R:P` |
| `--read_types` | `pe` | `pe`, `se`, `pb`, `ont`, comma-separated |
| `--long_read_bases` | `300000,1500000,6000000,30000000,150000000` | bases per long-read sample, one design point each |
| `--pb_setup` | `hifi:15000:3000:3` | `hifi` : length mean : length SD : SD of the reads' quality around their length's (`hifi_reads.py`); or a pbsim3 setup as `--ont_setup`'s |
| `--ont_setup` | `qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36` | pbsim3 method : model : length mean : length SD : accuracy mean (: error mix, for qshmm) |
| `--pbsim`, `--pbsim_models` | `pbsim`, its data folder | pbsim3 and the folder of its `.model` files (for ont, and pb with a pbsim3 setup) |
| `--archaea` | 0 | archaeal species per sample |
| `--congeners` | 0 | relatives that share a sample: `SHARE:MIN-MAX` (e.g. `0.25:2-5`), about SHARE of each sample's species in groups of MIN to MAX species of one genus, the genera drawn per sample (`simulate_metagenomes --congener_groups`); `N`, N species of one genus in every sample of a design point (the genus drawn per point); 0, none. Relatives share real samples at very different abundances, but hardly ever uniform draws from many genera; the relatives features need them in training ([below](#features)) |
| `--novel_species` | | species the database lacks (e.g. those a training database leaves out), optionally with the rank they were held out at and the clade (`heldout_species.txt`), for `meta_novel_*` |
| `--novel_clades` | 0 | species of held-out clades (ranks above species in `--novel_species`) in every sample, per rank: each design point takes one clade of each rank, in turn |
| `--taxonomy` | | the database's `internal_taxonomy.dmp`, for `meta_rep_genome`, `meta_relative_rank`, `meta_novel_level` and `meta_neighbour_rank` |
| `-t`, `--seed` | 4, 1 | threads of the protal run; seed |
| `--jobs` | `-t` | paired-end design points, and long-read samples, simulated at a time (ART simulates one genome at a time) |
| `--simulate_only` | | simulate and stop; a run without it profiles (the simulations need no database, so they can run while it is built; `build_gtdb_database.py` does so) |

The paired-end design points are simulated in parallel, then the long-read samples, then all the
samples are profiled in one protal run, which loads the database once. A rerun resumes: a design point already simulated from the same
inputs (simulator, genome table, seed and design) and profiled against the same database with the
same protal is skipped; one made from other inputs is simulated or profiled again (its folder holds
the key of what made it, `simulated.json` and `profiled.json`). Columns it adds start with
`meta_` and say where each row comes from: design point, sample, read type (`meta_read_type`), read
length (the mean for long reads), depth (`meta_read_pairs`: read pairs, reads for se, bases for pb
and ont), the taxon's
domain (`meta_domain`, from the genome table's lineages; `unknown` for species the table lacks),
how many species of `--novel_species` the sample holds (`meta_novel_species`), whether the taxon
shares a genus with one of them (`meta_novel_congener`: the taxa their reads land on), and whether a
present species was simulated from its representative, the database's reference, or another
genome (`meta_rep_genome` 1 or 0). With the ranks in `--novel_species` and `--taxonomy` also: the
sample's novel species by the rank they were held out at (`meta_novel_levels`, e.g.
`species:2,family:1`), the deepest rank a taxon shares with a species simulated in the sample
(`meta_relative_rank`; `species` for those species themselves), and, for an absent taxon whose
closest simulated species is a novel one, the rank that one was held out at (`meta_novel_level`):
its reads are the likely source of the taxon's. For a present taxon, `meta_neighbour_rank` is the
deepest rank it shares with another species of its sample: a congener's reads fit it nearly as
well, so it may be missed. The training report breaks its errors down by these.

The genome table should include species that the database lacks. Their reads land on relatives
the database has, and those species are the false positives the model must learn to reject: build
the training database with some species left out (`gtdb_to_protal_db.py --exclude_species`, which
keeps their taxids, so the model applies to the full database) and simulate from all of them.
`build_gtdb_database.py --holdout` does this, and `--holdout-clades` leaves whole families,
classes and phyla out too (below). Include archaea,
too (`--archaea`), since they have fewer marker genes than bacteria. Include genomes other than the
database's references (other strains of its species): real strains differ from the reference by up
to a few percent, and a model that has only seen reads of the reference itself may call them
absent. `scripts/mini_db/simulate_gtdb_release.py --strain_divergence 0.002-0.015
--species_divergence 0.015-0.04` makes a small world with such strains and close relatives.

## Features

The shipped model was trained on absolute counts: genes, k-mers and mates. These depend on the
database, on the domain (archaea have 52 marker genes, bacteria 119), on depth and on read length.
On simulated data it finds about a third of the archaea present, with probabilities pinned near 0.5.

`model_features.py` lists `NORMALIZED_FEATURES`, which do not depend on these:
- fractions of the hittable genes (`hit_gene_fraction`, `gene_presence_ratio`: genes hit against
  the number expected from the number of fragments);
- rates per aligned or covered kb (`lu_per_kb`, `lsu_per_kb`, `variant_sites_per_kb`,
  `multiallelic_sites_per_kb`);
- read identity with indels counted as differences (`identity`, `top_identity`);
- the share of low-identity reads, i.e. reads of relatives (`low_identity_share`);
- the depth's coefficient of variation across genes (`depth_cv`), uniqueness and allele frequency
  classes (`uniqueness`, `RAF0`-`RAF4`, the `*_gene_rate*` columns);
- `fragments` and `depth`, which say how much evidence there is;
- what the reads' other candidates and other records say: `mean_mapq` and `low_mapq_share` (the
  share of the taxon's reads with MAPQ below 10, i.e. nearly as good a second candidate);
  `congener_fit_share` and `other_genus_fit_share` (the share of its reads that another species of
  its genus, or a species of another genus, fits within one edit: the reads of a relative the
  database lacks fit several of its congeners about as well); `linked_share` (the share of its
  reads with two records on it: both mates of a pair, on one gene or two, or two genes of a long
  read). These count every read's best record, also those the profiler's MAPQ filter leaves out;
- how far the reads differ from the reference beyond what their base qualities explain:
  `excess_median` and `excess_high_share`, the median over the taxon's best records (also those the
  filters leave out) of a record's differences per aligned base (X, I and D over M, X, I and D)
  less the mean error probability of its bases (10^(-Q/10)), and the share of records above 0.02;
  both 0 for reads without qualities. A relative the database lacks gives reads that differ by
  several percent more than their errors, a present species' own reads by about one;
- which genes the reads hit, by how fast the genes diverge within species (the database's
  `gene_conservation.tsv`, [database-files.md](database-files.md)): `conserved_fast_depth_ratio`,
  log2 of the median depth of the taxon's hit genes with factor below 1 (conserved) over that of
  its other hit genes, each plus 0.001 (0 if either has none), and `conserved_hit_share`, the share
  of its hit genes that are conserved; 0 and 0.5 without the table. A relative the database lacks
  differs least from its congeners on the conserved genes, so its reads align there best; on the
  simulated benchmark world a taxon that only holds such a relative's reads has conserved genes
  deeper than its fast ones (median ratio +0.33, present taxa about 0), although many of those
  reads fit several congeners equally and fall below the profiler's MAPQ filter
  ([report](claude/2026-10-01-conservation-pattern/README.md); not yet checked on real genomes);
- `conserved_fast_record_ratio`: the same log2 ratio of depths over every best record of the taxon,
  before the MAPQ and length filters (MAPQ 0 included), each depth the records' reference bases
  over the summed length of the taxon's genes of that kind (+ 0.001); 0 without
  `gene_conservation.tsv`. A species' own reads cover both kinds alike; a relative's align as often as
  the species' own on the conserved genes and a fifth as often on the fastest, a drop that the MAPQ
  filter blurs, since it removes most of the relative's reads on the conserved genes;
- `conserved_fast_kept_ratio`: the same over the records the MAPQ and length filters keep (the
  taxon's genes' mapped bases); how far it falls below `conserved_fast_record_ratio` tells how many
  of the taxon's reads on conserved genes fitted several taxa equally.

  The four came after 0.7.1: refitted with them on the 0.7.1 benchmark's training tables, the forest
  gained 0.003 (paired-end) and 0.007 (ONT) of test F1, not significant on their own, and the
  cross-validated F1 rose too ([report](claude/2026-10-01-f1-opportunities/README.md)). A dump of
  an older protal lacks them.

The dump also has `adjacent_expected_share`, `adjacent_unlikely_share` and `adjacent_support`: of
the genes next to each other on a taxon's reads (a pair's mates on two genes, a long read's
consecutive genes within 3 kb of each other on the read; every read's best records, also those the
filters leave out), the shares whose ends are expected neighbours in the taxon's clades (facing each
other in 20% of their species or more, a sparse clade leaning on those above it) and unlikely ones
(5% or less), by the database's gene neighbours
([running.md](running.md#options-the-website-does-not-list); both 0 without them); and the mean
frequency of their pairings there, pulled towards 0.5 by one pairing ((sum of the frequencies + 0.5)
/ (pairings judged + 1); 0.5 without gene neighbours or reads across genes).
The models train on them by default (`ADJACENCY_FEATURES`, with `NORMALIZED_FEATURES` the set
`normalized+adjacency`, in the default `normalized+adjacency+distance+depth+divergence`). On a synthetic world they
changed the test F1 within noise (paired-end +0.0025, PacBio −0.0026, Nanopore +0.0015 over three seeds;
[report](claude/2026-10-01-gene-neighbours-run/README.md),
[2026-10-02](claude/2026-10-02-gene-neighbour-frequencies/README.md)): reads of a congener the
database lacks pair across the same genes as the species' own. Whether they help where gene order is
real, a build tells with `--features normalized` (without them) against the default, and every
`--evaluation full` training compares the feature sets with species held out ("Feature sets" in the
report). A model trained with them, used with `--no_gene_neighbours` or on a database without gene
neighbours, sees every taxon as one without reads across genes (0, 0 and 0.5). A database whose
`gene_neighbours.tsv` has species lines ([building-a-database.md](building-a-database.md#gene-neighbours))
judges a species' own unusual gene order expected, so its present taxa have more expected and fewer
unlikely pairings than with the clades' lines alone: train on the table the model will run with. The
`depth` feature leaves out foreign genes ([running.md](running.md#options-the-website-does-not-list)),
unless `--keep_foreign_genes`. Since 2026-10-02 it counts a fragment's bases on a gene once (overlapping
mates' overlap once, [report](claude/2026-10-02-fragment-depth/README.md)): with the collector's read
setups the depth of present taxa is 0-3% lower than before and no call of the build's test samples
changed, but a model is best trained on dumps of the protal it runs with.

The dump also has features from the other taxa of the taxon's sample, after the amplicon denoisers
(AmpliconNoise, DADA2, UNOISE), which judge a candidate sequence against the more abundant ones that
could have produced it ([report](claude/2026-10-02-amplicon-denoising/README.md)). protal computes them
once a sample's reads are in (`MicrobialProfile::ApplySampleContext`, `src/Profiling/SampleContext.h`),
from every taxon's fragments, the taxonomy and the references:
- `genus_skew`, log10((n + 1) / (m + 1)) for a taxon of n fragments whose genus's most abundant other
  species has m (UNOISE's skew); `family_skew`, the same against the species of most fragments of
  another genus of its family; `genus_share`, its share of its genus's fragments;
- `genus_spill`, log10((n + 0.5) / (0.001 × its congeners' fragments + 0.0001 × the family's other
  genera's + 0.5)): its fragments over those its relatives would spill onto it, by rank alone;
- by the distance of the references: of its congeners with more fragments (the four most), the one
  whose fragments times the spill rate 0.01 × 10^(−d / 0.05) are largest, d the median Mash distance
  (k = 12, from sketches of 64 hashes per gene) of the two references over the genes they share (10 or
  more, else d = 1); `relative_skew` against that congener, `relative_distance` its d (1 without one),
  `relative_spill`, log10((n + 0.5) / (the congeners' summed spill + 0.5)), and
  `relative_close_share`: of the taxon's fragments on the genes it shares with that congener, the share
  on the half where the two references are most alike, over that half's share of the genes' length.
  A species' own reads fall on its genes by their length (about 1); reads that spilled over from the
  congener land where the two are most alike (up to 2);
- `em_own_share` and `em_kept_own_share`: the share of the taxon's best records (all, and those the
  MAPQ and length filters keep) that an abundance-weighted assignment leaves to it, an EM over the
  reads' alternatives (`ZA`) in which a read fitting another taxon m edits worse goes to it in
  proportion to that taxon's records times (θ / (1 − θ))^m, θ the other taxon's own differences per
  aligned base in the sample (0.002 to 0.2): the error rate and strain divergence its reads show. A read
  that fits a congener a hundred times as abundant two edits worse is mostly the congener's.

They let a model call a thin taxon with no relative in the sample and reject one beside an abundant
congener whose reads it holds (`RELATIVE_FEATURES`, with the normalised and adjacency features the
set `normalized+adjacency+relatives`, opt-in). On the benchmark world they raised cross-validated
F1, AP, log loss and the F1 at each test set's best threshold, but changed the test F1 at the knobs protal
calls with by +0.003 to −0.005 and missed 1.6-2.3 times as many minor congeners (present species beside a
present congener of ten times their fragments); `normalized+adjacency+distance`, with only the four
`relative_*` features by the references' distance, was the best paired-end set at every knob (+0.001 to
+0.006), single-end within noise, still missing more minor congeners
([report](claude/2026-10-03-denoising-implementation/README.md)). At GTDB r226, where false positives are
the larger problem, the four distance features added 0.004 paired-end F1 on the test set at the knob curve
(0.007 at knob 0.5), single-end and long reads within noise; all the relatives features lowered the
paired-end log loss by 27% with species held out but did no better at the knob curve
([report](claude/2026-10-03-r226-v5-v6-training/README.md)). So the four distance features are in the
trainer's and the build's default set (`normalized+adjacency+distance+depth+divergence`,
[below](#the-samples-depth-and-the-divergence-features)); a table of a protal before them needs `--features
normalized+adjacency`.
They need training samples whose congeners share a sample (`--congeners SHARE:MIN-MAX`,
`build_gtdb_database.py`'s default): on the r226 tables, whose species were drawn uniformly, a model
with the rank features learnt that an abundant congener means absence and missed 69-73% of the present
species beside a congener ten times as abundant, against 21-29% without them. `genus_top_fragments`
(the fragments of the genus's most abundant other species) is in the dump for the singleton rule, not
a feature of these sets.

### The sample's depth and the divergence features

`SAMPLE_FEATURES` (`depth`): `sample_log_fragments`, log10 of the sample's fragments over all its
taxa, the number protal reads a knob curve at. No other feature says how deep a sample is, and what a
taxon of one perfect read is worth depends on nothing else: at GTDB r226 half of the paired-end false
positives were single reads at identity 0.97 or more, equal to true single-read species in every
feature and apart from them only by the sample's depth (true ones in samples of 50,000 pairs or fewer,
false ones in samples of 100,000 or more). As a feature it halved the false positives at knob 0.5
(pe 218 → 122, se 231 → 118) and raised the test F1 by 0.004 (pe) to 0.009 (se) over the distance set,
above its knob curve, PacBio +0.005, Nanopore the same ([report](claude/2026-10-03-false-positive-anatomy/README.md)). A knob curve fitted on top of
it corrects twice and loses (pe 0.9653 → 0.9555), so with the feature in the set the trainer fits none
even with `--depth-knobs`, and protal calls at `--knob` (0.5). The forest cannot extrapolate past the
deepest training sample any more than the curve could: train at the depths you profile.

`DIVERGENCE_FEATURES` (`divergence`), from the same report: `excess_scaled_median`, the median over
the reads of their divergence beyond their base qualities divided by their gene's conservation factor
(`gene_conservation.tsv`; 1 without), the genome's divergence from the reference as the species
definition (95% ANI) measures it rather than the marker genes', which are the most conserved part of
the genome and compress a 6-10% divergence into 2-4%; `excess_conserved_fast_ratio`, log2 of that
divergence on the genes of factor below 1 over the others (a species' own reads differ by the genes'
factors, a relative's reads that align only where a gene is conserved do not; 0 with fewer than 200
aligned bases on either kind); `third_position_share`, the share of the mismatches at third codon
positions of the reference, whose genes are coding sequences in frame (sequencing errors fall on the
three positions alike, a strain's differences mostly on the third, synonymous one; 1/3 without
mismatches); and `mate_lost_share`, of the paired fragments whose mate was expected on the taxon
(both mates with a record on it, or one kept record with room for the fragment inside its gene, the
room judged against the length within which 95% of the sample's fragments with both mates on one gene
lie), the share whose mate has no record on the taxon (a read of a relative that fits the reference
where it is conserved has a mate that fits nowhere on it; 0 for single-end and long reads). The
`linked_share` of the dump, both mates kept on the taxon whatever the room, stays out of the sets: it
gained 0.001 on the r226 tables. A dump of a protal before these features lacks the five columns:
train it with `--features normalized+adjacency+distance`. How they did on the benchmark world is in
[the implementation report](claude/2026-10-03-false-positive-fixes/README.md).

protal writes the alternatives as the `ZA` tag of a read's best record (`ZA:Z:<taxid>:<edits
more>,...`, the other taxa among the read's aligned candidates with at most 5 edits more, or `*`);
a SAM file of an older protal lacks it, protal warns, and the two fit shares are then 0, so such
files are aligned again (`--force`) for a model that uses them. A dump of an older protal lacks the
five columns, and the trainer stops with `--features normalized`; one without the gene neighbour
features stops with the default and trains with `--features normalized`.

## Training

    python3 scripts/random_forest_cmdline.py --truth-file training/training_data.tsv \
        --output-prefix training/model

The trainer needs Python 3 with numpy, pandas, joblib and scikit-learn; no Java. `model_pmml.py`
writes the forest as PMML itself, so that protal's probabilities equal scikit-learn's bit for bit:
scikit-learn compares a feature as float32 with its thresholds and cPMML as a double, so each
threshold is written as the largest double that rounds to a float32 at or below it; leaf counts are
written so that cPMML's count over total is scikit-learn's leaf probability; and the trees are
averaged in file order, as protal does. The trainer checks this on every training row and fails if
one differs.

| Option | Default | |
|---|---|---|
| `--truth-file`, `--output-prefix` | required | the training table; the prefix of the outputs |
| `--features` | `normalized+adjacency+distance+depth+divergence` | the feature groups joined by `+`, `normalized` among them: `normalized` (`NORMALIZED_FEATURES`), `adjacency` (the gene neighbours'), `relatives` (all `RELATIVE_FEATURES`) or `distance` (the four `relative_*` by the references' distance), `depth` (the sample's depth, `SAMPLE_FEATURES`; no knob curve is fitted with it), `divergence` (`DIVERGENCE_FEATURES`); `all`: every feature column of the dump. A table of an older protal lacks columns: `normalized+adjacency+distance` before the depth and divergence features, `normalized+adjacency` before the relatives features |
| `--reference-pmml` | | train on the input fields of an existing model instead |
| `--ntree`, `--maxnodes`, `--min-samples-leaf`, `--max-features` | 64, 256, 1, `sqrt` | the forest (`--maxnodes 0`: no limit on leaves; `build_gtdb_database.py` gives 512 for short reads and 128 for long reads: at GTDB r226 512 leaves gave short reads a lower log loss and fewer false positives, 128 long reads a lower log loss at the same F1) |
| `--knob` | 0.5 | the threshold protal will use; calls and their errors are counted at it |
| `--depth-knobs` | off | also fit a knob curve over the sample's depth and store it in the model ([below](#knobs-by-sample-depth)); `build_gtdb_database.py` passes it for every read type (`--depth-knob-read-types`). Not fitted when the sample's depth is a feature (the default set): the report says so |
| `--fdr-calls` | off | also calibrate the scores and choose a target share of false calls per sample, and store both in the model ([below](#calls-at-a-target-share-of-false-calls)); `build_gtdb_database.py --call-mode fdr` passes it |
| `--singleton-congener` | 0 | the singleton rule protal applies (its `--singleton_congener`; 0, the default: none): every call the trainer counts leaves out the rows it vetoes ([below](#the-singleton-rule)) |
| `--folds` | 5 | folds of the held-out evaluations |
| `--evaluation` | `full` | `basic`: the held-out evaluations only; `none`: fit and export only |
| `--previous-procedure` | off | also compare with the procedure this trainer used before (below), unless `--evaluation none`; `--no-previous-procedure` is the default |
| `--taxonomy` | | the database's `internal_taxonomy.dmp`: domains the table's `meta_domain` lacks, and the lineages for holding out whole clades |
| `--test-file` | | an independent test table (the collector with another design and seed), scored by the fitted forest: the report's section "Independent test set" (metrics, FN and FP rates by depth and by rank, the threshold with the highest F1 there), `<prefix>.test_predictions.tsv.gz`, and a warning when it scores clearly worse than cross-validation |
| `--seed`, `--threads` | 1, 4 | |

Rows of one sample share its reads, and rows of one species share its reference, so a random split
of rows scores a model on samples and species it was trained on. The trainer scores each row with
forests that saw neither its sample ("by sample") nor its species ("by species"), and out of bag
(scikit-learn draws each tree's rows by their class weight, so when absent taxa far outnumber
present ones, some present rows are drawn for every tree and have no out-of-bag score; they are left
out of that estimate and counted in the report). By species is the estimate that matters for a large database: of GTDB's ~130,000 species a
training set holds a few thousand, so most species protal meets in real samples were never in
training. With `--taxonomy` it also holds out whole genera, families, orders, classes and phyla (a forest
that saw no taxon of the row's clade): how far the model carries to parts of the tree the training
data barely cover. The report also compares with the model that profiled the training samples (the dump's
`probability`, e.g. the shipped model), and `--evaluation full` adds studies of whether the data
and settings suffice: the other feature set; the number of leaves and of trees; and a learning curve
with fewer training samples. `--previous-procedure` also compares with the procedure this trainer
used before: a grid search over `max_features`, then a 512-tree forest on only the top features,
judged on random rows and on species held out. Its grid search took 82% of the training time
([report](claude/2026-10-01-build-profiling/README.md)), so it is off by default and cheaper than
the procedure's own: every second value of `max_features` and then the two next to the best, 3
folds instead of 5, and at most 20,000 rows (whole samples, drawn at random); the forests it then
judges are the procedure's.

| Output | |
|---|---|
| `<prefix>.xml` | the model |
| `<prefix>.report.txt` | the evaluation (also printed): data summary with warnings, held-out results by domain, depth and evidence, the hardest taxa, the threshold, the studies |
| `<prefix>.metrics.json` | the report's numbers |
| `<prefix>.predictions.tsv.gz` | every taxon's probability out of bag, by sample, by species, by rows and by clade held out, with its main features |
| `<prefix>.thresholds.tsv` | precision, sensitivity and F1 by threshold, species held out |
| `<prefix>.varimp.tsv` | feature importances |
| `<prefix>.joblib` | the fitted scikit-learn forest |

On simulated worlds, 64 trees scored as well as 256 or 512 (the forest is 8 times smaller and
loads faster in protal), the leaf limit did not bind, and the grid search, which took most of the
old trainer's time, chose a few top features and did no better on species held out.

### Species and clades the database lacks

Real samples hold organisms the database lacks at every depth: a species of a genus it has, or of a
genus, family, order, class or phylum it has none of. Their reads land on the closest relatives the
database has, or nowhere, and those relatives are the false positives to avoid. `build_gtdb_database.py`
makes such samples: its training database lacks many single species (`--holdout`, 20%) and whole
clades of every rank (`--holdout-clades`, from phylum down to genus), every sample has species of
held-out clades (`--novel-clades-per-sample`), and the absent taxa their reads reach are negatives
the forest is trained on, like any other row.

The report's section "False positives and false negatives by taxonomic rank" gives rates in % at
the knob, for the model scored with species held out and for the collection model, in rows from
species to phylum:

- **False positives from what the training database lacks**: absent taxa whose closest species in
  their sample is one the database lacks (`meta_novel_level`), by the rank it was held out at: the
  novel species simulated, those taxa, how many are called (FP), the FP rate and FP per 100 novel
  species; and the other absent taxa for comparison. A second table splits the same taxa by the rank
  they share with that species (`meta_relative_rank`).
- **False negatives** of present taxa by the deepest rank they share with another species of their
  sample (`meta_neighbour_rank`): whether congeners, or relatives further up, cost sensitivity.
- **With the taxon's clade held out of training** (the cross-validation above, species to phylum):
  FN rate, FP rate, FP per sample and F1 when the forest has seen nothing of that species, genus,
  family, order, class or phylum.

The summary repeats the FP and FN rates by rank.

The section "The conservation features by class of taxon" (and the same for the test set) gives the
median and quartiles of `conserved_fast_record_ratio`, `conserved_fast_depth_ratio`,
`conserved_hit_share` and `excess_median` for present taxa beside a congener the database lacks
(`meta_novel_congener`) or not, for absent taxa whose closest species in the sample is a held-out
species of their genus (they hold its reads), and for the other absent taxa: whether the features
carry, in this training data, what a relative the database lacks does to its congeners' genes (on
a GTDB build, real genomes).

The section "Strains" (and the same for the test set) is about present species simulated from
another genome than the database's representative (the collector's `meta_rep_genome` 0: a strain, as
most species of real samples are): how often they are missed against those simulated from the
representative, by the taxon's fragments (1-10, 11-100, more), the identity of the strains found and
missed, and the conservation features of the missed strains next to the absent taxa that hold a
held-out congener's reads. A missed strain should look like a found one there, not like a missing
species' congener; at GTDB r226, 79% of the paired-end model's misses were strains, all those with
more than 10 fragments, at a median identity of 0.959
([report](claude/2026-10-02-r226-build-evaluation/README.md)). The summary gives their FN rate.

`check_model_parity.py` re-profiles saved training samples (`--profile_only` on the SAMs a
`collect_training_data.py` folder keeps) with a model and checks that protal's probabilities are
the model file's, and that protal computes the features as it did when the training data was
collected (another protal version may not). protal sums its features so that a sample gives the
same ones on any number of threads, alone or among others; a feature that differs in its last
digits only (relative difference up to 1e-12, as a sum added up in another order does) is noted in
`parity.txt`, and does not fail the check:

    python3 scripts/check_model_parity.py --db DB --model training/model.xml --training training

`--read_type se` (or `pb`, `ont`) checks a model of other reads, on the collection's samples of that
read type (protal is given it with `--model_se`, `--model_pb` or `--model_ont`).

`gradient_boosted_cmdline.py` and `hist_gradient_boosted_cmdline.py` train gradient-boosted trees
with the same inputs; they still export with sklearn2pmml, which needs Java. The histogram variant
is faster, but its PMML is not guaranteed to load in cPMML: check its scores before you use it. The
R scripts (`random_forest_cmdline.R`, `random_forest.Rmd`) are the older caret pipeline the Python
trainer was modelled on.

## Using a new model

Try it first without changing the database:

    protal --db DB --model training/model.xml --map test.map -t 16

or, cheaper, on existing alignments with `--profile_only`. To make it the database's default,
add it with `protal --add_model training/model.xml --read_type pe --db DB` (`se`, `pb` or `ont` for a
model of those reads), which checks it and replaces the database's model of that read type.

Train the production model on simulations from the database's own genomes (for example GTDB), with
held-out species as negatives. Choose `--knob` on held-out samples like the ones it will profile;
`<prefix>.thresholds.tsv` is a start.

### Knobs by sample depth

A taxon's features say how much evidence it has, but not how deep its sample is, and the threshold
that calls best differs with depth: a deep sample holds many more absent taxa with a few reads
(at GTDB r226 they grew as depth^0.84 while the present ones levelled off), so its best threshold
is higher. With `--depth-knobs` the trainer fits a knob curve over the sample's depth, log10 of its
fragments over all its taxa (rows): a point per half decade where the training samples are, placed
at the median depth of that half decade's samples, its knob the threshold (0.05 to 0.95, in steps of
0.01) with the highest F1 on species held out of the samples within half a decade of it (so that
neighbouring points share samples and the curve does not follow each bin's noise), where they have
more than 50 taxa and 10 present ones. A point needs 6 samples there: a half decade with fewer joins
the next deeper one, and one left at the deep end the point before, as protal keeps the last knob for
every deeper sample (at r226 the long-read points of 2 and 4 samples had their best knob anywhere from
0.33 to 0.91 over resamples of their samples, while the short reads' point of 6 at 10M read pairs held
at 0.68-0.87, below the 2M point's 0.94; [report](claude/2026-10-02-r226-v3-training/README.md)). A point keeps `--knob` unless its best knob gains
0.002 of F1 on its window (`DEPTH_KNOB_MIN_GAIN`): knobs that gained less varied between fits and lost
on test sets ([report](claude/2026-10-03-r226-v4-v2-rerun/README.md)). With the sample's depth among the
features (the default set, [above](#the-samples-depth-and-the-divergence-features)) no curve is fitted at
all. The report's section "Knobs by sample depth" lists the points
with their F1 at `--knob` and at their knob; on the test set (`--test-file`) it gives the F1 and the
errors at the curve, as protal calls by default, next to those at `--knob`. The curve goes into the
model's header as `<Extension name="protal_depth_knob_curve" value="1.300:0.12,4.320:0.92"/>`.

protal sums each sample's fragments over its taxa in the same way and reports the sample's taxa at
the curve's knob for that depth: linear in log10 of the fragments between the points, and the first
or last point's knob beyond them, so a sample deeper than any trained keeps the deepest knob (train
at the depths you profile: `build_gtdb_database.py` has points up to 10M read pairs and 6 Gb of long
reads). `--knob` on the command line applies to every sample instead. The log lists a model's curve
and each sample's knob (`Sample S: N fragments, knob K (the model's for that depth)`); the training
dump's `prediction`, the statistics and, unless `--msa_knob` is given, which samples enter the strain
MSAs follow the sample's knob. At GTDB r226 (0.7.0's features, design points up to 500,000 read
pairs) the best threshold went from ~0.1 at 1,000-5,000 read pairs to ~0.9 at 500,000, and thresholds
by depth chosen on the training rows raised the test sets' F1 out of sample by 0.006 (pe; 59% fewer
false positives), 0.009 (se), 0.033 (PacBio) and 0.013 (ONT)
([report](claude/2026-10-02-r226-build-evaluation/README.md)). `build_gtdb_database.py` gives every
model a curve (`--depth-knob-read-types`; `""` turns them off).

Models of 0.7.2 have knobs by whole decade instead (`protal_depth_knobs`, bins 2-6: the digits of the
fragments less one), which protal still reads: a sample in a bin without a knob keeps `--knob`'s
default (`Sample S: N fragments, knob K (the model's for depth bin B)`). On the v0.7.1 benchmark those
cost PacBio 0.007 to 0.016 and helped Nanopore not at all
([report](claude/2026-10-01-features-depth-knobs/README.md)): a bin's knob rested on a few training
samples, and samples near a bin edge switched between two knobs, which the curve avoids.

### Calls at a target share of false calls

A knob curve over depth stands in for what changes with depth: a deep sample has many more absent
candidates. DADA2 corrects its test for the number of sequences it tests instead. With `--fdr-calls`
the trainer calibrates the forest's scores, an isotonic fit of presence on the scores of species held
out, evaluated at 64 quantiles of the scores and at 0 and 1, and notes the share of present rows
among them (the prior). In a sample, protal makes each taxon's score a probability by that curve
(linear between its points), and adjusts the probabilities to the sample's own share of present
candidates (the fixed point of the EM of Saerens, Latinne and Decaestecker, 2002: the share is the mean
of the probabilities, each probability's odds times the shares' odds ratio, with 20 more candidates at
the training share, so that a sample of a few likely candidates does not reach a share of 1 and every
probability 1 with it). A deep sample's many unlikely candidates lower its share and every probability
with it. protal then reports the taxa by
score, the highest first, while the mean of their 1 − probability (the expected share of false calls
among them) stays at or below the target, and every taxon at the last one's score. The taxa the
singleton rule vetoes are left out. The target is the one of 0.005 to 0.5 with the highest F1 on
species held out (report section "Calls at a target share of false calls", which also gives the
calibration's log loss against the scores' and the range of the samples' adjusted shares). On the test
set the report gives F1 and errors at the target, at the knob curve and at `--knob`, by depth too.
The model's header holds `<Extension name="protal_calibration" value="0.0:0.0001,...,1.0:0.99"/>`,
`protal_prior` and `protal_fdr`.

protal uses them only with `--fdr F`: at GTDB r226 they called 0.001-0.007 F1 below the knob curve for
every read type ([report](claude/2026-10-03-r226-v5-v6-training/README.md)), so a model's calibrated
calls are off by default. With `--fdr` the log says per sample how many taxa it called,
their expected false calls, the adjusted share and the knob that came of it (`Sample S: N fragments, C
taxa at an expected share of false calls of at most F (E expected; prior adjusted to the sample P),
knob K`). `--fdr 0` and `--knob` leave them off. `build_gtdb_database.py --call-mode fdr` trains them; its default, `--call-mode curve`, does not: on the
benchmark world they called 0.0005-0.004 F1 below the curve for both read types; on test samples deeper
than any trained 0.003 below it for paired-end and 0.002 above it for single-end reads; and the chosen target (0.005-0.015) varied between models
([report](claude/2026-10-03-denoising-implementation/README.md)). At r226 they were 0.001-0.007 below the curve for all four read
types and missed most present species of the shallowest samples ([report](claude/2026-10-03-r226-v5-v6-training/README.md)).

### The singleton rule

A species with a single fragment beside a species of its genus with 100 fragments or more in the
same sample is never reported, whatever its score, if its read looks like that congener's: the
abundance-weighted assignment (`em_own_share`) leaves it less than half of it, or its identity is
below 0.95, further from the reference than a strain of the species would be (`protal
--singleton_congener N`; 0, the default, for none: on the r226 v5 training data the model called none of
the 5,961 rows the rule vetoed, so it changes nothing there, and where it would, it risks a minor
congener; [report](claude/2026-10-03-false-positive-anatomy/README.md)). Of the r226 v3 build's paired-end calls of single fragments beside
such a congener, all 40 (species held out) and all 19 (test set) were false; amplicon denoisers never
call a single read at all ([report](claude/2026-10-02-amplicon-denoising/README.md)). But on the
benchmark world's samples with congener groups the rule without the read's condition removed 2-11 true
calls per test set, of minor congeners whose read fits their own reference (EM share ~0.98, identity
~0.975), and with it none ([report](claude/2026-10-03-denoising-implementation/README.md)).
The training dump's `prediction`, the profiles, the statistics and the strain MSAs follow it; the
trainer counts it in every call (`--singleton-congener`, the same default of 0), from the dump's
`fragments`, `genus_top_fragments`, `em_own_share` and `identity` (a dump without them has no rule).
