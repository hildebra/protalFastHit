# Presence features that do not change with the read length or a sparse index, and training on several read lengths

The follow-up of [2026-10-14-sparse-index-calls](../2026-10-14-sparse-index-calls/README.md), which found that on a
grid's sparse index the reads align as on the whole index but count a quarter of their unique k-mers, and that 75 bp
reads under a 150 bp model call six times the false positives (MAPQ, relatives' short reads without a difference). The
user asked to use fractions instead of k-mer counts and to give the model the aligned read length, and to implement the
seven proposals of that report; during the work also: the read length's SD as a feature; build defaults of
paired-end reads of 75, 100, 150 and 250 bp (MiSeq for 250), single-end Illumina reads of 150 and 250 bp and Ultima reads
of 300 ± 50 bp, balanced so that the training does not grow much; Ultima reads that rarely exceed 320 bp, are poor over
their last 100 bases and whose errors are indels (substitutions at Q35-Q40); and random forests as the default model
again.

## What was built

On `2cfc9d1` (WSL `~/grid-build`, 4 cores), committed as `ce2940e`:

1. **`ZS`, the seeds of a record's anchor.** `ZU` and `ZT` count the unique seeds among the seeds of the anchor a record
   was aligned on; `ZS` (written where it is not 0) gives how many seeds that was (`ChainAlignmentAnchor::seeds`, through
   `AlignmentResult` to `SamEntry`). An anchor of a single grid seed (`AnchorLoneSeeds`) now counts its seed in `ZU`, `ZT`
   and `ZS` (it had `ZU` 0); a mate rescued on one seed (`RecoverAnchors`) has none (`ZS` absent).
2. **The unique seeds over their expectation.** Per gene the profiler sums its records' `ZS`, `ZU` and `ZT`
   (`profiler::Gene::m_seeds`, ...); per taxon `unique_seed_share` = ΣZU / ΣZS and `unique_seed_enrichment` = ΣZU / Σ(ZS
   x the gene's unique share x the sample's factor), the gene's unique share from `unique_kmers.tsv` ((long + short
   unique) / values; for ZT long at distance two + short). The factor (`SeedExpectation.h`) is 1 for a sample seeded
   from every core; for a sample seeded from a grid level alone (a sparse index, `--sparse_index emulate`, `--sparse_seeding`
   with the grid pick) it is that level's unique share of values over the whole index's, by the lookups' own rule
   (`Seedmap::UniqueSharesByLevel`: one parallel pass over the keys, direct or compact key map; a key with one value of the
   level counts unique at distance two whenever unique), computed once a run (`RunProtal::SeedExpectations`, logged
   "Seed expectation: ...") and written into the sample's SAM header (`@CO protal seed expectation: level=L unique=U
   two=T`), so that `--profile_only` finds it too.
3. **Ancestry sites per kb** (`ancestry_sites_per_kb`, `ancestry_indel_sites_per_kb`): per 1000 aligned bases of the
   records instead of per record.
4. **The best reads' identity bounds** (`top_identity_bound`, `allele_top_identity_bound`): the 98th percentile of the
   reads' identities, each the lower end of its 95% Wilson interval given its aligned bases (`IdentityLowerBound`: 150
   bases without a difference 0.975, 75 bases 0.951, 150 with one 0.963).
5. **No single-seed anchor among the failed candidates** (`ZF`): an attempt on a lone anchor (`FailedCandidate::lone`)
   is left out unless the taxon was also tried on a two-seed anchor.
6. **The sample's read length**: `sample_read_length`, `sample_read_length_sd` (mean and SD of the aligned bases of the
   sample's best records, all taxa) and `aligned_length_ratio` (the taxon's mean over the sample's).
   **Feature sets** (`model_features.py`): the group `length_free` puts the new features in place of the twelve they
   replace (the eight unique k-mer rates, `top_identity`, `allele_top_identity`, the two ancestry sites per record) where
   a set has those, and `read_length` adds the three; both are in the default set (90 features, 93 before).
   `unique_seed_share` and `unique2_seed_share` are in the dump, in no group (on a sparse index they are 0.79 of the whole
   index's). **Training on several read lengths**: `build_gtdb_database.py` defaults to paired-end setups of 75 (HiSeq
   2500), 100, 150 and 250 bp (MiSeq) at 9 samples a point (12 with three setups before: 36 a depth either way), the deep
   points `2000000:3,10000000:1,30000000:1` (deep paired-end pairs ~184 M instead of ~174 M), single-end Illumina reads
   of 150 and 250 bp (`--se-read-setups`, the first reads of those points' simulations) and Ultima reads
   (`--se-ultima-setup ultima:300:50:25:2`, one unit a depth of the depth's samples, replaying its communities:
   `collect_training_data.py --se_read_setups`, `--se_ultima_setup`): per depth 36 paired-end and 27 single-end samples
   (36 and 36 before), one more Illumina setup simulated. **Ultima reads** (`LongReadSimulator`): cut at 320 bases
   (`ultima:...:MAX`), base qualities falling over the last 100 bases to 15 below the read's at its end, errors indels
   by the qualities and substitutions apart at a rate of Q35-Q40 drawn per read. **Random forests**: the build's
   default `--model forest,gbm` and the trainer's `--model forest`.
7. **Validation**: below.

**Tests.** Unit tests 582 (new: `test_LengthFreeFeatures.cpp`, the failed candidates, the identity bound, the header
line, `ZS` written and read; `GridLevels.UniqueSharesCountTheValuesALevelsLookupsCountUnique`, direct and compact; the
Ultima model's tail, cut and substitutions in `test_LongReadSimulation.cpp`), e2e 182 (new:
`test_seed_tags_and_the_seed_expectation`), the scripts' 103 (the collector's split setups and Ultima units, the feature
sets, the forest by default) and the pipeline's 8 (its builds keep one setup for single-end reads and no Ultima units;
the scenario test no longer requires paired-end errors of its tiny world: its paired-end model made none) pass.

## Validation

**Setup** (WSL, 4 cores; `ce2940e`'s tree, its binaries and scripts copied to `~/lfv`):

- **Two builds of the small real world** (`~/realworld`: 6 families, 530 genomes; `lfv_build.sh`): paired-end reads
  only, seed 1, five depths (2,000 to 1,000,000 pairs), random forests, the new default features. `mix`: 15 samples a
  point at each of 75, 100, 150 and 250 bp (300 samples). `b150`: 60 samples a point at 150 bp (300 samples, as the
  strain study's model was trained).
- **Models** (`lfv_strain.sh`): each build's ungated model (`model_pe_c1.xml`, as the build trained it: `new`, 90
  features), and the same training table trained again (`machine_learning_cmdline.py`, the build's own command) as:
  - `prev`: the previous default set (93 features);
  - `nomapq` (`mix` only): `new` without `mean_mapq`, `uniqueness` and `low_mapq_share` (87);
  - `shares` (`mix` only): `new` with the plain `unique_seed_share` and `unique2_seed_share` in place of the two
    enrichments (90);
  - `both`: `new` with the eight unique k-mer rates it replaced added back (98);
  - `lf`: the previous set with `length_free` alone, without `read_length` (87).

  The variants by feature list were trained on the inputs of a stub model (`lfv_inputs.pl`, `--reference-pmml`).
  `nomapq` and `shares` were asked for; `both` and `lf` were added to explain the first results.
- **Test:** the strain study's world (`~/strainfn`; [2026-10-11-strain-false-negatives](../2026-10-11-strain-false-negatives/README.md)):
  532 held-out real strains of 30 species, 20 shallow and 20 deep samples, database `db_grid`. It was aligned four
  times by the new protal: ungated at 150 bp (`l150`), at 100 and 75 bp (`l100`, `l75`: the same reads cut, as in
  [2026-10-14-sparse-index-calls](../2026-10-14-sparse-index-calls/README.md)), and with the grid at C 4 on its
  sparse index at 150 bp (`g4`; seed expectation of level 2: unique 1.026, distance two 1.166). Each alignment was
  profiled with each model (`--profile_only`) and scored by the study's `evaluate.py` (`~/strainfn/eval_v.log`) and
  at equal false positives (`matched_fp.pl`).

**The builds' own evaluation** (c1 models, species held out): F1 0.9733 (`mix`) and 0.9743 (`b150`). Their
independent test sets differ in read lengths, so they are not compared. In `mix`'s model the enrichments rank 2nd and
3rd (importance 0.094 and 0.090, after `allele_identity`), `top_identity_bound` 6th, `allele_top_identity_bound` 9th
and the three read-length features 35th to 47th.

**The new features hold across read lengths and seeding** (`feature_shift.pl` of the earlier report on `mix_new`'s
dumps; present species' median ratio to `l150`):

| feature | `l100` | `l75` | `g4` |
|---|---|---|---|
| `unique_seed_enrichment`, `unique2_seed_enrichment` | 1.01, 1.01 | 1.02, 1.02 | 0.97, 0.99 |
| `unique_seed_share`, `unique2_seed_share` | 1.02, 1.02 | 1.02, 1.02 | 1.00, 1.15 |
| (replaced) `lu_per_kb`, `lsu_gene_rate3` | 0.91, 0.92 | 0.79, 0.76 | 0.26, 0.75 |
| `ancestry_sites_per_kb` (replaced: per record) | 0.99 (0.66) | 0.97 (0.49) | 1.00 (1.00) |
| `top_identity_bound` | 0.99 | 0.98 | 1.00 |

- **Identity bound.** `top_identity_bound` now shifts the other way: at 75 bp present species lose 0.65 of the
  present-absent gap, absent ones 0.26. A perfect 75 bp read is bounded at 0.951 and a perfect 150 bp read at 0.975.
  The model needs the read length to read it.
- **Plain shares.** `unique_seed_share` holds as well as the enrichment, except at distance two on the grid, where
  the seed expectation corrects ×1.15 to 0.99. On this database the grid level's unique shares are close to the
  whole index's.
- **Failed candidates.** Without the single-seed anchors' failures, `failed_candidate_rate` of absent species falls
  from 0.34 to 0.15 on `g4`. With them it had doubled, to 0.65. Dropping a lone anchor's failure also drops real
  failures, while its successes still count.

**At the knob (0.5).** Each cell gives the share of test strains called in shallow samples (hard strains in brackets),
then the false positives in the 20 shallow / 20 deep samples:

| model | `l150` | `l100` | `l75` | `g4` |
|---|---|---|---|---|
| `mix_new` | 0.918 (0.881) 3 / 12 | 0.903 (0.851) 9 / 17 | 0.922 (0.891) 12 / 29 | 0.915 (0.871) 6 / 15 |
| `mix_prev` | 0.902 (0.861) 2 / 9 | 0.908 (0.861) 7 / 11 | 0.907 (0.842) 11 / 10 | 0.847 (0.644) 1 / 5 |
| `mix_nomapq` | 0.905 (0.842) 7 / 20 | 0.888 (0.792) 11 / 26 | 0.913 (0.861) 14 / 45 | 0.888 (0.762) 8 / 23 |
| `mix_shares` | 0.913 (0.891) 6 / 14 | 0.913 (0.891) 9 / 24 | 0.923 (0.931) 15 / 32 | 0.913 (0.891) 4 / 15 |
| `mix_both` | 0.902 (0.851) 2 / 11 | 0.890 (0.822) 6 / 9 | 0.907 (0.861) 6 / 14 | 0.870 (0.693) 2 / 10 |
| `mix_lf` | 0.905 (0.842) 4 / 11 | 0.900 (0.832) 9 / 16 | 0.915 (0.891) 15 / 25 | 0.900 (0.832) 5 / 12 |
| `b150_new` | 0.910 (0.871) 6 / 21 | 0.895 (0.832) 9 / 23 | 0.898 (0.832) 10 / 16 | 0.913 (0.871) 6 / 24 |
| `b150_prev` | 0.898 (0.901) 2 / 7 | 0.897 (0.881) 5 / 9 | 0.905 (0.891) 7 / 11 | 0.862 (0.762) 1 / 4 |
| `b150_both` | 0.922 (0.941) 1 / 10 | 0.905 (0.911) 4 / 6 | 0.888 (0.851) 4 / 6 | 0.908 (0.901) 1 / 8 |
| `b150_lf` | 0.912 (0.911) 4 / 15 | 0.885 (0.842) 3 / 9 | 0.873 (0.802) 2 / 5 | 0.918 (0.901) 5 / 15 |

**At equal false positives** (`matched_fp.pl`): the share of test strains called when each model's knob lets through at
most 5 false positives in the shallow samples and at most 10 in the deep ones. This uses the probabilities, without
the singleton rule. Each cell is shallow / deep:

| model | `l150` | `l100` | `l75` | `g4` |
|---|---|---|---|---|
| `mix_new` | 0.930 / 0.982 | 0.892 / 0.973 | 0.902 / 0.962 | 0.913 / 0.975 |
| `mix_prev` | 0.927 / 0.983 | 0.872 / 0.975 | 0.872 / 0.973 | 0.907 / 0.983 |
| `mix_nomapq` | 0.870 / 0.982 | 0.855 / 0.972 | 0.875 / 0.957 | 0.857 / 0.963 |
| `mix_shares` | 0.910 / 0.968 | 0.887 / 0.970 | 0.907 / 0.962 | 0.918 / 0.962 |
| `mix_both` | 0.927 / 0.967 | 0.887 / 0.973 | 0.898 / 0.965 | 0.910 / 0.978 |
| `mix_lf` | 0.913 / 0.982 | 0.880 / 0.975 | 0.877 / 0.968 | 0.900 / 0.980 |
| `b150_new` | 0.897 / 0.972 | 0.867 / 0.967 | 0.882 / 0.987 | 0.898 / 0.972 |
| `b150_prev` | 0.937 / 0.975 | 0.898 / 0.963 | 0.887 / 0.965 | 0.925 / 0.965 |
| `b150_both` | 0.945 / 0.980 | 0.917 / 0.982 | 0.910 / 0.982 | 0.940 / 0.977 |
| `b150_lf` | 0.918 / 0.985 | 0.915 / 0.987 | 0.888 / 0.980 | 0.923 / 0.988 |

At 2 and 10 shallow false positives the order is about the same (`matched_fp.txt`). One test strain is 0.0017 of a
column. Differences below about 0.015, or of a few false positives, are noise: `mix_new` and `mix_shares` differ in two
nearly equivalent features, and they differ by 3 / 2 false positives at the knob.

**The false positives the new set adds** (`fp_species.pl`, `extra_fp.pl`):

- **Present in every model's calls:** *Weizmannia faecalis* is a false positive in 6 of the 20 deep samples under
  every model.
- **Added by the new set:** taxa with 1 to 5 fragments, a perfect best read and no unique seed (*Priestia
  megaterium* and *megaterium_C*, *Dorea_A longicatena_C*, single reads of *Actinobacillus* and *Lautropia*). The
  near-identical congeners are those of the earlier report. `b150_prev` gives them 0.25-0.45, `b150_new` 0.50-0.72.
- **The cause:** where the reads are 150 bp, the identity bound and the per-kb ancestry sites are only relabellings
  of what they replace (one read length). The difference is the eight unique k-mer rates that became two enrichments.
  All are 0 for these taxa. A forest drawing 9 of 90 features per split meets "no unique seed" in far fewer splits.
  With the rates back (`both`), the deep false positives fall to the previous set's (`b150` 21 → 10, `mix` at 75 bp
  29 → 14).

### What it means

1. **Sparse index: the probabilities keep their meaning.** At the knob, `mix_new` calls 0.915 on `g4` against 0.918
   ungated, and `mix_prev` 0.847 against 0.902 (hard strains 0.871 against 0.644). At equal false positives, however,
   `prev` on `g4` came within 0.02 of its `l150` and matched `new` there (`mix` 0.907 against 0.913, `b150` 0.925
   against 0.898). The sparse index lowered its probabilities without reordering them. The
   gain is calibration: one model, trained on fully seeded runs, serves a sparse index or `--sparse_seeding` at its own
   knob. Builds already avoided the problem for their gates by training a model per gate on gated runs.
2. **Read length: not solved.** At 75 bp eight of the ten models call 2 to 6 times the shallow false positives they call
   at 150 bp (`b150_new` 1.7 times, `b150_lf` fewer). Every `mix` model is among the eight, although its training had
   75 bp samples. At equal false positives `mix` is not better than `b150` at 75 or 100 bp: it is better with the new
   set and worse with each of the others. The features that should hold do hold. MAPQ, congener ties and relatives' short perfect reads do not, and
   15 samples a point of each length did not teach the forest to weigh them. The strain world's 75 bp false positives
   are near-identical congeners, which the training worlds may rarely hold.
3. **New against previous set at 150 bp:** about equal at the knob. At equal false positives the previous set leads
   at the strict end (2 shallow false positives: `mix` 0.910 against 0.887, `b150` 0.913 against 0.875). The cause is
   the deep false positives above.
4. **(a) The MAPQ features stay.** Without them the models call fewer strains in every run (−0.01 to −0.03 at the
   knob) with more false positives (deep 20-45), and they are last or nearly last at equal false positives.
5. **(b) Plain shares ≈ enrichments here,** because `db_grid`'s seed expectation is close to 1. The enrichment
   corrects the distance-two share on the grid (×1.15 → ×0.99), and a database whose grid levels keep fewer unique
   cores would need that. The enrichments stay.
6. **`both`** (the enrichments plus the eight rates) has the fewest false positives, and `b150_both` is best at equal
   shallow false positives in every run. `mix_both`, however, loses the sparse index at the knob again (0.870 on `g4`): the
   rates bring the seeding dependence back.

## Follow-ups (same day)

### Failed candidates over anchors of two seeds on both sides

The user asked for `failed_candidate_rate` (and `failed_gene_share`) to count only the reads anchored on the taxon by
two seeds or more, in the denominator as in the failures. This is built on `f98254a`, uncommitted (`Profiler.h`:
`LoneFragment`, `RecordEvidence::fragments_lone`, `GeneRecords::lone_records`; `test_LengthFreeFeatures.cpp`). The tree
was built from `f98254a` plus this patch only, because the shared checkout held another session's uncommitted index
work. Unit 583 and e2e 182 pass. The strain world was re-profiled from the same alignments (`~/strainfn/run_w_*`,
`eval_w.log`):

| absent species' mean | whole index (`l150`) | sparse index (`g4`), before | `g4`, after | `l75`, whole index |
|---|---|---|---|---|
| `failed_candidate_rate` | 0.340 | 0.155 | 0.163 | 0.147 |
| `failed_gene_share` | 0.333 | 0.153 | 0.163 | 0.141 |

- **Whole-index runs:** the same to the last digit. They have no single-seed anchors.
- **Calls:** the same on every run, with `mix_new` and `b150_new`.
- **Why so little:** few reads align through single-seed anchors alone, so the denominator hardly moves. The fall is
  in the failures. A relative's read meets the taxon on one grid seed where the whole index gives it about four, so
  its failure there is a single-seed failure and goes uncounted. 75 bp reads on the whole index lower the rate the
  same way, through fewer two-seed anchors.
- **A fix of the failures** would count a single-seed failure whose anchor extends as far as several seeds would.
  Not built.

### Can single-nucleotide sites stand in for the unique k-mers?

**The idea.** The user proposed marking the single nucleotides that are unique, and counting them as a read's
alignment covers them with the species' base, in place of the unique k-mers.

**Testing it with what exists.** The ancestry sites are such nucleotides already: where a species' copy carries a base
that 9 in 10 of its compared congeners do not, its known genomes sharing it. They are counted per read from the
alignment (`ancestry_agreement`, `ancestry_congener_share`, `ancestry_sites_per_kb`).
- **Models:** `mix` and `b150` trained again without the two enrichments (`nouniq`, 88 features). No unique k-mer
  feature is left; the ancestry features stay (`sites_test.sh`, `~/strainfn/eval_s.log`).
- **Result:** strains called at equal false positives (`matched_fp.pl`), at most 2 / 5 shallow and 10 deep false
  positives. Each cell is with / without the unique k-mers:

| run | `mix` ≤2 s | `mix` ≤5 s | `mix` ≤10 d | `b150` ≤2 s | `b150` ≤5 s | `b150` ≤10 d |
|---|---|---|---|---|---|---|
| `l150` | 0.887 / 0.645 | 0.930 / 0.830 | 0.982 / 0.893 | 0.875 / 0.748 | 0.897 / 0.857 | 0.972 / 0.908 |
| `l75` | 0.792 / 0.740 | 0.902 / 0.768 | 0.962 / 0.857 | 0.828 / 0.718 | 0.882 / 0.758 | 0.987 / 0.825 |
| `g4` | 0.878 / 0.657 | 0.913 / 0.780 | 0.975 / 0.880 | 0.890 / 0.727 | 0.898 / 0.883 | 0.972 / 0.897 |

At the knob (`evaluate.py`'s calls) the `nouniq` models call about as many strains (0.903-0.932, `new` 0.898-0.922)
with 1.5 to 6.7 times the false positives (`mix` on `l150`: 18 shallow / 80 deep against 3 / 12). The ancestry features rank 9th, 12th and 13th in `mix_nouniq`.

**What it means:**
- **Robustness:** counting nucleotides is more robust than k-mers in how it counts:
  - it does not depend on the seeding;
  - a site needs one covered base where a k-mer needs 31 inside the read (`ancestry_sites_per_kb` 0.97 at 75 bp,
    `lu_per_kb` 0.79);
  - a sequencing error costs one site where it costs every k-mer around it (at 1% errors 27% of 31-mers carry one).
- **Information:** the sites protal has today carry far less than the unique k-mers. They are defined against the
  congeners' consensus, while a k-mer is unique against every reference, and probably often through a combination of
  bases that are each shared with some other reference. How many unique k-mers hold a base no other reference carries
  at its column is not measured.
- **What a replacement would need:** a nucleotide design that replaces the k-mers would have to capture those
  combinations. For example: each site records which competing references its base excludes, and a read counts as
  unique where the sites it covers exclude every competitor (a unique k-mer is the same over a 31-base window).

**Memory at r226** (17.2 G gene bases from the 4.3 GB 2-bit gene store, ~2.9 G index entries; densities from the small
real world: 87 bytes of `ancestry_sites.bin` per kb of gene, 0.18 index entries per base):

| what is stored | per kb of gene | r226 |
|---|---|---|
| the index entries' two unique flags, kept by the sparse load (the k-mer reconstruction) | ~43 B | ~0.7 GB |
| unique-nucleotide positions, varint gaps, at 5% of positions | ~50 B | ~0.85 GB |
| the same, Golomb-coded gaps | ~36 B | ~0.6 GB |
| the same at 1% of positions | ~11 B | ~0.2 GB |
| a site's excluded competitors, 16 neighbours a site (2 bytes) at 10-20% of positions | ~300-600 B | ~5-10 GB |

Plain unique nucleotides are cheaper than the k-mer flags only if they are rarer than ~4% of positions; storing which
competitors each site excludes, which the combinations need, costs several times more. Not built.

### Unique k-mers counted from the alignment

The user asked for the k-mers to stay the measure of uniqueness, counted from each record's alignment so that the sparse
index does not matter. Built on `ab6bf25`, uncommitted, with the two-seed failed candidates above:

- **What a record counts:** the gene's cores whose whole 31-mer lies inside a run of matches of its alignment (`M`
  columns with no other column between them). A seed needs its 31-mer exact (16 of 16 flank bases), so these are the
  cores the whole index's seeds find on that stretch. `ZS` is their number, `ZU` and `ZT` the unique ones among them
  (`UniqueCores.h` `CountCovered`, called in `SimpleAlignmentHandler::AlignAnchor` once the alignment is valid).
- **Where the bits come from:** `unique_cores.bin`, made by `--build` (`WriteUniqueCores`). For every gene copy, in the
  order of its cores as the k-mer handler finds them in its 2-bit sequence, two bits: unique (flag 60) and unique at
  distance two (flag 61 for a core of several values, flag 60 for one of one, as seeds count them). It is read from the
  index once its flags are final (after `unique_kmers.tsv`); the genes are now released after it rather than before
  `unique_kmers.tsv`, at the same peak, because the uniqueness check before held them beside the index as long.
- **Why at the build:** the sparse load sees every entry, but in key order. Turning entries into per-gene arrays while
  loading would need a temporary table of several GB at r226.
- **At run time:** the table is read before the alignment, beside the genome preload. Each alignment thread finds a
  gene's cores with its own copy of the run's k-mer handler and keeps them for the last 4,096 genes. A gene whose cores
  do not fit its row (another database's table) keeps the anchor's counts.
- **What changes in a SAM:**
  - `ZS`, `ZU` and `ZT` are the alignment's counts wherever the database has the table;
  - the anchor moves to `ZL` (`ZL:i:1` a single seed, `ZL:i:0` a rescued mate, none for two seeds or more), which the
    failed-candidate rule now reads;
  - no seed expectation is written.

  Without the table, or with `--no_unique_cores`, the counts are the anchors' as before.
- **Cost at r226:** 2 bits for each of its 2.9 G entries (~0.73 GB) plus 4 bytes per gene row (~66 MB). On the strain
  world: 10.28 M cores of 53,833 copies, 68.2% unique, 2,631 without an index entry; 2.68 MB, built in 0.5 s.
- **Tests:** `test_UniqueCores.cpp` (the stretches, the table, the counter), `ZL` in the tag tests, and the e2e test's
  per-record check that a grid-seeded run counts what the whole index's run counts. Unit 590, e2e 182, scripts 103 pass
  (HEAD plus this patch, built in WSL).

**On the strain world** (`uc_strain.sh`: its database built again as `~/strainfn/db_cores`; the reads aligned at 150 and
75 bp and by the grid at C 4 on its sparse index; `uc_tags.pl` compares the records of the same alignment):

| runs | records of the same alignment | the same `ZS`, `ZU`, `ZT` | mean `ZU`, whole / sparse |
|---|---|---|---|
| counted from the alignment | 175,602 | 175,602 (100.00%) | 12.02 / 12.03 |
| counted from the anchors (`run_v_*`, before) | 175,602 | 163 (0.09%) | 12.02 / 3.13 |

At 75 bp a record covers fewer whole k-mers (mean `ZS` 7.16 against 18.22), and the unique share holds (0.651 against
0.660).

**Models trained on the new counts:** (running: `uc_models.sh`.)

### Ancestry sites in less memory

The user asked whether the ancestry sites really take ~1.5 GB, my estimate from the small world's bytes per kb of gene.
The r226 v27 build log gives **1.13 GB** for 11.94 M copies with 630 M sites (52.7 a copy). A run holds the whole file in
memory, read beside the alignment, so it counts toward the run's peak.

What the bytes hold (`anc_breakdown.pl`; reduced r226 database of the mock benchmark, 71,695 copies, 87 sites each):

| part | share |
|---|---|
| positions (a zigzag varint step each, ~1 byte a site) | 63% |
| bases (two bits each) | 16% |
| index (64-bit key and offset a copy) | 11% |
| identity (a float's bits as a varint, 5 bytes) | 3.5% |
| the rest of each copy's header and the indels | ~7% |

**Version 2** (`AncestrySites.h` `EncodeCompact`, `Table::Write`; a run reads both):
- positions as Rice-coded gaps, the parameter chosen per copy (0.57 of the varints);
- 4-byte keys (`taxid << 8 | gene`) and offsets where they fit;
- the identity in 16 bits.

The sites are the same. Measured:

| table | version 1 | version 2 |
|---|---|---|
| strain world (built twice, 45,926 copies) | 5.27 MB | 3.38 MB (0.64) |
| reduced r226 (computed from its version 1) | 10.2 MB | 6.6 MB (0.65) |
| small real world (computed) | 1.65 MB | 1.14 MB (0.69) |
| r226 v27 (estimated at 0.65) | 1.13 GB | ~0.73 GB |

**Thinning instead of compressing** loses evidence (`anc_density.pl`, reduced r226): a cap of one site per 8 bases of a
copy's span would thin 23.5% of the copies and keep 82.6% of the sites; one per 16 bases thins 54.7% and keeps 55.7%.
The median copy has a site every 13 bases. The few-read taxa the deep false positives come from need every site, so
nothing is thinned.

## Not done, and the website

- **The deep false positives without a unique seed.** One option is `both` as the default. It keeps the old rates and
  their seeding dependence, against the aim of no counts. Another is a fraction that tells "no unique seed in many
  reads" from "one read", for example the enrichment's lower bound by the seeds counted, as the identity bound does.
  Neither is tested.
- **The failed candidates on a sparse index.** The denominator is now over two-seed anchors too (above); the failures
  still fall. A single-seed failure whose anchor extends as far as several seeds would could count.
- **The ancestry table:** a dense index by gene row (4 bytes a row instead of 8 a copy) would save ~45 MB more at
  r226.
- **The read length.** Short reads need more than features: more training samples of short reads, or a model per
  read-length class. A test world of short reads with near-identical congeners would show which.
- **An r226 build with the new defaults** is the real test. Its `--test-grid` rows now see the grid's runs through
  features that do not change with the seeding.
- The website needs the `ZS` tag, the seed expectation header line, the model default (a forest again) and the new
  features; it is not in this repository.

## Files

- `lfv_build.sh`: the two builds.
- `lfv_strain.sh`: the models and the strain world's alignments, profiles and evaluation.
- `lfv_inputs.pl`: a stub model of a feature list, for `--reference-pmml`.
- `matched_fp.pl`: calls at equal false positives; `matched_fp.txt` holds its output.
- `fp_species.pl`: the false positives by species.
- `sites_test.sh`: the follow-up's models without unique k-mer features.
- `uc_strain.sh`, `uc_tags.pl`: the strain world's database built again with `unique_cores.bin`, its runs, and the
  records' seed tags compared; `uc_models.sh`: the small real world built again and its models on those runs.
- `anc_breakdown.pl`, `anc_density.pl`: what an `ancestry_sites.bin` holds, version 2's size, and what a density cap
  would keep.
- `feature_shift.pl` and `extra_fp.pl` are those of
  [2026-10-14-sparse-index-calls](../2026-10-14-sparse-index-calls/README.md).
- The runs stay in `~/strainfn/run_v_*` and the models in `~/lfv/models` (WSL).
