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

## Not done, and the website

- **The deep false positives without a unique seed.** One option is `both` as the default. It keeps the old rates and
  their seeding dependence, against the aim of no counts. Another is a fraction that tells "no unique seed in many
  reads" from "one read", for example the enrichment's lower bound by the seeds counted, as the identity bound does.
  Neither is tested.
- **The failed candidates on a sparse index.** They could be counted over two-seed anchors in both the numerator and
  the denominator (records with `ZS` ≥ 2), instead of dropping only the lone anchors' failures.
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
- `feature_shift.pl` and `extra_fp.pl` are those of
  [2026-10-14-sparse-index-calls](../2026-10-14-sparse-index-calls/README.md).
- The runs stay in `~/strainfn/run_v_*` and the models in `~/lfv/models` (WSL).
