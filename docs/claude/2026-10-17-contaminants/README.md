# Contaminants in the training simulations, and a slightly higher omega for the in-silico strains

The real mock communities showed what the simulated training data lacks
([2026-10-16-mock-regression](../2026-10-16-mock-regression/README.md)):

- **Low-level DNA of other organisms:** real libraries carry it (other samples on the flow cell, reagents, the lab),
  often of genera the community does not hold. 0.8.5 found 9-12 species the mock does not list per HiSeq library, 52-73
  per DNBSEQ one, at 1e-6 to 1.6e-5 of the profile. In the simulations every species of a sample came from its
  community's draw.
- **Strains further from their references than the in-silico ones:** real isolates at 98% identity were rejected.

This change adds contaminants to the simulated samples, on by default in GTDB builds. **Contaminants are true
positives:** they are real species of the sample, so its truth lists them. That is the user's call, made before this
was committed; a first draft labelled them absent. As the user asked, the change also raises the in-silico strains'
omega by a factor of 1.3. Nothing was rebuilt or retrained, so the effect on F1 is not measured yet (section 5).

## 1. Contaminants: what a sample gets

`simulate_metagenomes --contaminants MIN-MAX:LOW:HIGH` (also `N:LOW:HIGH`):

- **What is drawn:** each sample gets MIN to MAX (uniform) genomes of species its community does not hold, never an
  unclassified genome. Each is one genome of its species, with a share of the sample's read pairs drawn
  log-uniformly from LOW to HIGH and rounded at random, so a share below one pair still gives a pair now and then.
  The pairs come on top of the community's.
- **Its own random stream:** the draws, and the seeds of the contaminants' reads, come from a stream of their own
  (`MixSeed(seed, "cont" + sample)`). The communities' genomes, read pairs and seeds, and every read of theirs, are
  byte for byte those of a run without contaminants. The contaminants' reads follow the community's (before the
  host's), so the same samples can be compared with and without them. Only the relative abundances shift, by the
  contaminants' small share; the community's keep their proportions.
- **Where they go:**

  | Output | |
  |---|---|
  | reads | after the community's |
  | `protal_goldstd/*.profile_truth` | listed: present species of the sample |
  | `manifest.tsv`, `manifests/` | a row each, new last column `contaminant` 1 (0 for the community's genomes; a manifest without the column has none) |
  | relative abundance, `abundance_matrix.tsv` | as any species: all of a sample's add up to 1 |
  | `run_params.tsv` | `contaminants` (the spec) |
  | stdout | `Contaminants: N genomes in S samples, P read pairs (x% ...), present in the truth` |
  | `--from_manifest` | a replay keeps them, from their own seeds, and draws no number of the main stream for them |
  | `--test` | the same draws (a design run's manifest, which the long reads replay, has them) |

Files: `src/RandomForest/MetagenomeTypes.h` (`GenomeAssignment::contaminant`, the design options),
`CommunityProfileDesigner.{h,cpp}` (`draw_contaminants`), `MetagenomeSimulator.{h,cpp}` (`parse_contaminants`,
`add_contaminants`, the manifest column), `src/simulate_metagenomes_main.cpp` (the option, `run_params.tsv`, the summary
line).

## 2. Contaminants in the collector and the build

- **`collect_training_data.py --contaminants`:** passed to every paired-end simulation, the design's and the
  scenarios'. It is checked at the start (`contaminant_spec`) and is part of each simulation's key (the command), so
  points simulated without it are simulated again.
- **The long reads** (and the Ultima units) draw from the community's manifest rows, contaminants included, each by
  relative abundance times length: the same share of bases as of the paired-end reads.
- **The training table:**
  - contaminants are species of the sample like any other: present rows; an unseen species when they have reads on
    a marker gene and no row; relatives for `meta_relative_rank`;
  - a new last meta column, `meta_contaminant`, marks them with 1, so the trainer's reports and later analyses can
    tell how the models find them;
  - the trainer ignores every `meta_*` column (`model_features.feature_columns`).
- **Unchanged:** `error_reads.py`, `composition_accuracy.py` and `trace_relatives.py` read contaminants as the
  sample's species they are.
- **`build_gtdb_database.py --contaminants`**, default `0-30:1e-6:2e-5`:
  - passed to the training and test collections;
  - recorded in `build_metadata.tsv` (`classifier_training_design`: `contaminants ...`);
  - `--contaminants none` builds as before.

## 3. The default level

The unlisted species 0.8.5 found in the real mocks:

| Measure | Value |
|---|---|
| profile share | 1e-6 to 1.6e-5 |
| marker-gene fragments at 20M pairs | 2 to 25 |
| per Illumina HiSeq library | 9-12 |
| per DNBSEQ library | 52-73 |
| per Ion Torrent library | 116-150 |

`0-30:1e-6:2e-5` puts a median of 15 contaminants in a sample, each at 1e-6 to 2e-5 of its pairs. Assuming marker genes
take about 4-5% of a genome's reads:

| Sample depth | Marker-gene fragments per contaminant |
|---|---|
| 2M pairs and less | none mostly, one now and then |
| 10M pairs | 0.5-10 |
| 30-60M pairs (the deepest design points and the scenarios) | up to 50 |

That is between the HiSeq and the DNBSEQ libraries. The level is a judgement, which the next build's test set and the
mocks should calibrate (section 5).

## 4. Omega of the in-silico strains

`insilico_strains.py --omega auto` calibrates the dN/dS-like omega so that the in-silico strains put as many of their
marker substitutions on third codon positions as the real strains do (0.074-0.075 at r226 v17 and v18; 0.15 without
real strains).

- **The new options:** `--omega-scale FACTOR` multiplies the omega auto finds (or the default), and the run's log line
  says so (`omega 0.096 (0.074 x 1.3): ...`). `build_gtdb_database.py --insilico-omega-scale`, default 1.3, passes it.
  1 keeps the calibrated value.
- **What changes:** at the same divergence (each gene's target is unchanged), more substitutions change an amino acid.
  The share on third positions falls a little below the real strains'.
- **What it is for:** the Treichel misses: isolates at 98% identity to their references, MAPQ 21, 151 fragments, at
  P 0.04-0.69.
- **Caveat:** a spectrum further from the real strains' can also be something the models learn to separate. The next
  build's `insilico_strains` line and its strain F1 will tell.

## 5. Checks run, and what is not done

The binaries and every suite ran on `git archive HEAD` with the changed files laid over it (WSL `~/contam-build`,
`~/contam-tree`, `taskset -c 0-3`), so that no other session's uncommitted changes were in them.

**What these checks cover:** they tested the first draft, which labelled contaminants absent. The version committed
(present species) differs from it in a few places:

- the truth and abundance matrix list contaminants again;
- relative abundances are normalised over all of a sample's genomes again;
- the collector counts contaminants among a sample's species;
- `error_reads.py` and `composition_accuracy.py` are back to HEAD's.

The tests changed with them, but the commit was made at the user's request before the rerun of every suite on it
finished.

On the first draft:

- **C++ unit tests:** 634 pass (3 skipped: a debug-build check, the AVX-512 scan and a gene-sequence check, none of
  the simulator's), among them the new ones:
  - `CommunityDesign.ContaminantsParse`;
  - `CommunityDesign.ContaminantsAreOtherSpecies`: never a species of the community nor an unclassified genome, as
    many as asked for, shares in range, none where the sample holds every species;
  - `IlluminaSimulation.ContaminantsAddReadsBesideTheSameCommunity`: the community's genomes, pairs and seeds as
    without them, the contaminants' reads after the community's and named by their contigs, the manifest column, a
    replay of the same reads. Since the change to present species it also checks that all abundances add up to 1,
    with the community's in their old proportions.
- **Smoke test of the binary** (12 genomes, 3 samples of 20,000 pairs, `--contaminants 3:0.001:0.005`):
  - 3 contaminants per sample, 25-61 pairs each;
  - `run_params.tsv` records the spec, and a bad spec stops the run before any genome is read;
  - with the same seed and no contaminants, the manifests' member rows are identical, and the first reads files'
    first 80,000 lines are identical (the contaminants' 132 reads follow).
- **Python tests:**
  - `test_collector` (30, with `SIMULATE` set and without a skip);
  - `test_insilico_strains` (17, the omega scale in `test_spectrum_and_omega`);
  - `test_error_reads`, `test_composition_accuracy`, `test_trace_relatives`;
  - `test_gtdb_build` (46);
  - `test_profile_scripts`, `test_strain_scripts`, `test_foreign_rates`, `test_ancestry_sites` (25).
- **The mini GTDB pipeline** (`test_gtdb_pipeline.py`, 8 tests, 336 s) and the end-to-end tests (192, on a mini
  database built by the new protal) pass. At the mini world's shallow depths a contaminant seldom gets a pair, so the
  pipeline checks the plumbing, not the effect.

Not done:

- **No r226 build or training:** the effect on the models and on the mocks is unmeasured.
- **What to expect:** the models now see present species of a few fragments, of genera the community often lacks, at
  realistic levels, and the test set counts finding them. That makes calling such species likelier, as 0.8.5 already
  does on the mocks; there the truth tables list only the designed members, so these calls count as false positives.
  Scoring the mocks with that in mind (low-level unlisted species apart from false positives) is the natural next
  step.
- **What to check after the next build:**
  - recall on the `meta_contaminant` rows of the test set;
  - the false positives of the other absent rows (held-out relatives), against v31;
  - the thin present species' recall.
- **Not separated:** the two changes (contaminants, omega) land in one build. `--contaminants none` or
  `--insilico-omega-scale 1` takes either back.
