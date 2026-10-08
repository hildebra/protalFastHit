# The sample's composition: genome sizes, the share of the reads explained and the unknown share

2026-10-08. Code: `822c4cd` plus this change (uncommitted when the numbers were taken), built in WSL (gcc 13,
Release, 4 cores). Asked for: are genome sizes and the marker genes' share of a genome stored at database build
(no: now they are), and a run that tells how many of the scanned reads the detected genomes explain, the sample's
average genome size and its spread, and the species the rest would be, with the rest as a `?` line that fills
the profile to 100%; then a check in `build_gtdb_database.py` of how accurate `?` and the genome sizes are.

## What was there

`species_priors.tsv` held the representative's marker count, CheckM quality and GTDB's species cluster (for
the `+priors` features). GTDB's metadata, which the converter already read for CheckM, has every genome's
`genome_size` and its species representative, but no size was kept; `build_gtdb_database.py` measured the
FASTAs of the genomes it simulates from (`genomes.tsv`), not the database's species. A profile's abundances
were the called species' shares of their summed marker depth. The SAM did not say how many reads were read.

The converter called `read_quality` and `read_sp_clusters` with `args.release` instead of the release it had
detected: a conversion without `--release` (the tests', the mini databases') had no CheckM or cluster values.
`build_gtdb_database.py` always passes `--release`, so the r226 databases were not affected. Fixed.

## What changed

- **Database** (`scripts/mini_db/gtdb_to_protal_db.py`, `SpeciesPriors.h`): `species_priors.tsv` gains
  `genome_size` (the mean over the species' genomes of assembly size × (100 − contamination) / completeness,
  over its MIMAG high-quality genomes, else medium-quality, else the raw sizes), `sized_genomes`,
  `rep_genome_size`, `marker_bases` (the representative's marker genes, a copy each, of every marker) and
  `marker_share`. Runs read the size columns by name; older tables load as before (sizes unknown).
  `--priors_only` writes the table alone; `protal --add_tables species_priors.tsv` stores it in a built
  database and refuses a table whose taxids are not the database's species with the same representative.
- **Alignment** (`Classify.h`, `RunProtal.h`, `IO/ScannedReads.h`): the read loops count bases; the SAM
  header gets `@CO protal scanned reads: fragments=F reads=R bases=B` (not with `--full_sam_header`, whose
  header is written before the reads).
- **Profile** (`Profiling/Composition.h`, `Profiler.h`): per sample, the called species' depth × genome size
  (a species without a size at the average one) against the bases read times the aligned reads' fragment
  bases per read base (a pair's overlap counted once, as the depth counts it). The rest, at the called
  species' depth-weighted average genome size, is the unknown genome equivalents; their share of all is the
  profile's last line `?\t?\t<share>`, and the species' abundances are scaled by 1 − share
  (`--no_unknown_share`: as before). `<profile>.composition` (one line per sample) has the counts, shares,
  the average genome size, its SD and depth-weighted 10/25/50/75/90% quantiles, and the unknown genome
  equivalents over the median and the lowest called depth (`MissingSpeciesAt{Median,Lowest}Depth`).
  `.profile.log` gains `GenomeSize` and `GenomeFragments`. One console line per sample.
- **Scripts**: `protal_profile_utils merge` puts `?` last (NA where a profile has none), `composition`
  merges the `.composition` files; `prevalence_calls.py --profiles` keeps a sample's `?`;
  `simulate_gtdb_release.py --genome_length LOW-HIGH` draws each species' background length (log-uniform,
  its own random generator: single-length releases are byte-identical).
- **Build** (`scripts/composition_accuracy.py`, a step of `build_gtdb_database.py` after the parity
  check): every training and test sample's composition recomputed with the trained model's calls
  (`calls.tsv.gz`) from the taxa's depth and size in `.profile.log`, against the simulator's manifest;
  `model_logs/composition_accuracy*.tsv`, `composition_accuracy.txt`, `summary.txt`, a console line per read
  type.

## Check on a synthetic world

`scripts/run_validation.sh` (`evaluate.py`, `depth_check.py`): 40 species (`gtdb_like_lineages.py --seed 5`),
2 genomes each, background DNA drawn per species from 0.2-2 Mb (genomes 0.27-2.06 Mb, median 0.69 Mb, 168
markers), converted with the new converter; `db_all` (every species) and `db_held` (every 4th species left
out, 10); six samples of 200,000 pairs (2 × 150, HS25, fragments 350 ± 50), 12 species each, lognormal
abundances, 30% second strains (`simulate_metagenomes --seed 7`); protal with the mini database's shipped
model, default knob. The `--profile_only` rerun of `db_held`'s SAMs gave the same profiles byte for byte.

`results/evaluation.tsv`; shares are of read pairs (explained) and of cells (unknown):

| db | sample | absent | called (true) | explained: truth / protal | `?`: truth / protal | avg genome of the called: truth / protal | missing species: truth / at median / at lowest |
|---|---|---|---|---|---|---|---|
| all | s_1 | 0 | 12 (12) | 1.000 / 0.968 | 0 / 0.032 | 711,575 / 715,766 | 0 / 1.9 / 4.9 |
| all | s_2 | 0 | 12 (12) | 1.000 / 0.976 | 0 / 0.024 | 661,339 / 662,603 | 0 / 0.6 / 45 |
| all | s_3 | 0 | 12 (12) | 1.000 / 0.968 | 0 / 0.032 | 571,385 / 569,531 | 0 / 0.6 / 2.6 |
| all | s_4 | 0 | 12 (12) | 1.000 / 0.965 | 0 / 0.035 | 387,890 / 388,249 | 0 / 1.6 / 9.1 |
| all | s_5 | 0 | 12 (12) | 1.000 / 0.975 | 0 / 0.026 | 651,951 / 652,133 | 0 / 0.6 / 2.8 |
| all | s_6 | 0 | 12 (12) | 1.000 / 0.966 | 0 / 0.034 | 756,163 / 758,145 | 0 / 0.7 / 19 |
| held | s_1 | 3 | 9 (9) | 0.782 / 0.820 | 0.307 / 0.180 | 802,759 / 832,125 | 3 / 6.1 / 20 |
| held | s_2 | 3 | 9 (9) | 0.747 / 0.833 | 0.449 / 0.167 | 895,989 / 911,592 | 3 / 3.7 / 228 |
| held | s_3 | 5 | 7 (7) | 0.578 / 0.681 | 0.555 / 0.319 | 741,388 / 764,780 | 5 / 6.0 / 16 |
| held | s_4 | 2 | 10 (10) | 0.950 / 0.923 | 0.049 / 0.077 | 387,264 / 388,279 | 2 / 3.5 / 20 |
| held | s_5 | 2 | 10 (10) | 0.841 / 0.854 | 0.194 / 0.146 | 680,456 / 695,042 | 2 / 3.3 / 14 |
| held | s_6 | 3 | 9 (9) | 0.844 / 0.841 | 0.152 / 0.159 | 752,411 / 762,402 | 3 / 2.9 / 88 |

Where the error comes from (`results/depth_check.txt`, over the species called and present):

- **Genome sizes** are right: GenomeSize over the simulated genomes' length, median 1.0004 (10-90%
  0.992-1.008); the two genomes of a species differ by the markers each lost.
- **Depth** is 3% low with every species in the database: median 0.970 (10-90% 0.910-1.002), and so is
  the explained share (0.965-0.976), which puts 2.4-3.5% into `?`. It is not one deterministic loss:
  - By what the sample holds of a species (`results/depth_by_strain.txt`): its reference genome 0.979,
    another genome of it (0.5% divergent, its own 2% of markers lost) 0.953, two strains 0.957.
  - By gene length, among genes with 20 reads or more (`results/edge_correction.txt`): under 300 bp 0.904,
    300-600 bp 0.983, 600-1,200 bp 0.995, longer 1.004. Only short genes lose much; a constant loss per
    gene end fits none of it (fitted over all genes, 2.4 bases), and dividing each gene's depth by
    1 − c / length moves the species' depths by 0.4% (`edge_correction.py`).
  - A first reading of all genes (`results/depth_by_gene_length.txt`, 0.88 to 0.98-1.0 from short to long
    genes, which suggested 25-30 bases lost per gene) mixed in shallow genes: the median of a gene's few
    reads lies below their mean, and short genes have few. Among species with fewer than 10 reads per gene
    the median of the genes' depths is 0.967 and their length-weighted mean 1.024; with 30 or more, 0.973
    and 0.978 (`results/depth_estimators.txt`).
  The depth is a model feature, so it is not changed here. Relative abundances barely move with it (the
  species are 2-5% low alike); the absolute shares (`?`, the explained share) take the whole bias. The
  build's new step measures the species' depth over their reads' per read type at r226, the data a
  calibration of the composition would rest on.
- **Species the database lacks**: their reads land on present congeners within the depth's identity margin
  (0.08) and raise those congeners' depth (median 1.004, but 10-90% up to 1.35), so the explained share is
  too high where congeners are missing (s_2 +0.086, s_3 +0.103). The `?` is lower than the truth by more
  than that where the missing species' genomes are smaller than the called ones' (s_2: the missing cells'
  average genome 0.37 Mb against 0.90 Mb called): the unknown genome equivalents are counted at the called
  species' average size, an assumption no read can check. `AttributedShare` (reads) does not depend on it.
- **Average genome size** of the called species: within 0.6% with every species in the database, 0.1-3.7%
  high with species missing (the inflated congeners weigh more).
- **Missing species**: the median-depth estimate was 1.0-2.0 times the true number (3.7 for 3, 6.0 for 5,
  6.1 for 3), the lowest-depth one 4-75 times; with every species in the database the 3% bias reads as
  0.6-1.9 species.

The mini database's example (`examples/mini_db/run.sh`, 30,000 pairs of 3 species, all in the database)
still passes: abundances 0.192/0.480/0.291 against 0.2/0.5/0.3, `?` 0.037 (the same shortfall of the
depth, `results/mini.profile*`).

## Tests

- Unit tests: 468 pass (`ctest -j4`), among them `test_Composition.cpp` (the header line, the shares by hand,
  over-explained and empty samples, the priors' size columns by header, a profile on 1 and 3 threads that
  ends with `?` and `.profile.log`'s columns, `--no_unknown_share`).
- e2e (`tests/e2e/test_protal_e2e.py` on a mini database built with the new converter): 145 pass. Changed:
  `.profile` parsing skips `?`; reads of nothing in the database, a placeholder model, a header-only SAM
  and an unmet `--fdr` give `?\t?\t1`; the composition file, `--no_unknown_share`, and
  `--add_tables species_priors.tsv` (refused for another release, an old table without sizes, sizes that
  reach `.profile.log`) are new.
- Scripts: `test_mini_db.py` (18: the size columns against the release's metadata and reference,
  `--priors_only`, the size tiers, `--genome_length` ranges), `test_gtdb_build.py` (20), `test_profile_scripts.py` (11),
  `test_composition_accuracy.py` (4, numbers by hand, the same composition as `Composition.h`'s test).
  `test_gtdb_pipeline.py`: all 6 builds pass (260 s on 4 cores); the first now also checks the
  composition tables of the pe and se models (training and test samples, the species' genome sizes within
  5% of their genomes', every `?` between 0 and 1), `summary.txt` and the console line.

## Open

- Calibrate the depth's shortfall per read type on the next r226 build (`model_logs/composition_accuracy.txt`:
  species' depth over their reads', in the tables per sample and species) and decide whether the composition
  corrects it: a factor per read type, perhaps by the species' divergence from its reference (another strain
  read 2.5% lower than the reference itself here).
- Host reads count as unknown genomes of the called species' size: in host-rich samples the `?` is
  meaningless unless host reads are removed first; `AttributedShare` still holds.
- `--full_sam_header` SAMs carry no scanned counts (their header goes first).
- Seen in passing and fixed the same day: with `--map`, a relative `-o` (or `#OUTPUT_DIR`) put the folders
  the map did not name under `<o>/<o>/` (`Options::LoadFromMap` joined the output folder twice; new tests
  `SampleMap.ARelativeOutputFolderIsUsedOnce` and the e2e `test_a_relative_output_folder_is_used_once`).
- The website documents the profile format and does not know the `?` line, `.composition` or
  `--no_unknown_share` yet.

## Commands

```bash
SRC=~/protal-composition/src BIN=~/protal-composition/build \
    taskset -c 0-3 nice -n 5 bash docs/claude/2026-10-08-sample-composition/scripts/run_validation.sh ~/protal-composition/val
python3 docs/claude/2026-10-08-sample-composition/scripts/depth_check.py --manifest sims/manifest.tsv \
    --run all=run_db_all --run held=run_db_held
PROTAL=build/protal THREADS=4 bash examples/mini_db/run.sh --rebuild ~/protal-composition/example
```
