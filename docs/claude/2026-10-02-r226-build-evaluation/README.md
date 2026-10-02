# The GTDB r226 database build: results, costs, and what to tune next

2026-10-02. The first complete r226 build-and-train run with the scripts of `ff57266` (progress
lines and `--scratch` of `ccebade`, gene neighbours in the database) and a protal binary of that
time (it reports v0.7.0). SLURM job 23865669 on node q512n5: 16 threads, `--scratch` on the node's
SSD (157 GB free), default design (`--samples 12`, read pairs 1,000-500,000, three read setups,
20-200 species per sample, strains 0.3/0.1, 20% of species and 32 clades held out, four read
types, an independent test set). The logs the user copied are in `local/protal0.7_r226_v1/`
(git-ignored): the console output was pasted into the session, plus `sacct`. Analysis scripts here:
`errors.py` (false positives and negatives, thresholds by depth, out of sample on the test set) and
`depth_growth.py` (taxa and fragments per sample by depth), run with the protal-db-build env's
Python on `local/protal0.7_r226_v1/model_logs`; their output is in `errors_output.txt` and
`depth_growth_output.txt`.

What has changed since (current `audit-fixes`, to be 0.7.3), so as not to suggest it again:
0.7.1 and 0.7.2's model features (divergence beyond base qualities, conservation pattern before and
after the MAPQ filter, gene-neighbour support, `--features normalized+adjacency`), depth knobs
(`--depth-knob-read-types pb,ont` by default), parallel profiling, the collector drawing long reads
itself (one `pbsim --strategy templ` per sample, `c40cec3`), PacBio HiFi reads by
`scripts/hifi_reads.py` (`9037c8a`), simulations in the background during the builds, genome
lengths in the genome table, a 9.7x faster gene-neighbour search, the trainer's previous-procedure
study off by default (`d11381f`), and `--progress-every` off by default.

## Results

| read type | evaluated on | F1 | sensitivity | precision | FP per sample |
|---|---|---:|---:|---:|---:|
| pe | species held out | 0.952 | 0.963 | 0.941 | 3.5 |
| pe | independent test set | 0.948 | 0.958 | 0.938 | 4.0 |
| se | species held out | 0.947 | 0.955 | 0.939 | 3.5 |
| se | independent test set | 0.947 | 0.956 | 0.938 | 3.9 |
| pb | species held out | 0.917 | 0.889 | 0.947 | 3.4 |
| pb | independent test set | 0.886 | 0.828 | 0.954 | 4.2 |
| ont | species held out | 0.925 | 0.909 | 0.941 | 4.0 |
| ont | independent test set | 0.904 | 0.870 | 0.940 | 5.9 |

The short-read models generalise: holding out whole genera to phyla costs at most 0.005 F1, the
test set of another design agrees with cross-validation, and archaea are found as well as bacteria
(sensitivity 0.97). Parity with protal is exact for all four. At GTDB scale the paired-end model is
limited by **precision** (632 false positives against 390 false negatives), unlike on the 900-species
benchmark world, where 0.7.1 was recall-limited (`docs/claude/2026-10-01-f1-opportunities`). The
long-read models were trained on the old collector's reads (pbsim3 per contig: mostly ~100 bp reads
at the shallow points, PacBio qualities Q0), which 0.7.3 replaces; their numbers say little about
real long reads.

## Where the paired-end model errs

**False positives grow with depth.** Absent taxa per sample grow about as depth^0.84 (2 at 1,000
read pairs, 388 at 500,000), present ones level off near 90. The share of absent taxa called falls
from 8% to 2.2%, but the false positives per sample rise from 0.2 to 8.6 (12.5 at the test set's
1M pairs). The training design stops at 500,000 read pairs and the test set at 1M; real samples
are often 5-50M. At the deepest point's rate, a 10M-pair sample would get ~130 false positives
(extrapolated from the power law, not measured). This is the largest risk for real samples.

**The best threshold depends on depth.** Best F1 threshold on species held out: 0.12 at 1,000 pairs,
0.10 at 5,000, 0.40 at 20,000, 0.74 at 100,000, 0.92 at 500,000; a single knob gains nothing (best
0.60: F1 0.9526 against 0.9518). Chosen on the training rows and applied to the test set's other
depths (interpolated in log10 depth), thresholds by depth raise test F1 out of sample:

| read type | test F1 at 0.5 | with thresholds by depth | FP | FN |
|---|---:|---:|---|---|
| pe | 0.948 | 0.954 | 283 → 115 | 191 → 290 |
| se | 0.947 | 0.955 | 278 → 155 | 195 → 236 |
| pb | 0.886 | 0.919 | 83 → 156 | 358 → 180 |
| ont | 0.904 | 0.917 | 118 → 132 | 277 → 214 |

This is the opposite of the benchmark world, where 0.7.2's depth knobs cost PacBio F1. protal's
depth knobs bin a sample by the digits of its fragments: here 1,000-20,000 read pairs all fall in
bin 2 (median 20-764 fragments), where the best threshold goes from 0.10 to 0.40, and the deepest
point is bin 4 (21,000 fragments), so a real 10M-pair sample (bin 5 or 6) would get the last
trained bin's knob.

**Sources.** Of the 632 false positives, 69% are next to a species or clade the training database
lacks (64% a held-out single species, 327 of them its congeners), and 97 of those have more than
10 fragments: the missing species' reads, an information limit 0.7.2's conservation features aim
at. 31% are near a species the database has, 147 of those 193 with a single fragment: cross-mapping
singletons, which a higher threshold in deep samples removes. 55% of all false positives are
single-fragment taxa.

**Misses are mostly strains.** 79% of the 390 false negatives (308) are species simulated from
another genome than the representative, which are 53% of the present taxa; missed ones have median
identity 0.968, found ones 0.987. All 35 misses with more than 10 fragments are such strains
(median identity 0.959): to the model a divergent strain looks like a relative of a missing
species. The rest are 1-2 fragment taxa in shallow samples.

## Classifier settings

- Leaves: 512 instead of 128 lowers the log loss (pe 0.0970 → 0.0914, se 0.1155 → 0.1120) with
  fewer false positives (632 → 576) at the same F1; the report says the data support larger trees.
  Long reads: no gain. 64 trees are enough (128-256 change nothing).
- Training samples: halving them changes the short-read log loss little (pe 0.0998 → 0.0970);
  for Nanopore it still falls (0.2275 → 0.2139 from 30 to 60 samples).
- Feature set: all 71 dump columns against the 28 normalized ones, pe F1 0.9530 against 0.9518.
- The previous-procedure study took 220 of the pe trainer's 285 s (off by default since `d11381f`).
- The shipped model (`scripts/random_forest.xml`, the collection model here) finds 1.1% of the
  archaea; it is only used before training, or by databases built without it.

## Time and resources

8:54:24 wall, 46.3 CPU hours of 142 (16 threads), peak memory 88 GiB (both index builds at once),
2.4 TB read and 0.99 TB written (sacct), scratch at most 23.7 GB.

| stage | time | in 0.7.3 |
|---|---:|---|
| start-up to the first console line | 0:04:32 | the same (interpreter start and the tool check, probably the cold scikit-learn import) |
| genome table | 0:02:43 | also counts genome lengths once |
| conversion | 2:46:09 | the same |
| gene neighbours | 0:15:39 | 9.7x less CPU |
| training database files | 0:04:52 | the same |
| training database build (level 3) | 0:40:09 | the same; simulations run during it |
| finished database build (level 19, in the background) | 1:20:41 | the same |
| training collection: paired-end simulation | 1:52:18 | genome lengths from the table |
| training collection: long reads (pbsim3 per genome) | 1:23:00 | one templ run per sample, HiFi by hifi_reads.py |
| training collection: protal on 480 samples | 0:13:31 | parallel profiling |
| test set collection | 1:01:28 | as above |
| training (4 models in parallel) | 0:05:55 | no previous-procedure study |
| parity checks | 0:00:53 | |
| `--add_model` x 4 | 0:21:45 | the same |

All 15 paired-end design points took 1:46-1:52, the 1,000-pair ones as long as the 500,000-pair
ones: the simulator read every genome of the table (17,248 files on `/hpc-home`) for its length at
each point, 15 points at once; the test set's 18 points then took 9 minutes from the page cache.
`c40cec3` removes this. The final database's level 19 makes its index 18.05 GB against the training
database's 18.55 GB at level 3 (2.7% smaller) for 21.5 more minutes, and the bundle 20.45 against
21.95 GB for 1:20 against 0:40 of build. Writing `database.protal` took 20-39 minutes per build, and
each `--add_model` rewrites its 22 GB, on `/hpc-home`.

The conversion was the largest single stage. Its process group's memory in the status lines
suggests: ~35 minutes of parallel spooling of the representatives' genes (5-13 GB), ~45 minutes of
one process (2.4 GB), the all-genome workers (41 GB at 1:30), then ~65 minutes of one process
(2.6 GB) joining the 79.5M sequences through zstd, from chunks spooled uncompressed (~86 GB) to
`.convert_tmp` on `/hpc-home`. The converter logs no times of its own, so this is inferred.

## Suggestions, most important first

1. **Train and test at real depths.** Add paired-end points at 2M and 10M read pairs (fewer samples
   there, e.g. 6, to bound cost and scratch: ~45 GB of reads for 6 x 10M pairs x 3 setups), and a
   5M point to the test set; for long reads, points of 1.5 and 6 Gb (real HiFi and Nanopore
   metagenomes are 5-30 Gb; the deepest training point is 150 Mb). Without them the false positives
   of real deep samples are an extrapolation.
2. **Depth knobs for pe and se too** (`--depth-knob-read-types pe,se,pb,ont`): out of sample +0.006
   to +0.009 F1 for short reads and 59% fewer paired-end false positives here; the trainer reports
   the test set's F1 at its knobs, which tells whether to keep them. Finer bins than whole digits
   (e.g. half digits) would follow the threshold where it changes most (1,000-20,000 pairs); that
   needs protal's `DepthKnobBin` and the trainer to agree on them.
3. **The conversion** (2:46, untouched by 0.7.3): time its phases in `convert.log`; let the workers
   write their chunks of the full reference zstd-compressed (concatenated zstd frames are one valid
   file, so the main process only copies ~5 GB instead of piping ~86 GB through one zstd); put
   `.convert_tmp` on `--scratch` when given.
4. **One rewrite for all models**: `--add_model` for the four read types in one call, or the trained
   models packed by the build itself; ~16 of the 22 minutes.
5. **The finished database's compression level**: 19 gains 2.7% on the index for 21.5 minutes, and
   with 0.7.3's background simulations its 1:20 build may become the longest path; level 9-12
   (to measure) or the training database's 3.
6. **Strains**: check in the next run whether 0.7.3's conservation features (a relative's reads
   concentrate on conserved genes, a strain's do not) lower the misses of divergent strains (here
   79% of the misses); the trainer's table of the conservation features by class shows it.
7. **Larger trees** (`--maxnodes 256` or `512`): better calibrated probabilities, which matter more
   with knobs by depth; F1 the same here.
8. **Long reads**: retrain on 0.7.3's collector in any case; more samples per point (24) are cheap
   now and Nanopore's learning curve still falls.
9. **Run settings for 0.7.3**: `--progress-every 600` to keep the status lines (off by default
   since `c40cec3`); memory: ~95 GB with both builds at once, ~50 GB with `--one-build-at-a-time`;
   scratch: 23.7 GB here, ~70-110 GB with the deep points of 1.
10. Once a 0.7.3 r226 model exists, ship its paired-end model as `scripts/random_forest.xml` (the
    shipped one calls almost no archaea).

## Follow-up: suggestions 1-8 implemented (2026-10-02)

At the user's request, on `audit-fixes` after `aeb7bf4` (branch `r226-tuning`), with 256 leaves for 7:

1. **Real depths.** `--read-pairs` now ends `2000000:4,10000000:2`, `--test-read-pairs` `5000000:2`,
   `--long-read-bases` `1500000000:4,6000000000:2`, `--test-long-read-bases` `3000000000:2`
   (`DEPTH:SAMPLES`: a depth with other samples than `--samples`). So that the deep points do not
   take as long as their samples one after the other, `simulate_metagenomes` writes a point's samples
   on `-t` threads (their designs and ART seeds drawn first, in order, so the files are byte-identical
   for any `-t`), and BGZF-compresses each sample as each genome's reads are appended (no plain FASTQ
   on disk; the same bytes as compressing the whole file, `bgzf::Writer`); the collector gives the
   costliest points the most threads, in proportion to samples x read pairs, largest first.
2. **Knobs by depth for every read type**, as a curve: a point per half decade of log10 of the
   sample's fragments, its knob fitted on the samples within half a decade, read linearly between the
   points and at the end points beyond them (`protal_depth_knob_curve`; protal still reads 0.7.2's
   knobs by decade). `--depth-knob-read-types` defaults to `pe,se,pb,ont`.
3. **The conversion**: the full reference's chunks are compressed by the workers (zstd frames, joined
   by copying, so no 86 GB is written or piped uncompressed), every step's time goes to `convert.log`,
   and `--tmp` (with `--scratch`, the scratch folder) holds the spool.
4. **One rewrite**: `protal --add_model A,B,C --read_type pe,se,pb` checks every model, then writes
   `database.protal` once; the build script adds all models in one call (`final_package.log`).
5. **The finished database's level**: `--final-db-level`, default 9. `compress_levels.sh` on the
   tuning world's database (raw index 3.4 GB; WSL, 6 threads, shared and busy, so the times are
   rough; `compress_levels_output.txt`):

   | level | database.protal | time |
   |---|---:|---:|
   | 3 | 128.9 MB | 58 s |
   | 6 | 120.2 MB | 80 s |
   | 9 | 118.4 MB | 45 s |
   | 12 | 117.8 MB | 91 s |
   | 15 | 117.0 MB | 261 s |
   | 19 | 106.8 MB | 321 s |

   This index compresses ~30x (synthetic and sparse), GTDB r226's ~2x, where level 19 gained 2.7% over
   level 3 for 11x the time; levels between were not measured at r226. Level 9 takes about level 3's
   time and most of what levels up to 15 gain here.
6. **Strains**: the trainer's report has a section "Strains" (training and test set): misses of
   species simulated from another genome than the representative against the representatives', by
   fragments, the identity of strains found and missed, and the conservation features of missed
   strains next to the absent taxa holding a missing congener's reads; the summary gives their FN rate.
7. `--maxnodes` 256 by default (trainer and build script).
8. **Long reads**: `--long-read-samples` 24 per point; a long-read point now replays the communities
   of the paired-end points of its depth in every read setup (3 x 12 = 36 to draw from), and a
   long-read sample's templates are removed as soon as its reads are written.

Scratch with these defaults: estimated ~100 GB at the peak (the deep points' reads ~70 GB, the
converter's spool ~35 GB before the samples); give it 150 GB.
