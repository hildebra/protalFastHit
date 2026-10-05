# PacBio HiFi reads for training and benchmarks

- **Date**: 2026-10-02.
- **Code**: branch `audit-fixes` at `d11381f` with this work uncommitted: `scripts/hifi_reads.py` (new),
  `scripts/collect_training_data.py`, `scripts/build_gtdb_database.py`, `scripts/mini_db/test_mini_db.py`,
  `docs/databases.md`, `docs/databases.md`.
- **Data**: the operon world's build `~/opw/b_gn2` ([gene neighbour report](../2026-10-02-gene-neighbour-frequencies/README.md)):
  the communities of its first two paired-end test samples of 100,000 pairs (`rl100_p100000`, 69 and 25 genomes),
  30 Mb of long reads each, drawn by the collector; its training database and its PacBio model (trained on
  pbsim3's reads).
- **Machine**: WSL Ubuntu 24.04, Core Ultra 7 258V (6 cores), with another session's benchmark running (load 4-7),
  so wall times are inflated; CPU times are given beside them.
- **Commands**: [`compare_reads.sh`](scripts/compare_reads.sh) (templates by
  [`make_templates.py`](scripts/make_templates.py), reads by pbsim3 and by `hifi_reads.py`, profiles), run as
  `compare_reads.sh ~/opw/b_gn2 rl100_p100000 30000000 ~/hifi 2` and with `hifi:20000:10000:3` for a broad
  range of lengths; the unit test `MiniDbTest.test_hifi_reads`.

## Can pbsim3 simulate HiFi reads?

No, not directly. pbsim3 (3.0.5, the build environment's) simulates the subreads of multi-pass sequencing
(`--pass-num` 2 or more, also with `--strategy templ`), written as BAM (it needs `samtools`), and its README says
the HiFi reads are then made by PacBio's `ccs`, which is not in the environment. In one pass, as the collector
ran it for PacBio (`errhmm:ERRHMM-SEQUEL:15000:3000:0.999`), its error model writes every base quality as `!`:
the README says so, and every base of the operon world's PacBio training samples is Q0 (3,013,134 bases of one
sample; the Nanopore samples of its quality model: mean Q39.6, 60% Q30 or more).

That mattered: protal's `excess_*` features take a base's error probability from its quality (Q0: 1), so every
PacBio read's divergence beyond its qualities was about −1, and the PacBio models of every build and benchmark
until now were trained and tested on `excess_median` near −1, where real HiFi reads give about 0. A model trained
so meets real HiFi reads off its training range in that feature.

Alternatives looked at (not run here: neither is installed, and installing them is a download):

| | HiFi | ONT | quality values | speed |
|---|---|---|---|---|
| pbsim3 `--pass-num` + `ccs` (+ `samtools`) | the subreads of N passes, then PacBio's consensus | its own models | ccs's | pbsim3 simulates N times the bases, and ccs's consensus on top |
| [Badread](https://github.com/rrwick/Badread) ≥ 0.4.1 | error and quality models trained on HiFi reads (`pacbio2021`; its README: `--error_model pacbio2021 --qscore_model pacbio2021 --identity 30,3`) | `nanopore2018/2020/2023` (its default) | its models' | Python with alignments along each read; its README calls it slow and suggests running it in parallel; pbsim3 is C++ |
| SimLoRD (2016) | CCS reads by pass count | no | by pass count | RS II-era CCS, before HiFi chemistry |

## HiFi reads in Python

At the user's choice, the collector now makes its PacBio reads itself (`scripts/hifi_reads.py`, numpy), one read
of each template it draws (as before, so the read lengths, starts and the samples' bases are those of the pbsim3
runs); Nanopore reads stay with pbsim3's quality model.

- **Read quality by length** (the user's specification): a read's quality is Q50 up to 5 kb, Q30 at 25 kb,
  linear between, Q20 at 50 kb and lower beyond at that slope (`QUALITY_BY_LENGTH`), plus a normal deviation of
  the setup's SD (3), kept between Q10 and Q60.
- **Errors**: a base's error probability is the read's times its weight, scaled to the read's mean: a base of a
  homopolymer of n ≥ 3 weighs (n/2)², and each base has a lognormal factor of its own (σ 1); in a homopolymer an
  error is a base more or less of the run, elsewhere a substitution, insertion or deletion (30:35:35).
- **Base qualities**: each base's quality is its error probability (Phred, 1 to 93); an inserted base has that
  of the base it follows. So the qualities say how likely each base is wrong, as calibrated HiFi qualities do.
- **Setup**: `--pb_setup hifi:LENGTH_MEAN:LENGTH_SD:Q_SD` (collector) and `--pb-setup` (`build_gtdb_database.py`),
  default `hifi:15000:3000:3`; a pbsim3 setup still works. pbsim3 is needed only for ont (and pb with a pbsim3
  setup). A design point's key holds the setup and `hifi_reads.MODEL`, so points of the old reads are simulated
  again.

Calibration (unit-test templates with homopolymers, 13 Mb, one seed):

| read length | reads | drawn Q, mean | errors per kb | expected by the qualities per kb | ratio |
|---|---|---|---|---|---|
| 5 kb | 200 | 50.1 | 0.008 | 0.013 | 0.63 (8 errors in 1 Mb) |
| 15 kb | 200 | 40.1 | 0.130 | 0.122 | 1.07 |
| 25 kb | 200 | 29.7 | 1.375 | 1.410 | 0.98 |
| 50 kb | 80 | 19.5 | 13.26 | 13.33 | 0.99 |
| all | 680 | | | | 0.99 (60,319 errors, 60,732 expected) |

pbsim3 against `hifi_reads.py` on the same templates, both profiled with the build's PacBio model (trained on
pbsim3's reads), default lengths (15 ± 3 kb):

| sample | reads by | wall s, CPU s | Mb per CPU s | base Q, median | bases ≥ Q30 | present taxa: TP / FN | FP | present taxa's `excess_median`, median |
|---|---|---|---|---|---|---|---|---|
| s_1 (69 genomes) | pbsim3 | 16.6, 24.1 | 1.2 | 0 | 0 | 26 / 0 | 0 | −0.967 |
| s_1 | hifi_reads.py | 5.6, 5.9 | 5.1 | 43 | 0.977 | 27 / 0 | 0 | +0.014 |
| s_2 (25 genomes) | pbsim3 | 16.5, 24.7 | 1.2 | 0 | 0 | 3 / 0 | 0 | −0.959 |
| s_2 | hifi_reads.py | 4.2, 4.5 | 6.7 | 43 | 0.977 | 3 / 0 | 0 | +0.019 |

The reads' quality by length (from their qualities), default lengths and a broad range (20 ± 10 kb):

| read length | reads (15 ± 3 kb) | read Q, median (10-90%) | reads (20 ± 10 kb) | read Q, median (10-90%) |
|---|---|---|---|---|
| < 5 kb | 230 | 49.8 (46.4-53.9) | 273 | 50.0 (46.6-54.2) |
| 5-10 kb | 404 | 47.0 (42.5-50.6) | 587 | 47.3 (42.9-51.5) |
| 10-15 kb | 1,985 | 42.1 (37.9-46.2) | 714 | 42.5 (38.0-46.7) |
| 15-20 kb | 1,543 | 38.3 (34.1-42.3) | 706 | 37.7 (33.2-41.5) |
| 20-25 kb | 201 | 34.2 (29.7-38.3) | 490 | 32.6 (28.4-37.2) |
| ≥ 25 kb | 10 | 30.7 (27.7-33.1) | 688 | 27.2 (21.5-32.3) |

(The operon world's contigs are about 50 kb, so no read is longer.) `hifi_reads.py` takes 4-5 times less CPU than
pbsim3 for the same reads; present taxa's `excess_median` is about +0.015, a species' own reads exceeding their
errors by about one percent as on paired-end reads ([2026-10-01](../2026-10-01-f1-opportunities/README.md)),
instead of −0.96. On these deep samples the old model's calls hardly change (one taxon more found in s_1); the
PacBio models are to be trained again on the new reads all the same, and the PacBio results of the earlier
benchmarks ([v0.7.1](../2026-10-01-v071-benchmark/README.md), [v0.7.2](../2026-10-02-v072-benchmark/README.md),
which use the collector's default) rest on Q0 reads in training and test alike.

## Not done

- The `ccs` route and Badread's speed were not measured (both need installing).
- pbsim3's Nanopore reads (`QSHMM-ONT-HQ`) have 20% of their bases at Q93 (589,454 of 2,989,431), more than
  real R10.4 reads show; not looked into.
- The models trained on the new reads: the next build trains them (pb points are simulated again).
