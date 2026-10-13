# What `--test-c-gates` costs with its defaults: rows, CPU and disk

Date: 2026-10-13. Scripts: `build_gtdb_database.py` and `collect_training_data.py` at `c027c94` plus the stride rows
(`--test-stride`, added the same day). An estimate, not a measurement: no r226 build has run with these defaults yet.

Data:
- The test collection of r226 v25 (protal 0.8.1, ungated, 84 threads on a 1.5 TB node; `local/v25/protal0.8.1_r226_v25/logs/`):
  `protal_cpu_test.tsv` (each sample's `aligning` CPU and wall seconds) and `test_data.log` (which simulations were
  streamed).
- The alignment's time per gate relative to ungated, from run 4 of the sparse-first report (r226, EPYC 9634, 8 threads,
  compute-bound, `f66d520`): [r226_epyc8_runs.tsv and r226_epyc8_sparse.tsv](2026-10-12-sparse-first/results/), the
  `aligning_s` column.
- `--flex_indels`' cost from [the flex report](2026-10-13-flex-indels/README.md), section 6.

## The combinations

With the defaults (database gates `pe:1,se:1,pb:8,ont:1`; `--test-c-gates` `pe:0/1/2/4,se:0/1/2/4,pb:0/4/8/16,ont:0/1/4/8`;
`--test-flex-indels` `ont:0/1/4/8,ultima:0/1/2/4`; `--test-stride` auto) the test samples are profiled in these ways
(`s` = stride pick, `f` = `--flex_indels`; the database's gate is the test set itself):

| Read type | Database | Other gates | Stride | Flex | Ways | Besides the database's |
|---|---|---|---|---|---|---|
| pe | 1 | 0, 2, 4 | 2s, 4s | - | 6 | 5 |
| se (Illumina) | 1 | 0, 2, 4 | 2s, 4s | - | 6 | 5 |
| se (Ultima) | 1 | 0, 2, 4 | 2s, 4s | 0f, 1f, 2f, 4f | 10 | 9 |
| pb | 8 | 0, 4, 16 | 4s, 8s, 16s | - | 7 | 6 |
| ont | 1 | 0, 4, 8 | 4s, 8s | 0f, 1f, 4f, 8f | 10 | 9 |

That makes 33 ways in all, 29 of them besides the databases' own. Without the stride rows it is 24 and 20.

Only samples written to disk are profiled at other gates. A streamed sample is read once from its pipe. The r226
v25 test set had 312 samples: 98 pe, 90 se, 62 pb, 62 ont. Of those, 260 were written: 90 pe, 78 se, 46 pb, 46 ont.
The 12 Ultima samples (3 simulations) were streamed. `--test-flex-indels` writes them. So 272 samples get the
other ways: 1,638 more sample alignments, 590 of them for the stride rows. All of it runs in the test collection's
protal runs, with the index loaded once.

## The time

Relative alignment time by gate. Ungated is 1. Measured unless marked.

| | C 0 | C 1 | C 2 | C 4 | C 8 | C 16 |
|---|---|---|---|---|---|---|
| pe (also se, Ultima) | 1 | 0.649 | 0.477 | 0.377 | 0.323 | |
| pb | 1 | | | 0.421 | 0.382 | 0.360 |
| ont | 1 | ~0.5 (not measured) | | 0.339 | 0.253 | 0.225 |

Assumptions:
- `--flex_indels` adds 5%: +2 to 8% of the scan, plus the reads it keeps.
- The stride pick costs what the hash sample costs at the same C. It looks up ceil(n/C) of a read's n cores rather
  than ~n/C, and keeps a few more reads.

The written samples' ungated alignment at r226 v25, in CPU seconds: pe 5,929, se 1,310, Ultima 4,816, pb 4,221,
ont 3,912. The added alignment, in the same CPU seconds:

| Read type | Other gates | Stride | Flex | Sum | × its own gate's |
|---|---|---|---|---|---|
| pe | 10,991 | 5,062 | | 16,053 | 4.2 |
| se (Illumina) | 2,428 | 1,118 | | 3,547 | 4.2 |
| se (Ultima) | 8,928 | 4,112 | 12,655 | 25,694 | 8.2 |
| pb | 7,517 | 4,909 | | 12,426 | 7.7 |
| ont | 6,227 | 2,315 | 8,592 | 17,135 | 8.8 |
| **all** | **36,092** | **17,517** | **21,247** | **74,856** | **6.6** |

That is ~20.8 CPU hours, 6.6 times the written test samples' own alignment at the database's gates (11,389 CPU s).

Wall time:
- v25's run of the written test samples kept ~36 of its 84 cores busy (15,463 CPU s in 427 s; many small samples).
- At that rate the extra alignment takes ~35 minutes: other gates ~17, flex ~10, stride ~8.
- That is on a build of ~2 h 45 min, whose test collection's runs the trainers wait for: about +20%.
- The trainers score 29 more tables, seconds each. The tables at other gates are written in parallel worker
  processes.

Why this is an upper bound:
- v25 ran protal 0.8.1. Since 0.8.2, WFA2 takes 26% (pe) to 39-48% (long reads) less time.
- The relative factors come from 0.8.2 code.
- The next node may differ (v25's was the slower one).

## Cheaper defaults, if needed

| Change | Saves |
|---|---|
| Leave out the ungated rows (0 and 0f): known from earlier runs, the slowest of all | ~29,400 CPU s (39%) |
| The Ultima samples: they are deep, written only for `--flex_indels`, and a third of the added time | `--test-flex-indels ont:0/1/4/8` saves 12,655 + their gate and stride rows (25,694 in all, 34%), and their ~50 GB on disk |
| `--test-stride none` | 17,517 (23%) |

## Disk

Each extra row writes its own SAM, profile and dumps. The build puts these at 10% of a sample's reads
(`KEPT_SHARE`), and the stride rows add ~21 GB of the total below.
- The written test reads were 48.4 GB at v25, at ~5.9 extra rows a sample: ~29 GB.
- The Ultima simulations (~50 GB) at 9 rows: ~45 GB.

That is ~74 GB on the scratch disk that the build's estimate (`disk_needs`) leaves out: it counts each sample's SAM
once. v25 ended with 279 GB free, so it would have fit there; a smaller scratch disk might not.
