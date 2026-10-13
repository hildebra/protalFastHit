# A stride pick for the sparse-first gate: a read's sampled cores spread evenly

Date: 2026-10-13. Branch `audit-fixes`: the probe's stride rules in `574b5a4` (on `10a681b`), the stride gate
(`--sparse_first_pick`, [below](#the-stride-gate---sparse_first_pick-follow-up-same-day)) in the commit on `752f0e9`
that adds this sentence. Data: the r226 probe tables of `ccdbf52` (SLURM 24193130, `local/Perf0.8.2d/stage1_probe`, see
[the r226 probe report](../2026-10-12-stage1-probe-r226.md)) for the model; locally (WSL, 4 cores) the strain world of
[the strain false-negative study](../2026-10-11-strain-false-negatives/README.md) (`~/strainfn`) and the archaea world
of [the seeding profile](../2026-10-12-seeding-profile/README.md) (`~/seed`).

## The question

The user asked whether the sample of 1/C of a read's cores that `--sparse_first` looks up first, today chosen by the
cores' hash (`stage1_probe::Level`), could be chosen at database build time instead:

- (a) towards keys in conserved regions with few known alleles, and
- (b) evenly spaced along the gene.

Then, on adding stride rules to the probe: whether they are added at build time and stored in the database, so that
the C levels choose from precomputed scores.

## The answer in short

- **Spacing (b) needs no database.** The hash sample exists so that a read and an index agree, without a memory
  access, on which cores are sampled. That matters only for an index holding just those cores (the stage-1 idea,
  `--sparse_seeding`). `--sparse_first` looks its sample up in the full index, so a read may pick any 1/C of its
  cores: every C-th in read order (a stride) spaces them evenly. A stride depends on where the read starts, which only
  the read knows, so there is nothing to precompute or store. The probe measures every C from one run.
- **Conservation (a) does not fit the gate.**
  - A read cannot tell which of its cores are conserved without asking the index. The key map's 16-bit cells hold
    block offsets with no spare bit, so a flag would be a separate 4^15-bit (128 MiB) bitmap: one more cache miss
    per core, not per sampled core. At r226 a key-map probe is ~20 of the 175 ns of a lookup
    ([seeding profile](../2026-10-12-seeding-profile/README.md)), so the C = 8 gate's memory traffic would about
    double. That is the price of going to C = 4.
  - A key (a 15-base core) is shared by many genes and species, so "few alleles in species X" becomes a union over
    every copy holding it.
  - A sample concentrated in conserved stretches loses the reads of variable regions systematically. Those are the
    reads the alleles, polymorphic and ancestry features rest on.
  - It fits a real sparse index (the memory-saving form of `--sparse_seeding`). There the build chooses which entries
    to keep and a read probes all its cores. Per gene copy, one core per window of C cores (spacing) and, within the
    window, the one whose 31-mer covers the fewest sites of `strain_alleles.tsv` (conservation as a tie-break).
- **A read-side analogue of (a):** within each stride window, the core whose 31-mer has the best minimum base quality.
  It avoids sequencing errors rather than strain alleles, needs no memory access, and would matter most for ONT. Not
  built.

## What a stride could gain: a model of the r226 probe tables

With the hash sample, each core is in the sample on its own, so a pair with n good lookups (best flex score ≥ 14) has
none sampled with probability (1 − 1/c)^n: 7% at c = 8 for n = 20. A stride over a run of 20 good cores always
samples two or three. [scripts/model_fit.py](scripts/model_fit.py) fits the distribution of n among the reads with
records from all of a probe table's rules (T = 14, every c and m), then predicts the loss (the share of the reads with
records a gate at m = 1 gives up) for two cases: the good lookups in one run of the read (a read's overlap with a gene)
and scattered among bad ones (a fixed-size pick, hypergeometric). Real reads lie in between.

| Read type, c | Hash (observed) | Hash (fit) | Stride, good lookups in a run | Stride, scattered |
|---|---|---|---|---|
| pe, 2 | 0.37% | 0.56% | 0.13% | 0.53% |
| pe, 4 | 2.13% | 2.06% | 0.90% | 1.82% |
| pe, 8 | 7.10% | 6.59% | 2.17% | 5.00% |
| pe, 16 | 20.34% | 19.09% | 6.24% | 12.43% |
| HiFi, 8 | 0.83% | 0.83% | 0.43% | 0.82% |
| HiFi, 16 | 2.04% | 2.24% | 0.67% | 2.22% |
| HiFi, 32 | 5.29% | 4.73% | 2.59% | 4.61% |

- 2.5% of the pe pairs with records have 4 good lookups or fewer (of ~50), the HiFi reads 0.4%. No pick saves those
  reliably. At c = 4 they are most of the loss.
- The fit misses the observed hash loss by up to 2 points at c ≥ 16 (pe). The numbers are a model: the probe now
  measures the stride's rules directly (below).
- Outputs: [results/model_fit_r226_pe1.txt](results/model_fit_r226_pe1.txt),
  [results/model_fit_r226_pb2.txt](results/model_fit_r226_pb2.txt).

## What was built: stride rules in `--stage1_probe`

- `Stage1Probe.h`: `StrideOffset` and `StridePicked`, every c-th core of a mate in read order from offset
  ((cores − 1) mod c) / 2, so the cores left over are split between the read's ends (one more at the end when odd).
  A read (each mate, each chunk of a long read) with n cores has ceil(n/c) picks, at least one. Each `Mate` counts
  its picked lookups by c and score beside the levels. `Tally` keeps the stride's rules, and the taxa table gets the
  stride's columns `s<c>_t<T>_m<m>` after the hash sample's `c<c>_t<T>_m<m>`.
- Every rule, hash and stride, now has `lookups_left`: the share of the run's lookups (of cores in the index) a gate
  with the rule would make, the sample's of every read and all of those it keeps. It is counted per read, where the
  r226 report estimated it as 1/c + kept × (1 − 1/c). The stride takes ceil(n/c) of a read's n cores, more than n/c
  for short reads at large c, so recall must be compared at equal lookups left.
- `ChainAnchorFinder.h`: `LookupKmers` gives each lookup its k-mer's index in the read (`LookupPointer::probe_index`,
  in the struct's padding: still 48 bytes) and records the read's k-mer count; `FindSeeds` passes both to the probe.
  Only with `--stage1_probe`.
- `RunProtal.h`: `misc/<sample>_stage1_probe_stride.tsv` (the hash table's format), and the log line ends with the
  stride's c = 4, 8 and 16 at T = 14, m = 1.
- Tests: unit tests for the picks (ceil(n/c), c apart, ends balanced) and the stride's rules, taxa columns and lookups
  left; an e2e test (`test_the_probe_tallies_hash_and_stride_rules`): the records are those of a run without the probe,
  the c = 1 rules of both picks agree and make every lookup, and the taxa table's s1 and c1 columns agree.
- Scripts: `measure_performance.sh --stage1_probe` copies the stride table;
  [probe_summary.py](../2026-10-12-stage1-probe/scripts/probe_summary.py) takes `--pick stride`;
  [scripts/stride_compare.py](scripts/stride_compare.py) sets the two picks side by side, per c and at equal recall
  (the rule with the fewest lookups left at a recall of 99%, 99.5%, 99.8% and 99.9%).

## Local check

protal built in WSL from the working tree (`740de1f` plus these changes, and other sessions' then uncommitted
`--flex_indels` and zstd-frame hunks, since committed as `3cde4f4` and `c027c94`), 4 cores, niced. All unit tests pass (545, 3 skipped as usual; `Stage1Probe.*` 5); the e2e class
`SparseFirstGatesTest` passes (8, with the new test).
Runs: [scripts/run_local.sh](scripts/run_local.sh); outputs: [results/](results/) (`local_*_compare.md` from
`stride_compare.py`, `local_log_lines.txt`).

- **Strain world** (`~/strainfn`): 88,800 pairs of 532 test strains of 30 species (the 40 samples as one), against the
  representatives' database (`db_reps`), which has none of these strains; 88,597 pairs with records. The profile with
  the probe is byte-identical to the one without (33.3 against 33.4 s wall).
- **Archaea world** (`~/seed`): 100,000 simulated pairs, 2,887 with records. Most reads lie outside the marker genes,
  as in a real sample, so lookups left means something here. Only 43% of its cores are in the index, against 99% at
  r226 ([seeding profile](../2026-10-12-seeding-profile/README.md)), so it is not r226.
- **Long-read path**: the archaea world's first mates typed as PacBio reads (one chunk each), 2,451 with records.

Records lost (1 − recall) at T = 14, m = 1, hash → stride, and the lookups left in the archaea world's pairs:

| c | Strain world, pe | Archaea world, pe | Archaea pe, lookups left (hash / stride) | Long-read path |
|---|---|---|---|---|
| 2 | 0.00% → 0.00% | 0.80% → 0.00% | 53.1% / 54.1% | 0.78% → 0.00% |
| 4 | 0.02% → 0.01% | 2.91% → 0.00% | 29.6% / 31.3% | 3.47% → 0.20% |
| 8 | 0.62% → 0.02% | 7.97% → 0.62% | 17.7% / 19.7% | 13.38% → 1.10% |
| 16 | 7.67% → 0.11% | 20.26% → 1.39% | 11.1% / 14.4% | 33.86% → 2.08% |
| 32 | 27.72% → 1.17% | 40.80% → 9.59% | 6.8% / 10.0% | 56.59% → 12.53% |

- **At equal recall** (the rule with the fewest lookups left, over all c, T and m), archaea pairs: a recall of 99%
  takes the hash to c = 2 (53.1% of the lookups) and the stride to c = 8 (19.7%); 99.8% only the hash's c = 1 (100%),
  the stride's c = 4 (31.2%). The long-read path alike (99%: 52.9% against 19.7%).
- The stride makes a few more lookups at the same c, ceil(n/c) of a mate's n cores, most at large c (c = 16: 14.4%
  against 11.1%). It still wins at equal cost by a wide margin.
- Called species keeping less than 90% of their reads with records (strain world, 348 called): c = 8: 1 → 0;
  c = 16: 34 → 1; c = 32: 348 → 2.
- The archaea world's hash losses are close to r226's (7.97% against 7.10% at c = 8). The stride's are below even the
  model's "good lookups in a run" case (0.62% against 2.17%). The strain world's reads all lie inside marker genes
  with many good cores, so it loses little either way; its point is that strains missing from the database do not
  change the picture.

If r226 agrees, a stride gate at c = 4 to 8 would lose what the hash gate loses at c = 2 or less, with a third to a
half of the lookups.

## The stride gate: `--sparse_first_pick` (follow-up, same day)

The user then asked for the gate itself, with the names the build-script session had already wired up for
`--test-c-gates` rows (`752f0e9`: rows `<sample>__c<C>s`, added when `protal --full_help` lists `--sparse_first_pick`
and `SPARSE_FIRST_PICK`). The probe's rules above are in `574b5a4`.

- `--sparse_first_pick hash|stride`: one pick for every read type, read types and theirs (`pe:stride,pb:hash`), or
  `auto` (the default: the database's, hash where it records none). A map's `SPARSE_FIRST_PICK` column (`hash`,
  `stride`, `-`) sets one sample's. `--build` records the picks given in `build_settings.tsv` as
  `sparse_first_pick_<type>` (a protal from before skips these lines and gates by hash). Runs follow the database's
  and are warned when a gated sample picks otherwise; at C = 1 the pick makes no difference.
- The gate (`ChainAnchorFinder::SparsePicked`) takes `stage1_probe::StridePicked`, the probe's own pick: every C-th
  k-mer of a mate (of a long read's chunk) in read order. `PrefetchSparseKeys`, `SparsePrepare` and `SparseHits` use
  it; at C = 1 the gate's lookups are still reused for the seeding (`KeepSparseLookups`).
- `--sparse_seeding` stands for an index of the hash pick's cores, so a sample it would seed with the stride pick at
  C > 1 stops the run (exit 2).
- Logs: `sparse first: pe 8 stride, score 14, ...` (the pick named where it is stride) and
  `Sample <name> sparse first (1/8 of the cores by stride, score 14, 1 hit(s)): ...`.
- Tests: unit (`BuildSettings.SparseFirstPicksParsedWrittenAndReadBack`); e2e
  `test_each_pick_gives_up_what_the_probe_predicts` (for both picks the gate gives up exactly the pairs the probe's rule
  at C = 8, T = 14, m = 1 does not keep; at C = 1 the records of both picks are identical) and
  `test_runs_follow_the_databases_picks` (a database built with `--sparse_first_pick se:stride`: a run follows it,
  `--sparse_first_pick hash` and the map's column are warned, `--sparse_seeding` refused, a bad map value refused).

Local runs of the gate ([scripts/run_gate.sh](scripts/run_gate.sh), [results/local_gate.tsv](results/local_gate.tsv)),
the archaea world's 100k pairs and their first mates through the long-read path, against the ungated run's records
(2,887 and 2,484):

| C | Pairs: records, hash / stride | Sampled lookups, hash / stride | Long-read path: records, hash / stride |
|---|---|---|---|
| 2 | 2,864 / 2,887 | 879,030 / 897,758 | 2,465 / 2,484 |
| 8 | 2,657 / 2,869 | 221,317 / 253,627 | 2,146 / 2,457 |
| 16 | 2,302 / 2,847 | 110,487 / 155,184 | 1,637 / 2,433 |

The records lost are those the probe predicted (C = 8: 7.97% and 0.62%). The stride gate at C = 8 keeps more records
than the hash gate at C = 2, with 29% of its sampled lookups. The aligning times of these runs (0.1-0.5 s) are too short
to compare. The models of a database are trained at its gate's pick; the r226 builds still gate by hash, and
`--test-c-gates`' stride rows give the stride's F1 with those models.

## How to run it at r226

The same cluster run as the r226 probe, with a protal that has these rules:

```bash
bash scripts/measure_performance.sh --stage1_probe OUT_DIR DB pe:R1:R2 pb:HIFI
python3 docs/claude/2026-10-13-stride-pick/scripts/stride_compare.py --dir OUT_DIR/stage1_probe --sample pe1 --profile OUT_DIR/stage1_probe/pe1.plain.profile
python3 docs/claude/2026-10-13-stride-pick/scripts/stride_compare.py --dir OUT_DIR/stage1_probe --sample pb2 --profile OUT_DIR/stage1_probe/pb2.plain.profile
```

If the stride holds what the model says, the next step is a stride mode for `--sparse_first` (the gate's picks in
`SparsePrepare` and `SparseHits` by index instead of by level; `--sparse_seeding` keeps the hash, as it stands for an
index). Models trained at a gate would need retraining, since the reads given up change.
