# Sparse first at C = 1: seeding the reads kept from the gate's own lookups

- **Date**: 2026-10-13.
- **Question (the user's)**: with `--sparse_first`, are a read's keys all looked up again once the gate keeps it? They
  are: a pair kept is looked up in full by `FindSeeds`, so at C = 1, where the gate looks up every key, a kept pair's
  lookups are all made twice. Implement taking them from the gate, and check whether it saves time.
- **Code**: base `5334f2c`; reuse = `5334f2c` + the change (`ChainAnchorFinder::KeepSparseLookups`, called in
  `Classify.h`'s single-end and paired-end loops). Release builds with `-g`, built in WSL `~/sfreuse/{base,reuse}`
  (removed afterwards).
- **Data**: the archaea world of the [seeding profile](../2026-10-12-seeding-profile/README.md) (WSL `~/seed`: every
  r226 archaeal species' marker genes, 0.4 GB of values); its 1M simulated whole-genome pairs (`pe1m`) and the first
  100,000 of them (`pe100k`), as in the [sparse-first report](../2026-10-12-sparse-first/README.md).
- **Machine**: WSL2 on the Core Ultra 7 258V, cores 0-3, `nice -n 5`, 4 threads (callgrind: 1 thread).
- **Scripts**: [run_all.sh](scripts/run_all.sh) runs [sfr_exact.sh](scripts/sfr_exact.sh), the sparse-first report's
  [sf_callgrind.sh](../2026-10-12-sparse-first/scripts/sf_callgrind.sh) and [sfr_speed.sh](scripts/sfr_speed.sh)
  (`bash run_all.sh BASE_PROTAL REUSE_PROTAL WORK`). Per-run timings: [results/speed.tsv](results/speed.tsv).

## The change

A core's level is the trailing zero bits of its hash (`stage1_probe::Level`, at most 5). The gate at C looks up the
cores of level log2(C) or more, so C = 1 looks up all of them. Short reads go through the gate in two steps:

- `SparsePrepare` makes the pair's sampled lookups with `KmerLookupSM::Get` in k-mer order and fetches their values.
- `SparseScore` stops at the first hit.

Before, a kept pair then went through `FindSeeds` → `LookupKmers`, which calls `Get` again for every k-mer. `Get`
depends on nothing but the k-mer, so at C = 1 its list is exactly the gate's: the same lookups, in the same order.
Now, at C = 1, `KeepSparseLookups` swaps the gate's list in as the read's lookups and marks the read prepared. It
counts its k-mers and lookups in the seeding counts as `LookupKmers` does. `FindSeeds` then sorts and scans them as
before.

- At C ≥ 2 nothing changes: the gate's lookups are a part of the read's, and the read is looked up in full.
- Long reads are unchanged. `SparseHits` stops at the first hit, chunk by chunk, so a kept long read has made only a
  few lookups.

## Exact

[sfr_exact.sh](scripts/sfr_exact.sh): base and reuse on `pe1m`, 4 threads, paired-end and single-end (R1 alone), at
C = 1 (the changed path) and C = 2. In all 4 cases these were identical between the two binaries:

- the SAM records, as a sorted set;
- the profile;
- the alignment, seeding and sparse-first count lines.

| case | pairs or reads given up | k-mers the kept reads looked up |
|---|---|---|
| pe, C = 1 | 962,184 of 1M (96.22%) | 1,689,817 |
| se, C = 1 | 968,131 of 1M (96.81%) | 713,865 |
| pe, C = 2 | 963,220 of 1M (96.32%) | 1,643,215 |
| se, C = 2 | 968,862 of 1M (96.89%) | 697,283 |

**A new end-to-end test**: `SparseFirstGatesTest.test_every_core_sampled_keeps_every_record` in
`tests/e2e/test_protal_e2e.py`. On the operon release's 8,000 pairs it runs C = 1 with score 0, paired-end and
single-end. In that setting a read is given up only if none of its cores is in the index, and such a read has no seeds
anyway. The test asks for the same sorted records and the same header (minus `@PG`) as an ungated run. That header
holds the unaligned reads' failed candidates.

The whole `SparseFirstGatesTest` class (6 tests, ~11 s) was run three times:

| binary | result |
|---|---|
| reuse | passes |
| base | passes |
| broken variant | the new test fails on paired-end: no records against 6,819 |

The broken variant is reuse with the swap removed from `KeepSparseLookups`, so the read is marked prepared without the
gate's lookups.

## Instructions

callgrind on the paired-end loop, 1 thread, `pe100k`, `--sparse_first 1` (score 14). Of the 100,000 pairs, 96,977
were given up and 3,023 kept; the kept pairs looked up 135,409 k-mers.

| arm | instructions | `FindSeeds` (inclusive) | anchor finder's operator (inclusive) |
|---|---|---|---|
| plain (no gate) | 4,150.5 M | | |
| base, C = 1 | 3,353.1 M | 60.2 M | 147.2 M |
| reuse, C = 1 | 3,338.6 M (−14.5 M, −0.43%) | 45.1 M (−15.1 M) | 132.1 M (−15.1 M) |

The saving is exactly the second `LookupKmers` pass: ~111 instructions per k-mer of a kept pair. Everything else in
the profile is unchanged to the instruction (`SparsePrepare`, `SingleValueScore`, k-mer taking, reading).

## Time

[sfr_speed.sh](scripts/sfr_speed.sh): `pe1m`, 4 threads, aligning only (`--no_profile`), the arms alternated for 7
rounds. The stages are from `misc/<sample>_runtime.tsv`, in seconds summed over the 4 threads.

| arm | wall s | aligning s | "Seeding" s | "Sparse first" s | records |
|---|---|---|---|---|---|
| plain (base, no gate) | 3.31 | 2.350 | 4.251 | | 35,692 |
| base, C = 1 | 4.24 | 3.013 | 0.234 | 8.167 | 35,692 |
| reuse, C = 1 | 4.34 | 3.108 | 0.193 | 8.595 | 35,692 |

The table gives medians. Per round, reuse minus base:

- **Aligning**: −0.445 to +1.024 s, median +0.016 s. This is noise: the same arm varies by up to 35% between rounds
  on this machine.
- **"Seeding" timer**: lower in all 7 rounds, by 0.002 to 0.180 s, median −0.030 s (−18%). That timer covers only the
  kept pairs' `FindSeeds`.

So the change removes what it should, but 0.03 CPU-seconds of a 3 s, 4-thread stage cannot be seen in wall time
here.

In this world, C = 1 is still slower than no gate (wall 4.24 against 3.31 s), as the sparse-first report found. The
reasons: 45% of the world's keys have one value, which the gate scores from the gene (`SingleValueScore`, about 17% of the
loop's instructions), and lookups are cheap in a world whose index sits mostly in cache. At r226, C = 1 cut aligning
CPU by about a third ([sparse-first report](../2026-10-12-sparse-first/README.md#at-r226-cluster-runs)).

## What to expect at r226

The sparse-first report's run 5 (EPYC, 8 threads, BGZF input; [results](../2026-10-12-sparse-first/results/r226_epyc8_bgzf_sparse.tsv))
had these figures at C = 1:

- 90.1% of 50.6 M pairs given up;
- 214.4 M k-mers looked up by the kept pairs, each looked up twice before this change;
- 4,127 G instructions for the run;
- 503 CPU s of aligning, of which the kept pairs' "Seeding" was 30.0 s.

At ~111 instructions a k-mer the second pass is ~24 G instructions (0.6% of the run). Its memory should be warm:

- the key map blocks were touched by `SparsePrepare` one pair earlier;
- in between, only the next pair's sampled key blocks are fetched, and the kept pair's scoring stops at its first hit.

The local result fits this: the timer fell with the instructions. At r226 the pass would then cost ~3-4 CPU s of the
503, under 1% of the alignment stage. That estimate has not been measured.

## Conclusion

- **The change is exact, and it removes the redundant pass**: −15 M instructions per 100,000 pairs here, and the
  kept pairs' seeding timer is 18% lower.
- **It saves no measurable time locally.** The pass was 0.4% of the loop's instructions, below this machine's noise.
- **At r226 a gain under 1% of the alignment stage is expected**, for paired-end and single-end samples at the
  default build gate (`pe:1,se:1`). HiFi (C = 8) and ONT (long-read path) are unchanged.
- **It is kept** because it is a few lines and costs nothing. It does not need a cluster run of its own; the next
  `measure_performance.sh` run at r226 will show it as a lower "Seeding" stage at C = 1.
- **Not done**: at C ≥ 2 the kept read could take its sampled lookups too and look up only the other k-mers, merged in
  k-mer order. That would save 1/C of this pass, so at most ~0.3% at C = 2.
