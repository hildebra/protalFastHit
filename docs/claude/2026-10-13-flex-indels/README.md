# Flex cells scored with an indel per flank allowed (--flex_indels): ONT reads, the gate and sparse seeding

Date: 2026-10-13. Branch `audit-fixes`: the experiments on `ac4975e` with this change (`--flex_indels`) and another
session's `--sparse_seeding`, then uncommitted, since `78dce73` ([its report](../2026-10-13-sparse-seeding/README.md));
committed as `c027c94` on `3cde4f4`, built from `git archive` and tested there.

Data: the small real world `~/realworld_sp_s1` (six GTDB r226 families, 530 real genomes, 170 species; its
database and its paired-end test point `rl150_p100000`, 15 communities). ONT reads of those communities made by
`simulate_metagenomes` with pbsim3's qshmm model (QSHMM-ONT-HQ, the GTDB build's `--ont-setup` at four mean
accuracies), 40 Mb per sample, 75,100 reads per accuracy. The archaea world `~/seed` for the paired-end exactness
check. WSL on the Core Ultra 7 258V, cores 0-3 (AVX2; no AVX-512 here).

Commands: [`scripts/`](scripts/): `fi_ont_sim.sh` (the samples), `fi_ont.sh` (the arms and the read-level truth),
`fi_exact.sh`, `fi_speed.sh` (callgrind), `fi_kernels.sh` (the kernels' loops from the compiled code),
`fi_check_avx512.sh` (for a Zen 4 node). Tables: [`results/`](results/).

**Question.** Could the seeding of ONT reads tolerate errors in the flex part of the k-mer (the 8 bases either side
of the 15-base core), and at what cost? Then: implement it, with sparse seeding in mind, where a read has fewer
lookups.

## Summary

The flex part already tolerated a mismatch: a lookup takes the cells that share the most of the 16 flank bases with
the read, with no minimum. What it did not tolerate is an insertion or deletion, ONT's commonest error: the bases
beyond it are shifted and match by chance, so an indel next to the core cost ~5.8 of the flank's 8 bases (section 1).
`--flex_indels TYPES` (opt-in, per read type, a map's `FLEX_INDELS` column per sample) scores each flank as the most
of its bases equal in place and, from its first differing base on, the read's bases one place in or out. A single
indel then costs one base, as a mismatch does, and only an exact match still scores 16, so the unique seeds are the
same. It applies to the seeding, the `--sparse_first` gate and the stage-1 probe alike.

On ONT reads of the small real world (75,100 reads per error rate, every record checked against its read's source):

| | effect on the reads with a record of their own species |
|---|---|
| every lookup made (no gate, or the gate at C = 1, the GTDB builds' ONT default) | none (at most 11 reads of 9,500-12,900) |
| the gate alone at C = 8 to 32 | 12-31% of what the gate loses comes back, more the more errors |
| `--sparse_seeding` at C = 8 to 32 | 1-5% of its (much larger) loss comes back |
| any arm | the share of records on the read's species unchanged |

At the gate it does better than lowering `--sparse_first_score`: at 12.2% errors and C = 32 it keeps 742 more reads, 302
of them with an own-species record (0.41 per read kept), where score 12 keeps 777 for 268 (0.34) and score 11 1,152 for
361 (0.31).

Cost: the scan's loop takes 86 instead of 37 instructions per 8 cells with AVX2 and 51 instead of 24 per 16 cells with
AVX-512 (+6.1 and +1.7 instructions per cell). At r226's ONT soil sample (253 G flex cells, 18.4 T instructions on
Zen 4) that is about +2% of the run with AVX-512, about +8% on AVX2 nodes. In the small world (1.6 cells a lookup) the
seeding's instructions rose 2.6% and the run's did not move (−0.1%, from the candidates).
Without the option the outputs are those of `ac4975e`, byte for byte (paired-end and ONT, section 4).

What stays open: the AVX-512 kernel has not run (no AVX-512 here; `fi_check_avx512.sh` on a Zen 4 node), the small
world's index has 1.6 cells a lookup against r226's 160, where relatives compete in a block (the effect there is
untested), and the models are trained without it. Another session's `build_gtdb_database.py --test-c-gates` profiles
ONT and Ultima test samples with and without it at r226.

## 1. Why an indel costs a flank and a mismatch one base

A seed's k-mer is 31 bases: a 15-base core, looked up exactly, and 8 flank bases either side, the flex part
(`Seedmap.h`, `m_exact_k`, `m_flex_k`). A lookup scores each flex cell of the core's block by the bases it shares
with the read's flex part, position by position (`Seedmap::Similarity`), and takes the best cells as seeds, whatever
their score (`KmerLookupSM::GetFromLookup`). So a substitution costs one point and the true cell usually stays best.
An insertion or deletion shifts every base past it by one place, and those then match by chance (a quarter of them).
The unit test `FlexScan.OneErrorInAFlankCostsOneBase` measures the bases lost without the option, by the error's place
out from the core (20,000 cases, a fifth with a homopolymer next to the core):

| error at place | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| substitution | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| insertion | 5.78 | 5.11 | 4.39 | 3.74 | 3.02 | 2.23 | 1.52 | 0.73 |
| deletion | 5.76 | 5.18 | 4.51 | 3.72 | 2.96 | 2.28 | 1.49 | 0.76 |

With the option each of these scores 15 (16 where the read's flank is the gene's all the same).

Where it matters: a block with only the true genome's cell and chance cores elsewhere keeps the true cell best even
when it loses 6 bases. It matters where close relatives share the block (the shifted bases then decide among them),
and where a score is compared with a threshold: the gate keeps a read only if one of its sampled lookups scores
`--sparse_first_score` (14) or more.

## 2. The score and the implementation

Per flank (8 bases, place 1 next to the core), with `f` the first place where the read differs from the cell:

- equal: the bases equal in place (as before);
- deletion: the `f - 1` bases before `f`, and from `f + 1` on the read's base one place nearer the core (the cell base
  at `f` has none in the read);
- insertion: the `f - 1` bases before `f`, and from `f` on the read's base one place further out (the read's base at
  `f` is extra; the outermost cell base would need the read's base past the k-mer, which the key does not hold, and
  goes uncompared);
- the flank scores the most of the three; the cell scores the two flanks' sum.

A flank with any difference scores at most 7, so 16 still means an exact match, and `GetFromLookup`'s unique flags
(set only at 16) are unchanged. Both shifted read keys are the read's own flex key shifted by a base per flank
(`flex_scan::MakeIndelKeys`): no read sequence is needed. The first differing base is the left flank's lowest bit
(`x & -x`) and the right flank's highest (VPLZCNTD with AVX-512, a smear of the bit downwards with AVX2).

- `src/Hash/FlexScan.h`: `IndelKeys`, `MakeIndelKeys`, `IndelSimilarity` (the definition, scalar), and a template
  switch `kIndels` on the three kernels (`ScoreScalar`, `BestAvx2`, `ScanAvx512`); without it they compile as before.
  The AVX2 steps are functions of their own (a lambda inside a function compiled for AVX2 is not compiled for it).
- `src/Utilities/SimdLevel.h`: protal's AVX-512 level now also needs AVX-512 CD (VPLZCNTD) and BITALG (VPOPCNTW,
  each flank's bases counted in its 16-bit half). Every CPU with the level's VBMI2 has both (Ice Lake, Zen 4 and
  later), so no CPU should change level; `docs/installation.md` and the 2026-10-08 `check_avx512.sh` list them.
- `src/Hash/KmerLookup.h`: `KmerLookupSM::SetFlexIndels`, used by `GetFromLookup`, `BestScore` (the gate and the
  probe) and `FlexSimilarity` (the score of a key with one value, from the gene); copies keep it.
- `src/Core/ChainAnchorFinder.h`: `SetFlexIndels`; `src/RunProtal.h` sets it per sample (`Options::FlexIndels(index)`).
- `src/Options.h`: `--flex_indels` (`none`, `all`, or read types: `ont`, `pb,ont`), the map column `FLEX_INDELS`
  (`1`, `0` or `-`: the sample's read type's), the options summary line "flex indels:".
- `docs/running.md`: the option's row.

## 3. Tests

Unit (541 pass, WSL): `FlexScan.IndelSimilarityIsItsDefinition` (2 million cells and keys, a third of them a flank
shifted, against the base-by-base definition; never below the plain score; 16 only for equal),
`FlexScan.OneErrorInAFlankCostsOneBase` (section 1), the AVX2 and AVX-512 kernels against the scalar one with and
without the option (scores, best, ties; the AVX-512 tests skip here), `PackedIndex.WithFlexIndelsLookupsTakeTheCellsOfTheBestIndelScore`
(lookups through a packed index give the entries of the best indel score, unique only for the whole k-mer, a copy
keeps the setting), `PackedIndex.LookupsGiveTheSameSeedsAtEveryVectorLevel` (now with shifted keys and the option),
`Options.FlexIndelsByReadTypeAndSample` (values, the map column, a bad value). End to end: `OntTest.test_flex_indels`
(on the e2e ONT reads, 2% errors, 70% indels: every gene found once and at least as many as without; the map column
gives one sample the option's records and the other the plain ones; a bad value exits 2).
The commit `c027c94` (on `3cde4f4`), built from `git archive`: unit 543 and end to end 172 pass (WSL, 4 cores; the
AVX-512 tests skip).

## 4. Without the option nothing changes

`fi_exact.sh`, HEAD's protal (`ac4975e`) against the changed one without the option and with `--flex_indels none`:

| sample | SAM records (sorted, md5) | profile |
|---|---|---|
| paired-end, archaea world (`~/seed`, 100k pairs) | `9124952b06295b13`, all three | identical |
| ONT, small real world (4,931 reads at 7.8% errors) | `7be4956c4f186971`, all three | identical |

## 5. ONT reads of the small real world

`fi_ont.sh`: per error rate one protal run of all 15 samples per arm, every record (not secondary) checked against
its read's source genome (read names `g<k>x_<n>`, the community's k-th genome; the record's taxon's lineage from
`genome2tiid.tsv`). The database has no ONT model, so the profiles are not scored; the reads are.

Error rates measured from the simulator's statistics (`results/error_rates.tsv`; pbsim3's mean accuracy setting in
brackets, about a fifth of its reads error-free): 5.5% (0.97), 7.8% (0.93), 12.2% (0.88). The setting 0.99 gave 5.5%
again, so that run is a replicate of 0.97's; the model does not reach the 1-2% of today's R10.4.1 reads, where the
effect should be smaller, as it grows with the errors.

Reads with a record of their own species (of ~75,100 reads; `results/sweep1/reads_<acc>.tsv`,
`results/sweep3_0.99/reads.tsv`; arms: plain no gate, cC `--sparse_first C`, cCs with `--sparse_seeding`, f with
`--flex_indels ont`):

| arm | 5.5% errors (0.99) | 5.5% errors (0.97) | 7.8% errors (0.93) | 12.2% errors (0.88) |
|---|---|---|---|---|
| plain | 12910 | 12939 | 12353 | 9550 |
| flex | 12912 | 12938 | 12354 | 9556 |
| c1 | 12908 | 12937 | 12349 | 9544 |
| c1f | 12911 | 12938 | 12351 | 9555 |
| c8 | 12770 | 12798 | 12183 | 9340 |
| c8f | 12806 | 12829 | 12220 | 9405 |
| c16 | 12515 | 12548 | 11864 | 8904 |
| c16f | 12583 | 12619 | 11977 | 9079 |
| c32 | 11783 | 11764 | 11011 | 7984 |
| c32f | 11912 | 11905 | 11199 | 8286 |
| c8s | 12597 | 12628 | 11954 | 9049 |
| c8sf | 12608 | 12642 | 11969 | 9075 |
| c16s | 11958 | 11975 | 11205 | 8022 |
| c16sf | 11975 | 11991 | 11249 | 8083 |
| c32s | 10347 | 10335 | 9352 | 6095 |
| c32sf | 10377 | 10368 | 9401 | 6172 |

What the gate and sparse seeding lose against the ungated run, and the share the option wins back:

| reads lost against no gate: without / with the option (won back) | 5.5% | 5.5% | 7.8% | 12.2% |
|---|---|---|---|---|
| gate, C = 1 | 2 / -1 | 2 / 1 | 4 / 2 | 6 / -5 |
| gate, C = 8 | 140 / 104 (26%) | 141 / 110 (22%) | 170 / 133 (22%) | 210 / 145 (31%) |
| gate, C = 16 | 395 / 327 (17%) | 391 / 320 (18%) | 489 / 376 (23%) | 646 / 471 (27%) |
| gate, C = 32 | 1127 / 998 (11%) | 1175 / 1034 (12%) | 1342 / 1154 (14%) | 1566 / 1264 (19%) |
| sparse seeding, C = 8 | 313 / 302 (4%) | 311 / 297 (5%) | 399 / 384 (4%) | 501 / 475 (5%) |
| sparse seeding, C = 16 | 952 / 935 (2%) | 964 / 948 (2%) | 1148 / 1104 (4%) | 1528 / 1467 (4%) |
| sparse seeding, C = 32 | 2563 / 2533 (1%) | 2604 / 2571 (1%) | 3001 / 2952 (2%) | 3455 / 3378 (2%) |

- With every lookup made the option changes nothing that counts: in this world a lookup's block holds 1.6 cells on average,
  so the true cell is rarely in competition, and an ONT read has some 30-95 lookups with an intact core per kb of gene
  (2 in 9 of its places are syncmers, 14-43% of them with an intact core at these error rates).
- The gate keeps a read when one sampled lookup scores 14 or more. A lookup with an indel in a flank fell below that
  and now scores 15: the more errors and the fewer sampled lookups (larger C), the more reads come back.
- Sparse seeding loses many more reads (they lack seeds: an anchor needs two of one gene among 1/C of the cores), and
  the option wins back little of that.
- The records' precision does not move: the share of records on the read's species is 85-87% in every arm, with
  and without the option.

The gate's other knob, its score, as the alternative (12.2% errors; `results/sweep2_score/`, all four error rates
there; `scripts/fi_tables.sh` makes these tables):

| arm | reads kept by the gate | with an own-species record | more kept than at score 14 | more with a record | per read kept more |
|---|---|---|---|---|---|
| c16 (score 14) | 13616 | 8904 | | | |
| c16 (--flex_indels) | 14203 | 9079 | +587 | +175 | 0.30 |
| c16, score 13 | 13988 | 9001 | +372 | +97 | 0.26 |
| c16, score 12 | 14303 | 9053 | +687 | +149 | 0.22 |
| c16, score 11 | 14721 | 9118 | +1105 | +214 | 0.19 |
| c32 (score 14) | 11717 | 7984 | | | |
| c32 (--flex_indels) | 12459 | 8286 | +742 | +302 | 0.41 |
| c32, score 13 | 12122 | 8147 | +405 | +163 | 0.40 |
| c32, score 12 | 12494 | 8252 | +777 | +268 | 0.34 |
| c32, score 11 | 12869 | 8345 | +1152 | +361 | 0.31 |

For as many more reads kept, the option brings more of them with a record than a lower score does (C = 32: +302 for
742 against +268 for 777 at score 12); only score 11 gains more, for 55% more reads kept. At r226 a read
kept costs its whole seeding and its candidates' alignments, and most reads of a real sample are no marker gene's
(the soil ONT sample: 82% given up at C = 4), so the reads kept per read gained are what a gate is judged by. The
same holds at 5.5% and 7.8% errors (`results/sweep2_score/`).

## 6. Cost

The kernels' loops in the compiled code (`fi_kernels.sh`, GCC 13, `-O3`):

| kernel | cells a step | without | with | instructions per cell |
|---|---|---|---|---|
| `BestAvx2` | 8 | 37 | 86 | 4.6 → 10.8 (+6.1) |
| `ScanAvx512` | 16 | 24 | 51 | 1.5 → 3.2 (+1.7) |

The ties are found as before (their loops unchanged). At r226's real soil ONT sample (the sparse-first report's run 4:
253 G flex cells, 160 a lookup, 18.4 T instructions on the EPYC 9634) that is +0.43 T instructions (+2.3%) with
AVX-512, and +1.5 T (+8%) where the scan runs on AVX2 (Broadwell nodes). The seeding was 39 of the 173 s a thread of
that run's aligning, so a few seconds a thread, if the extra instructions are not hidden behind the scan's memory
waits (which they partly are).

Callgrind on the first 2,000 reads of an ONT sample at 7.8% errors, one thread (`fi_speed.sh`, `results/speed.txt`):

| arm | the run's instructions | the seeding's (`FindSeeds`, inclusive) | flex cells |
|---|---|---|---|
| without | 10.937 G | 250.6 M | 399,157 |
| `--flex_indels ont` | 10.922 G | 257.2 M (+2.6%) | 399,157 |

Here a block has 1.6 cells and takes one 8-cell step either way, so the option costs +49 instructions a block; at r226's
160 cells a block the per-cell figures above apply. The run's total fell by 0.1%: the candidates differ a little.

Wall times in the sweeps vary 2x between arms run back to back on this shared laptop and say nothing here.

## 7. Open

- The AVX-512 kernel: `scripts/fi_check_avx512.sh OUT DB ont:READS` on a Zen 4 node (unit tests with "this CPU:
  AVX-512", the option's outputs equal at both levels, the seeding's time with and without).
- r226: the effect where relatives fill the blocks, and its F1: another session's `--test-c-gates` rows profile ONT
  (gates 0, 1, 4, 8) and Ultima single-end test samples with the option, scored by models trained without it.
- If the gates C = 8 and above are used for ONT, the models would be trained with the option: `--build` does not
  record it yet (`build_settings.tsv` records the gates).
