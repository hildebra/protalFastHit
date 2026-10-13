# Sparse seeding: could a gated run hold only the gate's part of the index?

- **Date**: 2026-10-13.
- **Question (the user's)**: a run gated at C = 8 (the build default for PacBio) looks up 1/8 of each read's cores
  first. Does it then need the other keys of the map, or could that memory be saved?
- **Code**: `c448aeb` + the working tree's `--sparse_seeding` (this report). Release build with `-g` in WSL
  `~/sfseed`, removed afterwards.
- **Data**: two local worlds; no r226 run yet (the [command](#at-r226-the-run-that-decides) is below).
  - The archaea world of the [seeding profile](../2026-10-12-seeding-profile/README.md) (`~/seed`; its database has
    no working HiFi model, so its profiles are empty). Its 1M-pair community was read again as 300 Mb of HiFi reads:
    20,027 reads of 15,000 ± 3,000 bp, made by the sparse-first report's
    [sf_hifi_sim.sh](../2026-10-12-sparse-first/scripts/sf_hifi_sim.sh).
  - The small real world `~/realworld_sp_s1`: 6 GTDB families, 530 real genomes, its own HiFi models. Its 24 HiFi test
    samples were made again as in the [sparse-first report](../2026-10-12-sparse-first/README.md#hifi-small-real-world-accuracy):
    12 of 3 Mb and 12 of 30 Mb, with their truth.
- **Machine**: WSL2 on the Core Ultra 7 258V, cores 0-3, `nice -n 5`, 4 threads.
- **Scripts**: [run_local.sh](scripts/run_local.sh) runs the checks, [ss_compare.sh](scripts/ss_compare.sh) and
  [ss_accuracy.sh](scripts/ss_accuracy.sh). The read and species analyses are in [ss_read_truth.sh](scripts/ss_read_truth.sh)
  and [ss_lost_species.sh](scripts/ss_lost_species.sh). Tables: [results/](results/).

## Short answer

- **Today a gated run needs every key.** The gate is an early exit. A read it keeps is seeded from all its cores
  (`FindSeeds`), so every key is used for those reads and the whole index stays in memory.
- **Seeding the kept reads from the gate's cores alone barely changes the reads.** `--sparse_seeding` does exactly
  what an index of 1/C of the cores would do. On the archaea world's HiFi reads at C = 8 it lost 0.43% of the reads
  with a record against the plain C = 8 gate. The records' share on the read's true species was unchanged (99.3%).
  Aligning took 1.9 s instead of 4.3 s.
- **The calls change because the models do not know it.** On the small real world (models trained on fully seeded
  runs), 8 of 497 true species were lost (F1 0.7314 → 0.7239). Those species kept their reads and genes; their seed-count
  features fell ~8×, and their failed-candidate rate fell too, so the model scored them lower. PacBio models trained on
  sparse-seeded runs would be needed, as for the gate itself.
- **The memory at r226** (from the local level shares, which hashing makes 1/C):
  - an index of 1/8 of the cores would hold ~3.3 GB of values instead of 26.6 GB;
  - in today's layout the key map stays 3.2 GB whatever it holds;
  - with the genes (4.1 GB) and the rest, a HiFi run would need ~13 GB instead of ~36 GB, or ~10 GB with a compact
    key map.
- **Limits**: only runs whose samples are all gated so can use it, in practice HiFi-only runs. A paired-end sample at
  the default C = 1 needs the whole index. Saving the memory needs loader and build work
  ([below](#what-it-would-take-to-save-the-memory)). The decisive measurement is an r226 run.

## Why a gated run needs every key today

A core's level is the trailing zero bits of its hash (`stage1_probe::Level`, at most 5). The gate at C looks up the
cores of level log2(C) or more (1/C of them) and gives a read up unless one of those lookups is a near-exact hit.
A read it keeps goes through `FindSeeds` → `LookupKmers`, which looks up all its cores, so its seeds, anchors and
alignments are those of an ungated run. At C = 8 the HiFi reads kept at r226 (18%) are seeded from all of their cores.
So every key of the index can be needed by some read, and none can be dropped.

The sampled set is the same for every read and every gene, because the level depends only on the core. An index
holding only the cores of level ≥ 3 is therefore well defined:

- each of its keys holds exactly the values it holds now: entries, flex cells, the per-entry unique flags (a property
  of the k-mer in the whole database) and the ubiquity cap (per key);
- a run on it would look up exactly the gate's cores.

## The emulation: `--sparse_seeding`

With a gate, `--sparse_seeding` seeds the reads the gate keeps from the gate's cores alone:

- `ChainAnchorFinder::SetSeedLevel`; `LookupKmers` skips the cores below the level.
- That is what an index of only those cores would give, but the whole index is still loaded.
- For short reads, `KeepSparseLookups` (from [the C = 1 reuse](../2026-10-13-sparse-first-reuse/README.md)) now hands
  over the gate's lookups whenever the gate's level is the seed level, not only at C = 1.

After the index load it logs what such an index would hold for each C, from one pass over the key map
(`Seedmap::CountByLevel`, 0.88 s here). For the archaea world:

> Sparse seeding: an index of 1/C of the cores alone would hold, of the 48933202 entries in 67443358 slots (0.38 GB)
> of 23262723 keys: C 1: 48933202 entries (100.00%) … C 8: 6112781 entries (12.49%), 0.05 GB of values, 2907800 keys
> (12.50%) …; the key map stays 3.2 GB in this layout

Checks:

- **Paired-end at C = 1** (every core sampled, `pe100k`): with and without the option, the SAM records, the profile and
  the count lines are identical. The seeding line equals `c448aeb`'s (135,409 k-mers looked up, 110,635 in the index).
- **The counts add up**: the C = 1 row equals the index's own entry count, in both worlds.
- **The level shares are 1/C**: 12.49% of the entries and 12.50% of the keys at C = 8 in the archaea world; 12.48% and
  12.47% in the small world ([results](results/)).
- **New e2e test** `SparseFirstGatesTest.test_sparse_seeding`: the C = 1 records are the same with and without the
  option, and the log's C = 1 row is the index's entry count, falling at each larger C. The class's 7 tests pass;
  the unit tests: 540 pass, 3 skipped as always on this machine (debug-only, AVX-512).

## Archaea world, 300 Mb of HiFi: the reads

[ss_compare.sh](scripts/ss_compare.sh): the arms alternated for 3 rounds; medians below
([runs](results/archaea_hifi_runs.tsv)). "cCs" is `--sparse_first C --sparse_seeding`.

| arm | aligning s | wall s | reads given up | k-mers the kept reads looked up | candidate alignments tried | records |
|---|---|---|---|---|---|---|
| plain (no gate) | 11.01 | 13.29 | | 55.7 M | 750,891 | 11,093 |
| c8 (today's gate) | 4.33 | 7.28 | 14,493 (72.4%) | 15.8 M | 221,016 | 11,083 |
| c4s | 2.44 | 5.57 | 14,306 | 4.09 M | 34,809 | 11,076 |
| c8s | 1.86 | 5.10 | 14,493 | 1.97 M | 19,291 | 11,041 |
| c16s | 1.77 | 4.95 | 14,604 | 0.97 M | 14,993 | 10,936 |

The machine slowed over the rounds: plain took 6.8, 11.0 and 14.1 s. Within each round c8s took 0.36-0.60 of c8's
aligning time. With 1/8 of the seeds the long reads also have far fewer chance anchors, so 91% fewer candidate
alignments were tried.

**Read by read against c8** ([compare](results/archaea_hifi_compare.tsv)): records that are not secondary, as
(read, gene) pairs and as (read, taxon) pairs.

| arm | reads with a record | reads lost / gained | (read, gene) lost / gained | (read, taxon) lost / gained |
|---|---|---|---|---|
| c8 | 5,383 | | | |
| c4s | 5,383 | 5 / 5 | 18 / 11 of 11,082 | 10 / 10 of 5,408 |
| c8s | 5,360 | 23 (0.43%) / 0 | 52 / 10 | 33 / 8 |
| c16s | 5,286 | 97 (1.80%) / 0 | 153 / 6 | 111 / 2 |

**Against the truth** ([read truth](results/archaea_hifi_read_truth.tsv)): the simulator names a long read after its
source genome's place in the sample (`g<k>x_<n>`). That gives each record's read a true species, to compare with the
species of the record's taxon.

| arm | records | on the read's species | on its genus | reads with a record on their species |
|---|---|---|---|---|
| plain | 11,092 | 11,010 (99.26%) | 11,090 | 5,343 |
| c8 | 11,082 | 11,000 (99.26%) | 11,080 | 5,333 |
| c4s | 11,075 | 10,996 (99.29%) | 11,073 | 5,335 |
| c8s | 11,040 | 10,965 (99.32%) | 11,038 | 5,314 (−0.36%) |
| c16s | 10,935 | 10,870 (99.41%) | 10,935 | 5,240 (−1.7%) |

Sparse seeding loses some reads, but the records it keeps are as often on the right species, slightly more often in
fact. The reads lost are most likely those whose gene the read only grazes: one sampled core there is enough for the
gate, but an anchor needs two seeds on the gene.

## Small real world, 24 HiFi samples: the calls

[ss_accuracy.sh](scripts/ss_accuracy.sh): one run of all 24 samples per arm, the build's models
([f1](results/small_world_hifi_f1.tsv), [per sample](results/small_world_hifi_calls.tsv)).

| arm | TP | FP | FN | F1 | records |
|---|---|---|---|---|---|
| plain | 497 | 1 | 364 | 0.7314 | 20,044 |
| c8 | 497 | 1 | 364 | 0.7314 | 20,020 |
| c4s | 491 | 1 | 370 | 0.7258 | 19,993 |
| c8s | 489 | 1 | 372 | 0.7239 | 19,926 (−0.47%) |
| c16s | 487 | 1 | 374 | 0.7220 | 19,630 |

c8s gains no species and loses 8 that c8 calls ([lost](results/small_world_hifi_lost_c8s.tsv)). All 8 keep their genes
and their hits; only the model's probability falls:

| sample | species | probability c8 → c8s | genes | hits c8 → c8s |
|---|---|---|---|---|
| pb_b30000000_s_9 | Hydrotalea flava | 0.993 → 0.337 | 74 | 754 → 750 |
| pb_b30000000_s_12 | Hydrotalea flava | 0.998 → 0.441 | 1 | 1 |
| pb_b30000000_s_2 | Hydrotalea flava | 0.990 → 0.416 | 5 | 5 |
| pb_b3000000_s_8 | Hydrotalea flava | 0.985 → 0.458 | 2 | 2 |
| pb_b30000000_s_6 | Reichenbachiella sp033344935 | 0.895 → 0.276 | 2 | 2 |
| pb_b30000000_s_12 | UBA11889 sp002708205 | 0.795 → 0.107 | 3 | 3 |
| pb_b30000000_s_3 | Acidithrix sp021790495 | 0.676 → 0.297 | 1 | 2 |
| pb_b3000000_s_4 | Ferrimicrobium sp027409005 | 0.641 → 0.230 | 2 | 2 |

What moved, as means over the 502 true taxa of all samples:

| feature | c8 | c8s |
|---|---|---|
| `lu_per_kb` | 147.9 | 18.4 |
| `lsu_per_kb` | 137.2 | 17.1 |
| `failed_candidate_rate` (true taxa) | 0.539 | 0.105 |
| `failed_candidate_rate` (false taxa) | 0.717 | 0.262 |
| `mean_mapq` | 103.4 | 108.6 |

- The `lu`/`lsu` family counts a taxon's unique seeds, so it falls by C. Hydrotalea flava in sample 9: `lu` 73,102 →
  9,191, `lu_per_read` 97.0 → 12.3.
- With fewer chance anchors, fewer candidates fail, so `failed_candidate_rate` falls, as it did for the gate.
- MAPQ rises with fewer competitors.

These are the inputs of a model trained on fully seeded runs. They now come from a different distribution, which
explains a species whose evidence is unchanged dropping from 0.99 to 0.34. A database meant for sparse-seeded HiFi runs
would train its PacBio models on sparse-seeded simulations, as the build already does for the gates
(`--sparse-first`). How much F1 such models recover cannot be told here: it needs a build.

## What such an index would hold at r226

The r226 index has 2,903,557,366 entries: 26.6 GB of values at 48.69 bits a slot, plus a 3.2 GB key map, plus 4.1 GB of
packed genes. The run's peak is ~36.3 GB, so about 2.4 GB is everything else (tables, models, threads). Cores are
sampled by hash, and both local worlds hold 1/C of their entries at each C within 0.1 percentage point. So at r226:

| C | values | HiFi run with today's key map | with a compact key map |
|---|---|---|---|
| 1 (today) | 26.6 GB | ~36 GB | |
| 4 | ~6.7 GB | ~16 GB | ~13 GB |
| 8 | ~3.3 GB | ~13 GB | ~10 GB |
| 16 | ~1.7 GB | ~11 GB | ~8 GB |

- **Today's key map**: direct-addressed by core over all 2^30 cores, 2 bytes a key and 8 bytes per 8 keys. It is
  3.2 GB whatever it holds.
- **A compact key map** for 1/8 of the cores: at r226, 9.2% of the cores have values, so ~12 M sampled keys. A bitmap
  over the 2^30 cores (134 MB) with a rank directory, plus ~5 bytes a key, is ~0.2 GB. It costs a dependent cache line
  per lookup. The [memory audit](../2026-10-03-memory-audit/README.md) rejected this for the full index (~1.5 GB
  saved); for a sparse one it saves ~3 GB of ~13.
- **The genes stay whole**: a kept read can hit any gene.
- **The run's log will give the exact r226 figures**: the cluster run below prints them as the `Sparse seeding:` line.

## What it would take to save the memory

1. **A load that keeps one level.** When every sample of a run is gated at C or more with sparse seeding, the loader
   would keep only the cores of level ≥ log2(C):
   - count their slots in one pass over the key map;
   - allocate the packed values for those slots only;
   - give the dropped keys empty ranges in the key map while packing.

   The whole `index.prx` is still decompressed, because the levels interleave by core, so the load takes about as
   long. A build could instead write the sampled cores' values as a file of their own (more disk, a fast load).
2. **Models that match.** The build would train the PacBio models on sparse-seeded runs and record it in
   `build_settings.tsv`, beside the gates, so that runs follow it and are warned otherwise.
3. **Optionally the compact key map** above, for a further ~3 GB.
4. **Mixed runs** (a paired-end sample at C = 1 with HiFi samples) load the whole index, as now.

## At r226: the run that decides

`scripts/measure_performance.sh --test_sparse` now takes levels with an `s` (`--sparse_seeding`), and `0` for an
ungated run, which its help already asked for but its check refused. One run on the r226 HiFi sample gives, for c8
and the sparse arms against the timed runs:

- the reads given up, the records, the species lost and gained, and Bray-Curtis;
- the aligning time;
- the log's `Sparse seeding:` line with the r226 shares.

```
SPARSE_LONG="8 4s 8s 16s" sbatch --exclusive -c 32 --mem 100G --wrap \
    "bash scripts/measure_performance.sh --test_sparse OUT_DIR R226_DB pb:HIFI.fastq.gz"
```

The species there come from the r226 models, trained without sparse seeding, so the calls will move as they did
here. The read-level numbers are what to check: records lost against c8, and the reads given up. The F1 of
sparse-seeded models needs an r226 build with them (step 2).

## Conclusion

- **Today**: no memory can be saved. Every key serves the reads the gate keeps.
- **With a sparse index** at C = 8 for HiFi:
  - the reads barely change: 0.4% fewer reads with records against today's gate, their species as often right;
  - aligning is ~2× faster;
  - the index's values shrink 8× (~23 GB at r226), the run from ~36 to ~13 GB, or ~10 GB with a compact key map;
  - the PacBio models must be trained on sparse-seeded runs, because its seed-count features shrink by C, and the
    calls of models trained otherwise drop species that keep their reads.
- **Next**: the r226 run above. If its read-level loss stays near this, steps 1 and 2 above.
