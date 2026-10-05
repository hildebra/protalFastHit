# Outputs independent of record order and hash-map order

Date: 2026-10-05. Branch `audit-fixes`, base commit `48e8cd0`; the change described here sits on top of it (uncommitted
when this was written). Follows up two order dependences noted on 2026-10-04 in
[the GTDB-scale performance report](../2026-10-04-performance-gtdb-scale/README.md) (section on the fourth cluster
run and the binary gene table, and `results/followup_gene_table_bin.txt` there).

Data: the v0.7.3 benchmark world in WSL (`~/bench071/V073/protal_db/database.protal`, 765 species, text gene tables,
no `gene_table.bin`); 500k simulated read pairs (`~/bench071/samples/points/rl150_p500000/sim/reads/rl150_p500000_s_1_R{1,2}.fq.gz`)
and 90 Mb of simulated PacBio HiFi (`~/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz`).
Every run used `-t 6 --no_qcmsa`. Machine: WSL Ubuntu 24.04 on a Core Ultra 7 258V, 6 vCPUs, gcc 13.3, Release builds.

## Summary

- **The varying value** was `ANISum` (column 10) of `<prefix>.profile.gene.log`: the sum of the records' ANIs per gene,
  added up in doubles in record order. Multi-threaded alignment writes the reads in a different order each run, so the
  sum's last bits changed: 25.932187 or 25.932188 for taxon 176, gene 90, two runs each out of four at `48e8cd0`.
  Profiling the same SAM with its reads shuffled showed the same difference. More of the profile depended on record order
  the same way, below the printed digits: the taxon's ANI sum (the model feature `mean_ani`), the genes' identity-weighted
  bases (`BaseIdentity`), the gene map's iteration order (it decided the order of the floating-point sums behind
  `VCovStdDev`, `GeneVariance` and `DepthCV`), the sum behind the abundances, the row order of `misc/unreported_species.tsv`,
  and, on a tie, which allele a strain MSA called.
- **The genome map's order did not reach any output.** A `tsl::sparse_map` with `std::hash` (the identity on integers)
  iterates in taxid order whenever the taxids are dense, however the genomes were inserted and whatever was reserved; a
  rebuild in descending order with `reserve(4n)` iterated `1 2 3 4 5 ...` as before. With a scrambling hash forced on the
  genome map and on the profile's taxon and gene maps, `48e8cd0`'s outputs changed only by the `ANISum` digit above, and
  the change's not at all. The 13 PacBio alignments that moved on 2026-10-04 were most likely the long-read window change
  of `50ecb4b` (item 4a), as that day's follow-up comparisons (`det.sh`, `det2.sh`) concluded; the comparison build then
  was `c94bd0e`, from before 4a. The comments that said outputs break ties by the genome map's order are corrected.
- **After the change**, four runs of the same build give identical SAMs (as sorted sets) and identical other outputs
  (`diff -r`, without `*_runtime.tsv` and `*.statistics.tsv`), for paired-end and PacBio alike; so do the scrambled-hash
  builds and the profiles of read-shuffled SAMs. The outputs changed once, deterministically: `ANISum` settles on
  25.932188 (2 of the 4 runs of `48e8cd0` printed that value); nothing else printed changed. The profiling stage costs
  the same (alternated A/B below).

## The changes

| Where | What |
|---|---|
| `src/Utilities/ExactSum.h` (new) | `ExactSum`: doubles summed in fixed point (units of 2^-64) in a `__int128`, so the sum is the same in any order; `Value()` is the exact sum, rounded once. Exact for values of 2^-11 or more (identities, ANIs, bases times identity), finite values below 2^62. Follows the `kShareUnit` precedent (integer sums of `adjacent_support`) |
| `src/Profiling/Profiler.h` | `profiler::Gene::m_ani_sum`, `m_identity_bases` and `Taxon::m_ani_sum` are `ExactSum`; `BaseIdentity` sums the genes' `ExactSum`s. `GeneMap` is `std::map<uint32_t, Gene>`: a taxon's genes go in id order wherever they are iterated (depth vectors, their variance, the outputs). `Taxon::AddSam`/`AddHit` look the gene up once per record (`find`, `emplace`), which keeps the profiling stage's cost as it was. `PassingDepth` (the abundances' denominator) and the profile's `ToString` go over `SortedTaxa()` |
| `src/RunProtal.h` | `misc/unreported_species.tsv` rows in taxid order; `ExtractTaxa` (the order the strain MSAs take the taxa in) breaks ties in sample count by taxid instead of a `robin_map`'s order |
| `src/SequenceUtils/Variant.h`, `VariantHandler.h`, `src/Profiling/Strain.h` | `Variant::AlleleBefore`: a site's alleles ordered by what they are (bases, then deletions, then insertions; by base or by their bases). `RanksBefore` (the strain MSA's base call and IUPAC codes) and the bin sort in `PostProcessSNPBin` break ties with it, not with the order the reads first showed the alleles |
| `src/SequenceUtils/GenomeLoader.h` | `SortedKeys()`; `WriteSamHeader` (`--full_sam_header`) and `LoadAllGenomes` use it; comments on `m_genome_order` and the binary loader corrected (gene_table.bin keeps reference.map's order so that the same tables give the same bytes; no output depends on it) |
| `src/SequenceUtils/GeneTableFile.h`, `tests/test_GeneTableFile.cpp` | the same comment corrections |
| `src/Build.h` | the genera compared for the gene congener table list their members by taxid, not in the genome map's order (a build-time output; no effect on runs) |
| `tests/test_ExactSum.cpp` (new) | the sum is the same in shuffled and reversed orders where a double sum is not; parts add up to the whole; small sums exact |
| `tests/test_ProfileThreads.cpp` | `ProfileSam.TheProfileDoesNotDependOnTheOrderOfTheReads`: a 3000-read SAM and the same with its reads shuffled (each read's records together), profiled on 1 and 3 threads, give the same profile in an order-free dump (taxa by taxid, all model features, genes by id with their sums, sorted read identities, coverage, variant sites by position) |
| `docs/development.md` | one sentence on the new test |

Alignment itself needed nothing: each read's records depend only on the read (the candidates' order comes from the index
and the seeds), and the SAMs of repeated runs were the same sets before and after.

Unit tests: 349 pass (2 skipped benches). Sensitivity: with `ExactSum` patched back to a plain double accumulator
(`scripts/sensitivity.sh`), `ProfileSam.TheProfileDoesNotDependOnTheOrderOfTheReads` and the three `ExactSum` tests
fail ([results/sensitivity.txt](results/sensitivity.txt)).

## Results

Before the change ([results/experiment_base.txt](results/experiment_base.txt); `base` = `48e8cd0`). The `base-reorder`
lines in that file are from a first experiment build that only rebuilt the genome map in descending taxid order with
`reserve(4n)`; it iterated in taxid order as before, which is how the identity-hash behaviour was found.

| Comparison | pe | pb |
|---|---|---|
| run 1 vs runs 2, 3, 4 | SAM same; `.profile.gene.log` differs in 1 line in runs 2 and 4 | all same |
| profile of run 1's SAM vs of the same SAM, reads shuffled | `.profile.gene.log` differs in that line | same |
| run 1 vs the scrambled-hash build ([results/experiment_work.txt](results/experiment_work.txt)), 2 runs | SAM same; run 2 differs in that line | all same |

After the change ([results/final.txt](results/final.txt), [results/experiment_work.txt](results/experiment_work.txt)):

| Comparison | pe | pb |
|---|---|---|
| run 1 vs runs 2, 3, 4 | all same | all same |
| run 1 vs the scrambled-hash build, 2 runs | all same | all same |
| profile of run 1's SAM vs of the same SAM, reads shuffled | same | same |
| `48e8cd0` runs 1-4 vs the change's run 1 | SAM same; outputs same as runs 2 and 4, `ANISum` 25.932188 where runs 1 and 3 had 25.932187 | all same |

Profiling stage (`--profile_only` on `48e8cd0`'s run-1 SAM, 6 threads, 5 alternated rounds, "Profiling took"):

| | `48e8cd0` | first version of the change (two lookups per record) | the change |
|---|---|---|---|
| pe (456k records) | median 741 ms / 828 ms | median 778 ms (+5%) | median 841 ms (+1.6%, within the spread 822-950) |
| pb (36k records) | median 282 ms / 278 ms | median 286 ms | median 278 ms |

The two `48e8cd0` medians are from the two A/B sessions (the machine was slower in the second); compare within a column
pair. The first version looked each record's gene up twice in the `std::map` (`contains`, then `at`); one `find` per
record brought the cost back to that of the hash map. At GTDB scale profiling took 6 s on 32 threads (2026-10-04 third
cluster run); the gene map is per taxon (at most a few hundred genes), so its size does not grow with the database.

## Commands

The scripts are in [scripts/](scripts/); they build in `~/det-order/<name>/tree` and write runs to `~/det-order/runs`
(not in git; delete freely).

- `build.sh <name> <commit|tree> [reorder]`: a Release build with tests, from `git archive` or the working tree;
  `reorder` applies `scramble_maps.pl` (experiment only: a scrambling hash for the genome map, the profile's taxon map and,
  where it is a hash map, the gene map; the genomes inserted again in descending taxid order).
  `builds.sh` makes `base-reorder`, `work`, `work-reorder`.
- `runs.sh <name> [N]`: N runs (default 4) of pe and pb as above. `compare.sh A B`: SAMs as sorted record sets, the rest
  with `diff -r -x '*_runtime.tsv' -x '*.sam.zst' -x '*.sam' -x '*.statistics.tsv'`.
- `profile_shuffled.sh <name> <run dir> <pe|pb> <out>`: the SAM with its reads shuffled (`shuf` with a fixed random
  source; each read's records kept together), both profiled with `--profile_only`.
- `experiment_base.sh`, `experiment_work.sh`, `final.sh`: the comparisons above. `check_tests.sh`: the unit-test summary
  and the varying value per run. `sensitivity.sh`: the test check with plain double sums. `profile_ab.sh`: the A/B.

## What is left

- `GenomeLoader::m_genome_order` and the binary loader's insertion order are kept: they make `gene_table.bin` byte-for-byte
  reproducible and keep the two loaders alike (`TheGenomeMapIsBuiltInTheSameOrder`), but no output depends on them now.
- Other hash maps keyed by taxid (`SampleContext`, the record evidence, the profiler's per-taxon jobs) were checked: their
  users sort (taxid order, `std::map`), count or look up, so their order does not reach an output.
- At GTDB scale a profile's taxon map holds a few thousand of 143k taxids, so its buckets wrap and it iterates in insertion
  order; the loops over it that had an order-dependent result (the abundance denominator, `unreported_species.tsv`) now
  sort, and the scrambled-hash build exercised that case here.
