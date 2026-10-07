# Why the r226 build idles for an hour: one soil paired-end point on one thread

2026-10-07, scripts of `audit-fixes` at `9ecf7ee` (`collect_training_data.py`, `scenarios.py`; the simulator
`src/RandomForest/MetagenomeSimulator.cpp`, `ArtIlluminaWrapper.cpp`). The question: the r226 v15 build (SLURM job
24027266, 84 threads, `protal0.7.8_r226_v15`) goes through phases of very low CPU and I/O. At 10:05, 3:24 after
its simulations started (06:41:19), the node ran one `simulate_metagenomes` at ~18% CPU with a load average of 1.1
(49 fifteen minutes before). The test set's simulations had ended at 09:34. Is this I/O, and would compressing the
temporary files help?

Sources:

- **The r226 v14 build's logs** (`local/v14`, job 24005369, `-t 64`, scripts of `14632a8`/`cdce3c0`):
  `training_data_simulation.log`, `test_data_simulation.log`, `training_data.log`. v15 runs the same design and the
  same seed (1), so its scenario samples have the same depths. Its own logs were not available for this report.
- **A model of the collector's thread split** (`scripts/pe_threads_model.py`): it imports the collector, builds
  the build's design points (`units_of` with `build_gtdb_database.py`'s default arguments) and calls the
  collector's own `pe_threads`. It then estimates each point's time on those threads, splitting them as
  `write_all_reads` does. One estimate is the collector's own (`pe_point_seconds`). The other also counts each
  sample's genomes: species midpoint x 1.4 strains, at G seconds per genome and sample and P seconds per 150 bp pair.
- **The cost of one genome step on one core** (`scripts/art_call_cost.sh`), measured in WSL (Core Ultra 7 258V,
  `taskset -c 0`) with art_illumina 2.5.8 from `~/micromamba/envs/protal-db-build` on 10 synthetic genomes of real
  size (`~/colprof/genomes`, from the 2026-10-05 report: 3.4 Mb median, ~78 contigs).

## Summary

**This is a scheduling tail, not I/O.** The process left running is `sc_soil_shallow_pe_p5000000`: 6 shallow-soil
samples of ~10,000 species each (~14,000 genomes with strains) at 5M read pairs. The collector gave it **one
thread**, because its cost estimate (`pe_point_seconds`: 20 s per sample plus 96 µs per pair) assumes the design's
~200 genomes per sample. On one thread the simulator does every genome of every sample one after the other:

1. inflate the genome into a plain copy;
2. fork and exec `art_illumina`, which reads that reference and writes a few hundred pairs;
3. compress those pairs into the sample's reads.

That is ~84,000 ART calls in a row. The ~18% in `top` is the simulator's own share (inflating, compressing). The
rest is the ART children, which live ~0.1 s each and so rarely show up in `top`. The node does one core of work.

The model reproduces v14 (fitted on two points, checked on the others):

| v14 point | Collection | Threads (log = model) | Wall in the log | Model |
|---|---|---:|---:|---:|
| `sc_soil_shallow_pe_p5000000` | training | 1 | 4:22:09 (4.37 h) | 4.36 h (fitted) |
| `sc_gut_pe_p20000000` | training | 3 | 1:00:51 (1.01 h) | 1.01 h (fitted) |
| `sc_soil_pe_p20000000` | training | 3 | 2:14:21 (2.24 h) | 2.05 h |
| `sc_soil_shallow_pe_p5000000` | test | 1 | 2:18:01 (2.30 h) | 2.19 h |
| `sc_soil_pe_p20000000` | test | 4 | 1:24:48 (1.41 h) | 1.01 h |
| `sc_gut_pe_p20000000` | test | 4 | 0:44:37 (0.74 h) | 0.59 h |

At 64 slots the model gives exactly v14's logged threads (training: gut 3, soil 3, the others 1; test: gut 4,
soil 4, host 2). The costs fitted from the two training points are G = 0.166 s per genome and sample and
P = 53 µs per pair, on a node core while the protal runs and the other collection share the node. The test
points ran while the training collection did too, so they come out 5-40% slower than the model.

In v14 this one point ended the training collection **1:08 after everything else** (the last long-read point
ended at 3:13:56). The test collection's equivalent ended 1:16 after its long reads. **v15 at 84 slots gives it
1 thread again** (the model: soil 4, gut 4, host 2, the others 1), so it should end around 11:00 (06:41 + 4:22).
After that come the last protal run (~6-10 min) and the models.

**Compressing the temporary files would not help.** Per genome, the temporary files are a plain copy of the
genome (~4 MB) and ART's FASTQ of a few hundred pairs (~100 KB), removed as soon as they are appended. The node
has 171 GB of free memory and iowait was 0.0%. Linux writes dirty pages back after ~30 s by default
(`vm.dirty_expire_centisecs`), so files this short-lived most likely never reach the SSD. `art_illumina` reads
only plain FASTA, so a compressed copy would have to be inflated again: more CPU on the serial chain. The reads
themselves are zstd level 3 already. The 2026-10-05 measurements put level 12 at 11% fewer bytes for 13x slower
compression, and the disk was not the limit in v15 (527 GB free).

**Where a genome step's time goes** (one laptop core, 10 genomes):

| Step | Cost |
|---|---:|
| inflating a genome with `gzip -dc` (the simulator uses ISA-L, several times faster; not measured here) | 65 ms |
| `art_illumina`, 1 kb reference (its start alone) | 9.6 ms |
| `art_illumina`, 3.4 Mb genome, 100 / 400 pairs | 130 / 103 ms |
| `art_illumina`, 3.4 Mb genome, 20,000 pairs | 709 ms (≈ 31 µs per extra pair) |

**ART's cost is in reading its reference (~30 ns per base), not in the pairs.** At soil depths a genome gets
~100-400 pairs per sample, so each call is almost all reference.

## Fixes, ranked

| # | Change | Effect at r226 | Effort |
|---|---|---|---|
| 1 | **`pe_point_seconds` counts the genomes**: each sample's genomes (species midpoint x strains) x G + pairs x P, with G/P at the node's ~3,000 pairs per genome | v15 (84 slots): soil_shallow 1 → 5 threads, soil 4 → 7, gut 4 → 2, host 2 → 1; the longest paired-end point **4.36 h → ~1.4 h** in training, 2.19 h → ~0.7 h in test. The training collection then ends with its long reads (v14: 3:14 instead of 4:22): **~1 h off a ~5 h build**, and ~80 idle core-hours of the allocation | ~15 lines + a test of `pe_threads` |
| 2 | **One pool of threads for all of a point's samples** in `simulate_metagenomes` (`write_all_reads`): any thread takes the next genome of any sample, and each sample is appended in its own order as now, so the reads stay byte-identical | Now the threads are split per sample once: with 5 threads for 6 samples, five run with one thread each, then the sixth alone. The scenarios' depths differ up to 3.6x between samples (0.52-1.89). soil_shallow on 5 threads: 1.38 h → ~0.9 h; and any thread count works without rounding to the samples | ~80 lines of C++ + a byte-identity test |
| 3 | **One ART call per genome and point**, the pairs dealt to the samples that have the genome (by a seeded draw) | The soil points' ART calls (84,000 per training point, 42,000 per test point, ~11.6 core-hours in both collections) fall ~2-4x: a point's samples share most of their species (13,148 in the table for 9,000-11,000 per sample). Reads change (statistically the same; new simulation keys) | larger: per-sample streams fed from shared calls |
| 4 | The scenario's long reads profiled before their paired-end point is done (they wait for its `simulated.json`, because their truth is its samples'; the design run's identical truth could serve) | v14: soil_shallow's PacBio and Nanopore reads waited ~3 h on the scratch disk and went into the last protal run (18 samples, 6:21). After fix 1 this is ~1 h less | small-medium |

**Fix 1 was implemented the same day** (`1b76405`) with the long-read simulator in C++
([report](../2026-10-07-long-read-simulator/README.md)): `pe_point_seconds` counts `genomes_per_sample` x 0.17 s
besides the pairs x 53 µs. Fixes 2-4 are open.

Fix 1 alone removes the tail; fix 2 makes it robust to the next estimate error. Fix 3 saves CPU, not the critical
path, once 1 and 2 are in. After them, the next tail is the long reads: v14's last one ended at 3:14, with many
"waits for room" lines on 235 GB of scratch. v15 has 527 GB, so it should wait less; its
`training_data_simulation.log` will show.

Not worth doing: compressing the temporary files, or higher zstd levels for the reads (see the summary); a copy of
the genomes on scratch. Inflating is ~10% of a step even with gzip, and ART's reading of the plain copy, which is
already on scratch, is most of it.

## Checking the running v15

```bash
pgrep -af simulate_metagenomes        # expected: -o .../points/sc_soil_shallow_pe_p5000000/sim ... -t 1
grep -m1 "paired-end design points" /hpc-home/hildebra/DB/protal/protal0.7.8_r226_v15/training_data_simulation.log
L=/nbi/local/ssd/24027266/protalDBbuild_1/training/points/sc_soil_shallow_pe_p5000000/simulate.log
a=$(grep -c "ART command" $L); sleep 60; echo $(( $(grep -c "ART command" $L) - a )) ART calls a minute
```

The model's G = 0.166 s is ~360 ART calls a minute; with the node to itself it may be faster. `grep -c "ART command"
$L` against 6 x ~14,000 tells how far it is.

The website is not affected: nothing here changes how protal is run.

## Reproducing

```
python3 scripts/pe_threads_model.py SCRIPTS --slots 64 --genome_seconds 0.166 --pair_seconds 53e-6 [--test]   # v14
python3 scripts/pe_threads_model.py SCRIPTS --slots 84 --genome_seconds 0.166 --pair_seconds 53e-6 [--test]   # v15
bash scripts/art_call_cost.sh ~/colprof/genomes ~/idle_tail_art   # art_illumina on PATH; one core: taskset -c 0
```

`SCRIPTS` is a folder with `collect_training_data.py` and `scenarios.py` of `9ecf7ee` (`git archive 9ecf7ee scripts`;
later collectors count the genomes themselves, fix 1).
The model reads and writes no files; it runs in a second.
