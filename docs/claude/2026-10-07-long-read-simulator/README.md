# Long and Ultima reads made in simulate_metagenomes; streaming reads into protal; porting ART

2026-10-07, `audit-fixes` at `9ecf7ee` plus the changes described here, committed as `1b76405`. Three questions
followed the [idle-tail report](../2026-10-07-build-idle-tail/README.md):

1. Move the PacBio and Ultima read simulation (`scripts/hifi_reads.py`, with the collector's Python template drawing)
   into the `simulate_metagenomes` binary, so that no intermediate files are written and the reads come out as the
   compressed FASTQ that protal reads. Consider the Nanopore simulation (pbsim3) there too.
2. Could the simulated reads go straight into protal, so that the build needs no scratch disk?
3. Could ART's Illumina simulation be ported into `simulate_metagenomes`?

Sources: the code at `9ecf7ee`; pbsim3 3.0.x source (`src/pbsim.cpp`, GitHub `yukiteruono/pbsim3`, GPL-2.0, cloned to
the session scratchpad); ART MountRainier 2016.06.05 source (`artsrcmountrainier2016.06.05linux.tgz` from NIEHS, md5
matching the bioconda recipe); measurements in WSL (Core Ultra 7 258V, `taskset`, nice 5) with art_illumina 2.5.8,
pbsim3 3.0.5 and its models from `~/micromamba/envs/protal-db-build`, on the 400 synthetic genomes of real size of
the [2026-10-05 report](../2026-10-05-collector-profiling/README.md) (`~/colprof/genomes`). Build tree `~/lrsim`.

## Summary

- **Done: long and Ultima reads are made in `simulate_metagenomes`.** A new module (`src/RandomForest/LongReadSimulator`)
  draws the templates as the collector drew them and makes each into a read with a ported model: `hifi_reads.py`'s
  HiFi model, its flow model (Ultima), and pbsim3's qshmm model in template mode (Nanopore; it reads pbsim3's
  `.model` file, pbsim3 is not run). Each read is written as it is made, in zstd frames (or BGZF blocks) appended in a
  fixed order, so the files are the same for any number of threads. The collector runs it once per long-read design
  point (`--long_samples`). Gone: the template files (plain, gzipped per chunk), pbsim3's own files, the chunking
  (`--long_read_chunk`), the join, the renaming pass and the Python worker processes for long reads.
- **The models match their originals** on the same templates (`--long_templates`, below). Quality distributions and
  error rates agree within sampling noise. They are 2-25x cheaper per Mb: HiFi 0.033 against 0.086 CPU s per Mb, Ultima
  0.057 against 0.109, Nanopore 0.019 against pbsim3's 0.467 (before the collector's renaming pass).
- **Also done: the paired-end thread split counts genomes** (fix 1 of the idle-tail report): `pe_point_seconds`
  now costs a sample's genomes (0.17 s each) besides its read pairs, so the shallow soil point gets 5 threads instead of 1
  at 84 slots (its tail 4.4 h → ~1.4 h).
- **Streaming into protal: possible on protal's side, not worth it now.** protal already reads FIFOs (an e2e test
  does it), but the producers would need reworking: R1/R2 interleaved, samples in map order, the se unit reading the
  pe R1 a second time, the host reads appended later, the follower's file bookkeeping. With the long reads no longer
  leaving intermediate files, and `--profile-blocks` removing profiled reads, scratch capacity is no longer the
  constraint it was in v12-v14.
- **ART cannot be ported into protal as it is licensed.** art_illumina's source is GPL-3.0-or-later; protal is
  GPL-2.0-only, and GPLv3 code cannot be copied or translated into a GPLv2-only program. The options are relicensing
  protal (all copyright holders), a clean-room Illumina model with protal's own profiles, or a patched ART kept as a
  separate GPL-3 program (fork and exec is fine).

## 1. The long-read simulator in C++

### What changed

| File | Change |
|---|---|
| `src/RandomForest/LongReadSimulator.{h,cpp}` (new, ~1,350 lines) | `LongRng` (xoshiro256\*\*, splitmix64 seeds); `LongReadSetup::Parse` (the collector's setup texts); `hifi::Mutate` (HiFi and flow models of `hifi_reads.py`); `QshmmModel` (pbsim3's `simulate_by_qshmm_templ`, `set_qshmm`, `set_mut`); template drawing from genomes and a prepared host (`host.seq` by memory map); the engine: samples planned in rounds with their own random stream, work items (some reads of one genome in one round, at most 2 Mb) with their own streams, a thread pool across samples, each sample's items written in order; `MutateTemplates` (one read per FASTA record) |
| `src/simulate_metagenomes_main.cpp` | `--long_samples TSV --long_genomes TSV --long_setup TEXT [--long_model FILE] [--long_stats FILE] -t N`; `--long_templates FASTA --long_out FILE --seed N` |
| `tests/test_LongReadSimulation.cpp` (new, 4 tests) | RNG, setups, read lengths; HiFi calibration (errors made vs. expected by the qualities, Q30 at 25 kb) and flow model (mean base Q 25); qshmm tables (accuracy levels by exp(0.22 x level), the level-100 reads at Q93, read/template length, error rate) on a synthetic model; samples: names `g<i>x_<n>`, contig ends, weights, both strands, a host, zstd and BGZF, **the same files on 1 and 4 threads**, failures leave no file |
| `scripts/collect_training_data.py` | `long_unit_jobs`: writes the unit's samples and genomes TSVs, one job running the binary on `long_threads` threads (~2 min of work each, at most a quarter of the cores); removed `draw_templates`, `read_contigs`, `read_length`, `long_read_sample`, `draw_chunks`, `make_reads`, `long_read_chunks`, `join_chunks`, the draw group and `--long_read_chunk`; `parse_long_setup` keeps the text and refuses `errhmm` (pbsim3's one-pass errhmm reads have every quality at Q0); keys: `LONG_READS`, `LONG_MODELS`, the simulator's identity, the qshmm model file's; `pe_point_seconds` counts genomes (`genomes_per_sample`, `PE_GENOME_SECONDS`, `PE_PAIR_SECONDS`) |
| `scripts/build_gtdb_database.py` | checks that pbsim3's model file is found (`--pbsim-models` or pbsim3's data folder) instead of the pbsim3 executable; help texts |
| `scripts/mini_db/test_collector.py` | the long-read tests run the real binary (`$SIMULATE`, as CI sets it) with an error-free qshmm model, so the exact template checks still hold; a host share; the gzip and zstd files the same reads |
| `scripts/hifi_reads.py`, `scripts/error_reads.py`, `docs/databases.md` | notes: the model's reference, the read names' origin, how long reads are made |

The read names stay `g<i>x_<n>` (i: the genome's place in the sample's manifest list, the host after them; n = 1, 2,
... in the file's order): `error_reads.py` reads the genome from them. The samples are not byte-identical to the
Python ones (another random number generator); every long-read point is simulated again by the next build, since its
key changed.

### Faithfulness of the ported qshmm model

pbsim3 keeps its HMM in arrays sized for 50 states (`ip[101][51]`, `ep[101][51][94]`, `tp[101][51][51]`).
`QSHMM-ONT-HQ.model` has **56 states at accuracy level 92**. pbsim3 writes states 51-56 past the end of a row, into
accuracy 93's first rows of the same flat array, and transitions to them into the next state's first entries; then it
reads only states 1-50. The port uses the same flat layout and writes those records alike, refusing only what falls
outside the arrays. It therefore makes the reads pbsim3 makes, out-of-bounds writes included. The other two models
(`QSHMM-ONT`, 23 states; `QSHMM-RSII`, 37) stay within 50.

A property of pbsim3 worth knowing: at `--accuracy-mean 0.97` a read's accuracy level is drawn from 72 to 100 with
weights exp(0.22 x level). Level 100 has no HMM in the model and takes a uniform quality of Q93, so **about 20% of the
Nanopore training reads are error-free with Q93 throughout** (19.2% from pbsim3, 20.6% from the port, 19.8% expected).
Real Nanopore reads are not like that. This was so in every build since the qshmm setup; changing it (a lower accuracy
mean, or capping the level at 99) would change the Nanopore training data and is a separate decision.

### Validation on the same templates (`scripts/compare_models.py`, `results/compare_models.tsv`)

Templates drawn as the collector draws them from 40 of the synthetic genomes: 2,000 HiFi (15 kb ± 3 kb), 60,000 Ultima
(300 ± 40) and 4,000 Nanopore templates (8 kb ± 6 kb); one read of each by both programs, one core:

| Model | Program | Length / template | Mean base Q | Per-read Q 5/25/50/75/95% | Q93 reads | Errors by Q | Errors made | CPU s per Mb |
|---|---|---|---|---|---|---|---|---|
| HiFi | simulate_metagenomes | 1.0000 | 43.52 | 36.9/41.5/45.1/49.2/54.7 | 0 | 0.00014 | 0.00014 | 0.033 |
| HiFi | hifi_reads.py | 1.0000 | 43.49 | 37.1/41.6/45.1/49.2/54.8 | 0 | 0.00015 | 0.00015 | 0.086 |
| Ultima | simulate_metagenomes | 1.0000 (SD 0.0044) | 25.00 | 21.7/23.7/25.0/26.3/28.3 | 0 | 0.00654 | 0.00653 | 0.057 |
| Ultima | hifi_reads.py | 1.0000 (SD 0.0044) | 25.00 | 21.7/23.7/25.0/26.4/28.3 | 0 | 0.00654 | 0.00652 | 0.109 |
| Nanopore | simulate_metagenomes | 0.9947 (SD 0.0068) | 40.75 | 14.6/23.2/31.0/35.7/93.0 | 0.206 | 0.0427 | 0.0427 | 0.019 |
| Nanopore | pbsim3 3.0.5 | 0.9947 (SD 0.0066) | 39.67 | 14.4/23.0/31.0/35.5/93.0 | 0.192 | 0.0440 | 0.0440 | 0.467 |

Total variation distance of the base-quality histograms: HiFi 0.0046, Ultima 0.0007, Nanopore 0.0157. For Nanopore
that distance is mostly the share of Q93 reads, which differs by 1.4 points, within the sampling noise of 4,000 reads
(SD 0.6 points each around 19.8%). "Errors made" is each program's own count: `hifi_reads.mutate`'s events,
pbsim3's substitution, insertion and deletion rates per read base, the port's count per read base. pbsim3's CPU leaves
out the collector's pass that renamed and recompressed its reads (31 ms per Mb, 2026-10-05).

### Whole samples (`scripts/bench_community.py`, `results/bench_community.tsv`)

One sample of a community of the 400 synthetic genomes (3.4 Mb each, lognormal abundances, sigma 1.5): template
drawing, the read model and zstd compression together, at 1 and 4 threads (`taskset -c 0-3`, a shared laptop):

| Setup | Bases | Threads | Wall | CPU | CPU s per Mb |
|---|---:|---:|---:|---:|---:|
| HiFi 15 kb | 1 Gb | 1 | 50.1 s | 48.5 s | 0.049 |
| HiFi 15 kb | 1 Gb | 4 | 19.9 s | 73.9 s | 0.074 |
| Ultima 300 bp | 0.5 Gb | 1 | 37.2 s | 37.3 s | 0.075 |
| Ultima 300 bp | 0.5 Gb | 4 | 13.1 s | 49.7 s | 0.099 |
| Nanopore 8 kb | 1 Gb | 1 | 42.8 s | 43.9 s | 0.044 |
| Nanopore 8 kb | 1 Gb | 4 | 15.2 s | 54.3 s | 0.054 |

- **Genome reads.** A sample of the 400 genomes with one 1 kb read each on average took 3.0 s of CPU for the ~250
  genomes it read: **~12 ms per genome** (inflated with ISA-L and parsed). The Python path took 21 ms, 7 of them
  holding the interpreter lock.
- **Threads and memory.** 4 threads gave 2.5-2.8x. The CPU seconds rise with threads on this laptop's mix of fast and
  low-power cores, as in earlier reports. At most 107 MB of memory.
- **The collector's estimates** now use these numbers (`LONG_SECONDS_PER_MB`, `LONG_GENOME_SECONDS`) and give a
  long-read point threads for ~2 minutes of work, at most a quarter of the cores (`long_threads`).

Against the path it replaces (per Mb, laptop core; the 2026-10-05 unit costs and the table above):

| Reads | Before | Now | Factor |
|---|---|---:|---:|
| PacBio | templates drawn (27 µs each) and written gzipped, `hifi_reads.py` 86-112 ms, the reads written: ~0.12-0.15 s | 0.049 s | ~3x |
| Ultima | drawn (4 µs each, ~3,300 per Mb), flow model 109-223 ms: ~0.15-0.25 s | 0.075 s | ~2-3x |
| Nanopore | plain templates, pbsim3 467-591 ms, renamed and recompressed 31 ms: ~0.5-0.65 s | 0.044 s | ~11-15x |

At r226 the training collection makes about 134 Gb of PacBio, 134 Gb of Nanopore and 61 Gb of Ultima reads
(design and scenarios, at v14's depth factors), the test set about half that. The estimates:

- **Reads:** ~4.7 laptop core-hours, plus ~1.8 h for the soil samples' genomes (5 units x 6 samples x ~14,000
  genomes x 1.3 rounds x 12 ms).
- **Before:** the 2026-10-05 projection of the one-pass Python path was 12.7 h for the scenarios alone, at a quarter
  of today's soil communities.
- **Wall time** on the node is minutes per point, and long reads no longer wait for the GIL or for pbsim3.
- **Disk:** nothing besides the reads. Before, a chunk had its templates (0.35-1 byte per base) and pbsim3's files
  (1.4) on the disk with its reads.

### Follow-up (same day): no error-free Q93 Nanopore reads

The user asked that no error-free Q93 reads be made at all. `QshmmModel` now draws only the accuracy levels the model
has an HMM for (71-99 in `QSHMM-ONT-HQ`), with pbsim3's weights, exp(0.22 x level), in proportion; a model without an
HMM in the range keeps the uniform qualities below 100. The setup's key changed (`LONG_MODELS`), so the Nanopore
points are simulated again. On the same 4,000 templates (`compare_models.py`):

| Program | Length / template | Mean base Q | Per-read Q 5/25/50/75/95% | Q93 reads | Errors by Q | Errors made |
|---|---|---|---|---|---|---|
| simulate_metagenomes, now | 0.9934 (SD 0.0068) | 27.15 | 13.6/20.4/29.0/33.5/35.8 | 0 | 0.0546 | 0.0545 |
| pbsim3 3.0.5 | 0.9947 (SD 0.0066) | 39.67 | 14.4/23.0/31.0/35.5/93.0 | 0.192 | 0.0440 | 0.0440 |

The reads' error rate rises from 4.4% to 5.5%: the fifth of the reads that had none now has the others' mix of levels.
pbsim3's mean base Q of 39.7 was mostly its Q93 reads.

## 2. Streaming the reads into protal

Read-only investigation of protal's input path and the build (a subagent, with file and line references):

- **protal reads FIFOs already.** `tests/e2e/test_protal_e2e.py` (around line 1947) feeds R1 and R2 FIFOs, each from a
  writer thread. The startup check uses `fs::exists` and `access()` without opening (`Options.h` ~1996-2003,
  ~2290); compression is detected from the first bytes on a pipe (`ThreadedGzStream.h` ~384-397); read-type
  detection and the long-read check skip non-regular files and are not used when the map gives `READ_TYPE`, as the
  collector's maps do; nothing seeks, nor reads an input twice (except `--benchmark_alignment`).
- **The producers are the problem:**
  - **R1 and R2 deadlock.** protal reads R1 and R2 in lockstep, 32 records at a time (`SeqReader.h` ~167-221), with
    ~4 MiB of read-ahead per file. `simulate_metagenomes` appends a genome's whole R1, then its whole R2
    (`MetagenomeSimulator.cpp` ~741), and the host join does that for whole files. Either deadlocks on FIFOs.
    Interleaving in chunks of at most ~1 MiB would fix it.
  - **Sample order.** protal opens samples strictly in map order (`RunProtal.h` ~527-678), but the simulator writes
    several samples at once and opens its writers before ART runs. Writers for later samples would block, so
    simulation would be serialised to protal's pace. The same would hold for the new long-read engine, which keeps
    several samples open.
  - **The se unit reads the pe point's R1** again in the same map (`collect_training_data.py` `unit_map_rows`).
  - **The follower assumes regular files.** It uses `isfile` and `getsize` (`follow`, `reads_removed`,
    `remove_profiled_reads`), and re-simulates missing reads.
  - **Resuming or failing hangs.** A resumed run skips a sample whose SAM exists without opening its FIFO, so its
    writer blocks forever. A failed sample is not retried.
- **What it would save:** the reads on the scratch disk. The SAMs (29 GB in v13 with `--error-reads all`), profiles
  and ART's per-genome temporary files stay. The CPU and the network-FS genome reads stay the same. It would couple
  the simulations and protal into one pipeline, which then runs at the pace of the slower and can no longer be
  resumed.
- **Verdict.** Medium-to-large rework on the producer side, for disk capacity alone. v14 needed 253 of 263 GB. That
  peak shrinks now that long reads leave no templates or pbsim3 files (they took ~1.4-2.4 bytes per base besides
  the reads, [2026-10-05 report](../2026-10-05-collector-profiling/README.md), "Disk per read type").
  The follower already removes reads once profiled, and smaller `--profile-blocks` lowers the peak further. If the
  disk is still short, those two are the cheaper levers. A node's RAM (`/dev/shm`) is another option where the job's
  memory allows it, since tmpfs counts against the job.

## 3. Porting ART's Illumina simulation

From the ART source (a subagent's read and measurements, with file and line references):

- **License: no port into protal as licensed.**
  - Every `art_illumina_src` file says "License: GPL v3" and "either version 3 of the License, or (at your option)
    any later version" (e.g. `art_illumina.cpp:1-21`); the tarball ships `GPLv3.txt`.
  - The built-in profiles (`empdist.h`, 4.2 MB of string tables) say "Copyright(c) 2008-2011 Weichun Huang All
    Rights Reserved"; treat them as part of ART.
  - GPLv3 code cannot be copied or translated into a GPLv2-only program (GNU GPL FAQ, #v2v3Compatibility,
    #AllCompatibility, #TranslateCode). Running art_illumina as a separate program, as now, is fine.
  - pbsim3, by contrast, ships the GPL version 2 text (`COPYING`), which allows the qshmm port above.
- **Ways forward:**
  1. **Relicense protal as GPL-3.0-or-later.** Its vendored libraries (MIT, LGPL-2.1, zlib) allow it, but every
     protal copyright holder must agree. Then port: ~900-1,100 lines of C++, 1-2 weeks with validation.
  2. **Clean-room Illumina model.** Write it from a behavioural specification, not from ART's code, with protal's own
     quality profiles estimated from real FASTQ (or from art_illumina's output). It would not be ART, so the
     paired-end models would be retrained on its reads.
  3. **Keep ART as a separate GPL-3 program, but cheaper per genome.** A patched art_illumina, shipped separately,
    could take a genome once for all samples of a point and the pairs of each. Its fixed cost per genome is mostly
    its own preprocessing: 44% `ini_set` (`toupper` per base and a reverse-complement copy), 22% `mask_n_region`,
    15% reading the FASTA twice (once for a SAM header even with `-na`). Or protal's simulator could keep one plain
    copy of a genome per point instead of one per sample. Both keep ART's reads.
- **What a port would gain.** Per genome, ~3-10 ms in process instead of ~100-130 ms per genome per sample. Per
  pair, ~1-2 µs instead of 22-31 µs. For the soil scenarios that is minutes instead of hours of CPU.
- **Two quirks a port would have to copy to keep ART's reads:**
  - ART's indel count table is off by one (`seqRead.h` ~96-110), so a 150 bp read has an indel far less often
    than its rates say (measured 5e-5 reads with an insertion against ~1.5% nominal).
  - A pair with an N in the read's start window is dropped, not replaced.

  Both are in every paired-end training sample so far.

The fix that matters for the build's wall time needs none of this. The shallow-soil tail came from the thread
split, which now counts genomes (§ summary). ART's CPU (~11.6 core-hours for the soil scenarios at r226) is spread
over the node once the split is right.

## Tests run

In `~/lrsim` (a copy of the working tree, so including another session's uncommitted changes to `Options.h`,
`RunProtal.h` and the e2e tests), 4 cores:

- **Unit tests, Release:** 403 passed (the two debug-only tests skipped, as in every Release build), 5.5 s; the 4 new
  `LongReadSimulation` tests 0.5 s.
- **Unit tests, ASan and UBSan:** Debug with `PROTAL_NO_CLONES` and leak checking, as CI's sanitizer job: 403 passed.
- **`scripts/mini_db/test_collector.py`** with `$SIMULATE`: 14 tests OK.
- **`scripts/mini_db/test_*.py`:** with `$PROTAL`, `$SIMULATE` and the `protal-db-build` Python, 71 tests OK in
  173 s. These include the GTDB build tests, which simulate PacBio, Nanopore and Ultima reads of the design and the
  scenarios through `--long_samples`.
- **CI's other script suites** (`test_insilico_strains`, `test_trace_relatives`, `test_error_reads`,
  `test_profile_scripts`, `test_strain_scripts`, `test_model_pmml`): 67 tests OK. Without `$SIMULATE`,
  `test_collector.py` skips its two long-read tests (or fails them under `PROTAL_TESTS_REQUIRED=1`, as in CI).
- **Not run:** `tests/e2e` (protal itself is unchanged here) and a GTDB-scale build. The next r226 build re-simulates every long-read point; its `*_simulation.log` lines "simulated (N samples,
  R reads) in T on K threads" show the new path, and the paired-end line shows the soil points' threads.

## Reproducing

```
# build simulate_metagenomes, protal and protal_tests (Release, -DPROTAL_BUILD_TESTS=ON), then:
build/tests/protal_tests --gtest_filter='LongReadSimulation*'
python3 scripts/compare_models.py --scripts SCRIPTS --simulator build/simulate_metagenomes --pbsim PBSIM \
    --model .../QSHMM-ONT-HQ.model --genomes ~/colprof/genomes --out OUT          # one core: taskset -c 0
python3 scripts/bench_community.py --simulator build/simulate_metagenomes --model .../QSHMM-ONT-HQ.model \
    --genomes ~/colprof/genomes --out OUT --threads 1,4
SIMULATE=build/simulate_metagenomes python3 -m unittest scripts/mini_db/test_collector.py
```

The website is not affected: nothing here changes how protal itself is run.
