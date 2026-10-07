# Illumina reads without ART, and large samples streamed into protal

2026-10-07, `audit-fixes` at `c601cfc` plus the changes described here, committed as `c4cf9bc` (on `fe64d2d`, another
session's commit in between; the tree tested was `c601cfc` plus the change). This follows the
[long-read simulator report](../2026-10-07-long-read-simulator/README.md), which found two things. ART cannot be ported
into protal: it is GPL-3.0-or-later and protal is GPL-2.0-only. Streaming reads into protal needed the producers
reworked. The requests were:

1. Drop ART altogether. Write a reasonable model of Illumina reads from publications and expected error rates instead.
   An approximation is fine for training protal's models.
2. Stream a sample into protal through a pipe when it is large (for example over 2 GB of `.fq.zst`). Otherwise keep
   writing samples to the disk and profiling several in one run.

Sources:
- Published Illumina characteristics, collected by a literature subagent that did not open ART's code or profiles
  (list below).
- The code at `c601cfc`.
- Measurements in WSL (Core Ultra 7 258V, `taskset -c 0-3`, nice 5) on the 400 synthetic genomes of real size of the
  [2026-10-05 report](../2026-10-05-collector-profiling/README.md) (`~/colprof/genomes`).
- Scripts in [`scripts/`](scripts/). Build trees: `~/lrsim` (development) and `~/verify-illumina` (the committed tree,
  built from its `git archive`).

## Summary

- **ART is gone.** `simulate_metagenomes` makes the paired-end reads itself (`src/RandomForest/IlluminaSimulator`).
  There is no `art_illumina` process, no per-genome temporary FASTQ, and no join. ART is removed from the build
  environment (`envs/protal-db-build.yaml`), the conda recipe's run requirements, CI and the tests' prerequisites.
  `--art_path` and `--extra_art_args` are refused with an explanation. A read setup with ART's `file=R1.txt+R2.txt`
  quality profiles is refused when the design is made.
- **The model** is built from published numbers. Its main parts:
  - per-cycle mean qualities, R1 and R2, from InSilicoSeq's HiSeq, NovaSeq and MiSeq models;
  - run, cluster and read offsets and AR(1) noise along the read;
  - a low-quality state, and Q2 tails;
  - errors drawn from the hidden quality, more after GG, with a four-colour or two-colour substitution spectrum;
  - indels at their own rates, more in homopolymers;
  - Ns;
  - each instrument's quality bins.

  The written qualities are calibrated to the published curves (within about 0.5 Q over most of the read). The error
  rates are close to the published ones: HS20 a little above Schirmer 2016, NovaSeq's R1 a little above the measured
  share at Q30 or more, the rest within the published ranges (tables below). Five instruments are available:
  `HS20`, `HS25`, `HSXt`, `NovaSeq` and `MSv3`. The build's default read setups use `HS20` (100 bp), `HSXt` (150 bp)
  and `MSv3` (250 bp), as before.
- **Cheaper than ART.** A pair costs 18 µs of CPU and a genome 17.5 ms per sample. ART cost 31-96 µs per pair and
  85-166 ms per genome. Work items run on all threads (2.4x faster on 4 threads), and the files are byte-identical for
  any thread count.
- **Large samples are streamed.** A design point whose largest sample would take more than `--stream-above` GB of
  compressed reads (default 2) is never written to the disk. A protal run of its own reads the point's samples from
  named pipes while `simulate_metagenomes` makes them, one sample at a time, R1 and R2 in step. Smaller points are
  written and profiled in blocks as before. The mini GTDB build gives the same training tables when every point is
  streamed.
- **What this changes for builds.** Training reads now come from another simulator, so they differ statistically
  from the ART reads that existing models (r226 up to v15) were trained on. Retrain before comparing models across
  this change. Manifests from the ART days replay their communities with new reads.

## 1. The model

`IlluminaModel` (`src/RandomForest/IlluminaSimulator.{h,cpp}`) makes one pair per fragment.

- **Fragment.** Its length is normal (`--fragment_mean`, `--fragment_stdev`), between one read length and the longest
  contig; out-of-range draws are redrawn. Its position is uniform over the places in the genome's contigs where it
  fits (a position that runs off its contig is drawn again), on either strand. Read 1 is the fragment's first bases;
  read 2 is the start of its reverse complement.
- **Hidden quality of a base at cycle c of read R.** It is `μ_R(c) + run + cluster + read + noise_c`:
  - `μ_R(c)` is the instrument's published mean curve, calibrated so that the qualities *written* (after the low
    state, the tails and binning) average the curve. Six rounds of 2,000 reads of a random template run with a fixed
    seed when the model is built. Each round moves `μ` by the difference at each cycle, capped at 4 Q above the
    curve.
  - `run` ~ N(0, 1) per sample and `cluster` ~ N(0, 1.5) per pair (shared by R1 and R2). These are estimates: runs
    differ widely (Stoler & Nekrutenko 2021).
  - `read` ~ N(0, 1).
  - `noise` is AR(1) with ρ = 0.8. Its SD grows from 2 Q at the read's start to the profile's end value
    (MiSeq: 7 Q) as (c/L)². Neighbouring qualities are strongly correlated (Malysa 2015).
  - `--mean_quality Q` shifts the curves so that each read averages Q. Errors follow the shift.
- **Low state.** A read enters it at cycle c with probability `low_rate × exp(2 (c/L − 1))` (read 2 has its own rate),
  leaves it with 0.25 per cycle, and takes qualities N(15, 5) while in it. This gives the low-quality bases of binned
  instruments: the NovaSeq Q12 and Q23 bins hold 3.4% and 2.6% of bases (Illumina 770-2017-010). The rates were fitted
  to those shares and to the per-cycle curves (`--illumina_report`).
- **Q2 tails.** A share of the reads (`collapse`) fall to Q2 from a cycle in their last 30% to the end, with 5% errors
  there. Minoche 2011 saw Q2 ("B") tails in 26-40% of HiSeq 2000 reads. We use 10% for HS20, because more would push
  its substitution rate above Schirmer 2016's 0.26%/0.40%. The other instruments use 2-5%.
- **Errors.** A base is miscalled with probability `calibration × 10^(−q/10)` of its hidden quality. The calibration
  is 1 except on MiSeq, which is 1.5 because its high qualities overstate accuracy (Schirmer 2015). The probability is
  multiplied by `gg_factor` after GG in read direction:
  - ×6 on HiSeq and MiSeq: GGG is 1/64 of sites but holds about 10% of substitutions (Schirmer 2016; Nakamura 2011);
  - ×1.5 on HiSeq X and NovaSeq (Stoler & Nekrutenko 2021 found little context bias).

  Substitutions follow the instrument's spectrum:
  - four-colour (HiSeq, MiSeq): A>C is 66% of A's errors and C>A 58% of C's, with the complements likewise (Schirmer
    2015);
  - two-colour (NovaSeq): mostly transitions, with a G (no signal) bias.
- **Indels.** They are independent of quality: HiSeq has insertions 2.8e-6 and deletions 5.1e-6 per base (Schirmer
  2016), and MiSeq 4e-5 and 2e-5 (Schirmer 2015). They are ×10 in homopolymers of 5 or more. A read is always its full
  length. Past a short fragment's end the bases are random, with qualities at most Q10.
- **Ns.** An N at Q2 occurs at 1e-4 per base and 1e-3 at the first cycle, and wherever the genome has one.
- **Read 2 and long inserts.** Read 2's qualities drop by 0.005 Q per base of fragment above 500 bp (the shape is from
  Tan 2019).
- **Binning.** The written quality is the hidden one, clamped to the instrument's range and binned:
  - HS25: Illumina's 8 levels 6/15/22/27/33/37/40 (770-2012-058, 970-2012-013);
  - HSXt: 2/12/22/27/32/37/41 (as seen in HiSeq X and 4000 files; no Illumina document);
  - NovaSeq: 2/12/23/37 (RTA3, 770-2017-010);
  - HS20 and MSv3: unbinned (2-41).
- **Names.** Reads are named `<contig>-<n>/1` and `/2`, as ART named them, so `trace_relatives.py` and
  `error_reads.py` still find a read's contig. Host reads are named `h_<n>/1` and `/2`. `n` counts over the sample.

Not modelled:
- poly-G after short inserts on two-colour instruments (protal trims nothing, and training fragments are 300-550 bp);
- adapters;
- index hopping;
- GC coverage bias;
- per-tile effects.

### The instruments

| | HS20 | HS25 | HSXt | NovaSeq | MSv3 |
|---|---|---|---|---|---|
| instrument | HiSeq 2000, 2x100 | HiSeq 2500, 2x125 | HiSeq X Ten / 4000, 2x150 | NovaSeq 6000, 2x150 | MiSeq v3, 2x250 / 2x300 |
| mean curve (R1; R2) | InSilicoSeq HiSeq 2x126: 32.5 (cycles 1-5), 37 (8), 36 (70), 35 (110), 30 (126); R2 0.5 lower | as HS20 | InSilicoSeq NovaSeq 2x151: 36 (1), 37 (60), 36.5 (110), 35.5 (150); R2 36.5, 36, 34.5, 31 | as HSXt | InSilicoSeq MiSeq 2x301: 33.5 (1-9), 37.5 (10-100), 35.5 (160), 32 (220), 28 (264), 23 (300); R2 33, 37 (60), 33 (160), 27, 21, 17 |
| bins | none | 8-level | 7 values | 4 values | none |
| low state R1 / R2 | 0.020 / 0.032 | 0.020 / 0.032 | 0.022 / 0.045 | 0.016 / 0.040 | 0.012 / 0.022 |
| Q2 tails | 10% | 5% | 2% | 2% | 4% |
| calibration, GG | 1, ×6 | 1, ×6 | 1, ×1.5 | 1, ×1.5 | 1.5, ×6 |
| noise SD start → end | 2 → 3.5 | 2 → 3.5 | 2 → 3 | 2 → 3 | 2 → 7 |
| insertions, deletions | 2.8e-6, 5.1e-6 | 2.8e-6, 5.1e-6 | 5e-6, 5e-6 (none published; est.) | 5e-6, 5e-6 | 4e-5, 2e-5 |
| spectrum | four-colour | four-colour | four-colour | two-colour | four-colour |

## 2. What the profiles give

`simulate_metagenomes --illumina_report 20000 --sequencer X --read_length L --fragment_mean F --fragment_stdev 50`,
one core ([`scripts/profile_reports.sh`](scripts/profile_reports.sh); 0.5-1.2 s each):

| instrument, read setup | read | mean Q | ≥Q30 | substitutions | insertions | deletions | Ns | Q2 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| HS20, 2x100, fragments 300 | R1 | 35.9 | 94.3% | 0.34% | 3.5e-6 | 6.0e-6 | 1.1e-4 | 1.6% |
| | R2 | 35.3 | 92.4% | 0.46% | 3.5e-6 | 6.0e-6 | 1.1e-4 | 1.6% |
| HS25, 2x125, 300 | R1 | 35.5 | 93.8% | 0.30% | 4.0e-6 | 6.4e-6 | 1.0e-4 | 0.8% |
| | R2 | 35.1 | 92.1% | 0.44% | 4.4e-6 | 5.6e-6 | 1.1e-4 | 0.8% |
| HSXt, 2x150, 350 | R1 | 36.5 | 95.9% | 0.26% | 3.3e-6 | 7.3e-6 | 1.1e-4 | 0.3% |
| | R2 | 35.0 | 91.7% | 0.47% | 5.7e-6 | 5.0e-6 | 1.2e-4 | 0.4% |
| NovaSeq, 2x150, 350 | R1 | 36.2 | 96.4% | 0.20% | 6.7e-6 | 8.7e-6 | 1.2e-4 | 0.7% |
| | R2 | 34.9 | 90.3% | 0.43% | 4.3e-6 | 4.3e-6 | 1.2e-4 | 1.2% |
| MSv3, 2x250, 550 | R1 | 35.4 | 89.4% | 0.34% | 4.9e-5 | 2.4e-5 | 1.1e-4 | 0.6% |
| | R2 | 32.3 | 76.2% | 0.72% | 4.7e-5 | 2.4e-5 | 1.0e-4 | 0.7% |
| MSv3, 2x300, 550 | R1 | 33.9 | 81.0% | 0.52% | 4.6e-5 | 2.5e-5 | 1.1e-4 | 0.6% |
| | R2 | 30.1 | 64.3% | 1.42% | 4.5e-5 | 2.2e-5 | 1.0e-4 | 0.7% |

NovaSeq's written bins (R1 / R2): Q37 96.4% / 90.3%, Q23 1.2% / 4.4%, Q12 1.8% / 4.2%, Q2 0.7% / 1.2%.

Per-cycle means against the curve ("pub" is the curve at that cycle):

| instrument | R1: cycle 1 / 10 / middle / 80% / last | R2: the same |
|---|---|---|
| HS20 2x100 | 32.5 (32.5) / 37.0 (37.0) / 36.5 (36.3) / 35.8 (35.8) / **33.3 (35.2)** | 32.0 (32.0) / 36.6 (36.5) / 36.0 (35.8) / 35.6 (35.2) / **32.2 (34.8)** |
| HS25 2x125 | 32.3 (32.5) / 37.0 (37.0) / 36.1 (36.1) / 35.2 (35.2) / 30.4 (30.3) | 32.1 (32.0) / 36.6 (36.5) / 35.7 (35.6) / 34.8 (34.8) / 30.0 (29.8) |
| HSXt 2x150 | 36.0 (36.0) / 36.2 (36.2) / 36.9 (36.9) / 36.3 (36.2) / 35.8 (35.5) | 36.5 (36.5) / 36.5 (36.4) / 35.6 (35.5) / 33.7 (33.6) / 31.0 (31.1) |
| NovaSeq 2x150 | 36.1 (36.0) / 36.3 (36.2) / 36.5 (36.9) / 35.8 (36.2) / 34.9 (35.5) | 36.6 (36.5) / 36.3 (36.4) / 35.6 (35.5) / 33.6 (33.6) / 31.0 (31.1) |
| MSv3 2x250 | 33.5 (33.5) / 37.6 (37.5) / 36.7 (36.7) / 33.2 (33.2) / 29.3 (29.3) | 32.7 (33.0) / 33.3 (33.6) / 34.1 (34.4) / 28.6 (29.0) / 22.5 (22.9) |
| MSv3 2x300 | 33.6 (33.5) / 37.5 (37.5) / 35.8 (35.8) / 30.1 (30.2) / 22.8 (23.0) | 32.7 (33.0) / 33.2 (33.6) / 33.1 (33.4) / 24.0 (24.3) / 16.8 (17.0) |

The curves are met within 0.5 Q, except at two ends:
- **HS20's last cycles** fall about 2 Q below the curve. The curve is HiSeq 2500's, and HS20 adds Q2 tails to 10% of
  its reads, as HiSeq 2000 had them. Calibration may raise the hidden mean by at most 4 Q, which is not enough there.
- **NovaSeq's R1** ends about 0.6 Q low. With four bins, the written mean moves in steps, and the cap stops the
  calibration.

Against the published numbers:

| instrument | published | here |
|---|---|---|
| HiSeq 2000/2500 | ≥Q30: ≥80% at 2x100 and 2x125 (spec, 770-2012-041). Substitutions R1 0.26%, R2 0.40% (Schirmer 2016); 0.11-0.16% once Q2 tails are removed (Minoche 2011). Insertions 2.8e-6, deletions 5.1e-6 (Schirmer 2016) | 92-94%; HS25 0.30% / 0.44%, HS20 0.34% / 0.46% (its Q2 tails); 3.5-4.4e-6 and 5.6-6.4e-6 (homopolymers ×10) |
| HiSeq X / 4000 | ≥Q30 ≥75% (spec). PhiX R1 0.41%, R2 1.41% (770-2017-010); median over datasets 0.087% (Stoler & Nekrutenko 2021) | 95.9% / 91.7%; 0.26% / 0.47%, between the two |
| NovaSeq 6000 | ≥Q30 measured R1 92.2-94.2%, R2 89.6-91.3% (770-2020-013, M-US-00201). PhiX R1 0.14-0.35%, R2 0.20-0.61%. Bins Q37 94.0%, Q23 2.6%, Q12 3.4% (770-2017-010) | R1 96.4% (a little high), R2 90.3%; 0.20% / 0.43%; Q37 93.3%, Q23 2.8%, Q12 3.0% (and Q2 0.9%), averaged over R1 and R2 |
| MiSeq v3 | ≥Q30 >70% at 2x300 (spec). R1 0.41%, R2 0.99% (Schirmer 2016); amplicons 0.64% / 1.07% (Schirmer 2015). Insertions 4.0-4.3e-5, deletions 1.7-2.7e-5 (Schirmer 2015) | 72.6% over both reads at 2x300; 0.34% / 0.72% at 2x250, 0.52% / 1.42% at 2x300; 4.5-4.9e-5 and 2.2-2.5e-5 |

## 3. Cost

[`scripts/bench_pairs.sh`](scripts/bench_pairs.sh) on the 400 synthetic genomes (HSXt, 150 bp, zstd, seed 3):

| run | pairs | genomes | threads | wall | CPU | max RSS |
|---|---:|---:|---:|---:|---:|---:|
| deep | 2,000,000 | 400 | 1 | 57.4 s | 43.5 s | 22 MB |
| deep | 2,000,000 | 400 | 4 | 23.7 s | 45.5 s | 65 MB |
| shallow | 4,000 | 400 | 1 | 9.5 s | 7.0 s | 22 MB |

That is 17.5 ms of CPU per genome and sample (reading and packing it) and 18 µs per pair. ART cost 85-166 ms per genome
and 31-96 µs per pair ([idle-tail report](../2026-10-07-build-idle-tail/README.md): 31 µs and ~100 ms on a laptop core,
53 µs and 166 ms fitted on the cluster node). A shallow soil sample of ~14,000 genomes therefore costs ~4 CPU minutes
instead of ~20-40.

The collector's cost estimate follows (`PE_GENOME_SECONDS, PE_PAIR_SECONDS = 0.018, 18e-6`). The host's reads are now
made by the same run, after the community's (`--host_folder`, `--host_pairs`), instead of by a separate ART run per
sample and a join (`scenarios.host_pe_chunk` and the collector's `host_pe_jobs` and worker pool are gone).

## 4. Large samples streamed into protal

The long-read report listed what blocked streaming on the producers' side. Each was addressed as follows.

| Problem | Now |
|---|---|
| R1 and R2 deadlock: protal reads them in lockstep, the simulator wrote a genome's whole R1 then R2 | A FIFO output switches the pipeline (`ReadPipeline`) to streaming: one sample at a time, R1 and R2 written alternately in pieces of ~600 kB of FASTQ (a work item's pairs) |
| Samples out of map order | Streaming outputs open one sample at a time, in the map's order |
| The host's reads were joined in afterwards | They come from the same run, after the community's |
| The se unit read the pe point's R1 a second time | It gets its own run with `--first_reads_only`. Read 1s do not depend on whether read 2 is written, so they are the same reads |
| The follower assumed regular files | A streamed point gets a protal run of its own (`stream_run`), outside the block bookkeeping. The design is written first (`--test`: manifest, map, truth), which the point's long-read units need |
| Resuming or failing hangs | protal never skips a sample whose reads come through a pipe (`ReadsArePiped` in `RunProtal.h`), whatever SAM it has. `run_protal` watches the simulators beside protal and stops both when either fails |

Which points are streamed: those whose largest sample would exceed `--stream_above` GB, estimated at `PE_BYTES` per
paired-end base or `DRAWN_BYTES` per long-read base. The `--simulate_only` run and the follower must get the same value;
`build_gtdb_database.py --stream-above` (default 2) passes it to both. A streamed point leaves `streamed.json` and its
SAMs, profiles and dumps, but no reads and no pipes. Long-read units stream the same way: `--long_samples` writes into
FIFOs at its sample outputs.

What it saves: the disk for the largest samples (r226's 5M-pair soil and deep points, ~1.5-3 GB each). What it costs:
the point's simulation and protal run in step, at the pace of the slower, and a failed run must start the point again.
Smaller points keep the old path, which lets several samples share a protal run.

## 5. Tests

- `tests/test_IlluminaSimulation.cpp`:
  - each profile's written qualities match its curve, and error rates and bins are as specified;
  - pairs are their fragments (R1 forward, R2 reverse complement, errors as reported);
  - samples get their names, host reads and thread-independent files;
  - two FIFOs read in lockstep by a reader with a 60 s timeout.
- `test_collector.py`:
  - `test_streamed_simulations`: which points are streamed, and the follower's streamed run;
  - `test_host_genome`: host reads by the simulator;
  - instruments and `file=` setups refused.
- `test_gtdb_pipeline.py`:
  - `test_h_streamed`: a full mini GTDB build with every design point streamed (`--stream-above 1e-9`) gives the same
    training and test tables as the build that writes its reads, with no reads or pipes left;
  - `test_f_scenarios`: the host scenario's samples hold exactly the host's pairs, named `h_`.
- `tests/e2e/test_protal_e2e.py`: the simulator tests no longer need `art_illumina`.

**Verification.** The committed tree contains only this change: the staged index (tree `344bcd5`) was exported with
`git archive` and built in `~/verify-illumina`.
- ctest: 408 of 408 passed (2 debug-only tests skipped).
- The other script suites (`test_insilico_strains`, `test_trace_relatives`, `test_error_reads`): 16 passed.
- `scripts/mini_db/test_*.py` (72 tests, `$PROTAL` and `$SIMULATE` set, no ART on the path): 70 passed and 2 failed.
  - `test_h_streamed`: the streamed build passed (exit 0, the same training and test tables as the build that writes
    its reads, no reads or pipes left). The test then failed on its own glob: `stream_*` also matched the simulators'
    logs (`stream_pe.log`, `stream_se.log`), which are kept on purpose. The assertion now looks at folders only
    (`stream_*/`) and checks that the logs are there. The corrected glob was checked on a small tree. The test's
    later assertions (`streamed.json`, the SAMs, the log lines) were not reached and have not been rerun: the
    machine's cores were needed elsewhere.
  - `test_a_build_rerun_and_reduced_database`: not this change. Since `899d82f`, `--profile_only` refuses to
    overwrite an earlier run's profiles. `check_model_parity.py` profiles into the fixed folder `TRAINING/parity`, so
    the build's rerun fails its parity check (protal exits with code 30, "--profile_only would overwrite the results
    of an earlier run"). The fix is the parity script's own `--force`.
- The full e2e suite on the mini database was stopped before it ran, for the same reason. An earlier working-tree run
  passed the e2e simulator tests (3). The `RunProtal.h` change (piped samples never skipped) is covered here only by
  `test_h_streamed`'s streamed build.
- On the shared pipeline (`ReadPipeline`), the long reads are byte-identical to those of `acd4163` (md5 of HiFi,
  Ultima and qshmm output at 1 and 4 threads).
- Follow-up (the same day): `check_model_parity.py` passes `--force` to protal. With `--profile_only` that aligns
  nothing and only lets protal write the profiles again. On `d5ffc20` plus that change, built in `~/verify-illumina`
  on 4 cores, `test_a_build_rerun_and_reduced_database` and `test_h_streamed` both pass, `test_h_streamed` with all
  of its assertions.

## 6. Open

- The model has not been checked against real FASTQ files of each instrument. Fitting ρ, the offset SDs and the
  low-state rates to one public run per instrument would replace the estimates (per-cycle means, spread of per-read
  mean residuals, lag-1 autocorrelation).
- Models trained on these reads have not yet been compared with the ART-trained r226 v15 models. The next r226 build
  shows whether F1 moves, and by how much.

## Sources

1. Illumina 2014, *Understanding Illumina Quality Scores*, 770-2012-058.
2. Illumina 2014, *Reducing Whole-Genome Data Storage Footprint*, 970-2012-013.
3. Illumina 2017, *NovaSeq 6000 Quality Scores and RTA3*, 770-2017-010.
4. Illumina 2020, NovaSeq v1.5 reagents tech note, 770-2020-013.
5. Illumina specification sheets: HiSeq 2500 (770-2012-041), HiSeq 3000/4000 (770-2014-057), MiSeq, NovaSeq 6000.
6. Minoche, Dohm & Himmelbauer 2011, Genome Biol 12:R112, doi:10.1186/gb-2011-12-11-r112.
7. Schirmer et al. 2015, NAR 43:e37, doi:10.1093/nar/gku1341.
8. Schirmer et al. 2016, BMC Bioinformatics 17:125.
9. Stoler & Nekrutenko 2021, NAR Genom Bioinform, doi:10.1093/nargab/lqab019.
10. Ma et al. 2019, Genome Biol 20:50.
11. Nakamura et al. 2011, NAR 39:e90, doi:10.1093/nar/gkr344.
12. Tan et al. 2019, Sci Rep, doi:10.1038/s41598-019-39076-7.
13. Gourlé et al. 2019 (InSilicoSeq, MIT), Bioinformatics, doi:10.1093/bioinformatics/bty630. The per-cycle curves
    are read off its documentation's quality figure (±0.5 Q); no data file of it is used.
14. Malysa et al. 2015 (QVZ), Bioinformatics, doi:10.1093/bioinformatics/btv330.
15. Chen et al. 2018 (fastp), Bioinformatics, doi:10.1093/bioinformatics/bty560 (poly-G, not modelled).

Facts and rates from papers and Illumina's notes are used with citation. No code or data file of ART (GPL-3), or of
any other simulator, was used.
