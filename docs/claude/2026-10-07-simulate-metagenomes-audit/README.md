# simulate_metagenomes: logic errors, speed, and how the genomes are read

2026-10-07, `audit-fixes` at `c9b46f9`. The simulator's code last changed in `fb6c42d` (scenarios), `c4cf9bc` (the
Illumina model without ART, streaming) and `1b76405`/`acd4163` (long reads).

The user asked for an audit of `simulate_metagenomes` for (A) logic errors and (B) speed. During the audit they added
a question: how are the reference genomes read, and should the reading be merged across several simulations so that
the genomes are not read again and again?

**Sources.**
- Every line of the simulator:
  - [`simulate_metagenomes_main.cpp`](../../../src/simulate_metagenomes_main.cpp)
  - `src/RandomForest/`: [`MetagenomeSimulator`](../../../src/RandomForest/MetagenomeSimulator.cpp),
    [`CommunityProfileDesigner`](../../../src/RandomForest/CommunityProfileDesigner.cpp),
    [`IlluminaSimulator`](../../../src/RandomForest/IlluminaSimulator.cpp),
    [`LongReadSimulator`](../../../src/RandomForest/LongReadSimulator.cpp) and
    [`ReadPipeline`](../../../src/RandomForest/ReadPipeline.cpp)
  - the FASTA reader [`IO/ThreadedGzStream.h`](../../../src/IO/ThreadedGzStream.h)
- The collector's calls to the simulator (`scripts/collect_training_data.py`, `build_gtdb_database.py`,
  `scenarios.py`), mapped by a read-only subagent and spot-checked.
- Unit costs measured in earlier reports:
  - [Illumina model](../2026-10-07-illumina-model/README.md): 18 µs per pair, 17.5 ms per genome and sample.
  - [Long-read simulator](../2026-10-07-long-read-simulator/README.md): ~12 ms per genome.
  - [Collector profiling](../2026-10-05-collector-profiling/README.md): the synthetic genomes, median 3.4 Mb, ~78
    contigs, gzip level 6, ~1.1 MB each.

**Nothing was built or run for the audit (sections 1-6)**, because the shared machine's cores were to be spared.
Numbers marked "estimate" there come from counting operations against those measurements. Section 6 lists how to
check them.

**Follow-up, the same day ([section 7](#7-follow-up-implemented-and-measured)).** Implemented, tested and benchmarked
before/after:
- the exact fixes: a bulk FASTA reader, the header guard (L6), loading ahead instead of waiting, and the design's
  groups made once;
- the 2-bit genome store (`--genome_store`);
- plain FASTQ into named pipes (`--plain_pipes`);
- the three per-pair changes, with the homopolymer fix (L4).

The exact part writes the same bytes as `c9b46f9`.

**Follow-up 2 ([section 8](#8-follow-up-2-the-logic-findings-fixed)).** L1, L3, L5, L8, L9, L11 and L12 fixed:
- replays reproduce reads (seeds in the manifests);
- long-read templates placed where they fit;
- a continuous lognormal tail with a floor;
- no `_R2` named for first reads only;
- designs independent of the standard library;
- a failed stream cut off detectably;
- a note on genomes with contigs too short for a read.

## Summary

- **Logic errors.** Twelve findings in all. Five reach the build's training data:
  - long-read templates cut at contig ends, with the shortfall drawn again by weight (L3);
  - the genome table's header guessed from substrings (L6);
  - homopolymer detection for Illumina indels (L4);
  - the Poisson-lognormal abundances, which put ~40-45% of a sample's species at one floor abundance (L5);
  - read-count rounding (L10).

  The others concern options the build does not use:
  - `--from_manifest` does not reproduce the reads, although it says it does (L1);
  - strain-sharing floors use the first depth (L2);
  - the CLI default `HS25` extrapolates its quality curve at 150 bp (L7);
  - three smaller ones (L8, L9, L11).

  The pipeline itself checks out: ordering, determinism over threads, R1/R2 in step on pipes, and error cleanup
  (section 4).
- **Where the time goes.** A pair costs 18 µs. A genome costs 12-17.5 ms *each time a sample reads it*, as much as
  ~700-1,000 pairs.
  - A soil pe sample (~8,600 genomes, 20M pairs) spends ~150 s reading genomes against ~360 s making pairs.
  - A shallow-soil sample (5M pairs) spends ~150 s against ~90 s.
  - A shallow design sample spends nearly all its time reading genomes.
- **How often a genome is read.** Each sample of each simulation inflates and parses every genome it holds. Long
  reads do it again in every planning round. Nothing is shared between samples, read types, processes, or the
  training and test collectors.
  - A default r226 build therefore makes about **1.5M genome loads of ~54k distinct genomes: ~28 reads per genome**
    (estimate).
  - That is 5-7 core-hours and ~1.7 TB of gzip read from the network file system (page cache permitting).
- **The parse is most of a load.** Inflating takes ~3.4 ms of the 17.5 ms. The rest is a character-by-character copy
  with `isspace`/`toupper` calls (~10-12 ms) and a new thread with 4 MiB of zeroed buffers for every file.
- **Recommended for genome reading, in order:**
  1. a bulk parse: ~4× cheaper loads, the same reads;
  2. a decoded 2-bit genome store on `--scratch`, memory-mapped as `host.seq` already is. It removes ~95% of what is
     left, for every process and both collectors, without touching the scheduling. About 50 GB for ~59k genomes;
  3. an in-run genome cache across samples and rounds, if scratch is short;
  4. long-read placement that fits the template, so a long-read sample needs one round.

  **Merging the simulations themselves is not recommended.** pe, Ultima, PacBio and Nanopore of one community in one
  run would cut the loads most (scenarios 1.3M → ~0.35M). But nearly every scenario simulation streams into its own
  protal run, and several concurrent r226-sized protal processes cannot be afforded. The store gives the same saving.
- **Per pair**, an estimated 25-35% less work is possible:
  - a faster normal generator (one per base for the AR(1) noise);
  - skipping rare events geometrically instead of three uniform draws per base;
  - read 2 drawn on a stream of its own, so that `--first_reads_only` (the streamed se runs) skips it.

  These change the reads. That costs little now, because the next r226 build is the first on the in-process reads
  anyway. Writing **plain FASTQ into named pipes** saves the compression on both sides of every streamed simulation
  and leaves the reads unchanged.

## 1. Logic errors

| # | Finding | Where | Reaches the build | Severity |
|---|---|---|---|---|
| L1 | `--from_manifest` does not reproduce reads although it says so | `MetagenomeSimulator.cpp:579-580`, `IlluminaSimulator.cpp:403-404` | no (the collector never replays) | medium (claim is wrong) |
| L2 | strain-sharing `MIN_VCOV` floors computed at the first depth, per species not per strain | `MetagenomeSimulator.cpp:614-617`, `CommunityProfileDesigner.cpp:753-769` | no | medium for that option |
| L3 | long-read templates cut at contig ends, the shortfall drawn again by weight | `LongReadSimulator.cpp:497-541, 557-561` | yes (pb, ont, Ultima) | medium-low |
| L4 | homopolymer length for indels mis-measured | `IlluminaSimulator.cpp:278-280` | yes (pe, se) | low |
| L5 | Poisson-lognormal weights at `pln_mu` 0: ~40-45% of species at the floor | `CommunityProfileDesigner.cpp:283-294` | yes (every design) | a modelling decision |
| L6 | genome-table header guessed from substrings: a data row can be taken for it | `MetagenomeSimulator.cpp:137-174` | latent | medium if hit |
| L7 | default `HS25` at the default 150 bp extrapolates its 126-cycle curve to Q22 | `IlluminaSimulator.cpp:39-40, 135-143`, `simulate_metagenomes_main.cpp:269, 274` | no (the build uses HS20/HSXt/MSv3) | low |
| L8 | `--first_reads_only`: manifest and `protal.meta` name an R2 that is never written | `MetagenomeSimulator.cpp:562-566`, `simulate_metagenomes_main.cpp:118-123` | no (the collector writes its own map) | low |
| L9 | designs depend on libstdc++ (std distributions, hash-map order), reads do not | `CommunityProfileDesigner.cpp:306-333, 541-544` | yes | low |
| L10 | read-count rounding nudges random species by one | `CommunityProfileDesigner.cpp:21-43, 604-611` | yes | negligible at build depths |
| L11 | a streamed sample cut by a failure ends at a frame/block boundary and reads as complete | `ReadPipeline.cpp:313-323`, `ThreadedGzStream.h:82-83` | streamed points | low (the collector watches both) |
| L12 | contigs shorter than a read are left out of the reads but counted in `genome_length` | `ReadPipeline.cpp:123`, `MetagenomeSimulator.cpp:539-541` | yes | negligible |

### L1. A replay reproduces the community, not the reads

The help text (`simulate_metagenomes_main.cpp:234-238`) and the log (`:649-651`, "reads are reproduced exactly") say a
manifest with `art_seed` reproduces the reads. Under ART that was true: the seed was a genome's only randomness. Since
`c4cf9bc` two more seeds feed every read, and the manifest records neither:

- **The run offset.** Each sample's quality offset comes from `run_seed = MixSeed(seed_, "run" + i)`
  (`MetagenomeSimulator.cpp:580`, `IlluminaSimulator.cpp:403-404`). It moves every hidden quality, and so every error
  draw. An error draws one more uniform (`IlluminaSimulator.cpp:312-313`), so after the first changed outcome the
  item's random stream is out of step. Every later fragment, position and error differs.
- **The host reads.** They come from `host_seed = MixSeed(seed_, "host" + i)` (`:579`).

Both depend on the run's `--seed`, which a replay draws at random unless it is given, and on the sample's index `i`.
A per-sample manifest (`manifests/<sample>.tsv`) always replays as index 0. So a replay reproduces the reads only with
the original `--seed`, and only for the combined manifest. The comment at `simulate_metagenomes_main.cpp:689` ("the
seed is only consumed for rows without a recorded art_seed") is wrong too: `prepare_sample` draws from `rng_` for
every row, and the seed also makes the run and host seeds.

Fix: write `run_seed` and `host_seed` per sample in the manifest (two columns, appended on the right as the header
comment asks). Failing that, take `--seed` from `run_params.tsv` when it is not given, and write the sample's index.

### L2. Strain sharing's coverage floors

`assign_strains_across_samples` gets `profile_options.total_read_pairs`, the first value of `--total_read_pairs`
(`MetagenomeSimulator.cpp:614-617`). Sample `i` is then made at `ReadPairsForSample(i)` (`:629`). With several
depths, a sample at half the first depth gets half its `MIN_VCOV` floor, and a deeper one more than it needs.

Three more points:
- The floor is per species: the sum over its strains (`CommunityProfileDesigner.cpp:757-768`). The species' pairs
  are then split among its strains at random (`:621-635`). `MetagenomeTypes.h:44` documents a floor per strain.
- Species with a floor that already lay above it are rescaled with the rest (`:590-599`), so they can end below it.
- No collector run uses `--strain_sharing_file`.

### L3. Long-read templates cut at contig ends, the shortfall drawn again by weight

How it works now:
- A template's start is uniform over positions with at least 100 bases left on the contig, and the template ends
  where the contig does (`LongReadSimulator.cpp:557-561`). The Illumina fragments are instead placed only where they
  fit (`IlluminaSimulator.cpp:476-481`).
- `Plan` draws reads by weight until their planned lengths cover the bases still missing (`:497-515`). After a round,
  the bases lost to cuts are drawn again, by weight over **all** genomes.

The loss per read is about ℓ/(2·L), with L the base-weighted contig length. For HiFi at 15 kb:
- a MAG of 30 kb contigs loses ~25% of its planned bases;
- a 100 kb N50 loses ~7.5%;
- a closed genome loses ~0.2%.

The redrawn bases go to every genome, so fragmented genomes end ~10-25% below their weight's share and closed ones a
few percent above it. The data then disagree with the design's abundances by assembly quality. For presence training
this is a small depth shift; for abundance it is a bias.

Each round also re-reads every genome it touches. The pipeline's genome slot is keyed by (round, genome) and freed when
the round's parts are done (`ReadPipeline.cpp:296, 372-377, 463-466`). A round must finish before the next is planned.
- Rounds go on until the overshoot of the last read covers the cut losses. Their number grows with log(bases).
- Round 2 is a few percent of the bases, but spread over every genome by weight it touches most of them again.
- The collector assumes 1.3 loads per genome (`collect_training_data.py:1043`). By this arithmetic deep soil samples
  need ~1.5-2 (estimate). `--long_stats` writes `rounds` per sample: the next build's logs will tell.

The port reproduces the Python `draw_templates` faithfully, so this is inherited, not new. The fix fixes both problems:
place a template only where it fits (redraw the start, as the Illumina path does; cut only if no contig is long enough)
and stop at bases ≥ (1 − 10⁻³) × target. The bias goes, nearly every sample takes one round, and the read names and
numbering stay. The long reads change.

### L4. Homopolymers for indels

```cpp
std::size_t run = 1;
while (at + run < templ.size() && templ[at + run] == templ[at]) ++run;      // forward
while (run < 5 && at >= run && templ[at - run] == templ[at]) ++run;          // backward, from at - run
```

The backward scan starts `run` bases back, not one base back, so it skips the bases just before `at`:
- In a run of 6 (`AAAAAA`), positions 2 and 3 count 4 and get no boost.
- In a run of 5, only its first and last positions get the ×10.
- `AAxxAAA` at the second block's first A counts 5 and gets the boost although its run is 3.

Homopolymer indels are therefore about 40-60% of the intended rate. The rates themselves are tiny (3-5·10⁻⁶ per base,
×10), so the effect on the reads is small. Fix: count backward with its own offset (`at - 1 - k`).

### L5. The Poisson-lognormal abundances flatten the rare half

`draw_weights` makes each species' weight Poisson(λ) + 1 with λ ~ lognormal(`pln_mu`, σ)
(`CommunityProfileDesigner.cpp:283-294`). The build passes only σ: 1.3,2.0 in the design and 1.0-2.5 in the scenarios
(`collect_training_data.py:796`, `scenarios.py:87`). So μ = 0, the median λ is 1, and the Poisson draw puts
P(weight = 1) = E[e^−λ] ≈ 0.40-0.45 for every σ in use. A further 17-23% get weight 2 (numerical integration, ±3
points).

So about 40-45% of a sample's species share exactly the lowest abundance, 1/Σw, and two thirds sit at the two lowest
levels. σ changes only the top of the distribution: Σw grows with E[λ] = e^(σ²/2). Real rank-abundance curves have a
continuous tail. Here the rare tail is a block of ties, and nothing lies below the floor.

This may be intended (Poisson-lognormal is a model of observed counts). If a continuous tail is wanted, use the
lognormal λ itself as the weight, or a μ large enough (e.g. 5) that the Poisson step is negligible. Either changes the
training data and needs a retrain.

### L6. The genome-table header is guessed from substrings

The first non-comment row is taken for a header if any of its fields contains:
- "name", "genome" or "accession",
- and "tax",
- and "path", "fasta" or "file"

(`MetagenomeSimulator.cpp:142-170`).

A build's rows are data rows, but their FASTA paths often contain "genomes" and "files"
(`genomic_files_all/gtdb_genomes_all_r226/...`). If the path also contains "tax" anywhere (a project folder such as
`.../metataxonomics/...`), or the first genome's lineage does (e.g. a genus name with "tax" in it), the first data row
is taken for the header:
- that genome is dropped;
- the name column is mapped to the path column (the accession has none of the words), so every genome is named by
  its path and tables keyed by genome name no longer match;
- if the "tax" came from the path, the taxonomy column is mapped to the path too, so every genome becomes its own
  "species";
- the length column is lost (no header field says "length"), so `build_length_cache` reads every FASTA of the table
  for its length, on one thread, in every run: ~17 min per run at 59k genomes.

Nothing fails loudly. Fix: a header only if its fourth field is not a number, or match whole column names.

### L7-L12, briefly

- **L7.** The CLI defaults are `--sequencer HS25` and `--read_length 150`. HS25's curve ends at cycle 126 and
  `Interpolate` continues its last slope (−0.31 Q per cycle), so cycle 150 averages Q22.5 for R1 and Q22 for R2. The
  build's setups (HS20 100 bp, HSXt 150 bp, MSv3 250 bp) are within their curves. Make the default HSXt, or warn
  beyond the last knot.
- **L8.** `write_all_reads` names `_R2` whether or not it is written. The `--test` placeholders even create it. Both
  manifests and `protal.meta` then list a file that never exists.
- **L9.** The designs use `std::shuffle`, `std::uniform_int_distribution`, `std::lognormal_distribution` and
  `std::poisson_distribution`, all implementation-defined. Their random fill is shuffled from `unordered_map`
  iteration order (`CommunityProfileDesigner.cpp:541-544`, and `taxon_to_species` built in hash order). The same seed
  can give other communities with another standard library or version. The reads do not have this problem (`LongRng`).
  The congener groups sort before shuffling (`:513, 523`); the rest could do the same.
- **L10.** `adjust_counts_to_total` corrects the rounding one pair at a time at random entries (`:21-43`), so a
  species of 2 pairs loses one as readily as one of 100,000. At the build's depths the correction is a few pairs. A
  largest-remainder apportionment would be exact.
- **L11.** A failure closes the pipes (`ReadPipeline.cpp:313-323`). Pieces are whole zstd frames or BGZF blocks, and a
  pipe skips the BGZF end-of-file check (`ThreadedGzStream.h:82-83`, regular files only), so the reader sees a shorter
  but valid sample. `run_protal` stops protal when the simulator fails, which covers it in the build.
- **L12.** `LoadContigs` drops contigs shorter than the read (pe) or 100 bases (long reads), while
  `vertical_coverage` divides by the whole `genome_length`. This is negligible except for very fragmented genomes.

## 2. Speed: where the time goes

| Work | Cost (laptop core) | Source |
|---|---:|---|
| a read pair (model, FASTQ, zstd) | 18 µs | Illumina report |
| a genome, per sample and round (inflate + parse) | 12-17.5 ms | Illumina and long-read reports |
| HiFi / Ultima / Nanopore per Mb | 0.049 / 0.075 / 0.044 s | long-read report |

A genome costs as much as 700-1,000 pairs. The share of genome reading by sample (genomes per sample from the
scenarios' log-uniform species counts × 1.4 strains):

| Sample | Genomes | Pairs | Genome reads | Pairs' work | Genome share |
|---|---:|---:|---:|---:|---:|
| soil pe | ~8,600 | 20M | ~150 s | ~360 s | 30% |
| soil_shallow pe | ~8,600 | 5M | ~150 s | ~90 s | 63% |
| design pe, 200k pairs | ~190 | 200k | 3.3 s | 3.6 s | 48% |
| design pe, 2k pairs | ~190 | 2k | 3.3 s | 0.04 s | 99% |

### 2.1 How a genome is read

`LoadContigs` (`ReadPipeline.cpp:115-154`) does this for every (sample, round, genome):
1. It opens the FASTA through `ThreadedGzIstream`. That starts a thread and allocates 4 × 1 MiB of zero-filled
   blocks, plus an 85 KB inflate state and a 128 KB input buffer (`ThreadedGzStream.h:94, 103, 298`): ~0.5-1 ms per
   file before any data.
2. It reads line by line with `std::getline`: ~42k lines of 80 bases for 3.4 Mb.
3. It appends **base by base** with a `std::isspace` and a `std::toupper` call each: both glibc table lookups through
   a thread-local locale pointer.

The arithmetic for a 3.4 Mb genome:
- inflate ≈ 3.4 ms (ISA-L, ~1 GB/s, on the helper thread);
- the per-base loop ≈ 10-12 ms;
- lines ≈ 1-2 ms;
- the thread and buffers ≈ 0.5-1 ms.

That totals the measured 17.5 ms, so the parse is ~70% of a load.

**Fix (the same reads):**
- read the compressed file whole and inflate it in one call, sized from the gzip trailer's ISIZE or grown as needed;
  no thread, no 4 MiB;
- find lines with `memchr` and append whole lines;
- upper-case in place with a 256-byte table;
- scan for whitespace only in lines that have it.

A load becomes ~4-5 ms, about 4× cheaper and bound by the inflate. For the build that is 5-7 → ~1.5-2 core-hours
(estimate).

### 2.2 How often a genome is read, and what to merge

The collector's calls (default r226 build; the subagent's map, spot-checked):
- **Calls.** Only `collect_training_data.py` calls the simulator: ~100 runs in the training collector and ~70 in the
  test one. The two run side by side on the same `genomes_simulated.tsv` (~54k genomes, lengths in its fourth
  column, so `build_length_cache` reads nothing: `build_gtdb_database.py:375-432`).
- **The FASTAs.** They are NCBI's single-member gzip on the network file system. Nothing copies them to `--scratch`.
- **pe points.** One run each; design runs with `--test` read no genome.
- **se units.** Reuse the pe R1 files (`collect_training_data.py:1345-1356`). A streamed pe point gets a second, full
  run with `--first_reads_only` instead (`:1946-1950`).
- **Long-read and Ultima units.** One run per unit, replaying the same communities as the pe runs of their point or
  scenario.
- **Streaming.** Nearly every scenario simulation and the 10M/30M design points are streamed: one sample at a time,
  each simulation into its own protal run.

Genome loads per build (estimate):

| Part | Loads |
|---|---:|
| design, training (309 pe samples × ~190 genomes; 2 × 186 long samples × 1.3 rounds) | ~0.15M |
| design, test | ~0.05M |
| soil (14 samples × ~8,600 genomes × pe + Ultima + pb + ont at 1.3 rounds = 4.9) | ~0.59M |
| soil_shallow (× 3.6) | ~0.43M |
| moderate (14 × ~3,500 × 4.9) | ~0.24M |
| gut, host | ~0.04M |
| **total** | **~1.5M loads of ~54k genomes (~28 each)** |

At 12-17.5 ms that is 5-7 core-hours, and ~1.7 TB of gzip at 1.1 MB per genome. A soil-table genome is read ~30 times;
a design-only genome 3-4 times.

What could be merged or shared:

| Option | Loads saved | Reads change? | Effort | Conflicts |
|---|---|---|---|---|
| **a. bulk parse** (2.1) | none, but each ~4× cheaper | no | small (one function) | none |
| **b. decoded genome store**: each genome once per build as 2-bit sequence + N runs + contig index on `--scratch`, memory-mapped as `Host` maps `host.seq` | ~95% of the remaining cost; network reads once per genome; shared by every process and both collectors through the page cache | no (the same contigs; filtering by length at view time) | medium: a build step (~54k × ~5 ms ≈ 5 core-minutes, parallel), a `Contigs` view over the map, 2-bit decode of each fragment (~350 bases) | ~50 GB of scratch (~200 GB as plain bases); v12 ran out at 263 GB before streaming |
| **c. genome cache in a run** (by FASTA path, across samples and rounds, 2-bit, byte budget) | soil run: ~120k → ≤ ~18k loads (the scenario table) with ~16 GB; fewer hits with less | no | small-medium | memory beside protal; nothing across processes (pe, pb, ont of a community are separate runs) |
| **d. long-read placement that fits** (L3) | rounds → 1: ~25-50% of long-read loads | yes (long reads) | small | none |
| **e. read types of a community in one run** (pe + Ultima + pb + ont) | scenarios 1.3M → ~0.35M | no | large | each streamed simulation feeds its own protal run: a merged run needs four protal runs reading at once (4 × the r226 index, ~35 GB each), or a protal that profiles four read types in one pass |
| **f. samples of a run genome-major** (each genome once for all open samples) | soil point ~2.6× | order of reads in files | medium | streaming opens one sample at a time, so the streamed scenarios gain nothing |

Recommendation:
- **a now.** It is small, exact, and a ~4× cheaper load.
- **b next**, if `--scratch` can hold ~50 GB. It makes the remaining loads nearly free without touching the
  collector's scheduling, streaming or keys. The Host code is a template for it.
- **c** if scratch cannot.
- **d** with the next change that alters reads.
- **e and f do not pay** under streaming. b gives most of e's saving without merging anything.

### 2.3 Other speed findings

- **Threads wait for another thread's load** (`ReadPipeline.cpp:372-391`). Items are issued in order, and the parts of
  a genome follow each other. The first part loads the genome; every thread that takes another part meanwhile blocks
  in `loaded.get()`.
  - A 30M-pair design sample with ~190 genomes has ~95 items of 1,666 pairs per genome (~30 ms each).
  - At 84 threads, ~50 threads finish an item during a 17.5 ms load and all block on it: ~10-15% of the threads idle
    (estimate). At 4 threads it is ~2%.
  - Fix: start the next genome's load when a genome's first part is issued (load-ahead), or let a thread take the next
    item whose genome is ready. The pieces are written in order anyway, so the reads do not change. Fix a makes this
    3-4× smaller in any case.
- **Per pair.** `IlluminaModel::Read` draws per base:
  - one normal for the AR(1) noise, by Marsaglia's polar method: two uniforms, a log and a square root per two
    normals;
  - four uniforms: low state, indel, N, error.

  Estimate (operation counts against the measured 18 µs): the normals ~40-50% of `Read`, the uniforms ~20%, and zstd
  level 3 ~15% of a pair. Options:
  - a ziggurat normal (~3× faster);
  - geometric skipping for the rare events: the indels by thinning at their homopolymer maximum, the Ns, the
    low-state entries; leaving one uniform per base, for the error;
  - read 2 on a random stream of its own (`IlluminaSimulator.cpp:491` makes it "either way" to keep read 1 the same).
    `--first_reads_only` could then skip it. That run is the streamed se units' full second simulation, ~45% of its
    pairs' work.

  Together ~25-35% less per pair (estimate; a callgrind profile should rank them first). They change the reads, so
  bundle them before the next r226 build, which is the first on the in-process reads anyway.
- **Compressing into pipes.** A streamed sample's pieces are compressed (zstd level 3) only for protal to inflate
  them again. protal's reader takes plain FASTQ as it is (`ThreadedGzStream.h:396-397`, `Mode::Copy`). A raw packing
  for FIFO outputs saves the simulator's compression and protal's inflate (estimate ~10-15% of the simulator's pair
  work) with the reads unchanged. The ~600 kB pieces against protal's 4 MiB of read-ahead per file keep R1 and R2 in
  step as now.
- **The design is O(table) per sample, on one thread.** `design_profile` regroups the whole table for every sample
  (`group_by_species`: a copy of every `GenomeRecord`, the taxonomy split again). It rebuilds the genus and taxon maps
  (two keys per rank per species) and shuffles all species (`CommunityProfileDesigner.cpp:301-333, 541-544`).
  - At ~59k genomes that is ~0.1-0.2 s per sample (estimate), ~1-2 minutes per build. It runs before a run's reads
    start, with its threads idle.
  - Build the groups and maps once in the constructor; the designs stay the same.
- **Long reads per base.** `hifi::Mutate` takes an `exp` and a `log10` per base (`LongReadSimulator.cpp:178, 201`),
  and Ultima a `log10` and a `pow` more (`:181, 187-188`). The written quality follows from the weight's logarithm, so
  one transcendental per base can go: ~20-30% of the read model's time (estimate). It is minor at the build's volumes
  (~1.8 laptop core-hours of PacBio).
- **No problem:**
  - `Calibrate`: ~0.1-0.2 s once per run.
  - The in-flight cap (3 × threads).
  - `build_length_cache` with the table's lengths. It does print a line per 100 genomes: ~600 lines of
    "genomes processed" per run.

## 3. Recommended order

1. **Exact, no retraining needed:**
   - the bulk FASTA parse (2.1);
   - the header guard (L6);
   - load-ahead in the pipeline;
   - the design's groups built once;
   - raw FASTQ into pipes;
   - the homopolymer scan (L4) — this one changes reads slightly; bundle it with item 2 if preferred.
2. **Before the next r226 build, read-changing**, since that build retrains on new reads anyway:
   - long-read placement that fits, with the stop tolerance (L3);
   - R2 on its own stream and skipped with `--first_reads_only`;
   - the faster normal and rare-event skipping;
   - a decision on L5's flat tail.
3. **The decoded genome store** (2.2 b), or the in-run cache (c) if scratch is short.
4. **CLI-only:**
   - the replay seeds (L1);
   - strain-sharing floors per depth and strain (L2);
   - the default sequencer (L7);
   - R2 names with `--first_reads_only` (L8);
   - sorted keys before shuffles (L9).

## 4. Checked and found right

- **Pipeline ordering and determinism.** Each item has its own stream (`MixSeed(seed, part)`; for long reads
  `(seed, round, genome, part)`). Pieces are written in item order by one writer per sample. Plans for long reads come
  from a per-sample stream under the lock. So the files do not depend on the threads.
- **No deadlock.** Later samples open only when every open sample's items are issued, so the head item of each sample
  is always in flight.
- **R1/R2 in step on pipes.** protal takes 32 records per file per batch (`SeqReader.h:168, 219-220`). Its reader
  inflates up to 4 MiB ahead per file. The simulator leads R1 over R2 by one ~600 kB piece at most.
- **Fragments.** Placement is uniform over the positions where a fragment fits, by rejection that keeps the drawn
  length. Read 2 is the reverse complement's start, and the strand is chosen 50/50.
- **Random numbers.** Lemire's unbiased `Below`, Marsaglia-Tsang's gamma, splitmix-seeded xoshiro256**.
- **Error paths.** Partial files are removed, failures stop all threads, truncated FASTAs are refused
  (`read_failed`).
- **Host reads.** The host sequence is upper case (`scenarios.py:453`), so no soft-masked base turns into an N.
- **qshmm tables.** They follow pbsim3, including its >50-state overflow (long-read report).

## 5. What this report did not do

- It changed no code and ran nothing.
- The costs per pair and per genome are the earlier reports' measurements. Their breakdowns (the parse ~70% of a load;
  normals, uniforms and compression in a pair) and the build-wide counts (1.5M loads) are estimates from the code and
  the collector's defaults.

## 6. How to check

- **One core, about a minute:** a callgrind profile of a shallow sample (400 genomes, 4,000 pairs) and a deep one
  (2M pairs). It would confirm the load breakdown (`LoadContigs` against `ThreadedGzStreambuf::Inflate`) and the
  pair's (`LongRng::Normal`, `Uniform`, `ZSTD_compress2`).
- **The next r226 build's logs:** the `rounds` column of `--long_stats`, and `simulated (...) in T on K threads`
  against the genome counts.

## 7. Follow-up: implemented and measured

The user then asked for:
- the exact fixes;
- the 2-bit genome store;
- plain (uncompressed) pipes, behind a flag;
- the three per-pair changes;
- a benchmark before and after.

Committed as `746ea66` on top of `c9b46f9` (the code, scripts, tests and docs), this report in the next commit. That
commit, built from its `git archive`, passes the simulator's 17 unit tests. The work
was built and run in WSL (`~/simaudit`) with the scripts in [`scripts/`](scripts/) and the raw results in
[`results/`](results/).

### 7.1 What changed

| Part | Where | Reads |
|---|---|---|
| **Bulk FASTA reader.** `ReadWholeFile`: read whole, gzip (any members, BGZF) or zstd inflated in one go with ISA-L or libzstd; no thread and no 4 MiB of buffers per file; a cut or corrupt file refused. `ParseFasta`: lines found with `memchr`, whole lines copied and upper-cased in one pass, a line with white space character by character; the same records as the line-by-line reader, tested against it on quirky FASTA (CRLF, lower case, spaces, IUPAC, bytes above 127, text before the first header) | new `src/RandomForest/GenomeStore.{h,cpp}`; `LoadContigs` (`ReadPipeline.cpp`), `read_genome_length` (`MetagenomeSimulator.cpp`) | the same |
| **Genome store.** `--genome_store DIR`: a FASTA is parsed the first time a run needs it, and its file `DIR/<FNV hash of the path>.g2b` written through a temporary name and a rename. The file holds a header with the FASTA's path, size and modification time, then the contigs, the runs of non-ACGT characters, and the bases at 2 bits each. Later loads memory-map it (`MAP_POPULATE`) and decode only the fragments drawn (`Contigs::Extract`). A FASTA changed since is read again, and a cut or foreign file is ignored. Every write goes to a unique temporary name before the rename, so concurrent writers and readers are safe. A failed write warns once and the run goes on from the FASTA | `GenomeStore.cpp`, `Contigs` (`ReadPipeline.h`), `--genome_store` | the same |
| **Header guard (L6).** The first row is a header only if, besides the words, no field holds a lineage (`;`), a path (`/`) or a number | `read_genome_table` | the same |
| **Load-ahead.** A thread that would wait for another's load of the next item's genome loads a later genome of the round instead (at most 2 × threads genomes held) | `Engine::WouldWait`, `Preload` (`ReadPipeline.cpp`) | the same (any thread count) |
| **The design's groups once.** Species by genome, genus and taxon are made at the first design and kept. The random fill shuffles pointers to the species, which `std::shuffle` permutes exactly as it permutes the names | `CommunityProfileDesigner` | the same |
| **Plain pipes.** `--plain_pipes`: outputs that are named pipes get plain FASTQ (protal's reader takes it as it is; regular files stay compressed). The pipeline now packs the pieces itself after `Make` (`Packing::Plain`) | `ReadPipeline.cpp`, `--plain_pipes` | the same reads |
| **Ziggurat normals.** `LongRng::Gaussian`, Doornik's ZIGNOR with 128 blocks: one 64-bit draw a number in ~99% of cases. Used throughout the Illumina model; the long reads keep `Normal()` | `ReadPipeline.cpp`, `IlluminaSimulator.cpp` | new reads, same distribution |
| **Rare events skipped geometrically.** The low state's entries (thinned from the highest entry rate to each cycle's), indels (thinned from the homopolymer rate to each place's) and Ns are drawn as gaps between candidates, not as three uniform numbers per cycle | `IlluminaModel::Read` | new reads, same distribution |
| **Read 2 on its own stream.** Each item has a second random stream for read 2, so `--first_reads_only` skips read 2 altogether and its read 1s stay those of the paired run | `IlluminaModel::Pair`, `PairedJob::Make` | new reads |
| **Homopolymer fix (L4).** The backward scan counts from the base before | `IlluminaModel::Read` | a few more homopolymer indels |
| Collector and build | `collect_training_data.py --genome_store DIR`, `--compressed_pipes`; `build_gtdb_database.py --genome-store auto\|DIR\|none` (default `auto`, at the user's request: `SCRATCH/genome_store`, or `OUTDIR/genome_store` without `--scratch`, kept for the next build), `--compressed-pipes`. Streamed simulations pass `--plain_pipes` unless `--compressed_pipes` is given. The store's argument goes to the commands that run and never into a simulation's key | the same keys |
| Docs | `docs/development.md` (options, input), `docs/databases.md` (streaming, options) | |

### 7.2 Exactness and statistics

**Stage A** is everything except the per-pair changes and the homopolymer fix. Against `c9b46f9` it writes the same
bytes ([`verify_exact.sh`](scripts/verify_exact.sh)), with manifests compared after each run's folder is replaced in
their paths:
- paired-end reads: 3 samples with strains, species 30-60 and two depths, in zstd and in BGZF;
- the same with a genome store, both when the run writes it and when the next one reads it;
- first reads only;
- HiFi and Ultima samples (2 × 20 Mb), with and without the store;
- the designs of 20 samples from a 50,000-genome table (`--test`).

**Stage B** changes the Illumina reads by design. Its read statistics stay within sampling noise of the old ones on
every instrument ([`profile_stats.sh`](scripts/profile_stats.sh), `--illumina_report`, 50,000 pairs each;
[`results/profile_stats.tsv`](results/profile_stats.tsv)):
- mean written quality within 0.1;
- Q30 share within 0.003;
- substitutions within 2%;
- Ns ~1.1e-4.

Indels are up a few percent where homopolymers now get their factor everywhere (MSv3 R1 insertions 4.05 → 4.76e-5,
R2 4.45 → 4.41e-5; ~500 events each, so ±5%).

**Store files.** The table-driven writer writes the same files as the first, per-base one (400 of 400 the same).

**Tests.**
- New `tests/test_GenomeStore.cpp` (5 tests):
  - the parser against the old line-by-line reader;
  - inflating plain, gzip (several members, a full header), BGZF and zstd, and refusing cut, corrupt or trailing data;
  - every slice of a store file;
  - staleness;
  - contigs alike with and without the store;
  - the header guard.
- `test_IlluminaSimulation`:
  - the ziggurat's moments and tails (4M numbers);
  - rare events at exaggerated rates (indels and Ns per base, the low state's stationary share);
  - the store giving the same files;
  - plain pipes carrying exactly the FASTQ of the files.
- `test_LongReadSimulation`: the store gives the same files.
- `test_collector.py`: the store gives the same long reads, and the arguments.
- `test_gtdb_pipeline.py`'s `test_h_streamed` now also gives `--genome-store auto`: a mini GTDB build with every
  design point streamed (plain pipes) and the genomes from the store gives the same training tables as the build that
  writes its reads.

What ran:
- all 428 unit tests (426 pass; 2 skipped, in `test_PackedSequence` and `test_ExactSum`, which the change does not
  reach);
- the collector's 14;
- `test_h_streamed` (with the full build it compares against, 46 s);
- the end-to-end `SimulatorTest` (4).

With the store made the build's default:
- `test_gtdb_build.py` (16);
- `test_a` (full build, rerun and reduced database) with `test_h`: 73 s.

### 7.3 Benchmark

[`bench.sh`](scripts/bench.sh) ran three binaries on the 400 synthetic genomes of real size (median 3.4 Mb, gzip
level 6):
- `before`: `c9b46f9`;
- `A`: stage A;
- `B`: stage A plus the per-pair changes.

Setup:
- paired-end runs: HSXt, 150 bp, zstd, one sample of 400 species;
- every run on 4 pinned cores, niced;
- 3 repetitions, the variants alternated within each.

The table gives medians of CPU seconds (user + system); [`results/bench.tsv`](results/bench.tsv) has every run. The
first repetition ran on the laptop's low-power cores (up to 2.7× slower, every variant alike), the other two agree
within a few percent. Hence the medians, and the instruction counts below as the steadier comparison.

| Case | before | A | B | with the store | |
|---|---:|---:|---:|---:|---|
| shallow: 4,000 pairs of 400 genomes, 1 thread (genome reads) | 3.92 | 1.54 | 1.52 | 0.15 (A), 0.10 (B) | per genome ~9.4 → ~3.5 ms (2.7×), from the store < 0.4 ms (>20×) |
| the store written by that run (B2: table-driven writer) | | | | 2.10 against 1.47 without | writing a genome's file ~1.6 ms (the first writer ~12 ms) |
| deep: 2M pairs of 400 genomes, 1 thread | 27.6 | 24.4 | 19.0 | 18.7 (B) | per pair ~11.5 → ~8.8 µs (−23%); in all −31% |
| deep, 4 threads: wall (CPU) | 8.21 (29.7) | | 5.54 (21.3) | | −33% wall |
| first reads only, 2M pairs, 1 thread | 24.6 | | 9.2 | | −63%: read 2 no longer made |
| into named pipes, 2M pairs, 1 thread (the simulator's CPU) | | | 19.1 compressed | 15.2 plain | −20% (and protal's inflate) |
| HiFi, a 1 kb read a genome (genome reads), 1 thread | 1.79 | 0.57 | | 0.01 (A) | |
| HiFi, 250 Mb of 400 genomes, 1 thread | 17.5 | 12.0 | | 8.8 (A) | genome reads were half of it (rounds read genomes again) |
| design: 20 samples from 50,000 genomes, `--test` | 3.23 | 0.26 | | | 12× |

Memory:
- without the store, 30-33 MB (a whole decompressed FASTA in memory while it is parsed) against 22 MB before;
- with the store, 10-23 MB (pages mapped, not copied);
- 4 threads, B: 106-116 MB against 62 MB before (a FASTA per thread in memory while it is parsed, and the read-ahead
  genomes).

What it means at r226 (estimates from section 2.2's ~1.5M genome loads):
- **Bulk reader alone:** 5-7 → ~2 core-hours of genome reading.
- **With `--genome-store`:** ~54k FASTA reads (~2 min of writing on one core) plus 1.5M mappings at < 0.4 ms: under
  10 core-minutes.
- **Per-pair changes:** about a quarter off every Illumina pair; the streamed se runs about 60% cheaper.
- **Plain pipes:** a fifth off every streamed simulation, besides protal's inflating.

### 7.4 Instructions

[`profile_ir.sh`](scripts/profile_ir.sh) counts instructions with callgrind, on one thread and 40 genomes. Unlike the
timings above, these counts do not depend on which core a run lands on.

Billions of instructions; [`results/instructions.tsv`](results/instructions.tsv), with the functions in
[`results/instructions_by_function.txt`](results/instructions_by_function.txt):

| Case | before | A | B | B, store read |
|---|---:|---:|---:|---:|
| deep: 100,000 pairs | 23.85 | 17.17 (−28%) | 13.87 (−42%) | 11.47 (−52%) |
| shallow: 400 pairs (genome reads) | 10.68 | 4.00 | 3.64 | 1.14 (A: 1.51) |
| first reads only, 100,000 pairs | 22.53 | | 8.58 (−62%) | |

**Per pair** (deep minus shallow, over 99,600 pairs):

| Binary | instructions per pair |
|---|---:|
| before | 132.2 k |
| A | 132.2 k (the exact part leaves the pairs alone) |
| B | 102.7 k (−22%) |
| B, from the store | 103.7 k |

B's pair goes to:

| Function | Share |
|---|---:|
| `IlluminaModel::Read` | 21% |
| `Gaussian` | 16% |
| zstd | ~15% (what `--plain_pipes` saves) |
| `Uniform` | 7% |
| libm's `lround` | 7% |
| `ReverseComplement` | 5.5% |

In `before`, the normals, with their uniform numbers and logarithms, were ~35% of a pair, and the per-cycle uniform
numbers of the rare events another ~20%.

**Per genome read** (shallow minus the fixed costs that the store run leaves: calibration, design, 400 pairs; over 40
genomes):

| Binary | instructions per genome |
|---|---|
| before | ~229 M (~67 a base): the parse loop 3.6, `toupper` 1.6, `isspace` 1.0 and `getline` 0.2 billion, against 1.6 billion of inflating, of 40 genomes |
| A | ~62 M (3.7× fewer), two thirds of it ISA-L's inflate |
| from the store | too small to show beside the fixed costs |

**Not done, the next per-pair costs:**
- `lround`, a libm call per base (7%), could be an inline rounding;
- `ReverseComplement`'s switch per base (5.5%) could be a table;
- `Gaussian` could be inlined into `Read`;
- the inflate buffer's zero-fill is ~7% of a genome read.

None of these is measured.

### 7.5 Left open

- **`--genome-store auto` is the default of `build_gtdb_database.py`** (the user's choice after the benchmark). It
  needs ~50 GB at r226 beside the samples (v12 ran out at 263 GB before streaming and `--keep-free`), so
  `docs/databases.md` now asks for 200-250 GB of scratch. The store is kept for the next build, and the build's last
  lines give its size. `collect_training_data.py` on its own uses one only when given `--genome_store`.
- **Retraining.** The per-pair changes give new Illumina reads, with the same statistics. Models trained before
  (all of them on ART reads up to r226 v15) are retrained by the next build anyway.
- **Not fixed in this commit:** L1 (replay seeds), L2, L3 (long-read placement and rounds, which also cost genome
  reads: half of the HiFi time above), L5 (decision), L7-L12. L1, L3, L5, L8, L9, L11 and L12 followed: see
  [section 8](#8-follow-up-2-the-logic-findings-fixed).
- **Not run:** the rest of the end-to-end suite (protal on the mini database, which the simulator does not reach) and
  the other mini GTDB pipeline tests (`test_a`-`test_g`); see 7.2 for what ran.
- **The website:** none of this is on it. The simulator and the build are documented in `docs/` only.

## 8. Follow-up 2: the logic findings fixed

The user then asked to fix:
- L1 (replay), L3 (long-read placement) and L5 (the Poisson-lognormal tail): a long tail is normal in ecological
  species distributions, but its species should keep some abundance;
- L8, or its documentation;
- L9, L11 and L12.

L2, L4 (fixed in section 7), L6 (fixed in 7), L7 and L10 are not part of it. Committed as `436821c`, whose tree is
the one tested below.

| # | Fix | Where | Outputs |
|---|---|---|---|
| L1 | Each sample's `run_seed` and `host_seed` are drawn when it is designed (the values runs used before) and written to both manifests (new columns on the right). A replay takes them from there, so neither `--seed` nor the sample's place in the replay matters, and a per-sample manifest replays its sample's reads. A manifest without the columns falls back to the run's seed and says so. `run_params.tsv` records `--host_folder` and `--host_pairs` too, and a replay warns when they differ. The help text, the replay's messages and `docs/development.md` say what is reproduced | `MetagenomeTypes.h` (`SampleOutput::run_seed`, `host_seed`), `MetagenomeSimulator.cpp`, `simulate_metagenomes_main.cpp` | reads unchanged; manifests gain two columns |
| L3 | A long read's template is placed uniformly where it fits, a place that runs off its contig being drawn again (up to 1,000 times), as the Illumina fragments are. Only a template longer than every contig of its genome is still cut at its contig's end. A genome of short contigs now gets its weight's share, and a sample nearly always needs one round, so no genome is read again | `LongJob::Make` (`LongReadSimulator.cpp`) | new long and Ultima reads |
| L5 | A new distribution, `lognormal`, is the simulator's default and the collector's `--abundance lognormal:...` (the build's and the scenarios'). Species weights come from the lognormal itself, so the tail is continuous instead of 40-45% of the species tied at weight 1. No weight goes below `--abundance_floor` (0.001) times the median: at σ 2.5 that is the lowest ~0.3% of the species, below σ 2 hardly any. Every species still gets a read pair or more. The old model stays as `poisson_lognormal` (`--abundance poisson_lognormal:...`) | `CommunityProfileDesigner::draw_weights`, `MetagenomeTypes.h`, `collect_training_data.py`, `build_gtdb_database.py`, `scenarios.py`, docs | new designs |
| L8 | With `--first_reads_only` no `_R2` path is set: none is written (not even a `--test` placeholder), the manifests' `fastq_r2` is empty and the protal map's `SECOND` is `-` | `MetagenomeSimulator.cpp`, `simulate_metagenomes_main.cpp` | manifests and map |
| L9 | The designs no longer depend on the standard library. `PortableRandom.h` has the draws on `std::mt19937_64`'s numbers (whose sequence the standard fixes): uniform, bounded integers (Lemire), Fisher-Yates shuffles, normals (Box-Muller), gamma (Marsaglia-Tsang), Poisson (multiplication below 10, PTRS above, as numpy), negative binomial. Every list that a hash table ordered is sorted by name first: the species, each genus's and taxon's species, the genus and taxon requests and the capping of their quotas, forced strains, and the abundance matrix's rows. Only libm's `exp`, `log` and `lgamma` remain platform-dependent in their last bits | `PortableRandom.h`, `CommunityProfileDesigner.cpp`, `MetagenomeSimulator.cpp` | new designs |
| L11 | A failed run's open named pipes end, written without blocking, with the start of a zstd frame or gzip member, or a FASTQ header without its record (`--plain_pipes`). protal's reader reports each ("truncated file?"), so a cut sample can no longer pass for a whole one | `Engine::~Engine`, `Poison` (`ReadPipeline.cpp`) | failures only |
| L12 | `vertical_coverage` stays read bases over the genome's length, the genome-wide mean depth, which a `--test` design and its real run must agree on, and neither reads the genomes. A run now notes each genome with 1% or more of its bases in contigs shorter than any read (at most 20 a run), and `docs/development.md` says what the column is | `LoadContigs` (`NoteShortContigs`), docs | a note |

**Tests.**
- `CommunityDesign.LognormalTailIsContinuousAndFloored`: σ 2.5, 300 species, 20M pairs, 5 seeds. Under 3% of the
  species sit at the lowest count against over 30% with `poisson_lognormal`, and none falls below the floor; a floor
  of 0.1 floors the expected ~18%.
- `PortableRandom.Distributions`: the moments of Poisson (λ 0.3 to 1,000), gamma, negative binomial and normal
  numbers; the 6 orders of 3 equally often; the same numbers from the same seed.
- `IlluminaSimulation.ReplayFromManifestsMakesTheSameReads`: a run of 3 samples with host reads, replayed from its
  combined manifest with another seed and from one sample's manifest, gives the same bytes. First reads only gives
  the same `_R1`, no `_R2` and an empty `fastq_r2`.
- `IlluminaSimulation.AFailedStreamIsCutOff`: a genome missing after one that streamed, into zstd and into plain pipes.
  The zstd stream is refused as cut, and the plain one ends with the incomplete record.
- `LongReadSimulation.ShortContigsGetTheirShare`: a genome of 30 contigs of 2 kb and one of a single 60 kb contig,
  equal weights, 1 kb HiFi reads. Each gets half the bases (±3%), in at most 2 rounds.

**What ran.** The working tree held another session's unfinished changes to protal's core, which did not compile
then. So the tests ran on `55c5f35` plus this change alone (its patch on the commit's `git archive`, in
`~/simaudit/mine`):
- all 433 unit tests (431 pass, the same 2 skipped);
- `test_collector.py`, `test_gtdb_build.py`, the end-to-end `SimulatorTest`;
- `test_gtdb_pipeline.py`'s `test_a` (full build, rerun, reduced database) and `test_h_streamed`.

**L3 at benchmark scale** ([`long_shares.py`](scripts/long_shares.py), `results/long_shares.txt`). One HiFi sample of
250 Mb (15 kb reads) from the 400 synthetic genomes (median 3.4 Mb, ~78 contigs each), one thread. Templates cut at
contigs' ends (`746ea66`) against placed where they fit (this change). Ratios are each genome's share of the bases
over its weight's share, for the 135 genomes expected to get 300 kb or more, by thirds of their base-weighted contig
length:

| | cut | placed |
|---|---:|---:|
| rounds | 6 | 2 |
| reads | 20,577 | 16,656 (none cut short) |
| CPU (user) | 15.7 s | 13.6 s |
| genomes of the shortest contigs (~47 kb) | 0.845 | 0.991 |
| middle (~98 kb) | 1.044 | 0.978 |
| longest contigs (~221 kb) | 1.163 | 0.985 |
| ratio's SD over the 135 genomes | 0.228 | 0.161 (sampling noise of ~20+ reads each) |

The cut templates gave the most fragmented third of the genomes 15% less than their weight and the least fragmented
16% more, a 38% spread that is now gone. The sample needed 6 rounds (the collector's estimate assumes 1.3), each
reading again the genomes it drew; now 2. At r226, where most genomes are MAGs, the spread was likely larger.
