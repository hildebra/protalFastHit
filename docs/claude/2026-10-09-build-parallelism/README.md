# The GTDB build on 84-120 cores: what runs serially, and what can run in parallel

Date: 2026-10-09. Branch `audit-fixes` at `7ace9b8` (`scripts/build_gtdb_database.py`, `collect_training_data.py`,
`mini_db/gtdb_to_protal_db.py`, `mini_db/gene_neighbours.py`, `insilico_strains.py`, `trace_relatives.py`,
`error_reads.py`, `ancestry_sites.py`; protal's `src/RunProtal.h`, `src/Classify.h`, `src/Build.h`,
`src/SequenceUtils/CongenerGaps.h`, `src/Alignment/AlignmentScreen.h`).

Data: the r226 v19 build's share archive (`local/v19`; SLURM job 24099010, built at `a81cee3`, `-t 84`, a node with
502 GB, `--scratch` on the node's SSD, the inputs and OUTDIR on the `/hpc-home` network file system), and the
`console.log` of v17 and v18 for how much the steps vary.

Method: read-only. I read the code (four read-only agents went through the converter, the genome-level helpers, the
collector and the reports, and I checked their main claims in the code). The protal-run timings come from
[`scripts/run_timings.sh`](scripts/run_timings.sh) (awk over `logs/protal_runs_*.log`). Nothing was built or run: the
local cores are reserved, and the question is about the cluster's 84-120.

**Question.** `build_gtdb_database.py` runs on jobs of up to 120 cores. Which of its parts are serialised, and so
start to dominate as the core count grows? Can they be sensibly parallelised?

## Summary

v19 took 3:49 in all, and the database was ready after 3:31. The time splits into four parts:

| Part | v19 wall | What limits it |
|---|---|---|
| Setup: genome table to training database | 0:40 (v18 1:02) | A chain of steps, one after the other. Three of them read every genome from `/hpc-home` (I/O, not cores). The training database's files are derived on **one process** (5:17). The build itself is 14:41, of which "Congener gaps" is 10:24 |
| Collection (both collections) | 2:46 | Six protal runs, **strictly one after another** (the protal lock). Each aligns its samples one by one on all threads, then profiles them. ONT alignment is 70% of the aligning, and it is compute: the k-mer screen cannot refuse ONT candidates |
| Training to "Ready" | 0:04.5 | Already parallel (4 trainers × 21 threads). Parity is serial (27 s) |
| Reports after "Ready" | 0:18 | Serial waits: the contig names are read from `/hpc-home` before any report starts (7:33). The ancestry unpack starts only after every error-read report has ended. Each read type's ancestry pass repeats the same two reference passes |

**Truly serial on the critical path** (the cores idle or one process works): about **35-45 min** at 84 cores, plus
the 18 min of reports. This part does not shrink with more cores. The rest (ONT/PacBio alignment, congener gaps, the
index build, the trainers) is parallel and shrinks roughly with the core count. At 120 cores the serial share would grow
from about a quarter of the run to about a third, more if the short reads do not scale.

The serial parts, largest first:

1. **NFS-bound genome reads in setup** (genome table 8:47-22:32, gene neighbours 3:18-9:44 across v17-v19). Every
   genome is read from `/hpc-home` up to five times per build: lengths, gene neighbours, in-silico strains (one-genome
   species), the genome store, and contig names. More cores do not help; reading once and caching does.
2. **The derive of the training database's files** (5:17, one Python process; `gtdb_to_protal_db.py:642-735`). One
   `zstd -dc` stream feeds a Python loop over ~80M full-reference records, which writes through one `zstd -T8` per
   gene. The source has one zstd frame per gene and a seek table, so the work splits by frame.
3. **Profiling stages and start-ups inside the collection's protal runs** (~16 + ~5 min over six runs). After the last
   sample is aligned, each run profiles all its samples. `--profile_ahead` overlaps that stage with the alignment, but
   it is off by default.
4. **The setup chain itself.** In-silico strains (2:24) feed only the simulations, yet they sit on the path. The genome
   table does not depend on the converter. `gene_neighbours.py` could run beside the converter's last stage. This is
   change 4 of the [build-ordering report](../2026-10-07-build-ordering/README.md), still open.
5. **Reports** (18 min → ~5-6 min possible): contig names taken in the genome-length pass, the unpack started first,
   one reference pass for all read types, threads by SAM size.

Two costs are **not serialisation**, but they set the length of the collection:
- **ONT alignment: 100 of 143 min of aligning.** The k-mer screen is exact (q-gram lemma). At ONT's
  `max_score_ani` 0.85 the bound refuses next to nothing: 1,191 of 17.9M candidates on a host sample, against 17.2M
  of 17.8M on the PacBio sample of the same community. Every ONT candidate (~30 per read) goes through WFA2, and 99.8% of
  them fail on host reads. More cores help linearly here. A cheaper filter for ONT would help far more, but it would
  change outputs; that is a separate decision.
- **Short reads seem to use the cores poorly.** Paired-end alignment of 150 bp reads ran at ~1.0-1.4M pairs/s on 84
  threads. That is about 5× the thread time per pair of the 32-thread cluster benchmark on a real sample. This is
  either memory latency (the seeding is memory-bound and the node may span NUMA domains) or idle cores. The logs cannot tell which: **the build
  records no CPU time anywhere.**

Ranked changes (savings are v19-based estimates at 84 cores):

| # | Change | Saves on the critical path | Effort, risk |
|---|---|---|---|
| 0 | **Record CPU time**: per step (the `wait4` the script already calls returns it), per protal run, per sample in protal, and a cores-busy timeline from the job's cgroup | nothing by itself; it settles items 7-8 | small, none |
| 1 | Derive the training database's files on all cores: by zstd frame for the full reference, by gene for `reference.fna` | ~4 min | medium; same content, different compressed bytes |
| 2 | Setup as a graph: genome table ∥ converter; gene neighbours ∥ species clouds, beside the full reference's write; hold-out → derive → training build ∥ in-silico strains (background); the start-up checks concurrently | ~6-9 min (v19); more when NFS is slow | small-medium, low |
| 3 | Genome lengths and contig names in one pass (after the allele split, not before), cached per input folder by path, size and mtime; the reports' contig names from that cache or the genome store's headers | genome table 9-22 → ~1 min on a rebuild (~2-3 min less on a first build); contig names 7:33 → seconds always | small, low |
| 4 | Reports: start them right after training; unpack first (or keep the training database's `reference.fna`); one reference pass for all read types' ancestry; error-read threads by SAM bytes | ~12 min after "Ready" | small-medium, low |
| 5 | `--profile_ahead` on by default after one A/B on the cluster; the next batch chosen after the lock is taken, with streamed simulations merged into fewer runs | ~10-25 min | small, low (outputs identical by design; add a test) |
| 6 | Parity checks side by side; the collector's four tables in parallel | ~1-2 min | small |
| 7 | Two protal slots in the collection (two lock files, half the threads each, or one per NUMA node with `numactl`), when item 0 shows idle cores | 15-40 min if the short reads or tails leave cores idle; ~0 if not | medium; ~2 × 36 GB RSS beside the finished build's ~44 GB |
| 8 | Congener gaps (10:24 of the training build): already parallel; skip pairs whose sketch distance shows they are beyond 0.4 | perhaps 2-4 min | medium, changes no output if the cut is safe |

Items 1-6 are worth ~35-55 min of a 3:49 build at 84 cores, more on a rebuild from the same `--inputs` (item 3). Items 7-8 need item 0's numbers first. The
ONT screen is the largest single lever, but it is an alignment decision rather than parallelisation (section 4.2).

## 1. Where v19's time went

v19's `console.log` (offsets from the start at 21:18:27):

| Step | Starts | Wall | Parallel? |
|---|---|---|---|
| tool checks (`check_tools`: `--version`, git, Python imports from NFS) | 0:00 | 0:45 | serial |
| 1 genome table (`make_genome_table`) | 0:00:45 | 9:42 (v17 ~8:00, v18 ~21:45) | six recursive globs (serial), lengths in a process pool; NFS-bound |
| 2 converter | 0:10:28 | 4:04 | partly (section 2.2) |
| 2 `gene_neighbours.py` | | 3:18 (v17 6:28, v18 9:44) | NFS-bound (section 2.3) |
| 3 in-silico strains | 0:17:50 | 2:24 | feeds only the simulations |
| 4 species clouds, hold-out, scenarios | 0:20:14 | 0:12 | parallel, short |
| 5 derive the training database's files | 0:20:26 | 5:17 | **one process** |
| 5 training database build | 0:25:43 | 14:41 | parallel; "Congener gaps" 10:24 |
| 6+7 collections (simulations in the background: 17:40) | 0:40:24 | 2:46:18 | **one protal run at a time** |
| finished database, `SCHED_IDLE` in the background | 0:40:24 | 2:44:43 | ended 1:34 before the collection |
| 8 four trainers | 3:26:54 | 3:41 | 4 × 21 threads |
| 9 parity | 3:30:35 | 0:27 | read types one after another |
| 10 composition, 11 `--add_model` | 3:31:02 | 0:16 | |
| **Ready** | 3:31:18 | | |
| reports | 3:31:18 | 18:11 | section 5 |

In v18 the finished database's build ended 21 minutes after the collection (3 of them paused for the trainers). With
`--foreign-rates`, 16 more minutes followed. The idle-class build is on the path whenever the collection is fast. A
faster collection (items 5, 7) leaves more of it for after the collection, when it gets every core. The training
database, a similar build, took 14:41 alone.

## 2. Setup: a chain of steps, three of them reading `/hpc-home`

### 2.1 What each step needs

As written in `main()`, every step waits for the one before (`run()` blocks, `build_gtdb_database.py:534`). The data
dependencies are fewer:

- The genome table (lengths) is needed by `gene_neighbours.py`, the in-silico strains and the simulations. The
  converter does not need it.
- `gene_neighbours.py` needs `reference.fna/.map` (written by `gtdb_to_protal_db.py:970`). It does not need the full
  reference, whose write is the converter's longest stage (68 s + a 14 s join). This is P5 of the build audit, not
  done.
- The species clouds need the converted release. The hold-out needs the clouds (`--holdout-complex-distance`).
- The derive needs the hold-out and `gene_positions.tsv` (it derives the training database's `gene_neighbours.tsv`
  from the genomes kept).
- **The in-silico strains need the positions and the clouds; only the simulations read their output.** The hold-out
  reads `genomes.tsv`. Run them beside the derive and the training build: 2:24 off the path.

Proposed order (item 2):

1. At once: the start-up checks, the genome table (lengths) and the converter.
2. Once the converter has written `reference.fna/.map` (`:970`), and the genome table is done: `gene_neighbours.py`
   and the species clouds, side by side. Meanwhile the converter writes the full reference.
3. Hold-out (needs the clouds) → derive (needs the positions and the full reference) → training build → collection.
4. Beside step 3: the in-silico strains at idle priority, joined before the simulations start.

### 2.2 The converter (4:04) and the derive (5:17)

`convert.log` splits the 4:04 into stages. The derive logs no stages:

| Stage | Time | How it runs | Proposal |
|---|---|---|---|
| taxonomy + metadata | 6.6 s | serial (`gtdb_to_protal_db.py:815`) | read the metadata once: the species priors read it again (`:162-201`) |
| spool the representatives' genes | 52.7 s | `Pool(t)`, one task per marker (`:328-333`, `:850`). 173 tasks, of which 120 bac120 are heavy; the merge loop starts only after `pool.map` returns | `imap_unordered`, merge while the workers run |
| species priors | 12.6 s | serial (`:864-921`), the parent idle during the pools | in the parent during the spool |
| sort each gene's sequences | 32.5 s | parallel per marker (`:936`) | write sorted in the spool pass |
| write `reference.fna/.map` | 49.2 s | **one writer** (`:939-970`), 15 GB to OUTDIR on NFS, 14.5M map lines formatted in Python | in a thread beside the full-reference pool; or convert into SCRATCH |
| full reference's chunks | 68.4 s | parallel per marker, `zstd -T1` each (`:983`, `:398`) | 86 GB in 68 s ≈ 1.3 GB/s: probably reading the release from NFS |
| join | 14.4 s | **serial** copy to NFS (`:990-1001`) | append chunk i once chunks 0..i are done |
| **derive**: `reference.fna` filter | ~1-1.5 min (est.) | one text-mode loop over 14.5M records (`:697-710`), reading 15 GB from OUTDIR on NFS | split by gene (the map's offsets), workers, `pwrite` |
| **derive**: full reference | ~3-4 min (est.) | one `zstd -dc --long=31` stream (`:577`), a Python loop over ~80M records (159M lines), one `zstd -T8` per gene frame (`:511`, `:714-730`) | one worker per source frame: decompress, filter on the name's taxid, `zstd -T1`, then an ordered join with a new seek table. ~30-45 s |
| **derive**: `gene_neighbours.derive` | ~0.5-1 min (est.) | serial Python (`gene_neighbours.py:465-487`) | a process of its own beside the two filters |

Item 1 (the derive in parallel) is the safest large gain in setup: 5:17 → ~1 min. It sits directly before the
training build, the gate of the collection. The compressed bytes of `full_reference.fna.zst` would differ, but not
its content. protal reads its frames through the seek table, and `ancestry_sites.py` and `foreign_rates.py` stream
it whole.

A further step would make the derive unnecessary for the full reference: protal's build could skip the held-out taxa's
copies, as it already skips other genes' copies for `--build_gene_subset` (`src/Options.h:192`). But
`foreign_rates.py` and `ancestry_sites.py` read the training database's full reference too, and would need the same
filter. That moves the leak-proofing into three places, so I would not do it.

Converting into SCRATCH instead of OUTDIR would also help. OUTDIR's `reference.fna` (15 GB) is written once to NFS
and read about four times: gene neighbours, species clouds, the derive, and the finished build. Only
`database.protal` needs to end up in OUTDIR. The cost is ~22 GB more scratch, and a reconversion (2-3 min) on a new
node. Whether those reads really hit the network or the page cache is unknown.

### 2.3 Genome table, gene neighbours, in-silico strains

None of these steps logs its phases. The splits are estimates from the code and the sizes.

**The 45 s before the first line** are `check_tools` and `tool_versions` (`build_gtdb_database.py:2331-2332`), run
one after another:
- a cold `import joblib, numpy, pandas, sklearn` from the environment on NFS (`:1866`);
- each binary's `--version`, twice (`:1787`, `:1350`);
- `git diff` and `git status` on the checkout, also on NFS.

Run concurrently, with the `--version` output reused, they would take ~15-30 s less.

**Genome table (9:42; v17 ~8:00, v18 ~21:45).** The step's line is printed only after all its work.
- Serial parts, ~0.5-1 min together:
  - the taxonomy (`:544-557`);
  - six recursive globs over up to three roots (`:558-568`; one `os.scandir` walk would do, P4 of the build audit);
  - `read_representatives` (a gzip parse of ~730k metadata rows, `gtdb_to_protal_db.py:267-285`).
- Parallel: the lengths, a pool of 84 processes with chunksize 8 (`:590-592`). Each reads a whole genome, gunzips it
  in Python and counts the letters (`:614-622`). That is ~0.1 s of CPU per genome, about a minute on 84 cores. The
  rest of the ~8-9 min is reading ~49.6k files of ~1.2 MB from `/hpc-home`, ~100 MB/s.
- The swing from 8 to 22 min between builds is the network file system's load, not the build's.
- In the current code (`504de80` on, not in v19) the lengths are computed **before** `split_allele_genomes`
  (`:2437`). The split then rewrites `work/genomes.tsv` without the allele genomes. Their lengths, about a quarter
  of the files, are read for nothing. On a rerun into the same OUTDIR they are also missing from the table that serves
  as the cache (`:581-588`), so they are read again. Splitting first saves ~2-3 min.

**`gene_neighbours.py` (3:18; v17 6:28, v18 9:44).**
- Serial before: `read_map` parses all ~14.5M lines of `reference.map` into tuples, ~15-25 s
  (`gene_neighbours.py:146-154`, `:549`).
- Parallel: `Pool(84).imap`, one species per item (`:566-571`). Per genome: a gunzip and a Python k-mer dictionary.
  It reads the same genomes as the length pass, from the page cache as far as it still holds them. Its time follows
  the network file system too (3:18-9:44, slowest in v18, when the genome table was slowest as well).
- Serial after: `count_clades` and the summaries over ~4M gene ends (`:599-622`), ~1-1.5 min. The derive runs the same
  `count_clades` serially again for the training database (`gtdb_to_protal_db.py:676-679`).
- Proposals:
  - the species' ends and partners counted in the workers (each needs one species' genomes), only the clade merge in
    the parent;
  - the step started once the converter has written `reference.fna/.map` (P5);
  - the species clouds beside it, since they need only `reference.fna/.map` and the taxonomy.

**In-silico strains (2:24).**
- A pool of 84 processes over 9,030 strains, well balanced, ~1 s of CPU each. The cost is the per-ORF and per-frame
  Python loops and the gzip writes.
- `read_positions` reads the 5.2M position lines twice.
- Proposals:
  - vectorise the ORF scan and read the positions once, ~2:24 → ~1:30;
  - more simply, run the step in the background. Nothing reads its output before the simulations: the hold-out uses
    `genomes.tsv`; `pool_species(sim_table)` has the same species; `disk_needs` can count a strain at its
    representative's length.

**Every genome is read from `/hpc-home` up to five times per build**:
1. the lengths;
2. gene neighbours;
3. the in-silico strains' representatives;
4. the simulator's genome store on scratch (in the background, during the simulations);
5. the reports' contig names (section 5).

The page cache serves some of the repeats, until something pushes them out, but each build reads every genome over the
network at least once on the critical path. Caching what is derived from each genome in the input folder, keyed by
path, size and mtime, would spare these reads on every rebuild from the same `--inputs`:
- the length and the contig names (cheap: both come from the same decompressed bytes);
- the gene placements (`gene_neighbours.py`).

`download_gtdb.py` could write the first two as it checks each download.

## 3. The training database build: congener gaps

`training_db_index.log` (84 threads): Run build 14:21, of which:

| Stage | Time |
|---|---|
| gene conservation | 40 s |
| gene congeners, suspect copies, species neighbours | 40 s |
| **congener gaps** | **10:24** |
| index passes 1-2 | 15 s |
| uniqueness check | 33 s |
| unique k-mer statistics | 30 s |
| write the index, then `database.protal` | 0:42, then 1:13 in all |

Congener gaps (`Build.h:1367-1440`, `CongenerGaps.h:112-227`) go gene by gene (168). Within a gene, the alignments
are one `omp for schedule(dynamic, 64)` over every pair (85.6M pairs in all, 6.5M of them given up at 0.4). That part
uses every core. The serial parts per gene are short: the `holders` list, the per-genus maps, the sort of the pairs,
and the nearest/median pass. Two small `omp parallel` regions per sampled genus (75,507 groups) add fork/join
overhead, which is minor next to the WFA2 work.

So this stage is not serialised. It shrinks with more cores. It could get cheaper: a pair whose sketch distance shows
it is far beyond 0.4 could be given 0.4 without WFA2, which runs to its full budget on far pairs (item 8). The finished
database computes the gaps again (2:04:32 at idle priority, off the path).

From `504de80` on, the build adds **"Strain alleles"**: one pass over the full reference, aligning each allele
genome's copies. v19 did not have it. The next build's `training_db_index.log` will show its cost on the path.

## 4. The collection: one protal run at a time

### 4.1 The runs

The six runs end to end make 2:44:34 of the collection's 2:45:47 (`scripts/run_timings.sh`, times in s):

| Run | Samples | Aligning | Profiling stage | Index + preload + tables | Run |
|---|---|---|---|---|---|
| test 1 (7 streamed) | 28 | 1,472 | 164 | 8 | 1,783 |
| training 1 (6 streamed) | 44 | 2,583 | 111 | 77 | 2,751 |
| test 2 | 284 | 1,183 | 148 | 6 | 1,345 |
| training 2 | 1,062 | 1,874 | 248 | 80 | 2,183 |
| training 3 (7 streamed) | 54 | 1,343 | 209 | 80 | 1,616 |
| training 4 (1 streamed) | 10 | 132 | 49 | 6 | 192 |
| **sum** | 1,482 | **8,587 (143 min)** | **929 (15.5 min)** | ~4 min | **164.5 min** |

The simulations ended 17:40 into the collection. From then on, only one protal run was ever active. The lock is an exclusive `flock`
held around each run (`collect_training_data.py:1875-1887`), and each run gets the full `-t` (`:1487`). Each follower
picks its next batch **before** it waits for the lock (`:2146` vs `:2149`). So in v19, 165 GB of simulated training
reads sat on disk from ~0:17 to ~1:16, behind two streamed runs. In a serial chain that costs disk space, not time.

The index took 43 s to load in three runs and 2 s in the others. The 43 s runs presumably read `database.protal`
cold, after the background build and the reads had pushed it out of the page cache.

### 4.2 What the aligning time is

By read type, both collections (seconds):

| Read type | Samples | Aligning | Share |
|---|---|---|---|
| ONT | 298 | 6,001 | 69.9% |
| PacBio | 298 | 1,112 | 13.0% |
| paired-end | 457 | 902 | 10.5% |
| single-end | 429 | 571 | 6.7% |

**ONT is compute in WFA2, on candidates the screen cannot refuse.** Every deep ONT design point ran at 4,300-4,840
reads/s, whether all its reads came from bacteria or 90% of them from the human host:

| Sample | Reads | Candidates tried | Refused by the screen | Alignments made | Aligning |
|---|---|---|---|---|---|
| `sc_host_ont_b3000000000_s_1` | 622,132 | 17,926,550 | 1,191 | 39,711 | 2:06 |
| `sc_moderate_ont_b3000000000_s_1` | 191,586 | 5,692,086 | 95 | 100,125 | 0:40 |
| `sc_host_pb_b3000000000_s_1` | 332,066 | 17,846,585 | 17,211,101 | 42,467 | 0:17 |

That is ~29 candidates per read in both ONT samples. All of them go to WFA2, so the time per read is the same. The
screen (`AlignmentScreen.h`) is exact by the q-gram lemma. At `max_score_ani` 0.85 the score budget lets an alignment
touch every k-mer, so `Required` is ≤ 0 and nothing is refused. This is not a serial part, and it scales with the
cores. But it is 100 of the 143 min of aligning, and the host reads alone made 99.8% failing candidates. Two possible
fixes:
- an inexact pre-filter for ONT, such as an identity bound from the anchor chain;
- a cap on candidates per read segment.

Either would change outputs. That needs its own A/B on the ONT models, so it is not proposed here.

**Short reads probably leave cores idle, or wait on memory.** On 84 threads, paired-end ran at:
- 1.0M pairs/s for 150 bp (`rl150_p30000000`);
- 1.49M for 100 bp;
- 0.64M for 250 bp.

Those are ~690-715 MB of FASTQ per second for every read length. The 2026-10-06 cluster run aligned a real sample at
2.06M pairs/s on **32** threads, with similar counters per pair: 22.5 lookups per mate here against 22, and 3.8
candidates per pair against 3.5. 84 threads × 29.7 s for 30M pairs is ~83 µs of thread time per pair, against 15.6 µs
there. Flex cells per block are 115 here against 70, which explains part of it, not 5×.

The host paired-end sample settles the obvious cause. Its mates hit half as many flex cells per pair (2,123 against
4,208) and make a tenth of the alignments per pair, yet it ran only 40% faster than the community sample (1.40M
against 1.0M pairs/s). So the input path is not the whole limit either. That leaves two candidates:
- memory latency of 84 threads seeding a 25 GB index, perhaps across NUMA domains;
- cores idle at the reader lock (`SeqReader.h`: 32 records per acquisition, one lock for both files).

Short reads are 17% of the aligning (25 min), so this matters less than ONT. But it is where a second protal run
beside the first (item 7) could gain, and where 120 cores would help least.

### 4.3 Profiling stages, start-ups, streamed runs

- **Profiling stage**: 15.5 min over six runs, after each run's last alignment. The samples are profiled in parallel,
  largest SAM first, with threads in proportion to SAM bytes (`RunProtal.h:1560-1586`). The tail is one deep
  sample's SAM read: up to 82-95 s with a few threads.
- **`--profile_ahead`** (`RunProtal.h:1450-1535`, opt-in, `--profile-ahead` in the build): a worker with `-t/4`
  threads profiles each finished sample beside the alignment. The profiling stage then does only what the worker had
  not started, and the last sample. On a 6-core laptop it gained nothing, because the alignment already had every
  core. On the cluster, the worker's threads would come out of ONT's WFA2 (compute) but fill the short reads' idle
  cores. The profiling stage would shrink to roughly the last sample's 1-2 min per run. Estimate: 6-15 min, minus
  some alignment. One A/B on a streamed ONT+pe run decides it. Profiles are the same by design (a profile depends on
  the SAM and the database only); add a test that compares them.
- **Start-ups**: ~80 s for each cold run (index 43 s, tables 19-22 s, preload 15 s, partly concurrent). The streamed
  batches are capped like written ones, so training run 1, 3 and 4 were three runs. Choosing the batch after taking
  the lock, and dropping the cap for streamed batches (they write nothing), would make them one or two runs. That
  saves a start-up and a profiling stage each, ~4-10 min.
- **Streaming itself is not a limit**: streamed design points aligned at the same rates as written ones (ONT 4,776
  against 4,755 reads/s, PacBio 15,701 against 16,495).
- **Shallow samples are cheap**: training run 2's 1,062 samples spent 1,884 s in "Processing all samples" against
  1,874 s of their own alignment timers. Per-sample overhead outside the timer is ~10 ms.
- **After the last run**: the collector writes the four tables one after another in one process
  (`collect_training_data.py:2289-2290`). The estimate is 1-3 min at r226, not measured. Writing them in parallel, or
  per unit as runs end, saves most of it.

### 4.4 Two protal runs at once (item 7)

Nothing in the code needs one run at a time except the lock. Two slots (two lock files, `-t T/2` each) would let the
two followers run side by side, or a streamed long-read run beside a written paired-end block. The constraints:
- **Memory**: ~34-38 GB RSS per run at r226, plus the finished build's ~44 GB at idle priority, plus the collector.
  That fits on the 502 GB node, but no code checks memory before a run.
- **Cores**: each run at `T/2`. If ONT's WFA2 is the run's work, two halves take as long as one whole. The gain comes
  only from what one run leaves idle: profiling-stage tails, start-ups, short reads' idle threads, the deep sample's
  SAM read.
- **NUMA**: on a two-socket node, two runs pinned with `numactl --cpunodebind=N --membind=N` would each seed their own
  copy of the index from local memory. If the short-read slowdown is NUMA latency, this is where it would show.

Estimate: 15-40 min if 20-30% of the collection's core time is idle, near nothing if it is not. Item 0 says which.
A cheap test that does not need the build: `scripts/measure_performance.sh` with `THREADS=32,64,84,120` on one
simulated r226 paired-end sample and one ONT sample, once plain and once with `numactl` per socket.

## 5. After "Ready": the reports (18 min)

`reports()` (`build_gtdb_database.py:1008-1110`):

| Part | v19 | Why | Proposal |
|---|---|---|---|
| contig names of 56,684 genomes | 7:33, before any report starts | `trace_relatives.genome_contigs` (`trace_relatives.py:123-154`): a pool of 84 processes gunzips every genome from `/hpc-home` in Python for its header lines. ~25 ms CPU each (~20 s on 84 cores), so the 453 s are NFS reads. The cache is `SCRATCH/genome_contigs.tsv.gz`, gone with the job | collect the names in the genome-length pass of step 1 (the same files, already decompressed there; `write_genome_table`) into that cache; or read them from the genome store's `.g2b` headers (`GenomeStore.cpp:238-253`, which hold the names and the source's size and mtime). Seconds |
| `error_reads` × 4 + `trace_relatives` | 3:09 | side by side, 84 // 5 = 16 threads each (`:1053`). pe was the long pole (186 s); pb/ont end early. `trace_relatives`' SAM scan is one serial loop (`trace_relatives.py:196-253`) | threads and memory by SAM bytes |
| ancestry: unpack | 1:59 | `protal --unpack_db` runs only after every error-read report has ended (`:1107`, `:1121-1136`), though it needs none of them | start it first; or with `--share-logs` keep the training build's `reference.fna` |
| ancestry: per read type | ~5.5 min | four single-threaded processes side by side. Each makes the same two passes: `reference.fna` (15 GB, `ancestry_sites.py:231-251`) and `full_reference.fna.zst` (86 GB raw, through `zstd -dcq`, `:143-150`) | one pass for the four read types' wanted copies, written to scratch, then each read type counts its SAMs. ~3-5 min |

The reports need only the trainers' `calls.tsv.gz`, the SAMs and the training database. They could start right after
training, beside parity, composition and the wait for the finished build (~0.5-1 min more). They must not start
beside the trainers (`build_gtdb_database.py:3185`).

Together: ~18 → ~5-6 min.

## 6. Item 0: what the build should record

Nothing in v19's logs says how many cores a step used. Four small additions would make the next build answer items 7
and 8 directly:

1. **Per step**: `Job.poll()` already reaps with `os.wait4`. Its `usage.ru_utime + usage.ru_stime` is the CPU time of
   the command and every descendant it waited for. A line such as "built training_db in 0:14:41 on N cores on
   average, peak memory 40.7 GB".
2. **Per protal run** in the collector: the same from `wait4` on protal: "protal run 2: 36:24, N cores on average".
3. **Per sample** in protal: a `getrusage(RUSAGE_SELF)` difference around each sample's alignment and the profiling
   stage: "Aligning reads took 29.7 s (N of 84 cores busy on average)". This settles the short-read question.
4. **A timeline**: a thread in the build script that reads the job's cgroup `cpu.stat` (`usage_usec`; the path from
   `/proc/self/cgroup`, falling back to `/proc/stat`) every 15 s. It writes `logs/cpu_timeline.tsv` (time, cores busy,
   the labels of the running jobs), which the share archive carries.

For v19 now: `seff 24099010` (or `sacct -j 24099010 --format=Elapsed,TotalCPU,AllocCPUS`) gives the whole job's
average. Its TotalCPU / (3:49 × 84) is the job's efficiency.

## 7. How the steps scale from 84 to 120 cores

| | v19 at 84 | at 120, nothing changed (estimate) | with items 1-6 at 120 |
|---|---|---|---|
| setup | 40 min | ~35 (NFS, derive, chain unchanged; training build 14:41 → ~11) | ~20-25 (first build), ~12-15 (rebuild with the caches) |
| collection | 2:46 | ~2:05-2:15 (ONT and PacBio × 0.7; short reads, profiling stages, start-ups less) | ~1:50-2:00 |
| training to ready | 4.5 min | ~4 | ~3.5 |
| reports | 18 min | ~17 | ~5 |
| **total** | **3:49** | **~3:00-3:10** | **~2:20-2:40** |

The ONT screen (section 4.2) is the largest lever outside these. Halving the ONT candidates that reach WFA2 would take
~35-50 min off the collection at 84 cores.

## 8. Not affected

- The website: nothing here changes how protal is run. If items 1-5 are implemented, `docs/databases.md` (the build's
  stages and order, `--profile-ahead`) needs updating.
- Simulations: 17:40, in the background, off the critical path since `73607bb`.
- The trainers: 3:41 for four, 86-140 s of timed work each (the rest untimed, presumably start-up and imports from
  NFS).

## 9. Implemented (follow-up, 2026-10-09)

The same day the user asked for items 1-6, then item 0 written to log files only ("not reported to stdout") and carried
by the `--share-logs` archive, and then option 1 of the [ONT report](../2026-10-09-ont-wfa2-skipping.md). Nothing of it
has run at r226 yet.

| # | What was done | Where |
|---|---|---|
| 0 | **CPU logs, files only.** Every command the build runs: a row in `logs/cpu_jobs.tsv` when it ends (start, wall and CPU seconds from `os.wait4`, cores busy, nice or SCHED_IDLE, peak GB, exit, label, log). Every 15 s: `logs/cpu_timeline.tsv` (cores busy from the job's cgroup v2 `cpu.stat`, v1 `cpuacct.usage`, or `/proc/stat` for the node; iowait; the cores allowed; the build's step; the commands running). protal: `misc/cpu.tsv`, a row per stage (`start-up`, `aligning <sample>`, `profiling`, `rest of the run`; `src/Utilities/CpuLog.h`, `getrusage`). The collector: `protal_cpu.tsv`, per protal run the whole process (`wait4`) and protal's rows, copied to `logs/protal_cpu_<collection>.tsv`. All in the share archive | `build_gtdb_database.py` (`Job`, `CpuTimeline`), `RunProtal.h`, `collect_training_data.py` (`run_protal`, `record_cpu`) |
| 1 | **The derive on all cores.** `reference.fna` in pieces of whole records, each derived by a worker into its own file and map with local offsets, then placed at its offset; `full_reference.fna.zst` frame by frame through its seek table (one `zstd -dc`, filter, `zstd -6 --long=27` per frame, the converter's per-marker frames kept, so a species' copies of a gene stay in one frame and the build stays deterministic); the gene neighbours in a process beside them. A reference without a seek table takes the old serial path | `mini_db/gtdb_to_protal_db.py` (`_derive_reference_piece`, `_place_reference_piece`, `_derive_full_frame`, `read_seek_table`) |
| 2 | **Setup as a graph.** The converter starts in the background while the genome table is made (only when no finished database could be kept: a conversion into `protal_db` deletes `database.protal`); `gene_neighbours.py` starts as soon as the converter has written `reference.fna`/`.map` (its own stage key: an existing OUTDIR converts once more on its next run); the in-silico strains run in the background, joined before the converted files are removed and before the streaming decision (the genome store's room is planned from the table before them, `planned_bases`); the start-up checks side by side | `build_gtdb_database.py` `main` |
| 3 | **Lengths and contig names in one pass**, after the allele split, cached per genome by path, size and mtime (`--genome-cache`, by default in the `--inputs` folder); the names also fill the reports' contig cache, so the reports' 7:33 become a lookup | `genome_facts`, `GenomeCache`, `trace_relatives.add_to_contig_cache` |
| 4 | **Reports right after training**, beside parity, packaging and the wait for the finished database; every report on all threads; the training database's `reference.fna` kept by a hard link for the ancestry report (no unpack; its ~15 GB in the scratch estimate); one `ancestry_sites.py` run for every read type, which reads both references once and analyses the read types side by side, each read type's outputs byte-identical to a run of its own (gzip headers now without a time stamp); one log, `logs/ancestry_sites.log` | `reports`, `ancestry_reports`, `ancestry_sites.py` |
| 5 | **`--profile-ahead` on by default** (`--no-profile-ahead`); the collector chooses a run's batch after it has taken the protal lock, never holding the lock while it waits for simulations; the streamed simulations whose communities are there go into one run, uncapped (their reads are on no disk, and protal reads one sample at a time, so the simulators it has not reached wait at their pipes) | `collect_training_data.py` (`follow_state`, `next_batch`, `follow`) |
| 6 | **Parity checks side by side** (one per read type, `--profile_only`, logs concatenated in order); **the collector's tables in parallel** (forked workers, their output printed in read-type order, failures re-raised in the serial order) | `build_gtdb_database.py`, `collect_training_data.py` (`write_tables`) |
| ONT | **The exact indel-distance bound** before WFA2 for ONT reads, in place of the k-mer screen (option 1 of the ONT report): its own counter in the per-sample counts line and an "Indel bound" timer; `--no_indel_bound` gives the k-mer screen back | `AlignmentScreen.h` (`IndelBound`), `AlignmentStrategy.h`, [ONT report](../2026-10-09-ont-wfa2-skipping.md#7-implemented-follow-up-2026-10-09) |

Tests (WSL on 4 cores, shared with another session's work, the mini database, protal built at the change):

| Suite | Result | New or changed |
|---|---|---|
| protal unit tests | 479 passed, 3 skipped | `IndelBound.*` (3) |
| end to end (`tests/e2e`) | 149 OK, 136 s | `test_the_run_cpu_table`, `test_profiles_ahead_are_the_same` (with `--profile_ahead` the profiles, sorted SAM records and strain MSAs are the baseline's), `test_the_indel_bound_changes_no_alignment` |
| `test_gtdb_build.py`, `test_mini_db.py`, `test_error_reads.py`, `test_trace_relatives.py` | 61 OK | one pass and the cache of genome facts, the folders walked once, the bases planned before the in-silico strains, the CPU logs; the derive at `-t 3` byte for byte the derive at `-t 1` (files, the full reference's content and frames). `test_cpu_logs` first missed a CPU threshold on the loaded machine; its busy loop now counts CPU time, not wall time |
| `test_collector.py`, `test_ancestry_sites.py` | 18 and 8 OK | the batch chosen after the lock, streamed runs merged, the tables in parallel (byte-identical to serial), `protal_cpu.tsv`; several read types in one ancestry run, each byte-identical to its own run |
| `test_gtdb_pipeline.py` (whole mini GTDB builds) | 6 OK, 304 s | one assertion changed: the reports' header now comes after "trained in", not after "Ready" (item 4), and the ancestry report writes one log |

What the next r226 build should show: `logs/cpu_timeline.tsv` over the collection (idle cores in the paired-end
runs settle item 7), `logs/cpu_jobs.tsv` for the derive and the training build (item 8), the setup's wall time against
v19's 40 min, the reports against 18 min, and the ONT rows of `logs/protal_cpu_*.tsv` against v19's aligning.

## Reproducing

```bash
bash docs/claude/2026-10-09-build-parallelism/scripts/run_timings.sh local/v19/protal0.7.9_r226_v19/logs
```

`local/v19` is the unpacked `protal0.7.9_r226_v19_share.tar.gz` (not in git). The step times come from its
`console.log`, and the build's stage timers from `logs/training_db_index.log`, `logs/convert.log` and
`logs/ancestry_unpack.log`.
