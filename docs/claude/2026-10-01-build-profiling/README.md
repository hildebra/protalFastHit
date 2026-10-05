# Where the database build spends its time

2026-10-01, branch `audit-fixes` at `25d457e` (protal, `simulate_metagenomes` and the scripts built and
run from `git archive 25d457e`; the trainer of `a32b60f` for the scaling runs, see below). Question: where
does `scripts/build_gtdb_database.py` lose most of its time, and does protal profile the simulated
metagenomes one run per sample, or several samples per run?

Machine: the WSL Ubuntu 24.04 laptop (Core Ultra 7 258V, 6 vCPUs, 23 GB, on mains power), shared with
other sessions' jobs the whole time (load 9-16): wall-clock times are inflated and noisy, CPU times less
so. pbsim3 3.0.5, ART 2.5.8 (2016-06-05), scikit-learn 1.9.1.

Data: the GTDB-like tuning world (765 species, `~/tune/release_p`; its genome table of 2,295 genomes,
765 representatives and 1,530 other strains, `~/tune/V3/genomes.tsv`). Its genomes have ~270 kb in 3
contigs; GTDB's have ~3.4 Mb in 1 to hundreds of contigs, and the simulators cost per base and per
contig. So the genomes to simulate from were padded to real sizes (`scripts/pad_genomes.sh`): each keeps
its own contigs (and marker genes) and gains random contigs up to a size drawn lognormal (median
3.4 Mb, 1.5-8 Mb) in a number of contigs drawn lognormal (median 40; mean 78.6). The databases are the
tuning world's (`~/tune/V3/training_db`, 83 MB), so everything protal does is at the tuning world's
size, ~1/190 of GTDB r226. Work data: WSL `~/bprof`.

## Summary

**protal already profiles all samples of a collection in one run.** `collect_training_data.py` writes
one map of every sample of every read type (`profile_all/samples.map`) and runs protal once on it, so
the database is loaded once per collection: twice per build (training data, test set), plus four
`--profile_only` parity runs and four `--add_model` rewrites. Inside the run, samples are aligned one
after the other with all threads each (nothing database-sized is set up per sample) and profiled in
parallel, the largest first. At r226 a load takes 1-1.5 min from the page cache
([load-and-output](../2026-09-30-load-and-output/README.md)), so one run per sample would cost 664 loads
(about 14 h); merging the two collections' runs would save one load. protal is not where the time goes.

**Most of the time goes to simulating the long reads.** pbsim3 recomputes its read-length table (100 bp
to 1 Mb, three `pow`/`exp` per length) for every sequence of a genome, and the collector calls it once
per genome of every sample: 45 ms per contig for Nanopore, 8 ms for PacBio, whether the contig gets a read
or not. With 78 contigs per genome that is 3.5 CPU-s per Nanopore genome. In the whole build run
here at a reduced scale (one sample per design point, 1:48 h), simulating was 92% of the time, the
long reads 71%, and pbsim3 69% of the CPU. At r226's design (30,000 genome calls) it is ~18 of the
~32 CPU hours estimated for the collections and the gene neighbours (below).

**The shallow long-read samples are not what the design asks for.** pbsim3 makes reads for a sequence
until their bases reach depth × its length and cuts the last read to what is left (at least 100 bp). At a
shallow point every contig's quota is a few bases, so every contig gets one ~100 bp read: the 300 kb
Nanopore sample of real-sized genomes has **4,471 reads with a median of 100 bp** (and 2.0-4.7 times the
bases asked for) instead of ~40 reads of ~8 kb, and even the deepest point's (150 Mb) reads have
medians of 1.4-2.5 kb (Nanopore) and 1.7-4.9 kb (PacBio, 15 kb asked). The tuning world's own builds
have it too, less (3 contigs per genome): its Nanopore 300 kb samples' reads have a median of 217 bp,
its 6 Mb samples' 1.8 kb, the 60 Mb samples' 5.9 kb. The pb and ont models, and their depth knobs,
were trained on these.

**One pbsim3 call per sample fixes both.** pbsim3's `--strategy templ` turns each template of a FASTA
into one read of its length with the model's errors. Drawing the reads in the collector (genome by
abundance × length, length from the setup's gamma, start uniform, cut at the contig's end) and calling
pbsim3 once per sample took 2-3.5 s (the templates 1.8-3 s in Python, pbsim3 0.5 s) instead of ~450
CPU-s for a shallow Nanopore sample, and gives 41 reads with a median of 4.8 kb
(`scripts/long_read_templates.py`, a prototype). Deep samples cost about the same per base either way
(0.44-0.75 s per Mb).

Where the time is lost, ranked by what it costs at r226 (details and the fixes below; none
implemented):

| # | Where | Cost at r226 (estimated) | Fix | Saves |
|---|---|---|---|---|
| 1 | pbsim3 once per genome, its length table per contig | ~18 CPU h; ~45-55 min wall at 32 threads, ~2.3 h at 8 | one `--strategy templ` call per sample | ~17 CPU h; and the reads as designed |
| 2 | collection waits for the training database's build | the simulations' wall time (~1-1.5 h at 32 threads now, ~15-20 min with 1 and 3) after a 45-80 min build | simulate during the builds (the simulator needs no database) | the simulations' wall time, up to the build's |
| 3 | `simulate_metagenomes` reads every genome of the table for its length, every design point | 33 × 20,000 genomes × 21 ms = 3.8 CPU h; ~7 min before each point's first sample | lengths in the genome table (it takes a fourth column) | 3.8 CPU h, ~7 min per collection |
| 4 | trainer: the previous procedure's grid search (`--evaluation full`) | 82-83% of training; trainer time grows as rows^1.3, so ~2.5 CPU h per paired-end model at ~70k rows, ~7 CPU h for the four (~20 min wall at 32 threads, ~1.2 h at 8) | leave it out after the first r226 comparison, or 3 folds and fewer values | most of the training step |
| 5 | `gene_neighbours.py`: up to three whole-genome `find`s per gene | 1.45 CPU-s per representative, 3.2 CPU h, before the builds can start | an index of each genome's k-mers at a stride | ~3 CPU h, 5-25 min wall |
| 6 | the collector's Python copies each genome for pbsim3 under the GIL | 90 ms per genome, ~45 min serial | gone with #1 | |
| 7 | four `--add_model` runs, each rewriting the 21 GB `database.protal` | 4 rewrites | one run with all models | 3 rewrites |

Follow-up (2026-10-02, [below](#follow-up-2026-10-02-fixes-1-2-3-and-5-implemented)): 1, 2, 3 and 5
implemented. The same whole build took ~15:30 instead of 1:48:20 (pbsim3 643 instead of 10,070
CPU-s), the long-read samples have the asked lengths and bases, and the genome lengths and gene
neighbours leave the outputs byte-identical. Estimated at r226: ~8 instead of ~32 CPU hours for the
collections and gene neighbours, and the simulations off the critical path.

## What was run

- `scripts/collect_real_size.sh` (`c1`): `collect_training_data.py` with the build's training design
  (20-200 species, strains 0.3,0.1, two archaea, one species of a held-out clade per rank) on the padded
  genomes, two samples at 1,000 and 500,000 pairs (150 bp), 300 kb and 150 Mb of PacBio and Nanopore
  reads; ART and pbsim3 behind `scripts/timed_tool.sh`, which logs every call's wall and CPU time.
- `scripts/build_real_size.sh` (`b1`): the whole `build_gtdb_database.py` on the padded genomes with the
  default design but one sample per design point (training and test).
- `scripts/length_scan.sh`, `scripts/pbsim_fixed_cost.sh`, `scripts/templ_prototype.sh` (with
  `long_read_templates.py`), `scripts/gene_neighbours_cost.sh`, `scripts/trainer_scaling.sh` (with
  `enlarge_table.sh`): the single measurements below; `scripts/tool_summary.sh`,
  `scripts/protal_log_summary.sh`, `scripts/build_summary.sh`, `scripts/long_read_bases.sh`,
  `scripts/genomes_per_sample.sh`: the summaries; `scripts/r226_estimate.py`: the extrapolation.
  `scripts/make_wrappers.sh` makes the timing stand-ins.

## The whole build, measured (`b1`)

`build_gtdb_database.py` of `25d457e` with 4 threads on the tuning world's release, the padded genomes
to simulate from and the default design with one sample per design point: 15 paired-end training
points (15 pe, 15 se samples), 5 + 5 long-read training points, 18 + 5 + 5 test points. 1:48:20 wall,
14,675 CPU-s (13,047 user, 1,629 system), at most 4.5 GB.

| Step | Wall | Share | What its time is |
|---|---:|---:|---|
| convert the release | 0:00:03 | 0% | |
| gene neighbours (765 representatives) | 0:04:05 | 4% | ~1,100 CPU-s |
| the training database's files, its build | 0:00:46 | 1% | the finished database built meanwhile (6:40) |
| training data: paired-end simulation | 0:11:15 | 10% | 15 points, 4 at a time, each ~2 min: a ~48 s length scan, then one sample (ART: 2,504 calls, 330 CPU-s) |
| training data: long reads | 0:36:58 | 34% | 10 points one after another; pbsim3 1,746 calls, 4,400 CPU-s (Nanopore 3,592) |
| training data: protal (40 samples, one run) | 0:00:38 | 1% | index 1.1 s, aligning 34.5 s (pe 8.3, se 4.3, PacBio 8.5, Nanopore 13.3), profiling 1.0 s |
| test set: paired-end simulation | 0:11:51 | 11% | 18 points (ART: 4,698 calls, 556 CPU-s) |
| test set: long reads | 0:39:58 | 37% | pbsim3 2,802 calls, 5,670 CPU-s (Nanopore 4,479) |
| test set: protal (46 samples, one run) | 0:00:49 | 1% | |
| training the four models | 0:01:41 | 2% | the previous procedure 37-77 s of each model's 49-99 s |
| parity, `--add_model` | 0:00:11 | 0% | |

So 92% of the run is simulation, 71% long reads, and pbsim3 alone is 69% of the CPU (10,070 of
14,675 s); the 33 length scans are ~1,600 CPU-s (11%) and ART 886 CPU-s (6%). With one sample per
point the length scan is a large part of each paired-end point; with the default 12 (training) or 4
(test) it is less, but it grows with the genome table (2,295 here, ~20,000 at r226). The databases,
protal and the trainer are minor at this scale; at r226 the builds (45-80 min each) and the trainer
grow, the simulations do not (their cost is per genome of a sample, and per base).

## Unit costs

Each measured on the padded genomes (mean 3.58 Mb, 78.6 contigs; 1.13 MB gzipped), CPU seconds unless
said otherwise.

| What | Cost | From |
|---|---|---|
| `simulate_metagenomes`: a genome read for its length (`build_length_cache`, every genome of the table, every run) | 20.9 ms (2,295 genomes: 47.9 s; with a length column: 0.02 s) | `length_scan.sh` |
| ART per call (a genome of a sample) | 0.08-0.10 s fixed + 32 s per million 150 bp pairs | `c1` |
| the simulator's plain copy of the genome for ART | ~20 ms (as the length scan) | |
| pbsim3 PacBio (`errhmm`) | 8 ms per contig (3, 40, 300 contigs: 0.06-0.11, 0.32-0.41, 2.3-2.7 s) + 0.42-0.75 s per Mb | `pbsim_fixed_cost.sh`, `c1` |
| pbsim3 Nanopore (`qshmm`) | 45-60 ms per contig (0.19-0.26, 1.7-2.5, 13-21 s) + 0.44 s per Mb; 3.5 s per genome of 78 contigs | same |
| the collector's `long_contigs` copy of a genome for pbsim3 (Python, holds the GIL) | 90 ms | |
| `gene_neighbours.py` per representative genome | 1.45 s (100 genomes, one thread: 145 s) | `gene_neighbours_cost.sh` |
| protal, training database of the tuning world, 4 threads | pe 150k pairs/s, se 380k reads/s, PacBio 11 Mb/s, Nanopore 7 Mb/s; index load 2.9 s; profiling 16 samples 1.0 s | `c1` |

Where pbsim3's time per contig goes (callgrind, Nanopore, 40 contigs, depth 0.001: 15.2 G
instructions): 79% in `pow` and `exp` (libm), 12% in `simulate_by_qshmm`. Its source
(`src/pbsim.cpp`, `main`) calls `simulate_by_qshmm` once per sequence of the genome file, and that
fills the read-length table for every length from `--length-min` (100) to `--length-max`
(1,000,000) with a gamma density (`pow(i, kappa-1) * exp(-i/theta) / pow(theta, kappa) / gamma`)
before it makes a read. Reads are made `while (len_total < sim.len_quota)`, `len_quota = depth ×
sequence length`, and a read that would pass the quota is cut to it (`mut.len = sim.len_quota -
len_total`, at least `len_min`).

What that does to the long-read samples (`long_read_bases.sh`):

| Sample | Asked | Got | Reads | Median read |
|---|---:|---:|---:|---:|
| PacBio 300 kb, real-sized genomes (`c1`, 2 samples) | 300 kb | 607 kb, 1.40 Mb | 4,318, 13,293 | 100, 100 bp |
| Nanopore 300 kb, real-sized genomes | 300 kb | 615 kb, 1.40 Mb | 4,471, 13,396 | 100, 100 bp |
| PacBio 150 Mb, real-sized genomes | 150 Mb | 150.7, 150.5 Mb | 30,983, 21,199 | 1,720, 4,933 bp |
| Nanopore 150 Mb, real-sized genomes | 150 Mb | 150 Mb | 48,715, 35,892 | 1,377, 2,521 bp |
| PacBio, tuning world (V3, 3 contigs per genome; sample 1 of each point) | 300 kb, 6 Mb, 60 Mb | 306 kb, 6.03 Mb, 60.0 Mb | 241, 1,209, 4,428 | 340, 2,577, 14,191 bp |
| Nanopore, tuning world | 300 kb, 6 Mb, 60 Mb | 310 kb, 6.04 Mb, 60.0 Mb | 314, 1,892, 8,275 | 217, 1,805, 5,892 bp |
| PacBio, one `templ` call per sample (prototype, the `c1` communities) | 300 kb, 150 Mb | 271 kb, 124 Mb | 22, 10,018 | 13,436, 13,260 bp |
| Nanopore, one `templ` call per sample | 300 kb | 264 kb | 41 | 4,908 bp |

The setups ask for PacBio reads of 15 kb (SD 3 kb) and Nanopore reads of 8 kb (SD 6 kb). With
genomes of GTDB's size and fragmentation even the deepest point's reads are a fifth to a third of
that; the shallow points are reads of 100 bp, one per contig of every genome of the sample, so that
every species of a shallow sample has reads (the pb_b300000 point of `c1` has 74 present taxa in 2
samples, the paired-end point of the same communities and bases, 1,000 pairs of 150 bp, 20). The
tuning world's builds (V3 above) had the same at the shallow points, less at the deep ones: pb and
ont models, and depth knobs, trained on such builds learnt shallow long-read samples unlike real
ones. That fits what [features-depth-knobs](../2026-10-01-features-depth-knobs/README.md) found: the
PacBio depth knob of a shallow bin was fitted on samples with few absent taxa, and did not carry
over. (The prototype's 150 Mb sample came out short of its quota as it counted the reads'
lengths before the cut at the contig ends; the script now counts the templates' bases.)

## At GTDB r226 with the defaults (estimated)

`scripts/r226_estimate.py`: the unit costs above times the default design's calls and bases. The
training data are 15 paired-end design points of 12 samples (180 samples, 22.5M pairs) and 5 + 5
long-read points of 12 (2.25 Gb each type); the test set 18 points of 4 (72 samples, 15.2M pairs) and
5 + 5 of 4 (1.12 Gb each). Measured in `c1`, `V3` and `b1`: 170 genomes per training sample, 240 per
test sample, so 47,900 ART and 30,000 pbsim3 calls. Assumed: a genome table of 20,000 genomes
(download_gtdb.py's defaults allow 6,000 species × 3 + 2,000; fewer if NCBI has fewer strains), 8,000
representatives, 78 contigs per genome as in the padded genomes (`genomes.tsv` of download_gtdb.py
has each strain's contig count; GTDB's representatives are often MAGs with more), and protal 3×
slower per read than on the tuning world's database (unknown; the r226 run's log has "Aligning reads
took" per sample).

| stage | part | now, CPU h | proposed, CPU h |
|---|---|---:|---:|
| training | simulator: genome lengths | 1.7 | 0.0 |
| training | ART and its genome copies | 1.3 | 1.3 |
| training | pbsim3 PacBio | 2.0 | 0.3 |
| training | pbsim3 Nanopore | 10.2 | 0.3 |
| training | collector: genome copies for pbsim3 (GIL) | 0.5 | 0.2 |
| training | protal alignment (×3 at GTDB size) | 2.5 | 2.5 |
| test | simulator: genome lengths | 2.1 | 0.0 |
| test | ART and its genome copies | 0.7 | 0.7 |
| test | pbsim3 PacBio | 1.0 | 0.1 |
| test | pbsim3 Nanopore | 4.8 | 0.1 |
| test | collector: genome copies for pbsim3 (GIL) | 0.2 | 0.1 |
| test | protal alignment (×3 at GTDB size) | 1.3 | 1.3 |
| before the builds | gene_neighbours.py | 3.2 | 0.3 |
| all | | **31.7** | **7.2** |

Not in the table: the two index builds (45-80 min each at 16 cores, estimated in
[index-build-gains](../2026-09-30-index-build-gains/README.md); the r226 logs will tell), the
conversion (not measured at r226), the trainer (below) and protal's start-up (1-1.5 min per run from
the page cache, more from network storage).

How the wall time adds up: every step waits for the one before. The simulations need neither
database, but start after the training database's build; the test set's start after the training
data's; within a collection the paired-end points run in parallel (`--threads` at a time) with their
samples one after another, the long-read points one at a time with their genomes in parallel, and the
collector's genome copies for pbsim3 hold the GIL (~45 min of the two collections at r226, whatever
the threads). Rough wall times from the CPU hours, the chains and the GIL:

| step | now, 32 threads | proposed, 32 threads | now, 8 threads |
|---|---|---|---|
| gene neighbours | ~6 min | <1 min | ~25 min |
| training database build | 45-80 min | 45-80 min, the simulations meanwhile | 45-80+ min |
| training data: paired-end simulation | ~18 min (the 500k-pair points' 12 samples one after another, after a 7 min length scan) | in the build | ~25 min |
| training data: long reads | ~30-40 min (GIL ~30 min) | ~3 min, in the build | ~95 min |
| training data: protal | ~7 min | one run with the test set, ~10 min | ~20 min |
| test set: paired-end, long reads, protal | ~13 + ~15 + ~4 min | in the build; in the run above | ~20 + ~45 + ~10 min |
| trainer, parity, `--add_model` | below; parity 4 runs, 4 rewrites | | |

That is about 1.5 hours of the run at 32 threads that the proposals remove (and the long reads come
out as designed), and about 3 hours at the default 8 threads (where the simulations compete with the
build for the cores and hide less well), before anything in the builds or the trainer.

## The trainer

`scripts/trainer_scaling.sh`: the V3 build's paired-end table (7,366 rows, 36 samples) and the same
three times over (`enlarge_table.sh`: each copy's samples renamed, its features ×(1 ± 2%) noise),
`--evaluation full`, one thread, the trainer of `a32b60f` (later ones want the features of `e680d91`,
which V3's table lacks):

| Table | fit | evaluation | feature sets | previous procedure | capacity | learning curve | total | CPU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 7,366 rows (V3's own run, light load) | 1.1 s | 31 s | 24 s | 510 s | 39 s | 7 s | 615 s | ~600 s |
| 22,098 rows (load 10-19) | 8.9 s | 259 s | 142 s | 2,739 s | 142 s | 25 s | 3,321 s | 2,471 s |

The previous procedure (a grid search over 20 values of `max_features` × 5 folds with 128 trees,
then a forest of 512 trees on the chosen features) is 82-83% of the time either way, and the
whole grows ~4.1× for 3× the rows (rows^1.3). How many rows r226 gives is not known: the tuning
world's database gave 116-205 per paired-end sample; GTDB's, with 190 times the species, puts reads
of a sample on more relatives. At 400 rows per sample (70k for 180 samples) a paired-end model would
take ~3 CPU h, 2.5 of them the previous procedure; with the long-read models' smaller tables (60
samples) the four come to ~7 CPU h, run in parallel with a quarter of `--threads` each.

## protal in the build: one run per collection

`collect_training_data.profile` joins the rows of every unprofiled design point of every read type
into one map (`OUT/profile_all/samples.map`, with a `READ_TYPE` column) and runs
`protal --db TRAINING_DB --map ... --no_strains --no_qcmsa` once; `check_model_parity.py` re-profiles
two points' SAMs per read type with `--profile_only` (one run per read type), and the build ends with
one `--add_model` per read type. In the run (`RunProtal.h`): the index is loaded once
(`Load Index`, 2.9 s here, 1-1.5 min at r226); `RunWrapper` aligns the samples one after another,
each on all threads, and sets up only per-read buffers for a sample (`ChainAnchorFinder`,
`SimpleAlignmentHandler`: nothing of the database's size); the SAM header lists only the genes
aligned to; `ProfileWrapper` profiles the samples in parallel, the largest SAM first, each on
threads in proportion to its size. In `c1` (16 samples, 4 threads): 80.9 s aligning, of which the
two 150 Mb Nanopore samples 43.5 s and the two 150 Mb PacBio samples 27.5 s, the 500k-pair samples
6.7 s (pe) and 2.6 s (se), the shallow samples 0.01-0.2 s each; 1.0 s profiling; 0.1 s for the SAM
headers; 85.9 s in all.

So the samples already share one run; per sample there is nothing left to amortise. What a run
per collection still costs: the start-up twice (training data, test set) and the four parity
runs' gene loads; the training data and the test set could be profiled in one run once both are
simulated (fix 2). At r226 the alignment itself (2.5 + 1.3 CPU h above, if 3× the tuning world's
cost per read) is the part to watch in the first run's `profile_all/protal.log` ("Aligning reads
took" per sample).

## The fixes, in detail (as proposed; 1, 2, 3 and 5 implemented the next day, see the follow-up)

1. **Long reads by templates.** In `collect_training_data.simulate_long`, per long-read sample: draw
   the reads (genome by relative abundance × genome length, as now; length from a gamma with the
   setup's mean and SD; start uniform over the genome's contigs; cut at the contig's end; either
   strand) until their bases reach the sample's, write them as one template FASTA and run
   `pbsim --strategy templ --method {errhmm,qshmm} ... --template T --accuracy-mean A
   [--difference-ratio R] --seed S` once. pbsim3 writes one `.fq.gz` and one `.maf.gz` per call,
   instead of three files per contig of every genome. The samples change (as they should: these are
   the reads the design asks for), so the long-read models must be retrained and their tests
   compared; the paired-end and single-end data are not touched. Reading the genomes for the
   templates should be done in processes, not threads, or with bytes operations (`gzip.decompress`,
   `split(b">")`), so that the GIL does not serialise it as it does `long_contigs` now.
   Cheaper alternative that keeps per-genome calls: `--length-max` a few times the setup's mean
   (the table then has 50-100k instead of 1M lengths, ~10× less per contig); it does not fix the
   truncated reads.
2. **Simulate during the builds.** The simulator needs the genome table and, for the held-out clades
   in every sample, `heldout_species.txt`, which needs only the conversion's taxonomy. A collector
   option to simulate only (both collections at once, `--threads` shared with the build) started as
   a background job right after the holdout is chosen, and the profiling after the training
   database's build, removes the simulations from the critical path. The collector's resume keys
   (`simulated.json`) do not depend on the database, so the profiling run takes the samples as
   they are. Profiling both collections in one protal run saves one start-up.
3. **Genome lengths in the table.** `make_genome_table` writes a fourth column, the genome's length
   as the simulator counts it (letters outside header lines), computed once in parallel; every
   reader of the table takes fields 0-2 or the lineage, and `read_genome_table` takes a fourth
   column as the length. Simulations stay byte-identical; the table's hash, in the stage keys,
   changes once, so the first rerun converts and simulates again.
4. **The trainer's previous-procedure study.** It exists to compare the new training procedure
   with the old one on r226; once that comparison is in hand, `--evaluation basic` (or a cheaper
   grid: 3 folds, every second `max_features`) saves most of the training step (below).
5. **Gene placement by an index.** In `gene_neighbours.place`, index the genome's k-mers (k ≈ 24)
   at a stride s (≈ 16) and look up each gene's first s k-mers (both strands); verify the hits.
   Every exact occurrence is found, so the first one and the repeat flag are the same as with
   `find`.
6. **One `--add_model` call for all models**, or packing the models in the build: three fewer
   rewrites of `database.protal` (21 GB at r226).

## The console output (changed)

Asked for alongside: the status lines every 10 minutes (`--progress-every 600`) made the console
hard to read. `build_gtdb_database.py` now prints numbered steps, one line when a step starts (with
its log) and indented lines when it ends: time, peak memory and a few numbers of what it made
(sizes of the databases, present/absent taxa per read type of each collection, each model's F1 with
species held out and on the test set, the parity result). The genome table and what the training
database leaves out are one line each (the details stay in `genome_table.txt` and
`model_logs/holdout.txt`); the removal of a `full_reference.fna.zst` is part of its build's line; the
finished database's background build gets its line when it ends. `--progress-every` stays, default
0 (off). `scripts/mini_db/test_mini_db.py` (`GtdbBuildTest`) checks the new lines: 3 tests OK on
`25d457e` with the two files changed (471 s), and again on `a212559` (287 s), after another session
had committed its own changes to the script (`trace_relatives`, whose line follows the same indented
style). The example in
[building-a-database.md](../../databases.md#what-it-prints) is from a small build of the
tuning world with the changed scripts (first `demo`, 22:45 min on the loaded machine; since the
follow-up `demo2`, with the background simulations, 1:42 on an idle one).

## Follow-up (2026-10-02): fixes 1, 2, 3 and 5 implemented

On `audit-fixes` at `5c7d64c` (uncommitted when measured): `scripts/collect_training_data.py`,
`scripts/build_gtdb_database.py`, `scripts/mini_db/gene_neighbours.py`, their tests in
`scripts/mini_db/test_mini_db.py`, and `docs/databases.md`, `docs/databases.md`,
`docs/development.md`.

1. **Long reads by templates.** `collect_training_data.simulate_long` draws each long-read sample's
   reads (`long_read_templates`: a read's genome by relative abundance × genome length, its length
   from the setup's gamma between 100 bp and 1 Mb as pbsim3's table, its start uniform over the
   contigs of 100 bases or more, cut at the contig's end, either strand; until the reads' bases reach
   the sample's, each genome read once per round) and runs `pbsim --strategy templ` once per sample
   (`long_read_sample`); the reads, which pbsim3 names `r_<n>` in template order, are renamed
   `g<genome>x_<n>` as before. The samples of all long-read points run `--jobs` at a time. A point's
   key has `reads` (`LONG_READS`), so points of the old collector are simulated again. `long_contigs`
   and `long_read_genome` are gone.
2. **Simulating during the builds.** `collect_training_data.py --simulate_only` simulates and stops
   (it reads no database). `build_gtdb_database.py` starts it for both collections as background jobs
   at `nice` 10 as soon as the species to leave out are chosen (`training_data_simulation.log`,
   `test_data_simulation.log`), before the builds start; each collection's step waits for its
   simulations if they are still running, then runs the collector, which finds every point simulated
   and profiles. A failing simulation stops the run within seconds, as a failing build does.
3. **Genome lengths in the table.** `make_genome_table` (and, for `--genome-table` without lengths,
   `with_lengths`, which writes `OUTDIR/genomes.tsv`) adds each genome's length as the fourth column
   (`genome_length`: the letters outside header lines, as `simulate_metagenomes` counts them),
   counted in `--threads` processes and kept on a rerun for the FASTAs not changed since.
5. **Gene placement by an index.** `gene_neighbours.GenomeIndex`: the genome's 32-mers at every 32nd
   position; a gene's 32 first 32-mers looked up, hits checked against the whole gene; shorter genes
   by `find`. Every occurrence is found, so placements and repeats are those of the scans.

**Same outputs where they should be** (`scripts/check_lengths.sh`, `scripts/check_gene_neighbours.sh`):

- The Python lengths of all 2,295 padded genomes equal the padding plan's sizes and the simulator's
  own counts (`simulate_metagenomes --test` on the table without lengths, every species once in 3
  samples: 2,295 genomes, 0 differ). A sample simulated with the same seed from the table with and
  without lengths has byte-identical reads and the same manifest (but for its output folder).
- `gene_neighbours.py` on the 765 real-sized representatives: `gene_neighbours.tsv` (659,550 lines)
  and the positions (85,930) byte-identical to those of the version before, in 70 instead of 680
  CPU-s (25 instead of 175 s wall on 4 threads), 9.7× less.

**The long reads** (`collect_real_size.sh` again, `c2`, the same communities and seeds as `c1`;
`long_read_stats.sh`):

| Sample | Before (`c1`): reads, median | Now (`c2`): bases, reads, median |
|---|---|---|
| PacBio 300 kb (2 samples) | 4,318 and 13,293 reads, 100 bp | 310 and 309 kb, 25 and 28 reads, 14,198 and 12,040 bp |
| Nanopore 300 kb | 4,471 and 13,396, 100 bp | 301 and 299 kb, 47 and 40, 5,326 and 5,217 bp |
| PacBio 150 Mb | 30,983 and 21,199, 1,720 and 4,933 bp | 150.4 Mb each, 12,275 and 11,804, 13,185 and 13,624 bp |
| Nanopore 150 Mb | 48,715 and 35,892, 1,377 and 2,521 bp | 149.2 Mb each, 21,954 and 21,488, 5,475 and 5,555 bp |

The setups ask for 15 kb (PacBio) and a gamma of mean 8 kb, SD 6 kb (Nanopore, median ~6.5 kb);
reads cut at contig ends make the medians a little shorter. The 300 kb points now see what so few
reads can see: 2 present taxa in 2 samples (PacBio), against 74 before; the paired-end point of the
same communities and bases sees 20. (In the padded genomes every marker gene lies in the genome's
own 270 kb, 8% of it, so a long read hits one less often than in a real genome, whose markers are
spread over it.) The pb and ont models need to be trained again on such data, and their test sets
compared.

pbsim3: 8 runs and 469 CPU-s instead of 1,252 runs and 2,999 CPU-s; the whole collection 905 CPU-s
and 4:50 instead of 3,551 CPU-s and 26:44 (the load was lower than during `c1`, 5-9 against ~15, so
the CPU figures compare better than the wall times); the long-read points 2:01 instead of 22:42.

**The whole build** (`b2`: `build_real_size.sh` as `b1`, 4 threads, the changed scripts on `5c7d64c`,
the genome table without lengths, so that the run writes `OUTDIR/genomes.tsv` with them; `run_b2.sh`).
The first run used `a212559`'s binaries and stopped at the trainer, which at `5c7d64c` wants the
feature that commit added (`conserved_fast_kept_ratio`), which protal of `a212559` does not write;
rerun into the same OUTDIR with `5c7d64c`'s protal, it kept the samples
(their simulations took 5 s each, finding every point done), built both databases again for the
other protal, profiled again and finished (2:57). Load 5-10, as during `b1` (9-16).

| Step | `b1` (before) | `b2` (after) |
|---|---:|---:|
| genome table | 0:00:00 | 0:00:59 (lengths of 2,295 genomes; 16 s on an idle machine with the process pool) |
| converting, gene neighbours | 0:04:08 | 0:00:21 |
| the training database's files and build | 0:00:46 | 0:00:29, the simulations started before it |
| training data: simulation | 0:48:14 (pe 11:15, long reads 36:58) | 0:09:51 in the background (pe 5:46, long reads 4:04) |
| training data: protal (40 samples) | 0:00:38 | 0:00:31 |
| test set: simulation | 0:51:51 (pe 11:51, long reads 39:58) | 0:12:21 in the background, beside the training data's (pe 10:06, long reads 2:14) |
| test set: protal (46 samples) | 0:00:49 | 0:00:37 |
| **up to the trainer** | **1:46:26, ~14,300 CPU-s** | **0:14:21, 2,795 CPU-s** |
| trainer, parity, `--add_model` | 0:01:52 | 0:01:08 (the rerun) |

pbsim3: 10,070 CPU-s in 4,548 runs before, 643 CPU-s in 20 runs after; the 33 length scans
(~1,600 CPU-s) are one count of 63 CPU-s; ART 886 and 807 CPU-s. The whole build takes ~15:30
instead of 1:48:20 here. The paired-end simulations are now the longest part (they run beside each
other and the builds, at `nice` 10); at r226 the builds (45-80 min each) are, and the simulations
run during them.

`r226_estimate.py` with the measured costs of the changed collection (pbsim3 in template mode 0.78
CPU-s per Mb for both read types, more per base than per genome before; the lengths 27 ms per genome
once): the collections and the gene neighbours ~8.0 CPU h instead of ~31.7, of which ART 2.0, protal
3.8 (×3 the tuning world's cost per read, unknown), pbsim3 1.4.

**Tests.** `scripts/mini_db/test_mini_db.py`: the long-read tests rewritten (`test_long_read_replay`:
a stand-in pbsim3 for `--strategy templ`, reads from either strand of their genome, shares by
abundance × length, the same reads on a rerun; `test_long_read_templates`: contigs under 100 bases
left out, the read lengths, pbsim3 failing or making fewer reads than templates), new
`test_genome_table_lengths` and `GeneNeighboursTest.test_genome_index_finds_what_find_finds`;
`GtdbBuildTest` checks the background simulations (before the training database's build line, the
collector's log with nothing simulated), the genome table's lengths against the simulator's own
counts, and stops the run by SIGTERM while the simulations run. On `5c7d64c` with the changed
scripts and its binaries: the mini-database suite, 36 tests OK (80 s), `GtdbBuildTest` among them
(3 tests, 67 s).

## Follow-up 2 (2026-10-02): fix 4, the previous procedure cheaper and off by default

`random_forest_cmdline.py` and `build_gtdb_database.py` take `--previous-procedure` (and
`--no-previous-procedure`, the default): the comparison with the trainer's previous procedure runs only when
asked, with any `--evaluation` but `none`, and `build_metadata.tsv` records it (`classifier_previous_procedure`:
`compared` or `not compared`). For the first builds of a release, `build_gtdb_database.py ...
--previous-procedure`.

When it runs it is cheaper. Its grid search takes every second value of `max_features` (6, 8, ..., 24 and 25),
then the two next to the best (the same folds and forest seeds), with 3 folds instead of 5, on at most 20,000 rows
(whole samples drawn at random, `PREVIOUS_GRID_ROWS`); and the current procedure's side of the comparison is the
evaluation's own forests with species held out instead of 5 more. What it judges (a forest of 512 trees on the top
features, on a random 20% of rows and with species held out) is the procedure's as before.

Measured with [`previous_procedure_cost.sh`](scripts/previous_procedure_cost.sh) on the paired-end training table
of an operon-world build (`~/opw/b_gn2`, 72 samples, the gene neighbour report of the same day) and on that table
five times over (`enlarge_table.sh`), `--evaluation full`, 4 threads, the trainer before against after; the study's
CPU is the training's less that of a training without it. Another session's benchmark ran meanwhile (load 6-7 on 6
cores during the before and after runs, ~2 during those without the study), so the seconds are inflated and the
ratios are what counts:

| rows | trainer | grid | study, wall s | training, CPU s | the study, CPU s | its `max_features` | its F1, species held out | the current procedure's F1 |
|---|---|---|---:|---:|---:|---:|---:|---:|
| 4,883 | before | 20 values × 5 folds, all rows | 62 | 167 | 146 | 7 | 0.9954 | 0.9979 |
| 4,883 | after | 12 values × 3 folds, all rows | 35 | 90 | 69 | 6 | 0.9954 | 0.9979 |
| 24,415 | before | 20 values × 5 folds, all rows | 194 | 833 | 722 | 7 | 0.9970 | 0.9976 |
| 24,415 | after | 12 values × 3 folds, 19,998 rows | 51 | 291 | 181 | 7 | 0.9972 | 0.9976 |

The study costs 2.1 times less CPU at 4,883 rows and 4.0 times less at 24,415; above 20,000 rows its grid stays
at the cost of 20,000, so at the ~70,000 rows estimated for an r226 paired-end table (above) it is about 12 times
less on the grid (100 forests on 56,000 rows against 36 on 13,300, and more for the trainer's rows^1.3). The
result is the same: the previous procedure chose 7 (or 6) features and did no better with species held out than
the current one. Every second value alone (11 values, before the step to the neighbours) chose 6 at 24,415 rows,
where every value chooses 7, and lost 0.003 F1 there (0.9940); every value with 3 folds on 20,000 rows chose 7
(0.9972) for about three times the grid time of every second value. Tests: `test_model_pmml.py` (the study only
when asked, both flags, the grid's values, folds and rows, the reused predictions; whole samples up to the row
limit), the mini-database `GtdbBuildTest` (`classifier_previous_procedure` `not compared` by default).

## Reproducing

All scripts take their inputs from WSL `~/tune` (the tuning world and its V3 build) and write to
`~/bprof`; binaries from `git archive 25d457e` built in `~/build-prof`.

```
bash scripts/pad_genomes.sh ~/tune/V3/genomes.tsv ~/bprof/world 4      # real-sized genomes
bash scripts/length_scan.sh ~/bprof/world/genomes.tsv ~/bprof/lenscan
bash scripts/pbsim_fixed_cost.sh ~/bprof/pbsim_fixed
bash scripts/collect_real_size.sh ~/bprof/c1 4                          # c1
bash scripts/tool_summary.sh ~/bprof/c1/tool_times.tsv ~/bprof/world/plan.tsv
bash scripts/long_read_bases.sh ~/bprof/c1/collect ~/tune/V3/training
bash scripts/templ_prototype.sh ~/bprof/c1/collect ~/bprof/templ
bash scripts/gene_neighbours_cost.sh ~/bprof/gn 100
NICE=0 LEVELS=full bash scripts/trainer_scaling.sh ~/bprof/trainer 3 1
bash scripts/build_real_size.sh ~/bprof/b1 4                            # b1
bash scripts/build_summary.sh ~/bprof/b1
python3 scripts/r226_estimate.py
# the follow-up, with the changed scripts (SCRIPTS: a copy of HEAD's scripts with them; BIN: HEAD's binaries)
bash scripts/check_lengths.sh REPO/scripts ~/bprof/lengths2
bash scripts/check_gene_neighbours.sh REPO/scripts/mini_db/gene_neighbours.py ~/bprof/gncheck 4
bash scripts/collect_real_size.sh ~/bprof/c2 4 && bash scripts/long_read_stats.sh ~/bprof/c2/collect
BIN=... SCRIPTS=... bash scripts/build_real_size.sh ~/bprof/b2 4 && bash scripts/build_summary.sh ~/bprof/b2
# follow-up 2 (BEFORE: a copy of the scripts with the trainer before the change)
bash scripts/previous_procedure_cost.sh ~/opw/b_gn2/training/training_data.tsv BEFORE REPO/scripts ~/bprof/prev/run
```
