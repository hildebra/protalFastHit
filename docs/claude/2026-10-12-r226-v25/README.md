# The r226 v25 build: slower on a slower node, and the soft fragment count does not help

**Data.**
- **v25:** SLURM 24201031, 2026-10-12 14:25:54 to 17:07:53, 84 threads, `-t 84 --scratch` on the node's SSD, built
  with `--fragment-count soft`. protal and the simulator at `5222386` (`build_metadata.tsv`); the scripts at
  `5222386`, and `18a95b1` by the end of the run (the only script that changed between them is
  `measure_performance.sh`, which the build does not run; the protal binary stayed at `5222386`). The weighted site shift is off (protal's default since `1aadfd2`). The node has 1.5 TB of
  memory (protal's runs: "Total available memory is 1510GB"); the archive does not name it. The run is a restart:
  an earlier job (24198933, from 12:15) stopped after simulating, and v25 reused its converted release. The share
  archive `local/protal0.8.1_r226_v25_share.tar.gz` is unpacked to `local/v25`.
- **Compared with v24b** ([2026-10-12-r226-v24b](../2026-10-12-r226-v24b/README.md)): SLURM 24171430, q512n11 (one
  AMD EPYC 9634 with 84 cores, 502 GB), protal `4a42f3b`, the shift on in the training runs, the hard (MAPQ-filtered)
  fragment count. The two builds share the download, seed, design and hold-outs, and **they profiled the same reads**:
  every one of the 1,482 samples looked up the same k-mers, scanned the same index cells and wrote the same records
  in both ([run_log_work.txt](run_log_work.txt)).
- **Commands** (WSL, cores 0-3, `nice 5`; the outputs are in this folder):
  - [run_log_work.py](run_log_work.py) `local/v24b/... local/v25/...` → [run_log_work.txt](run_log_work.txt): each
    sample's alignment work and its aligning and profiling timers, matched across the builds.
  - `logs/cpu_jobs.tsv`, `logs/cpu_timeline.tsv` and `logs/protal_cpu_{training,test}.tsv` of both builds, summed with
    short perl one-liners (section 1).
  - [soft_pairs.py](soft_pairs.py) → [soft_pairs.txt](soft_pairs.txt): both builds' tables and held-out and test
    predictions, matched by (sample, taxon); which features differ, how `fragments` moved, the calls that flip, and F1
    with a paired bootstrap over samples (1,000 resamples).
  - [approx_soft.py](approx_soft.py) → [approx_soft.txt](approx_soft.txt): the soft count rebuilt from a table's own
    columns, checked on pb and ont, where both counts are known.
  - [refit.sh](refit.sh) (v24b's trainer, `git archive 4a42f3b scripts`, species folds only) on v24b's pe and se
    tables, as they are and with `fragments` replaced by the rebuilt soft count; pe also with `--seed 2`.
    [refit_pairs.py](refit_pairs.py) compares two refits on the same rows → [refits.txt](refits.txt).

**Questions (the user's, 2026-10-12).** a) Why is this build slower? b) Does `--fragment-count soft` help?

## Summary

- **a) The node is slower, not the code.** v25 took 2:41:59 against v24b's 1:51:26 (+50 min). The same alignment work
  took 1.40× (se), 1.53× (pb) and 1.54× (ont) as long (aligning timers of samples with identical work), with
  ~77 of the 84 cores busy in both builds. The trainer's single-threaded table load took 2.2-2.4× as long, and the
  simulator 1.5× the CPU. The node had 1.5 TB of memory (v24b's 502 GB, the EPYC 9634 nodes), like q1536n4, where v18
  and v21 ran. Of the +50 minutes: the collection +23.5, waiting for `protal_db`'s background build +15.3, the setup
  +10.5 (the training database +8.7; ~11 minutes of cold reads from `/hpc-home` with the cores idle took the place
  of v24b's 10 minutes of converting), the training +1.4.
- **The soft count costs no time.** It is a few multiplications per taxon (`MicrobialProfile::SoftFragments`), and
  pb/ont profiling slowed by the same factor as their alignment (1.58/1.52 against 1.55/1.55). pe's aligning even
  slowed less (1.14×) than the node's factor, because the weighted shift is off and whole index blocks are prefetched
  (`e04f4a1`).
- **b) The soft count does not help, in any read type.** For pb and ont the two builds differ only in the `fragments`
  feature (same rows, same values otherwise): F1 −0.0005 and +0.0002, and the thin present taxa it was meant to
  rescue stay missed (pb 142 → 145, ont 161 → 168 false negatives). The feature lifts 12-19× as many absent taxa as
  present ones from ≤ 2 to ≥ 3 fragments: the MAPQ-filtered reads are mostly those of absent taxa. Refits of v24b's
  pe and se tables with the soft count rebuilt in place of the hard one change F1 by no more than another random
  seed does (pe −0.0001 at knob 0.5, −0.0002 at the best knob, seed 2 alone: the same; se +0.0001 and +0.0002).
- **v25 is still better for short reads** on the rows of both builds (pe test set +0.006, held out +0.0008 at the
  knob), but that is the weighted shift being off: with it on, some present taxa lost all their records.
- **The console's marker-truth F1 is wrong** (pe 0.954 → 0.633): the collector counts held-out species, which the
  training database lacks, as missed (section 3; fixed with this report).

## 1. Why the build is slower

### Where the 50 minutes went

| | v24b | v25 | difference |
|---|---|---|---|
| start to the collections (genome table, species clouds, training database) | 0:34:46 | 0:45:14 | +10:28 |
| both collections (to the later one's end) | 1:03:25 | 1:26:53 | +23:28 |
| training the four models (step 8/11 to 9/11) | 0:02:33 | 0:03:54 | +1:21 |
| training's end to the end of the run (`protal_db` built in the background, reports) | 0:10:42 | 0:25:58 | +15:16 |
| **whole run** | **1:51:26** | **2:41:59** | **+50:33** |

The steps, from `logs/cpu_jobs.tsv` (wall, CPU seconds, mean cores busy):

| step | v24b | v25 |
|---|---|---|
| genome table (no CPU log; console) | 3:43 | 6:43, ~0 cores busy (`cpu_timeline.tsv`) |
| converting the release | 6:22 (beside the genome table) | reused from the stopped job |
| species clouds | 10 s, 552 CPU s, 54 cores | 245 s, 782 CPU s, 3.2 cores (33 cores in iowait) |
| leaving the species out | 134 s, 1,235 CPU s | 289 s, 2,668 CPU s |
| training database index | 1,318 s, 73,983 CPU s, 56 cores | 1,680 s, 86,593 CPU s, 52 cores |
| simulating the training data / test set | 697 / 553 s, 7,043 / 3,180 CPU s | 836 / 631 s, 10,662 / 4,840 CPU s |
| profiling the training data / test set | 3,654 / 3,805 s, 107,274 / 43,105 CPU s | 4,982 / 5,212 s, 174,702 / 69,785 CPU s |
| `protal_db` in the background (idle class) | 4,527 s, 150,591 CPU s | 7,049 s, 185,230 CPU s |
| training pe / se / pb / ont | 122 / 152 / 100 / 108 s | 231 / 158 / 108 / 112 s |

- **The setup waited on `/hpc-home`.** v25 skipped the conversion (the stopped job had done it), so its genome table
  and species clouds read the release and genomes cold from the network file system: 6.5 minutes with no core busy,
  then 4 minutes with 33 cores in iowait. In v24b the conversion had just written those files, and they were in
  the page cache. Up to the training database, both builds took about as long (11:34 against 10:00, v24b waiting for
  its conversion and gene neighbours); the stopped job's own conversion took 28:29 against v24b's 6:22. The training
  database then took 32:51 against 24:12: its index +6:02, leaving the species out and writing its files +2:35.
- **Everything compute-bound took 1.2-2.4× as long for the same work.** Identical work: the training database
  index (+17% CPU, +27% wall), the finished database's index (+23% CPU), the simulator (+51% CPU), the trainer's
  table load (2.2-2.4× wall, single-threaded).

### The collection: the same reads, aligned 1.4-1.5× slower

The two builds aligned exactly the same work: every sample's k-mers looked up, flex cells scanned, alignments made and
records written are equal ([run_log_work.txt](run_log_work.txt); pe's "candidate alignments tried" differ per sample
but not in sum). The time it took, v25 / v24b:

| | pe | se | pb | ont |
|---|---|---|---|---|
| aligning, all samples | 1.14 | 1.40 | 1.55 | 1.55 |
| aligning, samples with identical work counts, ≥ 1 s | — (3 samples) | 1.40 | 1.53 | 1.54 |
| profiling, all samples | 1.53 | 1.85 | 1.58 | 1.52 |
| profiling, per-sample median (≥ 1 s) | 1.53 | 1.55 | 1.54 | 1.55 |
| protal process CPU during the aligning stages, training design samples (`protal_cpu_training.tsv`) | 1.74 | 1.77 | 1.72 | 1.52 |

- Both builds kept the machine busy during the collection: 78.3 cores in v24b, 76.8 in v25
  (`cpu_timeline.tsv`). So v25 used 1.34× the core-seconds (382,000 against 285,000) for the same reads.
- Long reads, which neither the shift nor the prefetch change much, show the node's own factor: ~1.55. pe aligned
  relatively faster (1.14): the weighted shift is off in v25 (in v24b it re-scored short reads' alignments with the
  column weights), and `e04f4a1` prefetches each lookup's whole block. se (1.40) gains less of both.
- Profiling slowed by the same ~1.55 (medians). The sums for se (1.85) are inflated by two scenario samples (one in each
  collection's streamed run) whose "record evidence and sample context" took 43.2 s each against 4-5 s in v24b:
  the same length in two processes points to a stall of the node, not to the code. 7 of 166 samples with at least
  0.5 s there slowed more than 3×.
- The protal runs took 25-45% longer: training run 1 26:59 → 34:26, the streamed run 2 27:59 → 40:31; test runs
  14:37 → 18:16, 31:53 → 41:30, 3:28 → 4:34. Loading the index did not slow (43 s in both).

### The end of the run: the background build

`protal_db` is built at the idle scheduling class on the cores the collections leave free. In v25 it needed 23% more
CPU for the same index, and the longer, busier collection left it less. It finished 26 minutes after the models were
trained, against 9 minutes in v24b.

### Which node

The archive does not name v25's node. Its 1.5 TB of memory differs from the q512 nodes (EPYC 9634, 502 GB) that ran
v24 and v24b. v18 and v21 ran on q1536n4, also with 1.5 TB, and v21's CPU report already found that short reads gained
little there from more cores ([2026-10-09-r226-v21-cpu.md](../2026-10-09-r226-v21-cpu.md)). A single-threaded step
2.2-2.4× slower and the threaded ones 1.2-1.6× slower fit an older, lower-clocked CPU with more cores, or other jobs on
the same node. Only `lscpu` in the job would tell which.

## 2. Does `--fragment-count soft` help? No

### The headline numbers mix it with the shift

F1 from each build's `console.log` (species held out at knob 0.5; the test set as protal calls):

| | pe held out | pe test | se held out | se test | pb held out | pb test | ont held out | ont test |
|---|---|---|---|---|---|---|---|---|
| v24b (hard count, shift on) | 0.9513 | 0.9652 | 0.9591 | 0.9630 | 0.9714 | 0.9842 | 0.9649 | 0.9732 |
| v25 (soft count, shift off) | 0.9492 | 0.9666 | 0.9576 | 0.9619 | 0.9710 | 0.9832 | 0.9650 | 0.9736 |

These are not comparable for short reads: each build scores only its own rows (taxa with any evidence), and v25's
pe table has 29,831 rows v24b lacks, 435 of them present taxa, while v24b has 6,261 (29 present) that v25 lacks (se:
14,400 and 254 against 775 and 2). With the shift on, some present taxa lost all their records to other taxa and got
no row, so v24b's F1 never counted them as missed. The pb and ont tables have exactly the same rows in both builds.

### Long reads: only `fragments` changes, and nothing improves

For pb and ont the shift does not apply (protal applies it to short reads only), and on the shared rows the only
feature that differs is `fragments` ([soft_pairs.txt](soft_pairs.txt)). That makes the builds a clean test of the
soft count, at the same seed and trainer settings:

| species held out | pb | ont |
|---|---|---|
| rows whose `fragments` moved: present / absent | 2.7% / 27.6% | 4.0% / 29.7% |
| rows lifted from ≤ 2 to ≥ 3 fragments: present / absent | 288 / 3,849 | 405 / 4,817 |
| F1 v25 − v24b at knob 0.5 (95% interval, paired over samples) | −0.0005 [−0.0008, −0.0001] | +0.0002 [−0.0002, +0.0006] |
| at each build's best knob | −0.0004 (0.46 / 0.48) | +0.0001 (0.53 / 0.53) |
| test set at knob 0.5 | −0.0010 [−0.0041, +0.0013] | +0.0004 [−0.0015, +0.0024] |
| the intended class (present, ≤ 2 kept, ≥ 3 before the filters): FN v24b → v25 | 142 → 145 of 374 | 161 → 168 of 511 |
| the same pattern among absent taxa: FP v24b → v25 | 39 → 40 of 5,305 | 54 → 50 of 7,075 |

The soft count moves absent taxa 7-10 times as often as present ones. The reads the MAPQ filter drops are mostly reads
of absent taxa (a held-out congener's reads tying with a database species), and the read EM leaves them a large
enough share. So 13× (pb) and 12× (ont) as many absent as present taxa cross from ≤ 2 to ≥ 3 fragments, and the
model cannot use the feature to rescue the thin present ones: they stay missed.

### Short reads: the same, isolated by a refit

The builds cannot separate the soft count from the shift for pe and se. So v24b's own tables (shift on) were refitted
with v24b's trainer twice: as they are, and with `fragments` replaced by the soft count. The tables lack the record
counts protal uses, so the soft count is rebuilt with records taken in proportion to fragments
([approx_soft.py](approx_soft.py)). On pb and ont, where the true soft count is known, this matches it closely: rank
correlation 0.9991 on the rows it moves, the mean increase 11.87 against 11.91 (pb), 97% of moved rows within 0.5
([approx_soft.txt](approx_soft.txt)). A third refit with `--seed 2` measures the trainer's own noise
([refits.txt](refits.txt)).

| pe, species held out | hard → soft | seed 1 → seed 2 (hard) |
|---|---|---|
| F1 at knob 0.5 | 0.9514 → 0.9513 (−0.0001 [−0.0005, +0.0002]) | 0.9514 → 0.9513 (−0.0002) |
| F1 at the chosen knob (0.72 for all three) | 0.9560 → 0.9561 (+0.0001) | 0.9560 → 0.9561 (+0.0001) |
| F1 at the best knob | 0.9563 → 0.9561 (−0.0002) | 0.9563 → 0.9561 (−0.0002) |
| calls that flip at knob 0.5 | 1,343 | 1,552 |
| test set at the chosen knob | +0.0003 [−0.0012, +0.0019] | −0.0009 [−0.0028, +0.0010] |
| intended class (1,616 present rows lifted to ≥ 3): FN, best knob | 464 → 445 | — |
| the same pattern among absent taxa (31,243 rows): FP, best knob | 137 → 167 | — |

The soft count changes pe's model less than another random seed does. It rescues 19 of the 464 thin present taxa at
the best knob and adds 30 false positives among the 19× as many absent taxa that it lifts the same way.

se, the same refits ([refits.txt](refits.txt)): hard → soft F1 0.9590 → 0.9591 at knob 0.5 (+0.0001 [−0.0004,
+0.0006]), 0.9610 → 0.9613 at the best knob (0.68 in both; +0.0002 [−0.0004, +0.0008]); the test set −0.0014
[−0.0033, +0.0004] at knob 0.5. The 1,030 thin present rows the soft count lifts: 245 → 251 missed at the best knob;
the 16,701 absent ones: 82 → 83 false positives. (At the refits' own chosen knobs the difference is +0.0022, but only
because the hard refit's knob search kept 0.5 and the soft one's chose 0.67; v24b's build chose 0.7 on the same
table. The best-knob comparison takes the knob search out.)

### What made v25's short reads better

On the union of both builds' rows (a row a build lacks counts as not called), v25 is ahead for short reads: pe
+0.0008 [+0.0002, +0.0014] at the knob both call at (0.72), +0.0058 [+0.0034, +0.0082] on the test set; se +0.0007
[−0.0002, +0.0018] held out at the builds' knobs, +0.0061 [+0.0037, +0.0087] on the test set ([soft_pairs.txt](soft_pairs.txt)).
These intervals resample samples only; the trainer's seed adds about ±0.0002 held out and ±0.001 on the test set.
The refits give the soft count no measurable part in this, so it is most likely the shift being off, largely through
the present taxa it no longer strips of all their records (pe test set: 37 of v24b's 180 false negatives at knob 0.5
are taxa without a row in v24b). That agrees with turning the shift off (`1aadfd2`), which was
decided on depth alone ([2026-10-12-r226-v24b](../2026-10-12-r226-v24b/README.md), section 4). The composition
check agrees too: pe species' depth −4.8% → −4.0%, explained share −0.011 → −0.002; se −4.8% → −3.4%, −0.018 →
−0.007.

## 3. The marker-truth F1 counts held-out species as missed

The console's new line reads "pe 0.954 -> 0.633 (… 78388 unseen counted as missed)": with only the species whose
marker genes got reads, F1 should rise, since species protal cannot see are left out. It falls because of the
`unseen` rows (`<table>.unseen.tsv`, present species with marker reads but no row, counted as false negatives).

- `collect_training_data.py` writes a species there if `taxon in db_lineages` (line 1759). `db_lineages` comes from
  `--taxonomy`, the converted release's `internal_taxonomy.dmp`. It names the species the training database leaves
  out: 19 of 20 sampled from `model_logs/heldout_species.txt` are in it. So every held-out species with a read on a
  marker gene counts as missed, though protal profiled against a database without it.
- The counts show it: soil's 10 hold-in pe samples have 22,681 present rows and 34,081 unseen species, while 60% of
  soil's species are held out (the scenario's table: 14,034 of 23,390). Over all pe samples with species held out,
  78,388 unseen against 76,888 present rows.
- **Fixed in the commit that adds this report** (the user's request, 2026-10-12): `unseen_species()` also skips the
  species of `--novel_species`, the held-out list; `test_collector.py`'s `test_unseen_species` checks it, and
  `docs/databases.md` says so. The marker-truth rows of this build cannot be corrected from the archive: the unseen
  tables are not in it. The next build's will be right.

## 4. What to do

1. **Pin builds to one node type, and log it.** Build times are only comparable on the same hardware: this run's
   +50 minutes is the node. Ask SLURM for the q512 (EPYC 9634) nodes that v24 and v24b ran on, and have
   `build_gtdb_database.py` write the host name and CPU model (`lscpu`) into `console.log` and `build_metadata.tsv`.
   Today only the memory size in protal's log hints at the node.
2. **Keep `--fragment-count mapq` (the default); do not build with soft.** It costs nothing, but gains nothing in any
   read type: long reads within ±0.0005, pe within the trainer's seed noise. The feature lifts 12-19× as many absent
   taxa as present ones over the 2-fragment line. If thin strains are to be rescued, the evidence has to separate a
   strain's own MAPQ-0 reads from a held-out congener's ([2026-10-12-r226-v24b-errors](../2026-10-12-r226-v24b-errors/README.md)
   found no such read signal).
3. **v25 is the better database for short reads**: the shift off, at no cost to long reads. The soft count in it is
   harmless (the models were trained with it, and protal follows `build_settings.tsv`). The next build can go back
   to `mapq` without loss.
4. **The marker-truth unseen rows** (section 3): fixed; the next build's marker-truth F1 means what it says.
5. **Restarting a build on another node costs its page cache**: the reused conversion was read cold from `/hpc-home`
   (genome table and species clouds 11 minutes with the cores idle). Not worth changing; worth knowing when timing.
