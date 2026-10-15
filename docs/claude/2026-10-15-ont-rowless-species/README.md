# Why ONT loses 3.7% of present species entirely: the identity floor counts indels twice

Date: 2026-10-15. Branch `audit-fixes` at `55f850d` (protal 0.8.4) for every run; the recommendation's code
references are to `d6ad7b9` (0.8.5). Nothing in protal was changed; the simulator was patched locally to name reads by
their origin (below).

**The question.** In the r226 v29 build ([2026-10-15-r226-v29](../2026-10-15-r226-v29/README.md), section 7) ONT gets
no row at all for 3.7% of the present species whose genomes gave reads on a marker gene (species held out; 5.2% of
the test set), against 1.4-1.6% for pe, se and pb. These "unseen" species are half of ONT's marker-truth misses
(1,946 of 3,887). v27 had the same rate, so the sparse-first gates are not the cause. Where are their reads lost?

**Data and commands** (WSL, cores 0-3, `nice 5`; scripts in this folder):
- **r226 v29's error reads** (`local/v29/.../model_logs/error_reads/<type>/taxa.tsv.gz`): for every rowless
  species, its own reads that got a record and where they went. They hold no SAMs or manifests of these species.
- **A small world built at `55f850d`**: the real small world (`~/realworld`: 6 GTDB families, 170 species, 352
  genomes; [real small world](../2026-10-10-real-small-world/README.md)), ONT only
  ([rw_build.sh](rw_build.sh): `build_gtdb_database.py --read-types ont --model gbm --model-gates none --scenarios
  none`, 60 training and 20 test samples of 0.3-150 Mb, the default ONT setup
  `qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36`). Its training database leaves 51 species out, as r226's does.
- **Reads named by their origin.** [read_origin_names.patch](read_origin_names.patch) makes `simulate_metagenomes`
  (with `PROTAL_SIM_READ_ORIGIN=1`) name each long read `g<i>x_<n>|<contig>|<start>|<end>|<errors>|<strand>`. The
  reads stay the same (the patch draws nothing): re-simulating the build's 150 Mb training point gave the same marker
  read counts per genome, and protal the same rows for all 12 samples ([trace_b150.sh](trace_b150.sh)).
- **The trace** ([trace_marker_reads.py](trace_marker_reads.py), [aggregate_marker_reads.py](aggregate_marker_reads.py)):
  each read's overlap with the marker genes of its genome (`gene_positions.tsv`; a marker read shares 50 bases or more
  with one, the simulator's `--marker_min_overlap`), its true accuracy (1 − the simulator's errors / template length),
  and its fate in protal's SAM (written with unmapped records for reads that seeded but aligned nowhere): its primary
  record on its own species at MAPQ ≥ 4 (`own`, as the profiler counts it), on its own species below MAPQ 4, on
  another species (`elsewhere`), seeded but aligned nowhere (`unaligned`), or no record at all (`none`). Rows: the
  taxa of the sample's `.truth_annotated` profile, as the collector takes them. The species held out of the training
  database are left out.
- **The floor experiment** ([exp_floor.sh](exp_floor.sh)): the 150 Mb point's 12 communities simulated at 6, 30 and
  150 Mb with origins, profiled on the small build's training database as the collector profiles
  (`--no_strains --no_qcmsa --sparse_index emulate`) at protal's ONT identity floor `-a` 0.85 (the default), 0.80 and
  0.75. The three depths share the communities and seeds, so many of a shallow sample's reads recur in the deeper
  ones: they are not independent samples.
- **Three builds** ([exp_builds.sh](exp_builds.sh) `0.85 0.80 0.75`) of the same world, seed and samples as the
  trace build but twice as many (120 training samples of 0.3-150 Mb, 40 test samples of 0.15-250 Mb; GBM, no gate
  models, no scenarios) that differ only in protal's ONT floor. protal at `55f850d` does not record `-a` in the
  database (0.8.5 does, in `build_settings.tsv`), so a wrapper adds `-a` to every protal call of the build. A first
  attempt gave `-a` to `--build` only, and its runs used 0.85: its 0.80 and 0.75 training tables were identical to
  the baseline's. The baseline is from that attempt (0.85 is ONT's default either way); the 0.80 and 0.75 builds were
  rerun with the wrapper, and their logs' `max score ani` lines and their tables' checksums differ. Compared with
  [compare_floors.py](compare_floors.py): F1 per set with the trainer's calls (`trained_model_ont.calls.tsv.gz`) and
  with the marker truth (unseen species counted as missed), with a paired bootstrap over samples (1,000 resamples).
  Its inputs (the builds' `model_logs` and unseen tables, 16 MB) are kept in WSL `~/ontu_floor_builds`; the reads,
  builds and traces were removed.

## Summary

- **Where the reads are lost: in the alignment's identity floor.** Mostly not seeding, the sparse-first gates or
  MAPQ. A rowless ONT species is one whose one or two marker reads protal fails to align. These reads seed but
  no alignment reaches protal's ONT identity floor (`-a` 0.85).
- **Why the floor drops them.** The identity counts a one-base indel as two differences, so the floor drops most
  reads below ~0.88 true accuracy. pbsim3's ONT model gives 12% of marker reads such accuracy. Of the marker reads,
  7.4% seed and align nowhere: 65-76% of those below 0.85 accuracy, 36% of those at 0.85-0.88, under 1% above 0.92.
- **Lowering the floor brings most of them back.** Three small builds differ only in the floor (0.85, 0.80, 0.75):
  - Rowless species: 126 → 72 → 58 with species held out (3.0% → 1.7% → 1.4% of the present species), 46 → 23 → 20
    on the test set (4.0% → 2.0% → 1.8%). r226 pe, se and pb lose 1.4-1.6%.
  - Marker-truth F1 at 0.80: +0.0036 with species held out, +0.0075 on the test set (the test gain is not
    significant).
  - At 0.75: +0.0050 and +0.0128, both significant.
- **The cost:**
  - more false positives, nearly all absent species with a congener left out of the database: +23/+6 at 0.80,
    +26/+2 at 0.75;
  - 17% (0.80) and 45% (0.75) more profiling time, and 0.11-0.15 GiB more peak memory at 4 threads.
  - F1 as protal calls falls, because the lower floors make rows for present species without marker reads. It
    cannot compare floors; the marker truth can.
- **Recommendation, implemented:** make 0.80 ONT's floor for databases built from now on, and keep 0.85 for
  databases that record none. 0.75 is the next step if GTDB-scale time allows.

## 1. What r226 v29's error reads show

The rowless species' own reads that got any record, by where their best record went (`error_reads/*/taxa.tsv.gz`;
the counts include reads off the marker genes, which the table cannot tell apart):

| read type | rowless species | with a record of their own reads | own reads with a record | best on itself (MAPQ ≥ 4) | elsewhere | seeded, aligned nowhere |
|---|---:|---:|---:|---:|---:|---:|
| pe | 24,901 | 3,856 | 15,606 | 1,938 (219) | 8,589 | 5,079 |
| se | 14,545 | 2,103 | 10,511 | 1,083 (128) | 6,063 | 3,365 |
| pb | 44,968 | 5,401 | 10,523 | 846 (148) | 2,000 | 7,677 |
| ont | 43,891 | 7,171 | 14,068 | 983 (131) | 2,011 | 11,074 |

(The "unseen" taxa of error_reads.py are every species of the sample without a row, most of them without a single
read; the marker truth's 1,946 are those with a marker read.) For ONT, 79% of the rowless species' own reads with a
record seeded and aligned nowhere; for pe the commonest fate is a record on another species. That points at the
alignment, not at the seeding or the gate.

## 2. The small world reproduces it

The small build's marker truth (`trained_model_ont.report.txt`): 65 rowless species of 2,371 present with species
held out (2.7%), 23 of 530 in the test set (4.3%), against r226's 3.7% and 5.2%. 61 of the 65 had a single marker
read, 4 had two: a species loses its row when its one or two marker reads are lost.

## 3. Which marker reads protal loses

Over the 12 communities at 150 Mb, at the default floor 0.85: 35,395 marker reads of species the training database
has. 87.9% are on their own species at MAPQ ≥ 4, 1.9% on it below MAPQ 4, 1.1% elsewhere, **7.4% seeded but aligned
nowhere** and 1.8% without a record.

**By the read's true accuracy** (the simulator's errors):

| accuracy | marker reads | own | own, MAPQ < 4 | elsewhere | unaligned | no record |
|---|---:|---:|---:|---:|---:|---:|
| < 0.80 | 462 | 0.0% | 0.0% | 0.0% | 64.5% | 35.5% |
| 0.80-0.85 | 1,117 | 11.0% | 0.4% | 0.1% | 75.6% | 13.0% |
| 0.85-0.88 | 2,606 | 57.6% | 2.3% | 0.4% | 36.4% | 3.3% |
| 0.88-0.90 | 3,131 | 81.4% | 4.0% | 0.9% | 11.8% | 1.9% |
| 0.90-0.92 | 2,525 | 89.3% | 4.5% | 0.6% | 3.7% | 1.8% |
| 0.92-0.94 | 2,953 | 93.8% | 3.1% | 1.3% | 0.9% | 0.9% |
| 0.94-0.96 | 2,944 | 95.7% | 1.9% | 1.2% | 0.4% | 0.8% |
| 0.96-0.98 | 7,672 | 96.7% | 1.4% | 1.3% | 0.2% | 0.5% |
| ≥ 0.98 | 11,985 | 97.6% | 0.8% | 1.3% | 0.1% | 0.3% |

The reads holding a whole gene behave the same (30 Mb: 82% of the whole-gene reads at 0.80-0.85 unaligned, 36% at
0.85-0.88, 11% at 0.88-0.90, ≤ 2.5% above 0.90). How much of a gene a read holds matters little: of the reads that
share only 50-99 bases with a gene, 45% have no record (too little to seed), but they are 1% of the marker reads.

**Why.** protal keeps an alignment only if its identity (`AlignmentInfo::GetProxyANI`) reaches the read type's
floor, `-a` (`max_score_ani`): 0.9, and 0.85 for ONT (`ReadType.h`). The identity is a score: 1 + score / (4 ×
length), where a mismatch costs 4 and a gap 6 + 2 per base (`AlignmentInfo::Score`), so a one-base indel counts as
two differences. The WFA budget is set from the same floor (`MaxScore`: (1 − floor) × overlap mismatches). The
simulated ONT errors are 61% indels (the setup's substitution/insertion/deletion ratios 39/24/36), and ONT's floor of
0.85 was set lower for that reason. In practice it still drops most reads below ~0.88 true accuracy, and pbsim3's ONT
model draws each read's accuracy level between 0.75 and 1.05 times the mean (0.97), so 12% of the marker reads fall
below 0.88 (0.85: 4.5%).

**The rowless species** (6 and 30 Mb, floor 0.85): their marker reads are mostly whole-gene reads that seeded and
aligned nowhere (14 of 19 whole-gene reads at 6 Mb, 13 of 26 at 30 Mb). Most of the rest aligned best to another
species (2 and 10 whole-gene reads; why was not looked at), or to their own species below MAPQ 4 (3 and 3 reads that
share 1,000 bases or more with a gene but do not hold it). A species with one marker read went rowless 11-13% of the
time at 6-30 Mb (9% at 150 Mb), about the share of marker reads lost.

## 4. Lowering the ONT floor

### Rows (the floor experiment)

The same reads profiled at three floors. Rowless: species with marker reads and no row, of the (sample, species)
pairs with marker reads; rows of present and absent taxa from the `.truth_annotated` profiles; time and peak memory
of the protal run (12 samples, 4 threads).

| depth | floor | marker reads own / unaligned | rowless species | present rows | absent rows | protal run |
|---|---|---|---:|---:|---:|---|
| 6 Mb | 0.85 | 88.4% / 7.3% | 24 of 418 (5.7%) | 429 | 12 | 1.0 s, 0.25 GB |
| | 0.80 | 92.8% / 2.2% | 14 (3.3%) | 451 | 15 | 1.3 s, 0.26 GB |
| | 0.75 | 94.2% / 0.4% | 12 (2.9%) | 451 | 17 | 1.1 s, 0.28 GB |
| 30 Mb | 0.85 | 88.5% / 7.1% | 28 of 719 (3.9%) | 749 | 25 | 5.0 s, 0.27 GB |
| | 0.80 | 93.1% / 2.1% | 20 (2.8%) | 769 | 36 | 4.5 s, 0.28 GB |
| | 0.75 | 94.6% / 0.3% | 16 (2.2%) | 773 | 39 | 3.5 s, 0.31 GB |
| 150 Mb | 0.85 | 87.9% / 7.4% | 9 of 944 (1.0%) | 973 | 52 | 12.7 s, 0.39 GB |
| | 0.80 | 92.8% / 2.1% | 6 (0.6%) | 984 | 73 | 13.8 s, 0.51 GB |
| | 0.75 | 94.2% / 0.4% | 4 (0.4%) | 986 | 77 | 16.0 s, 0.54 GB |

- 0.80 aligns 70% of the reads 0.85 loses and has 29-42% fewer rowless species; 0.75 aligns nearly all (43-56%
  fewer rowless species).
- What stays rowless at 0.75 is not the floor's. Of these species' 40 marker reads, at 6, 30 and 150 Mb:
  - 18 align best to another species, at 0.85-0.99 accuracy, 10 of them at MAPQ 26-115. Probably their species'
    copy of the gene is missing from the database or farther than a congener's; not checked.
  - 11 are on their own species below MAPQ 4, reads that hold only part of a gene.
  - 7 have no record.
  - 2 are unaligned.
  - 2 are on their own species at MAPQ 118 and still give no row. Perhaps a gene copy the profiler leaves out of the
    evidence (suspect or foreign copies); not checked.
- The cost: more rows of absent taxa (the held-out species' poor reads now align to their congeners): +3/+5 at 6 Mb,
  +11/+14 at 30 Mb, +21/+25 at 150 Mb. Whether they become false positives is the model's call (next section).
- The table's times are single runs while other jobs loaded the machine; [Cost](#cost) has the 150 Mb runs re-timed.

### Calls (three builds)

The trainer's calls at its knob (0.5 in all three builds), as protal calls (truth: every simulated species with a
row) and with the marker truth (present: the species with marker reads; unseen ones count as missed). Change against
0.85 with a paired bootstrap over samples ([compare.txt](compare.txt)):

| set | floor | as protal calls TP/FP/FN | F1 | marker truth TP/FP/FN | unseen | F1 | change (95% interval) |
|---|---|---|---:|---|---:|---:|---|
| species held out (120 samples) | 0.85 | 4219/82/94 | 0.9796 | 4057/82/145 | 126 | 0.9728 | |
| | 0.80 | 4321/105/133 | 0.9732 | 4109/105/94 | 72 | 0.9764 | +0.0036 (+0.0006, +0.0067) |
| | 0.75 | 4342/108/143 | 0.9719 | 4124/108/79 | 58 | 0.9778 | +0.0050 (+0.0023, +0.0081) |
| test (39 samples) | 0.85 | 1136/28/50 | 0.9668 | 1084/28/58 | 46 | 0.9618 | |
| | 0.80 | 1174/34/72 | 0.9568 | 1106/34/36 | 23 | 0.9693 | +0.0075 (−0.0007, +0.0163) |
| | 0.75 | 1188/30/62 | 0.9627 | 1114/30/28 | 20 | 0.9746 | +0.0128 (+0.0056, +0.0208) |

(The test set has 39 samples: one of the 40, at 0.15 Mb, has no row and no unseen species in any build.) 0.75
against 0.80: +0.0015 (−0.0008, +0.0037) with species held out, +0.0053 (−0.0036, +0.0150) on the test set.

By sample depth (marker truth TP/FP/FN, unseen in brackets, F1):

| set | depth | 0.85 | 0.80 | 0.75 |
|---|---|---|---|---|
| species held out | 0.3 Mb | 107/1/12 (12) 0.9427 | 111/2/8 (8) 0.9569 | 114/6/5 (5) 0.9540 |
| | 1.5 Mb | 275/10/18 (17) 0.9516 | 284/12/9 (8) 0.9643 | 286/14/7 (6) 0.9646 |
| | 6 Mb | 677/20/37 (34) 0.9596 | 691/24/24 (21) 0.9664 | 692/27/23 (20) 0.9651 |
| | 30 Mb | 1317/24/52 (48) 0.9720 | 1338/40/31 (26) 0.9742 | 1344/33/25 (21) 0.9789 |
| | 150 Mb | 1681/27/26 (15) 0.9845 | 1685/27/22 (9) 0.9857 | 1688/28/19 (6) 0.9863 |
| test | 0.15 Mb | 11/2/0 (0) 0.9167 | 11/2/0 (0) 0.9167 | 11/2/0 (0) 0.9167 |
| | 1 Mb | 73/0/11 (11) 0.9299 | 78/0/6 (6) 0.9630 | 78/0/6 (6) 0.9630 |
| | 5 Mb | 191/5/16 (16) 0.9479 | 197/5/10 (9) 0.9633 | 200/2/7 (6) 0.9780 |
| | 25 Mb | 248/9/18 (11) 0.9484 | 257/12/9 (6) 0.9607 | 258/8/8 (6) 0.9699 |
| | 250 Mb | 561/12/13 (8) 0.9782 | 563/15/11 (2) 0.9774 | 567/18/7 (2) 0.9784 |

- The marker truth gains at every depth but the shallowest and the deepest test points, most at 1-30 Mb, where a
  species has few marker reads. The species found again are the rowless ones: the misses with a row barely change
  (19, 22, 21 with species held out).
- F1 as protal calls falls by 0.004-0.010. Lower floors make rows for present species the simulator counts no marker
  read for (none shares 50 bases with a marker gene), and the model leaves most of them uncalled: misses as
  protal calls rise 94 → 133 → 143 while the marker truth's fall 145 → 94 → 79. This F1's truth changes with the
  floor; the marker truth is the same reads' in all three builds, so only it compares floors.
- The new false positives ([new_false_positives.py](new_false_positives.py)): 45 at 0.80 and 53 at 0.75 with species
  held out (22 and 27 of 0.85's gone), 10 and 9 on the test set (4 and 7 gone). Of the new ones, 42 of 45, 48 of 53,
  9 of 10 and 9 of 9 are in samples with a congener left out of the training database: a novel species' poor reads
  now reach its congener's gene. Half to most of them had a row at 0.85 already (20, 31, 8, 7), which the model now
  calls.

### Cost

protal on the floor experiment's 150 Mb samples (12 samples, 1.8 GB of reads, the small training database, 4
threads), three rounds in rotating order ([time_floors.sh](time_floors.sh)). Other jobs ran on the machine (load 4-5):

| floor | wall (3 rounds) | median | user CPU | peak memory |
|---|---|---:|---|---:|
| 0.85 | 13.7, 13.8, 13.7 s | 13.7 s | 49-52 s | 0.37 GiB |
| 0.80 | 16.1, 15.9, 16.5 s | 16.1 s (+17%) | 61-63 s (+22%) | 0.48 GiB |
| 0.75 | 17.5, 19.9, 22.9 s | 19.9 s (+45%) | 67-84 s (+52%) | 0.52 GiB |

The floor also sets the WFA budget (`MaxScore`: (1 − floor) × overlap mismatches), so every alignment that fails
runs longer, not only the reads that now pass. The memory grows too, probably WFA's buffers (not checked).
At GTDB scale, where seeding costs more, the share may be smaller; not measured.

## Recommendation

1. **Make 0.80 ONT's identity floor for databases built from now on** (`ReadType.h`). 0.8.5 records each read type's
   floor in `build_settings.tsv`, and runs follow it, so a database's model sees reads at the floor it was trained
   at. 0.80 gives most of the gain for +17% ONT profiling time: rowless species −43% with species held out and
   −50% on the test set, F1 +0.004 and +0.008. 0.75 gives a little more for +45%: −54% and −57%, +0.005 and +0.013
   (not significant against 0.80). It is the next step if ONT time at GTDB scale allows.
2. **Keep 0.85 where a database records no floor** (built before 0.8.5). `Options::GetMaxScoreAni` falls back to
   `ReadType.h`'s floor when the database records none. Without a fallback, an older database's model would see
   reads at a floor it was not trained at. `Options::FollowDatabaseRunSettings` should set the ONT floor to 0.85 when
   the database records none, as it sets mate guidance to `guide`. 0.8.5 databases record 0.85 already.
3. **Check it in the next r226 build:** ONT's unseen share (3.7% in v29; ~1.5-2% expected) and ONT's profiling time.
4. **Not tried:** an identity that counts a one-base indel once for ONT, so that the floor means the read's accuracy.
   It changes the identity the model's features see, and the WFA budget has to grow the same way, so the time cost
   would be about the same.
5. **The simulation sets the size of the effect.** pbsim3's ONT model draws read accuracies down to 0.73 (0.75 ×
   the mean 0.97), and 12% of marker reads are below 0.88. The 0.85 floor was set
   ([2026-09-29-pacbio-long-reads](../2026-09-29-pacbio-long-reads/README.md)) on reads with one error rate for all
   (1.6-5% errors, at most 5.4% differences per aligned base with the strains' divergence), not on such a spread.
   Real runs' share of such reads depends on chemistry, basecaller and read filtering (a mean-Q filter removes many
   of them); it was not measured here.

**Implemented** (1 and 2) in the commit that adds this report, on `audit-fixes` after `0deffcf`:
- `ReadType.h` gives ONT the floor 0.8, which `--build` records.
- `build_settings::DatabaseMaxScoreAni` gives a database that records no floors 0.9 and ONT 0.85. Its runs follow it
  (`Options::FollowDatabaseRunSettings`), the run-settings block shows it as the database's (no note), and
  `--profile_only` checks SAMs against it.
- Tests: `test_RunSettings.cpp`, `test_Parsing.cpp`, and the e2e ONT and run-settings tests, with a new one: a
  database without floors aligns ONT reads at 0.85, and a run told 0.8 is warned.
- Docs: `-a`'s help, [running.md](../../running.md), [databases.md](../../databases.md) and
  [version_changes.md](../../version_changes.md).

The website is out of date where it gives `-a`'s default as 0.85 for ONT.
