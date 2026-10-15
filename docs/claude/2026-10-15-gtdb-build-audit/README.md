# The GTDB build script after the r226 v28 crash: an audit for bugs

Date: 2026-10-15. Branch `audit-fixes` at `b42cbc3`. That commit fixes the bug the r226 v28 build crashed on: the parity
check took the first and last design point by name, which since `ce2940e` may have no single-end samples. Nothing else
was changed for this audit.

**Question.** After v28 crashed at step 9 (`no samples with a SAM and a training dump in .../rl100_p1000,
.../rl75_p500000`), the user asked for an audit of `scripts/build_gtdb_database.py` to find other bugs.

**Scope.** `scripts/build_gtdb_database.py` (4905 lines) and what it runs, as far as the build depends on it:
- `collect_training_data.py`, `check_model_parity.py`, `composition_accuracy.py`, `error_reads.py`,
  `trace_relatives.py` and `ancestry_sites.py`;
- what `machine_learning_cmdline.py` writes that the summaries read.

Not in scope: protal itself, the simulator, and the trainer's own logic. An earlier audit,
[2026-10-07-database-build-audit](../2026-10-07-database-build-audit.md), covered the conversion and protal's builds;
its fixes (L1-L3 in `54b1aa1`) still hold.

**Method.**
- **Four parts read in parallel**, each followed into the helper scripts:
  - training to the end, with the reports;
  - the collections, disk and streaming;
  - the early steps and the resume logic;
  - the summaries and gates.
- **Each finding re-checked** against the code path before it went in here. CONFIRMED means read end to end or
  reproduced; PLAUSIBLE means the mechanism is in the code but its trigger was not reproduced.
- **v28's console** (SLURM 24248584, 120 threads, `/hpc-home/hildebra/DB/protal/protal0.8.4_r226_v28`), as the user
  pasted it, gave the configuration and the numbers quoted below.
- **Probes:**
  - full mini builds on the pipeline test's synthetic release (60 species), WSL, 4 cores;
  - protal and the simulator built from a `git archive` of `b42cbc3`;
  - scripts in [scripts/](scripts/) (see Probes).

## Summary

**The new single-end layout no longer crashes the build.** With `b42cbc3`, a mini build with v28's layout runs every
step to the end:
- paired-end-only, both and single-end-only points;
- the design's and the scenarios' Ultima units;
- pb and ont;
- both model kinds and the models at other gates;
- `--test-c-gates`, streamed simulations, error reads and the share archive.

The parity check now picks `rl100`/`rl150` for pe and `rl150`/`rl250` for se. Rerunning into the same OUTDIR after a crash at parity, with a new empty scratch
folder, keeps `protal_db` and finishes.

**Two bugs change what the next r226 build reports.** Both are small to fix:

1. **G1, the rerun after v28 loses the marker-truth evaluation.** Once `protal_db` is kept, the collections get no
   `--gene_positions`. The tables then lack `meta_marker_reads`, there are no `.unseen.tsv` tables and no "Marker
   truth" sections, and the line "F1 as protal calls, with every species simulated -> with only those whose marker
   genes got reads" disappears. The probe of exactly this rerun has 0 of 244 rows with marker reads, against every row
   in the first run. The models are not affected.
2. **G2, `c_gates.tsv` and its console lines compare the gates unfairly.** The tables at other gates never get their
   `.unseen.tsv`: it stays `.unseen.tsv.partial`. So every row except the database gate's leaves out the present
   species a gate drops entirely, which was the fix of 2026-10-13. This is why the database's gate is the worst
   row of every read type in v28:
   - pe: 4 stride 0.9419, against 0.9501 at 2 and 0.9511 at 8;
   - ont: 8 stride 0.9250, against 0.9597 at 16;
   - pb: 16 hash 0.9632 against 16 stride 0.9537, with the same model and samples.

   Do not choose gates from v28's c_gates lines.

**G3 makes the CPU column of c_gates misleading.** The rows at the models' gates also cover the streamed (deepest)
samples, so their F1 and CPU do not compare with the other rows. That explains v28's pb 8 stride 6326 CPU s against
2420 ungated, and ont 1 9158 against 2636.

**The rest are a possible stop on a tight disk (G4), stale results after an interrupted run (G5), signal and logging
issues, and small things.** Table below, details after.

## Findings

| ID | Severity | Where | What | Status |
|---|---|---|---|---|
| G1 | wrong result (evaluation missing) | `build_gtdb_database.py:4125-4127`, `4184`, `4193-4212` | A rerun that keeps `protal_db`, with a new scratch, collects without `--gene_positions`: no marker truth, no unseen tables, no warning | CONFIRMED, probe |
| G2 | wrong result | `collect_training_data.py:2116-2128` | A gate's table returns before `os.replace(unseen_partial, unseen_path)`: every `c_gates/c<C>/*.unseen.tsv` stays `.partial`, the trainer counts 0 unseen at every gate but the database's | CONFIRMED, probe |
| G3 | misleading log | `collect_training_data.py:1222-1225`; `build_gtdb_database.py:2157-2174` | The model gates' rows include the streamed samples (profiled again only at `--model_levels`); every other row and the database row have the written samples only | CONFIRMED, reading and v28's console |
| G4 | stop on a tight disk | `build_gtdb_database.py:4471` against `4486-4487` | The check that keeps the test set's Ultima reads on disk ignores `streamed_levels`, which the threshold subtracts next: a room between the two stops the run "even with every simulation streamed" | CONFIRMED, reading |
| G5 | stale result after an interruption | `build_gtdb_database.py:3670-3673`, `3683-3686`, `3763-3769`, `3792-3798` | A new conversion into a folder starts without `stages.forget("convert")`; killed midway, the next run takes the folder as converted (no full reference, no neighbours, uniform in-silico ANI) | CONFIRMED, reading; needs a kill in a minutes-long window |
| G6 | stop outside a batch job | `build_gtdb_database.py:4900-4901` | The SIGHUP handler replaces nohup's ignore: a `nohup ... &` run stops when the terminal goes | CONFIRMED, probe |
| G7 | misleading log | `build_gtdb_database.py:2375` | `contextlib.redirect_stdout` is process-wide: the reports thread's lines said meanwhile go into `composition_accuracy.txt` and `summary.txt` | CONFIRMED, probe |
| G8 | minor | `build_gtdb_database.py:1574`, `1632`, `1681` | The reports' commands are plain `Popen`s, not `Job`s: a build that stops leaves them running | CONFIRMED, reading |
| G9 | misleading log | `build_gtdb_database.py:4742-4753` | The marker-truth, c_gates, feature-choice and scenario lines describe the first kind's models (forest) and follow the gbm lines unlabelled | CONFIRMED, reading |
| G10 | wrong result | `build_gtdb_database.py:3395` against `3740` | `--insilico-strains 0` given explicitly also drops the `--insilico-multi` strains (default 0.25); `auto` resolving to 0 keeps them | CONFIRMED, reading |
| G11 | minor | `build_gtdb_database.py:4460-4462` | The ancestry report's reference (a hard link already on the disk) is counted again in the disk needs: ~14.7 GB too little room at r226, more streamed than needed | CONFIRMED, reading |
| G12 | minor | `collect_training_data.py:2269-2273` | A scenario point's disk cost for the scheduler leaves out the host reads, which the same simulator run writes (host scenario: ~10x) | CONFIRMED, reading |
| G13 | minor | `build_gtdb_database.py:1674`; `ancestry_sites.py:91` | `ancestry_sites.py` requires `--heldout`, which the build leaves out without a hold-out: with `--share-logs`, every read type's ancestry report "failed (2)" | CONFIRMED, probe |
| G14 | minor | `error_reads.py:637-704` | Its output folders are never emptied: a rerun keeps an earlier run's per-sample SAMs and `.partial` files, which the ancestry report and the archive read | PLAUSIBLE |
| G15 | minor | `build_gtdb_database.py:2085-2102`; `src/Utilities/CpuLog.h` | CPU seconds looked up by sample name, which training and test scenario samples share; and `aligning_cpu_s` includes the profiling of the previous sample (`--profile-ahead`) | name clash CONFIRMED, wrong value PLAUSIBLE |
| G16 | latent | `composition_accuracy.py:161-163`, `244-248` | A drawn sample whose community point is single-end only has no pe composition: host share NaN, under "host share not known". Not with v28's setups (its se setups are pe setups too) | CONFIRMED, probe |
| G17 | latent | `collect_training_data.py:888-926`, `2920-2921` | Points that no unit profiles (e.g. `--read_types se,pb,ont`) are simulated for their communities and their reads never removed nor counted on the disk | CONFIRMED, probe |
| G18 | reruns on a kept samples disk | `build_gtdb_database.py:4342`, `4457`; `collect_training_data.py:2119-2122`, `2820-2825` | The outputs already on the disk counted again; `stream_design` replaces a profiled point's `sim/` (marker counts lost); a gate without samples keeps an earlier run's table | PLAUSIBLE |
| G19 | hang | `collect_training_data.py:2792-2803`, `1626-1632`, `2428-2433` | With a fixed `--stream-above`, the follower waits for 20 GB of ready reads while the simulations wait for room | PLAUSIBLE; not under `auto` |
| G20 | stale result | `build_gtdb_database.py:4377-4389` | The finished database's foreign-rates stage is keyed on `final_key`, its table stored only by `--add_tables`: a rebuild under an old key skips the scan and the table | CONFIRMED, reading; `--foreign-rates` only |
| G21 | coverage | `check_model_parity.py:95-105`; `build_gtdb_database.py:4774-4776` | Parity never re-profiles an Ultima sample (no `sim/protal.meta`) nor checks a model of another gate | CONFIRMED, probe |
| G22 | minor | see below | Step count with `--foreign-rates`; argument checks; the docstring's `--model` default; the database row of c_gates depends on `levels[0]` | CONFIRMED, reading |

## Details

### G1: a rerun that keeps protal_db loses the marker truth

How the collections get `gene_positions.tsv`:
- `collect_command` passes `--gene_positions` only if `marker_positions()` finds `gene_positions.tsv` in `converted`,
  `db` or `full`, or already at `SCRATCH/gene_positions.tsv` (4193-4212).
- On a first run, `protal_db` holds the converted release until its build in the background packs it; the collections
  start long before that.
- On a rerun that keeps `protal_db`, its copy is packed. The release is converted again into `work/converted` (L2 of
  the earlier audit, for the in-silico strains). That folder is removed at 4127, before the first `collect_command`.
  A new `--scratch` has no earlier copy, so `marker_positions()` returns None. Nothing says so.

Consequences:
- The collectors run without `--gene_positions`: empty `meta_marker_reads`, no `<table>.unseen.tsv`.
- The trainers write no "Marker truth" section, and `marker_scores` prints nothing.
- The database row of c_gates, which G2 leaves as the only one counting unseen species, also loses them.

Probe (`crash_rerun.sh`): the rerun's console has 0 "Marker truth" lines (the first run 24), and its training table
has 0 of 244 rows with marker reads.

Fix: take `positions = marker_positions()` before `shutil.rmtree(converted)` (or keep a copy in `work/`).

Without a code change, a copy of the release's `gene_positions.tsv` at `SCRATCH/gene_positions.tsv` before the run is
used. It can be unpacked from `protal_db/database.protal` with `protal --unpack_db`, or taken from a conversion.

### G2: the gate tables' unseen species

What happens:
- `write_table(..., gate=C)` writes `c_gates/c<C>/<table>.unseen.tsv.partial` beside the gate's table. Since 2026-10-13
  the comment at 2019-2020 says it should.
- It then returns at 2128 without the rename that the main table gets at 2131. A gate without samples also leaves the
  `.partial` (2119-2122).
- The trainer's `--unseen-file auto` looks for `<table>.unseen.tsv` (`machine_learning_cmdline.py:1561`), finds none,
  and counts 0.

Which numbers are affected:
- Every c_gates row other than the database's, and the marker-truth sections of the gate models. They leave out the
  present species a gate drops entirely.
- The database row and every `baseline_F1` read `test/<table>.unseen.tsv`, which exists, and count them.

In the v28-like probe every gate folder holds only `*.unseen.tsv.partial`, and `c_gates.tsv`'s `unseen` column is 0 on
every row.

v28's numbers this explains:
- the database gate is the lowest of each read type;
- the same model on the same samples gives pe 4 hash 0.9506 against 4 stride 0.9419, and pb 16 hash 0.9632 against
  0.9537;
- the gap is about the marker-truth drop of the test set (pe 0.963 -> 0.953).

The 2026-10-13 report ([2026-10-13-r226-v27-c-gates](../2026-10-13-r226-v27-c-gates.md)) described the flaw that
unseen counting was to fix; G2 means the fix never took effect at the other gates.

Fix: `os.replace(unseen_partial, unseen_path)` before the gate branch returns, and remove the `.partial` where it
returns for want of samples. `scripts/mini_db/test_gtdb_pipeline.py` test_j could check that `unseen` is not 0 on a gate
row.

### G3: rows over different samples in c_gates

Why the sample sets differ:
- Streamed simulations are streamed again only for `--model_levels` (`unit_streamed_level_ids`).
- So `c_gates/c<C>/samples_<t>.txt` at a model gate (pe and se 1 and 8, pb 8 and 32, ont 1 and 16) lists written and
  streamed samples.
- At 0, 2, 16, hash, grid and `--flex_indels` it lists the written samples only, and so does the database row, which
  takes `levels[0]`.

The streamed samples are the deepest (v28: the 3-6 Gbp long-read samples, the 30 M-pair design points and the 20
M-pair scenario samples). The model gates'
rows therefore have more CPU and other F1 than their neighbours:
- pb: 8 stride 6326 CPU s against 2420 ungated;
- ont: 1 9158 against 2636;
- pe: 1 4941 and 8 2882 against 4 2077.

The console line prints neither the baseline nor the sample count. A grid row at a model gate "(Y by stride)" compares
written samples with all samples.

Fix: report the model gates' rows over the written samples too (or give the streamed ones a row of their own), and
put the sample count in the console line.

### G4: a spurious "Not enough room" stop

- At 4471, the test set's Ultima simulations (`--write_simulations`, for `--test-flex-indels`) stay on disk unless
  `room - pinned - pinned_levels < 0`.
- The threshold at 4486-4487 then gets `room - pinned - pinned_levels - streamed_levels`.
- When the room falls between the two, the Ultima reads stay on disk, `stream_threshold` returns None, and the run stops.
  Streaming the Ultima simulations instead would have let it run.

v28 had room to spare (766 GB free against ~483 GB needed besides the reads). A smaller node disk could hit this.

Fix: subtract `streamed_levels` in the check at 4471 too.

### G5: an interrupted reconversion taken as converted

How it happens:
- `converted_before(folder)` is the convert stage's key plus a few files.
- A new conversion into the same folder (after a change to `final_base` or the genome table, 3683-3686 and
  3792-3798; or into `work/converted`, 3763-3769) does not forget the stage first.
- The converter writes `internal_taxonomy.dmp`, `reference.fna` and `reference.map` in place; only the full reference
  is a `.partial`.

Scenario: a run killed during that reconversion (Ctrl-C, a time limit, a full disk) leaves the old stage key with
half-written files. The next run takes them:
- `clear_build_outputs` drops the full reference's `.partial`;
- the stale `gene_neighbours` stage skips the neighbours;
- both databases are built without a full reference (no strain alleles, no gene conservation) and without neighbours;
- the in-silico strains fall back to the uniform `97-99.5` ANI.

Fix: `stages.forget("convert")` before every `start_convert`/`convert`, and `converted_before` requiring the full
reference when the release has `genomic_files_all`.

### G6: SIGHUP and nohup

`signal.signal(SIGHUP, on_signal)` at 4900 replaces the ignore that `nohup` set; a probe in WSL showed the disposition
ignore before and the handler after. A build started with `nohup ... &` on a workstation or login node stops ("Stopped
by SIGHUP", with everything it started) when the terminal or SSH session goes. Batch jobs are not affected. Fix:
install the SIGHUP handler only when `signal.getsignal(SIGHUP)` is not `SIG_IGN`.

### G7: composition output takes the reports' lines

`composition_report` runs `composition_accuracy.main` under `contextlib.redirect_stdout(StringIO)`, which swaps
`sys.stdout` for the whole process. The reports thread (`report_all`, started before the parity check) prints through
`say()` meanwhile, so its lines land in the StringIO:
- the line that starts "Reports of what the models' errors rest on";
- the contig names;
- the error reads' results.

That text becomes `model_logs/composition_accuracy.txt`, is appended to `summary.txt`, and reaches the console only
inside the final summary. `console.log` has them in place. Fix: have `composition_accuracy.main` write to a given
stream, or run it as a subprocess.

### G8: the reports outlive a stopped build

`reports()` starts `error_reads.py` (one per read type), `trace_relatives.py`, `protal --unpack_db` and
`ancestry_sites.py` with `subprocess.Popen`, and the thread is a daemon.

When the main thread stops (a failed Job, as v28 at step 9, after the thread had started):
- `stop_jobs()` does not know these commands;
- the thread dies at exit without its `finally`.

Under SLURM the job's cgroup ends them. On a workstation they keep writing `model_logs/error_reads/` and
`work/ancestry_files` while a rerun starts. Fix: register them as `Job`s, or have `stop()` kill a list the thread
keeps.

### G9: lines of one kind, unlabelled

With `--model forest,gbm`, the lines after training are printed after the "gbm:" lines with no label, but they are
the forest's (the first kind, `prefixes`):
- "F1 as protal calls, with every species simulated -> ...";
- "sparse-first gates (c_gates.tsv)";
- the feature choices;
- the scenario F1.

The gbm models' extra tests are computed but never summarized. Fix: label them with the kind, or print both.

### G10 to G22, briefly

- **G10:** `insilico_wanted` at 3395 ignores `--insilico-multi` when `--insilico-strains 0` is given; at 3740 `auto`
  resolving to 0 keeps the multi strains. The step total (3577) and `make_clouds` (3576) use the value before `auto`
  is resolved, and the provenance note says "none".
- **G11:** `ancestry_reference.fna` is linked at 4124, before `free` is measured at 4456, and added to `needs` at
  4462 as well (v28: "the ancestry report's reference.fna ~14.7 GB").
- **G12:** `pe_bytes` says a host's reads "are added by jobs of their own", but `scenario_command` writes them in the
  same run (1374-1375). `sample_bytes` and `simulation_bytes` count them, so `auto` is right. The scheduler's room
  check per job is not.
- **G13:** `ancestry_sites.py --heldout` is `required=True`. The build passes it only with a hold-out file, so
  `--share-logs` builds without species held out report every ancestry run as failed (exit 2). They have also
  unpacked the training database for nothing.
- **G14:** `error_reads.py` writes into `model_logs/error_reads/<type>/` without emptying it. `*.sam*` globs match
  `.partial` files left by killed workers, and `ancestry_sites.open_text` reads a name without `.zst` as text.
- **G15:**
  - `aligning_times` keys CPU by sample name over both collections' tables. The scenario samples
    `sc_<scenario>_..._s_<n>` exist in both. The test table is read last and wins, but a test sample without a row
    would take the training sample's CPU.
  - `aligning_cpu_s` is the process's CPU since the previous mark, so with `--profile-ahead` it includes the previous
    sample's profiling. The first sample of each level block carries the previous level's last profile.
- **G16:** `composition_accuracy.sample_files` takes a drawn sample's host share from its community's paired-end
  `.composition`. With a single-end setup that `--read-setups` lacks, the long-read and Ultima units replay
  communities whose pe samples are never profiled. A fallback to `protal_se/profiles/<community>_se.profile.composition`
  would serve.
- **G17:** `design_points` ignores `--read_types`, and long-read and Ultima units take every point of a depth as a
  community. With `--read_types se,pb,ont`, the rl75/rl100 points are simulated with reads (probe: 22 points) and
  never profiled, removed or counted on the disk.
- **G18:** on a rerun with the same samples disk:
  - `disk_needs` counts the SAMs, profiles and gate outputs again, although `free` already excludes them;
  - `follow()` redesigns every pe point newly streamed (`stream_design`, `rmtree` of its `sim/`) even when profiled,
    which loses its marker counts;
  - a gate table whose samples are all streamed keeps an earlier run's table and `samples_<t>.txt`.

  A new node-local scratch, as on the cluster, avoids all three.
- **G19:** the follower starts a protal run only once `--profile_block` GB of reads are ready or the simulations
  have ended. Long-read and Ultima reads are not ready until every community point of their depth is simulated, and the
  simulations wait for room. With a fixed `--stream-above` and too little room both can wait forever.
- **G20:** run A with `--foreign-rates`, run B with other options dies while rebuilding `protal_db`, run C with A's
  options rebuilds it. The scan counts as done ("kept"), and `database.protal` has no `foreign_rates.tsv`.
- **G21:** single-end points need `sim/protal.meta`, so `se_ultima_r*` is never chosen. The se model is never parity
  checked on an Ultima sample, and the models at other gates are never checked.
- **G22:**
  - `Steps.total` counts one foreign-rates step, but `foreign_rates()` starts one per database (13/12 on the console).
  - Arguments that are not checked:
    - `--read-types pe,se,pe` passes, and two trainers then write the same prefix;
    - `--holdout` above 1, or below 0 with clades, raises in `rng.sample`;
    - a typo in `--depth-knob-read-types` gives no knob curve silently.
  - The module docstring says `--model` defaults to "gbm,forest" and "Four scenarios"; the default is `forest,gbm`
    and five scenarios.
  - The c_gates database row takes the samples of `levels[0]`. If that is a model gate, it covers the streamed
    samples too; if it is se `0f`, the Ultima samples only.

## What the rerun of v28 will do

Same command, same OUTDIR, a new cluster job and so an empty `--scratch`. From the code (the review of the resume
logic) and the crash-and-rerun probe:

**Kept:**
- `protal_db` ("built by an earlier run from the same release and protal; kept"), as long as the protal binary is
  the same file: `file_identity` is its path, size and mtime;
- the species clouds;
- the hold-out, which is the same: `choose_holdout` is deterministic for a seed.

**Done again:**
- the release, converted again into `work/converted` for its gene positions (a few minutes);
- the in-silico strains, whose stage key names the scratch folder;
- the training database, whose stamp is on the old scratch;
- both collections and all the training (no stage).

**Added at the end:** `--add_model` puts the models into the kept database in place, since its placeholders match the
same key.

**Lost, unless fixed or worked around:** the marker truth (G1).

The c_gates lines stay biased (G2, G3), as in v28.

## Checked and found fine

- **The helpers find every sample of the new layout:**
  - `error_reads.sample_files` and `composition_accuracy.sample_files` cover pe-only points, points with both,
    `<point>_se`, `se_ultima_r*`, `sc_*_se_ultima_r*`, pb and ont;
  - `trace_relatives` reads `<point>/protal/alignments` only and skips `sc_*`;
  - no other code takes a read type from a point name's prefix.
- **Units:** in a v28-like `units_of` probe:
  - no unit names, profile folders or level folders collide;
  - Ultima units get their depth's sample count, and `community_of` resolves through `_s_`.
- **Streamed samples** keep their SAMs, profiles, dumps, truth files, `sim/protal.meta` and manifests. Nothing after
  training needs reads.
- **Reads are removed** only once every level has profiled them.
- **The threshold arithmetic** `r + c - s`, less `streamed_levels`, matches the collector's comparison.
- **Turn-taking** holds the lock only around protal runs.
- **Kinds and gates:**
  - the trainers, `gate_model_scores`, `summary_lines`, parity (every kind x read type) and the `--add_model` list
    (every kind x read type, plus every extra gate) are complete;
  - `nearest_gate` matches protal's `NearestGate`, ties included;
  - the database row and the levels are scored by the first kind's own models.
- **Stages** are marked only on success. `stop_jobs` records only jobs that exited 0, and `kill()` continues a
  paused build before it is ended.
- **The earlier audit's L1-L3** still hold.
- **The full references:**
  - the training database's full reference stays until after `reporting.join()`;
  - `ANCESTRY_REFERENCE` is removed only after it.

## Probes

Setup (WSL Ubuntu 24.04, `taskset -c 0-3 nice -n 5`):

```bash
git archive -o head.tar b42cbc3                    # from the checkout
mkdir -p ~/bgaudit/src && tar -xf head.tar -C ~/bgaudit/src
cmake -S ~/bgaudit/src -B ~/bgaudit/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF
ninja -C ~/bgaudit/build -j4 protal simulate_metagenomes
```

**[scripts/probe_build.py](scripts/probe_build.py)** runs `build_gtdb_database.py` as the pipeline test does. It uses
the test's synthetic release and its small scenarios, with long reads added to gut and host. The options follow v28:
- `--read-types pe,se,pb,ont`;
- pe setups 100 and 150 bp, se setups 150 and 250 bp, so the points are pe only, both, and se only;
- Ultima units;
- `--model forest,gbm`, `--model-gates pe:1/8,se:1/8,pb:8/32,ont:1/16` and `--test-c-gates`;
- the scenarios, with `--profile-blocks 0.00001` and `--share-logs`.

[scripts/run_probe.sh](scripts/run_probe.sh) runs it.

Results:
- **`run_probe.sh v28like`:** exit 0 in 1 min.
  - Parity passes for all eight models, on `rl100_p1000`/`rl150_p4000` (pe) and `rl150_p1000`/`rl250_p4000` (se).
  - Every gate folder has only `*.unseen.tsv.partial` (G2).
- **[scripts/crash_rerun.sh](scripts/crash_rerun.sh) `<checkout>`:** the first run uses the pre-fix parity check,
  delayed 90 s so that `protal_db` finishes first, with setups 75/100 bp pe and 150 bp se and `--stream-above 0.0005`
  (6 of 40 simulations streamed). It stops at `parity_se` as v28 did.
- **The rerun** with the fixed scripts and a new scratch folder:
  - exits 0;
  - keeps `protal_db`;
  - converts the release again;
  - has no marker truth (G1).
- **The other agents' probes** (G6 SIGHUP, G13 `--heldout`, G7 redirect, the units of a v28-like configuration, and a
  fake collection for the helpers) were small Python and shell scripts. They were not kept; each is described in its
  finding.

## Not done

- No fix was made beyond `b42cbc3`.
- Nothing was timed.
- No r226-scale run was made.
- The website is not affected by any finding.
