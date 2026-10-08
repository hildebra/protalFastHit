# The GTDB build's outputs: sorted, each file once

2026-10-08, branch `audit-fixes` at `e24dcff` plus this change. What `scripts/build_gtdb_database.py`
leaves in `--outdir`, taken from its code, from the r226 v17 build (SLURM job 24060119, `--scratch`,
`--share-logs`, `-t 84`, the default design: 4 read types, 5 scenarios; its share archive unpacked in
`local/v17`, its `console.log`), and from the mini pipeline test (`scripts/mini_db/test_gtdb_pipeline.py`).

## Summary

A build left about 11,900 files in `--outdir`. Three things made them:

1. **9,030 in-silico strain FASTAs** (`insilico_strains/`, one per species with a single genome), about 9 GB
   by estimate, on `--outdir`'s network file system. The simulations read them, and the error reports read
   their contig names.
2. **2,759 files of error reads** (`model_logs/error_reads/`): 1,480 tables of one sample each (19 MB in all),
   1,275 SAMs of one sample each (with `--share-logs`) and 4 summaries.
3. **About 80 files beside `protal_db/`**: 32 step logs, 40 outputs of the trainer (32 of them also copied
   into `model_logs/`), and the genome tables, taxonomy, gene ids and held-out species (some of them copied
   into `model_logs/` too). The build's own reports (`gene_congeners.tsv`, `gene_incongruence.tsv` of
   51 MB) sat in `protal_db/` and were copied into `model_logs/`.

Now `--outdir` holds four folders and the console. Each file is there once:

```
OUTDIR/
  protal_db/     database.protal, build_metadata.tsv, gene2geneid.tsv, genome2tiid.tsv      4 files
  console.log
  model_logs/    the evaluation and what the run chose                                       ~63 files
  logs/          every step's log                                                            ~30 files
  work/          what a rerun reuses: stages/, genome tables, taxonomy, gene ids, the
                 finished database's foreign-rates table, training/ and test/ tables         ~20 files
```

At r226 v17's settings that is about 118 files without `--share-logs`, and about 1,400 with it (the
per-sample SAMs of the error reads, 1,275 at v17, and the share archive). Without `--scratch`, `work/`
also holds the samples, the training database, the genome store and the in-silico strains, as
`--outdir` did before.

## What a build wrote, and where it goes now

Counts are r226 v17's (its `console.log`, the share archive, the code). "Copied" means the file was in two
places.

| Before (path in `--outdir`) | Files | Read by | Now |
|---|---|---|---|
| `*.log` (32) | 32 | people; `build_metadata.tsv` reads `index_and_package.log` | `logs/`. The foreign scan's `_add.log` is merged into the scan's log, the four `parity*.log` into `parity.log` (an appended command follows a `---` line) |
| `console.log` | 1 | people | stays at the top |
| `trained_model{,_se,_pb,_ont}.*`: `.xml`, `.joblib`, `.report.txt`, `.metrics.json`, `.thresholds.tsv`, `.varimp.tsv`, `.predictions`, `.test_predictions`, `.scenario_predictions`, `.calls.tsv.gz` | 40 | the parity check, `--add_model`, `error_reads.py`, `composition_accuracy.py`, the summary | written once, into `model_logs/` (the trainer's `--output-prefix`); the 32 copies are gone |
| `model_logs/` copies of `training_data*.log`, `test_data*.log`, `genome_table.txt`, `heldout_species.txt`, `species_clouds.tsv`, `build_metadata.tsv`, `gene_ranking.tsv`, `gene_subset.txt` | 4-10 | people, the share archive | the logs only in `logs/`; `build_metadata.tsv` only in `protal_db/`; the others written into `model_logs/` once, and read there by the steps that need them |
| `parity.txt`, `parity_{se,pb,ont}.txt` in `model_logs/` | 4 | people | one `model_logs/parity.txt`; each check's text names its model and read type |
| `genomes.tsv`, `genomes_simulated.tsv`, `internal_taxonomy.dmp`, `gene2geneid.tsv` | 4 | the steps; a rerun (the genome lengths, about 8 minutes at r226) | `work/` |
| `.stages/` | ~6 | a rerun | `work/stages/` |
| `.converted/` (a gene subset's run, removed at its end) | | the run | `work/converted/` |
| `insilico_strains/` | 9,031 | the simulations; the error reports (contig names) | the samples' disk: `SCRATCH/insilico_strains`, `work/insilico_strains` without `--scratch`. A new `--scratch` makes them again, 2.5 minutes at r226. The strains are seeded, so they are the same strains, and the collector's sample keys hash the genome table's content |
| `training/`, `test/`: 4 tables each, and with `--scratch` `protal_runs/*.log` (4 and 2 runs at v17) | 14 | the trainer; people (stage timers) | `work/training/`, `work/test/`; the runs' logs one after the other in `logs/protal_runs_<collection>.log`, each after a `==> <run> <==` line |
| `protal_db/gene_congeners.tsv`, `gene_incongruence.tsv` | 2 (+2 copies) | people, `make_gene_rates.py` | moved into `model_logs/` at the end |
| `protal_db/foreign_rates.tsv` (the finished database's scan, also stored in it) | 1 | a rerun's stage check | `work/foreign_rates/foreign_rates.tsv` (`--add_tables` takes a table by its file name). The training database's table stays beside it on the samples' disk |
| `model_logs/error_reads/<type>/<set>/<point>/<sample>.taxa.tsv` | 1,480 | people (the r226 analyses) | one `model_logs/error_reads/<type>/taxa.tsv.gz` per read type: the rows of every sample, with the sample's design point and scenario after its set (`error_reads.py`) |
| `model_logs/error_reads/.../<sample>.FP.sam.zst`, `.FN.sam.zst` (`--share-logs`) | 1,275 | `ancestry_sites.py`, people | unchanged: one per sample and kind, and only with `--share-logs` |
| `<name>_share.tar.gz` (`--share-logs`) | 1 | people, off the cluster | laid out as `--outdir`: `console.log`, `logs/`, `model_logs/` without the models (`.xml`, `.joblib`), `protal_db/build_metadata.tsv`, `work/internal_taxonomy.dmp` (v17's analysis had to borrow v15's taxonomy), the shortened tables under `work/training/`, `work/test/` |

The run's last console line counts the files in each folder, for example `In OUT: protal_db/ the
database (4 files), model_logs/ the evaluation (63), logs/ and console.log every step's log (30),
work/ what a rerun reuses (20; not needed to use the database)`.

## Choices

- **Each file once.** The trainer writes into `model_logs/` directly. Inputs that also tell what the run
  chose (`heldout_species.txt`, `species_clouds.tsv`, `gene_subset.txt`, `gene_ranking.tsv`,
  `genome_table.txt`) are written there and read from there. A file given on the command line that is
  already that file (`--holdout-species OUTDIR/model_logs/heldout_species.txt` on a rerun) is not copied
  onto itself (`copy_file`).
- **`work/` holds the training and test tables**, about 1.2 GB at r226, which retraining needs. Nothing
  in it is needed to use the database. Removing it costs a rerun its conversion, its genome lengths and
  its tables.
- **The in-silico strains moved to the samples' disk.** They were 9,030 small files written to the
  network file system, the kind of small-file churn that stalled r226's collection there. They are
  inputs to the simulations, like the genome store beside them. Their stage key now includes their
  folder, and they count as made only while their `insilico_strains.tsv` is there.
- **The error reads' tables are merged.** `error_reads.py` still writes a table per sample in its
  workers, which run in parallel and in separate processes. It then merges them into `taxa.tsv.gz` and
  removes them, along with the folders left empty. Without `--share-logs` an `error_reads/<type>/` folder
  holds two files.
- **Not changed.** With `--share-logs` the SAMs stay one file per sample and kind (1,275 at v17):
  `ancestry_sites.py` and the r226 analyses read them by sample. Merging them into one SAM per read
  type and kind would need a sample tag on every record and a merged header. With `--scratch` the
  samples stay on the scratch disk for a rerun, as before.
- **No migration.** A rerun into an `--outdir` of the old layout does its steps again (with `--scratch` on
  a new node it simulated and profiled again anyway). Moving an old folder's files would not be enough:
  `genomes_simulated.tsv`, the collector's manifests and maps hold absolute paths into it.
  `build_gtdb_releases.py` reads a build's taxonomy and gene ranking at either layout.

## Checks

- Unit tests (WSL, Python 3.14, the `protal-db-build` environment, 2 cores), in the checkout and again
  on the committed tree: `python3 -m unittest scripts/test_error_reads.py` (11 tests, `taxa.tsv.gz` and
  no folder per sample with `--sams none`) and `test_gtdb_build.BuildOptionsTest`, `BuildFunctionsTest`,
  `CladeHoldoutTest`, `BuildReleasesTest` (15, among them a step's commands in one log, `copy_file` and
  the share archive's members) pass.
- End to end (`scripts/mini_db/test_gtdb_pipeline.py`, `GtdbBuildTest`, PROTAL and SIMULATE a build of
  `e24dcff`'s `src/`, 4 cores): see [below](#end-to-end).

## End to end

`scripts/mini_db/test_gtdb_pipeline.py` (`GtdbBuildTest`, 6 tests: the build with every gene, its rerun
that builds nothing, a reduced database, a failing background build, another seed stopped by SIGTERM,
scenarios with the error reads and the share archive, samples profiled as simulated, streamed samples)
passes in 198 s on 4 cores. It ran on a `git archive` of the committed tree in WSL, against a
snapshot of a protal build of `e24dcff`'s `src/`.

A first run from the shared checkout failed `test_a` (the rerun built everything again) and `test_c`
(the finished database rebuilt). Another session edited `gtdb_to_protal_db.py` there during the run,
which changed the conversion's stage key. Run on the committed tree alone, both pass.

The mini build's `--outdir` (60 species, pe and se, `--scratch`, `--rank-genes`) holds `console.log`,
`protal_db/` (4 files), `model_logs/` (32), `logs/` (21) and `work/` (17). The test now checks:

- the top level and `protal_db/` exactly;
- the trainer's files in `model_logs/`, with no logs or metadata copies there;
- one `parity.txt` covering every read type;
- the in-silico strains on the scratch disk, not in `work/`;
- `--add_tables` in the scan's log;
- `protal_runs_training.log`;
- the error reads' `taxa.tsv.gz` with each sample's point and scenario, and no per-sample tables;
- the share archive's members in the new layout, without the models or the database.
