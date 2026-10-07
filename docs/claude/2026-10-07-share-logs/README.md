# What a build's logs hold, and a lean archive to share (--share-logs)

**Data.** The r226 v15 build's two archives as the user packed them (`local/protal_r226_v15_logs.tgz`, 2.4 GB and cut
off; `local/protal_r226_v15_training.tgz`, 443 MB), extracted to `local/v15`, and the Nanopore error reads that arrived
(WSL `~/v15/model_logs/error_reads/ont`, 264 samples). The v15 report
([2026-10-07-r226-v15](../2026-10-07-r226-v15/README.md)) has the sizes of the parts that did not arrive. Code: the
working tree on `2e16c68`. The measurements ran in WSL on one niced core (`nice -n 19 taskset -c 0`), because the
machine was busy.

## Where the bytes were

| Part | v15 | What it is |
|---|---|---|
| `model_logs/error_reads/` | se 4.3 GB, PacBio 2.4 GB, Nanopore 2.2 GB, pe failed | the SAM records of the reads behind each model's errors in every sample (`652ed53`) |
| `training/`, `test/` tables | 1.3 GB, 390 MB gzipped | the training data, numbers written as Python's repr (17 digits) |
| `trained_model*.predictions.tsv.gz` and the other trainer outputs | 45 MB | packed twice: they are copied into `model_logs/` too |
| `internal_taxonomy.dmp` | 11 MB | no report has used it |
| the rest of `model_logs/`, the logs | ~15 MB | the reports, metrics, calls, predictions; the logs are 0.5 MB at most |

The v14 logs archive was 76 MB; the error reads made v15's 2.4 GB, and the full archive would have been about 9 GB.

**The error reads are mostly reads of species nobody looks at.** Of the Nanopore samples' error-read records (summary
tables of `error_reads.py`, sum per set and scenario):

| Set, scenario | Samples | FP | FN | Unseen | Fragments | MB |
|---|---|---|---|---|---|---|
| test design | 42 | 52 | 67 | 4,951 | 545,894 | 94 |
| test soil | 3 | 530 | 599 | 18,127 | 1,995,577 | 453 |
| test soil_shallow | 3 | 432 | 503 | 20,950 | 550,608 | 138 |
| training design | 186 | 190 | 169 | 13,828 | 1,337,416 | 223 |
| training soil | 6 | 1,124 | 1,105 | 38,831 | 4,139,856 | 963 |
| training soil_shallow | 6 | 895 | 867 | 43,582 | 1,091,837 | 264 |

"Unseen" (species in the training database that protal profiled no reads to) outnumber the errors 20 to 1, and most
are species held out of the training database (95% in soil, 82% in shallow soil, half in the design):
`genome2tiid.tsv` keeps them without genes, so no read of theirs can align to them (v15 report). Per kind of reason, in one deep soil sample
(`test/sc_soil_ont_b6000000000_s_2`, 181 MB of `.sam.zst`; [caps.pl](caps.pl), uncompressed bytes):

| Kind | Fragments, all | MB, all | Fragments, 20 per taxon and reason | MB | Fragments, 50 | MB |
|---|---|---|---|---|---|---|
| FP (on the taxon) | 906 | 3.9 | 780 | 3.4 | 852 | 3.7 |
| FN (on it, seeded on it, from its genomes) | 60,205 | 66.2 | 9,897 | 12.3 | 19,304 | 22.5 |
| unseen (seeded on it, from its genomes) | 752,729 | 436.5 | 114,579 | 64.3 | 198,358 | 114.5 |

By fragment, the same sample: 610,774 fragments from an FN or unseen genome that aligned nowhere (160 MB), 42,544 that
aligned elsewhere (174 MB), 118,129 that seeded on an error taxon (107 MB) and 4,426 on an FP or FN taxon (21 MB). In a
PacBio soil sample, 80% of the records are unmapped and 89% were taken only for their source. The one analysis of the
records in the v15 report (`own_seed_check.py`: did an unseen species' own reads seed on its genes?) needs a count,
which the taxa tables now give (`own_seeded_not_aligned`).

**The tables.** On 30 MB of v15's pe training table: gzip -6 8.32 MB (3.6 s), gzip -9 8.16 MB (10.9 s), zstd -19
5.97 MB (54 s), xz -6 5.92 MB (44 s); numbers to 9 significant digits, then gzip -6: 6.02 MB; to 7 digits 5.35 MB. Nine
digits keep every float32 (what a forest compares) and far more than the features' precision, at gzip's speed.

## What changed

- `error_reads.py`: one SAM per kind of error (`<sample>.FP.sam.zst`, `<sample>.FN.sam.zst`; `<sample>.unseen.sam.zst`
  only with `--sams FP,FN,unseen`; `--sams none` for none); at most `--max-fragments` (20) fragments per taxon and reason,
  the lowest CRC-32 of the read name (a smaller cap keeps a subset of a larger one's); QUAL `*` unless `--qualities`;
  `--heldout` takes the held-out species out of "unseen"; the taxa tables count every fragment as before, with
  `own_seeded_not_aligned` added; the summary has `sam_fragments`, `sam_records` and `sams` (the files) instead of `sam`.
- `build_gtdb_database.py --share-logs` (off by default): the error reads' FP and FN SAMs are written (without it, the
  taxa tables and summaries only), and the run ends by packing `OUTDIR/<name>_share.tar.gz`: every `*.log` of OUTDIR,
  `console.log` (new: every line the build says, for every run, after a line with the time and command),
  `model_logs/`, and the training and test tables to 9 digits, all under `<name>/`. Each file once; the database, the
  models and `internal_taxonomy.dmp` are left out. The build passes `heldout_species.txt` to `error_reads.py`.
- Tests: `test_error_reads.py` (7, from 4: the files per kind, the cap, `--sams`, `--qualities`, `--heldout`);
  `test_gtdb_pipeline.py`'s `test_f_scenarios` builds with `--share-logs` and checks the SAMs per kind and the
  archive (its members, the tables' rows and values within 1e-8).

## Expected sizes

| Part | v15 | With `--share-logs` | Without |
|---|---|---|---|
| error reads | ~9 GB | the Nanopore soil sample above 15.7 MB uncompressed with QUAL, about 2% of its records' bytes once QUAL goes; all read types an estimated 100-200 MB | the taxa tables, an estimated 1 MB gzipped per read type (Nanopore 4.6 MB with the held-out species) |
| tables | 390 MB gzipped | ~280 MB | (not packed) |
| `model_logs/` rest and the logs | ~60 MB (twice) | ~60 MB | |
| total to copy | 2.4 GB + 443 MB (and 9 GB uncut) | about 0.45-0.55 GB, one file | |

The error-read estimate is extrapolated from one Nanopore sample; the next build's console line gives the size. On a
30 MB table slice the archive step took 7 s on one core (35 MB of tables to 26 MB, 7 MB packed), so the 1.3 GB of a
GTDB build's tables take one to two minutes at its threads.

## Verification and what is open

- `test_error_reads.py`: 7 tests pass (WSL, one core). `share_archive` ran on a copy of 30 MB of v15's pe table and
  5 MB of its se test table: text columns identical, numbers shortened, `.partial` files and the work folder left out.
- Not run, to spare the machine: `test_gtdb_pipeline.py` (`test_f_scenarios` builds a mini GTDB with scenarios).
- The held-out species stay in "unseen" of earlier builds' tables; `error_taxa_summary.py --heldout` of the v15
  report splits them off.

Commands (WSL):

    zstd -dc ~/v15/model_logs/error_reads/ont/test/sc_soil_ont_b6000000000/sc_soil_ont_b6000000000_s_2.sam.zst \
        | perl caps.pl ~/v15/model_logs/error_reads/ont/test/sc_soil_ont_b6000000000/sc_soil_ont_b6000000000_s_2.taxa.tsv 20,50
    PROTAL_TESTS_REQUIRED=1 python3 -m unittest scripts/test_error_reads.py
