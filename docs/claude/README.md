# Audits and benchmarks written with Claude

This folder keeps the reports that Claude (Claude Code) wrote while working on protal: code
audits, reviews, benchmarks and experiments. They are records of what was checked, with which
commit, data and numbers. They are not user documentation; that is in [`docs/`](../README.md)
and on [protal.earlham.ac.uk](https://protal.earlham.ac.uk/main.php?site=documentation).

| Date | Report | Scope |
|---|---|---|
| 2026-09-30 | [Clade holdouts for training](2026-09-30-clade-holdouts.md) | the training database without whole families, classes and phyla as well as single species; attribution of false positives by rank; cross-validation by clade; results on the tuning world |
| 2026-09-30 | [Faster index builds](2026-09-30-index-build-gains/README.md) | what else costs time in `protal --build` at GTDB scale: an unused pass, the unique k-mer statistics, parts that grow faster than the database; measured on the tuning world at three sizes, with exact fixes and estimates for GTDB r226 |
| 2026-09-29 | [Parallelising the index build](2026-09-29-index-build-parallel.md) | where `protal --build` spends its time (measured, and extrapolated to GTDB r226), why the index build is serial, and a parallel design that keeps the index byte-identical; follow-up: the uniqueness check without its lock, 5x faster on 8 threads |
| 2026-09-29 | [Training the presence model](2026-09-29-model-training-tuning/README.md) | the Java-free trainer against the old procedure; what simulated training data need (strains, species held out) on a GTDB-like world, and the settings for GTDB r226; placeholder models; plan for se/PacBio/ONT training data |
| 2026-09-29 | [Single-end and long reads](2026-09-29-single-end-and-long-reads.md) | what single-end, HiFi and ONT reads need (evaluation), and single-end support as implemented and tested |
| 2026-09-29 | [Website documentation review](2026-09-29-website-review.md) | protal.earlham.ac.uk against branch `audit-fixes`: what to update |
| 2026-09-26 | [Code audit](2026-09-26-code-audit.md) | three rounds: read mapping and SAM output, integration contracts, the profiling stage; fixes and measurements |

## Conventions

- One Markdown file per report, named `YYYY-MM-DD-<topic>.md` (the date the work was done), e.g.
  `2026-10-02-index-load-benchmark.md`. A report with its own scripts or tables gets a folder of
  the same name instead, with the report as its `README.md`.
- Start with the date, the branch and commit it describes, and the data it used (database,
  simulated or real samples, thread counts, machine). Numbers without their setup cannot be
  compared later.
- Record what was run: the commands, or the scripts next to the report. Keep generated data out
  of git unless it is small; say where it lives and how to regenerate it.
- A later follow-up adds a section or a new file that links to the old one. Do not rewrite old
  results; they describe the code as it was.
- Add each new report to the table above, newest first.
