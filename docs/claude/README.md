# Audits and benchmarks written with Claude

This folder keeps the reports that Claude (Claude Code) wrote while working on protal: code
audits, reviews, benchmarks and experiments. They are records of what was checked, with which
commit, data and numbers. They are not user documentation; that is in [`docs/`](../README.md)
and on [protal.earlham.ac.uk](https://protal.earlham.ac.uk/main.php?site=documentation).

| Date | Report | Scope |
|---|---|---|
| 2026-09-29 | [Strain MSAs and trees](2026-09-29-strain-audit/README.md) | audit round 4: how the strain MSAs are built, what qcmsa filters out, and the hand-off to IQ-TREE, measured against true genotypes and known trees; findings and fix plan |
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
