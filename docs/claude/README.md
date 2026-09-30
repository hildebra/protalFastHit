# Audits and benchmarks written with Claude

This folder keeps the reports that Claude (Claude Code) wrote while working on protal: code
audits, reviews, benchmarks and experiments. They are records of what was checked, with which
commit, data and numbers. They are not user documentation; that is in [`docs/`](../README.md)
and on [protal.earlham.ac.uk](https://protal.earlham.ac.uk/main.php?site=documentation).

| Date | Report | Scope |
|---|---|---|
| 2026-09-30 | [Loading the database and reads, writing the output](2026-09-30-load-and-output/README.md) | start-up, FASTQ input and SAM output after the performance work: per-gene and per-value load costs measured and extrapolated to GTDB r226 size, the reader lock and zlib inflate, pigz against in-process libdeflate, the all-genes `@SQ` header, profiling's SAM parsing; opportunities ranked. Follow-up: `.sam.zst` and in-process BGZF `.sam.gz` output, headers listing only the genes aligned to, zstd by default, libdeflate as the gzip handler, implemented and measured; gzip or zstd for the strain MSAs, evaluated; the gene tables parsed in parallel, the gene arena and the FASTQ reader lock (#3–#5), implemented and measured |
| 2026-09-29 | [Performance profiling](2026-09-29-performance-profiling/README.md) | where time and instructions go on simulated paired reads (64- and 900-species worlds, realistic 5% mix): gzip reader lock, syncmers, WFA, huge pages, start-up; cPMML exception patches, a threaded gzip reader, huge pages and fixed stage timers, `--x_drop`, a branch-free (and AVX2) syncmer scan, and alignment from the anchors' exact matches, implemented and measured |
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
